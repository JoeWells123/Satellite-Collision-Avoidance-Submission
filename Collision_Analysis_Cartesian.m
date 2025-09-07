startTime = datetime(2020,5,11,17,35,38);
stopTime = startTime + days(2);
sampleTime = 30; % seconds
timeVector = startTime:seconds(sampleTime):stopTime;
numTimeSteps = length(timeVector);

numSamples = 1000; 
fprintf('Running collision risk Monte Carlo with %d samples...\n', numSamples);

collisionThreshold = 10.0; % km
closeApproachThreshold = 20.0; % km

earthRadius = 6.378137e6; % meters
mu_earth = 3.986004418e14; % m³/s²

%% Load ISS Data
ISS_omm_file = 'ISS.xml';
ISS_data = readOMMFile(ISS_omm_file);
fprintf('  Loaded: %s (ID: %s)\n', ISS_data.objectName, ISS_data.objectId);

[ISS_nominal_cartesian, ISS_uncertainty_cartesian, P] = convertOMMToCartesian(ISS_data);

mu = [ISS_nominal_cartesian.position(1);
    ISS_nominal_cartesian.position(2);
    ISS_nominal_cartesian.position(3);
    ISS_nominal_cartesian.velocity(1);
    ISS_nominal_cartesian.velocity(2);
    ISS_nominal_cartesian.velocity(3)];

%% Initialize Parallel Pool
if isempty(gcp('nocreate'))
    fprintf('Starting parallel pool...\n');
    parpool();
else
    fprintf('Using existing parallel pool...\n');
end

%% Generate all samples beforehand (outside parfor)
fprintf('Generating Monte Carlo samples...\n');
ISS_ic_samples = mvnrnd(mu, P, numSamples)'; % 6 x numSamples

% Initialize result arrays for parfor compatibility
ISS_positions_par = zeros(numSamples, 3, numTimeSteps);
ISS_velocities_par = zeros(numSamples, 3, numTimeSteps);
storedICs_ISS = cell(numSamples, 1);

% Progress tracking
fprintf('Starting parallel Monte Carlo simulation with %d samples...\n', numSamples);

%% Parallel Monte Carlo Simulation Loop
tic
parfor sample = 1:numSamples
    try
        % Get the pre-generated sample for this worker
        ISS_ic_sample = ISS_ic_samples(:, sample);
        storedICs_ISS{sample} = ISS_ic_sample;
        
        % Each worker creates its own satellite scenario
        sc_local = satelliteScenario(startTime, stopTime, sampleTime);
        
        % Convert Cartesian to Keplerian elements
        [a, ecc, incl, RAAN, argp, nu, ~, ~, ~] = ijk2keplerian(...
            ISS_ic_sample(1:3), ISS_ic_sample(4:6));
        
        % Create satellite with local scenario
        ISS_temp = satellite(sc_local, a, ecc, incl, RAAN, argp, nu, ...
            "OrbitPropagator", "sgp4");
        
        % Get states
        [ISS_pos_sample, ISS_vel_sample, ~] = states(ISS_temp);
        
        % Store results
        ISS_positions_par(sample, :, :) = ISS_pos_sample;
        ISS_velocities_par(sample, :, :) = ISS_vel_sample;
        
    catch ME
        % Handle any errors that occur in individual workers
        fprintf('Error in sample %d: %s\n', sample, ME.message);
        % Fill with NaN to indicate failure
        ISS_positions_par(sample, :, :) = NaN;
        ISS_velocities_par(sample, :, :) = NaN;
    end
end

simulation_time = toc;
fprintf('Parallel simulation completed in %.2f seconds\n', simulation_time);

%% Check for any failed samples
failed_samples = any(any(isnan(ISS_positions_par), 2), 3);
num_failed = sum(failed_samples);
if num_failed > 0
    fprintf('Warning: %d samples failed and contain NaN values\n', num_failed);
end

%% Rearrange data for statistical analysis
% Concatenate positions and velocities: [numSamples, 6, numTimeSteps]
next_Particle_ECI = cat(2, ISS_positions_par, ISS_velocities_par); 

% Rearrange dimensions to [6, numTimeSteps, numSamples]
next_Particle_ECI = permute(next_Particle_ECI, [2, 3, 1]);

%% Statistical Analysis
fprintf('Computing statistical moments...\n');

Particle_weighted_mean = zeros(6, numTimeSteps);
Particle_weighted_covariance = zeros(6, 6, numTimeSteps);
Particle_weighted_skewness = zeros(6, numTimeSteps);
Particle_weighted_kurtosis = zeros(6, numTimeSteps);

tic
for t = 1:numTimeSteps
    % Extract particles at time step t: 6 x numParticles
    particles_t = squeeze(next_Particle_ECI(:, t, :)); % 6 x numParticles
    
    % Remove failed samples (containing NaN)
    valid_mask = ~any(isnan(particles_t), 1);
    particles_t = particles_t(:, valid_mask);
    num_valid = sum(valid_mask);
    
    if num_valid > 1
        % Equal weights for all valid samples
        w = ones(1, num_valid) / num_valid;
        
        % Weighted mean: 6 x 1 vector
        Particle_weighted_mean(:, t) = particles_t * w';
        
        % Center the particles
        particles_centered = particles_t - Particle_weighted_mean(:, t);
        
        % Weighted covariance: 6 x 6 matrix
        if num_valid > 1
            Particle_weighted_covariance(:, :, t) = (particles_centered .* w) * particles_centered' / (1 - sum(w.^2));
        end
        
        % Weighted skewness and kurtosis for each component
        for j = 1:6
            x = particles_centered(j, :); % 1 x num_valid - already centered!
            
            % Calculate weighted moments correctly
            mu2 = sum(w .* (x.^2));       % 2nd central moment (variance)
            mu3 = sum(w .* (x.^3));       % 3rd central moment
            mu4 = sum(w .* (x.^4));       % 4th central moment
            
            % Check for near-zero variance to avoid division by zero
            if mu2 < eps
                Particle_weighted_skewness(j, t) = 0;
                Particle_weighted_kurtosis(j, t) = 3; % Normal distribution kurtosis
            else
                % Correct skewness calculation
                Particle_weighted_skewness(j, t) = mu3 / (mu2^(3/2));
                
                % Correct kurtosis calculation (raw kurtosis)
                Particle_weighted_kurtosis(j, t) = mu4 / (mu2^2);
            end
        end
    else
        % Not enough valid samples for statistics
        Particle_weighted_mean(:, t) = NaN;
        Particle_weighted_covariance(:, :, t) = NaN;
        Particle_weighted_skewness(:, t) = NaN;
        Particle_weighted_kurtosis(:, t) = NaN;
    end
end

stats_time = toc;
fprintf('Statistical analysis completed in %.2f seconds\n', stats_time);
fprintf('Total processing time: %.2f seconds\n', simulation_time + stats_time);

%% Summary
fprintf('\n=== Simulation Summary ===\n');
fprintf('Total samples: %d\n', numSamples);
fprintf('Failed samples: %d\n', num_failed);
fprintf('Successful samples: %d\n', numSamples - num_failed);
fprintf('Time steps: %d\n', numTimeSteps);
fprintf('Simulation duration: %.1f hours\n', hours(stopTime - startTime));


fprintf('Generating skewness and kurtosis plots...\n');

%% Convert time vector to duration for plotting
time_hours = hours(timeVector - startTime);

%% Create figure with subplots for skewness
figure('Position', [100, 100, 1200, 800]);

% Skewness plots
subplot(2, 2, 1);
plot(time_hours, Particle_weighted_skewness(1,:), 'b-', 'LineWidth', 1.5); hold on;
plot(time_hours, Particle_weighted_skewness(2,:), 'r-', 'LineWidth', 1.5);
plot(time_hours, Particle_weighted_skewness(3,:), 'g-', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Skewness');
title('Position Components - Skewness vs Time');
legend('X Position', 'Y Position', 'Z Position', 'Location', 'best');
grid on;

subplot(2, 2, 2);
plot(time_hours, Particle_weighted_skewness(4,:), 'b-', 'LineWidth', 1.5); hold on;
plot(time_hours, Particle_weighted_skewness(5,:), 'r-', 'LineWidth', 1.5);
plot(time_hours, Particle_weighted_skewness(6,:), 'g-', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Skewness');
title('Velocity Components - Skewness vs Time');
legend('X Velocity', 'Y Velocity', 'Z Velocity', 'Location', 'best');
grid on;

% Kurtosis plots
subplot(2, 2, 3);
plot(time_hours, Particle_weighted_kurtosis(1,:), 'b-', 'LineWidth', 1.5); hold on;
plot(time_hours, Particle_weighted_kurtosis(2,:), 'r-', 'LineWidth', 1.5);
plot(time_hours, Particle_weighted_kurtosis(3,:), 'g-', 'LineWidth', 1.5);
yline(3, 'k--', 'Normal Distribution', 'LineWidth', 1); % Reference line for normal distribution
xlabel('Time (hours)');
ylabel('Kurtosis');
title('Position Components - Kurtosis vs Time');
legend('X Position', 'Y Position', 'Z Position', 'Normal Ref', 'Location', 'best');
grid on;

subplot(2, 2, 4);
plot(time_hours, Particle_weighted_kurtosis(4,:), 'b-', 'LineWidth', 1.5); hold on;
plot(time_hours, Particle_weighted_kurtosis(5,:), 'r-', 'LineWidth', 1.5);
plot(time_hours, Particle_weighted_kurtosis(6,:), 'g-', 'LineWidth', 1.5);
yline(3, 'k--', 'Normal Distribution', 'LineWidth', 1); % Reference line for normal distribution
xlabel('Time (hours)');
ylabel('Kurtosis');
title('Velocity Components - Kurtosis vs Time');
legend('X Velocity', 'Y Velocity', 'Z Velocity', 'Normal Ref', 'Location', 'best');
grid on;

sgtitle('ISS Orbit Uncertainty: Skewness and Kurtosis Evolution', 'FontSize', 14, 'FontWeight', 'bold');

%% Create a combined overview plot
figure('Position', [150, 150, 1200, 600]);

subplot(1, 2, 1);
% Plot all skewness components
plot(time_hours, Particle_weighted_skewness(1,:), 'b-', 'LineWidth', 1.2); hold on;
plot(time_hours, Particle_weighted_skewness(2,:), 'r-', 'LineWidth', 1.2);
plot(time_hours, Particle_weighted_skewness(3,:), 'g-', 'LineWidth', 1.2);
plot(time_hours, Particle_weighted_skewness(4,:), 'c-', 'LineWidth', 1.2);
plot(time_hours, Particle_weighted_skewness(5,:), 'm-', 'LineWidth', 1.2);
plot(time_hours, Particle_weighted_skewness(6,:), 'y-', 'LineWidth', 1.2);
yline(0, 'k--', 'Symmetric', 'LineWidth', 1); % Reference line for symmetric distribution
xlabel('Time (hours)');
ylabel('Skewness');
title('All Components - Skewness Evolution');
legend('X Pos', 'Y Pos', 'Z Pos', 'X Vel', 'Y Vel', 'Z Vel', 'Symmetric Ref', 'Location', 'best');
grid on;

subplot(1, 2, 2);
% Plot all kurtosis components
plot(time_hours, Particle_weighted_kurtosis(1,:), 'b-', 'LineWidth', 1.2); hold on;
plot(time_hours, Particle_weighted_kurtosis(2,:), 'r-', 'LineWidth', 1.2);
plot(time_hours, Particle_weighted_kurtosis(3,:), 'g-', 'LineWidth', 1.2);
plot(time_hours, Particle_weighted_kurtosis(4,:), 'c-', 'LineWidth', 1.2);
plot(time_hours, Particle_weighted_kurtosis(5,:), 'm-', 'LineWidth', 1.2);
plot(time_hours, Particle_weighted_kurtosis(6,:), 'y-', 'LineWidth', 1.2);
yline(3, 'k--', 'Normal Distribution', 'LineWidth', 1); % Reference line for normal distribution
xlabel('Time (hours)');
ylabel('Kurtosis');
title('All Components - Kurtosis Evolution');
legend('X Pos', 'Y Pos', 'Z Pos', 'X Vel', 'Y Vel', 'Z Vel', 'Normal Ref', 'Location', 'best');
grid on;

sgtitle('ISS Orbit Uncertainty: Complete Statistical Evolution', 'FontSize', 14, 'FontWeight', 'bold');

%% Optional: Create statistical summary
fprintf('\n=== Statistical Summary ===\n');
fprintf('Final skewness values:\n');
fprintf('  Position (X,Y,Z): %.4f, %.4f, %.4f\n', Particle_weighted_skewness(1,end), Particle_weighted_skewness(2,end), Particle_weighted_skewness(3,end));
fprintf('  Velocity (X,Y,Z): %.4f, %.4f, %.4f\n', Particle_weighted_skewness(4,end), Particle_weighted_skewness(5,end), Particle_weighted_skewness(6,end));

fprintf('Final kurtosis values:\n');
fprintf('  Position (X,Y,Z): %.4f, %.4f, %.4f\n', Particle_weighted_kurtosis(1,end), Particle_weighted_kurtosis(2,end), Particle_weighted_kurtosis(3,end));
fprintf('  Velocity (X,Y,Z): %.4f, %.4f, %.4f\n', Particle_weighted_kurtosis(4,end), Particle_weighted_kurtosis(5,end), Particle_weighted_kurtosis(6,end));

fprintf('\nSkewness interpretation:\n');
fprintf('  > 0: Right-skewed (tail extends toward positive values)\n');
fprintf('  < 0: Left-skewed (tail extends toward negative values)\n');
fprintf('  ≈ 0: Approximately symmetric\n');

fprintf('\nKurtosis interpretation:\n');
fprintf('  > 3: Heavy-tailed (more extreme values than normal distribution)\n');
fprintf('  < 3: Light-tailed (fewer extreme values than normal distribution)\n');
fprintf('  ≈ 3: Similar to normal distribution\n');

fprintf('\nPlots generated successfully!\n');

%{
%% Times (in hours) to plot
times_hours = [0, 1.75, 3.5, 5.25, 7.02, 8.8, 10.6];
orbits_labels = [2, 3, 4, 5, 6, 7, 8]; % orbit counts corresponding to times_hours
numSamplesToPlot = 100;
padding = 10; % km padding around point cloud

% Select samples once for all times to keep consistent samples
sampleIndices = randperm(numSamples, numSamplesToPlot);

all_points = []; % to collect all points for zooming

% Colors for each time
colors = lines(length(times_hours));

figure('Name','ISS Monte Carlo Samples at Multiple Times (Single Plot)','Position',[200 200 900 700]);
hold on; grid on; axis equal;

for idx = 1:length(times_hours)
    % Calculate target index for each time
    targetTime_sec = times_hours(idx) * 3600;
    targetIndex = round(targetTime_sec / sampleTime) + 1;
    
    % Extract sample points at this time (km)
    ISS_sample_points = squeeze(ISS_positions_par(sampleIndices, :, targetIndex)) / 1000;
    all_points = [all_points; ISS_sample_points]; % append for zoom limits
    
    % Extract mean position at this time
    mean_at_target = Particle_weighted_mean(:, targetIndex) / 1000;
    
    % Plot samples with orbit label
    scatter3(ISS_sample_points(:,1), ISS_sample_points(:,2), ISS_sample_points(:,3), ...
        20, colors(idx,:), 'filled', 'MarkerFaceAlpha', 0.3, ...
        'DisplayName', sprintf('Samples @ %d orbits', orbits_labels(idx)));
    
end

% Set zoom limits around all points with padding
x_min = min(all_points(:,1)) - padding;
x_max = max(all_points(:,1)) + padding;
y_min = min(all_points(:,2)) - padding;
y_max = max(all_points(:,2)) + padding;
z_min = min(all_points(:,3)) - padding;
z_max = max(all_points(:,3)) + padding;

xlim([x_min x_max]);
ylim([y_min y_max]);
zlim([z_min z_max]);

% Labels and title
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Orbital Object Cartesian Distribution');
legend('Location', 'best');
view(3);
hold off;

%}
%% Helper Functions

function ommData = readOMMFile(filename)
    if ~exist(filename, 'file')
        error('OMM file not found: %s', filename);
    end
    
    try
        xmlDoc = xmlread(filename);
        ommData = struct();
    
        objectNameNodes = xmlDoc.getElementsByTagName('OBJECT_NAME');
        if objectNameNodes.getLength() > 0
            ommData.objectName = char(objectNameNodes.item(0).getFirstChild.getNodeValue);
        else
            ommData.objectName = 'Unknown';
        end
        
        objectIdNodes = xmlDoc.getElementsByTagName('OBJECT_ID');
        if objectIdNodes.getLength() > 0
            ommData.objectId = char(objectIdNodes.item(0).getFirstChild.getNodeValue);
        else
            ommData.objectId = 'Unknown';
        end
        
        epochNodes = xmlDoc.getElementsByTagName('EPOCH');
        if epochNodes.getLength() > 0
            epochStr = char(epochNodes.item(0).getFirstChild.getNodeValue);
            ommData.epoch = datetime(epochStr, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');
        else
            ommData.epoch = datetime('now');
        end
        ommData.meanMotion = getXMLNodeValue(xmlDoc, 'MEAN_MOTION', 0);
        ommData.eccentricity = getXMLNodeValue(xmlDoc, 'ECCENTRICITY', 0);
        ommData.inclination = getXMLNodeValue(xmlDoc, 'INCLINATION', 0);
        ommData.raOfAscNode = getXMLNodeValue(xmlDoc, 'RA_OF_ASC_NODE', 0);
        ommData.argOfPericenter = getXMLNodeValue(xmlDoc, 'ARG_OF_PERICENTER', 0);
        ommData.meanAnomaly = getXMLNodeValue(xmlDoc, 'MEAN_ANOMALY', 0);
        
        semiMajorAxisNodes = xmlDoc.getElementsByTagName('USER_DEFINED');
        ommData.semiMajorAxis = 0;
        for i = 0:semiMajorAxisNodes.getLength()-1
            node = semiMajorAxisNodes.item(i);
            if strcmp(char(node.getAttribute('parameter')), 'SEMIMAJOR_AXIS')
                ommData.semiMajorAxis = str2double(char(node.getFirstChild.getNodeValue));
                break;
            end
        end
        
        if ommData.semiMajorAxis == 0 && ommData.meanMotion > 0
            mu = 3.986004418e14; % m³/s²
            n = ommData.meanMotion * 2 * pi / 86400; % rad/s
            ommData.semiMajorAxis = (mu / n^2)^(1/3) / 1000; % km
        end
        
    catch ME
        error('Failed to parse OMM file %s: %s', filename, ME.message);
    end
end

function value = getXMLNodeValue(xmlDoc, tagName, defaultValue)
    nodes = xmlDoc.getElementsByTagName(tagName);
    if nodes.getLength() > 0
        value = str2double(char(nodes.item(0).getFirstChild.getNodeValue));
    else
        value = defaultValue;
    end
end

function [nominal_cartesian, uncertainty_cartesian, P] = convertOMMToCartesian(ommData)


    mu = 3.986004418e14; % m³/s²
    

    semiMajorAxis = ommData.semiMajorAxis * 1000; 
    eccentricity = ommData.eccentricity;
    disp(eccentricity)
    inclination = deg2rad(ommData.inclination);
    raan = deg2rad(ommData.raOfAscNode);
    argPeriapsis = deg2rad(ommData.argOfPericenter);
    
    
    M = deg2rad(ommData.meanAnomaly);
    E = M; 
    for i = 1:10
        E = M + eccentricity * sin(E);
    end
    trueAnomaly = 2 * atan(sqrt((1 + eccentricity)/(1 - eccentricity)) * tan(E/2));

   
    [r_ijk, v_ijk] = keplerian2ijk(semiMajorAxis, eccentricity, inclination, raan, argPeriapsis, trueAnomaly, mu);
    
   
    nominal_cartesian = struct();
    nominal_cartesian.position = r_ijk;
    nominal_cartesian.velocity = v_ijk;
    
    uncertainty_cartesian = struct();
    
    pos_uncertainty = 0.00005;
    vel_uncertainty = 0.00005;

    uncertainty_cartesian.position = [pos_uncertainty * r_ijk(1); pos_uncertainty* r_ijk(2); pos_uncertainty*r_ijk(3)];
    uncertainty_cartesian.velocity = [vel_uncertainty*v_ijk(1); vel_uncertainty*v_ijk(2); vel_uncertainty*v_ijk(3)];

    P = diag([uncertainty_cartesian.position(1)^2, ...
          uncertainty_cartesian.position(2)^2, ...
          uncertainty_cartesian.position(3)^2, ...
          uncertainty_cartesian.velocity(1)^2, ...
          uncertainty_cartesian.velocity(2)^2, ...
          uncertainty_cartesian.velocity(3)^2]);
end

function cartesian_ic = generateRandomCartesian(nominal_cartesian, uncertainty_cartesian)
    % Generate random Cartesian initial conditions with uncertainties
    cartesian_ic = struct();
    
    % Add Gaussian noise to position and velocity
    cartesian_ic.position = nominal_cartesian.position + randn(3,1) .* uncertainty_cartesian.position;
    cartesian_ic.velocity = nominal_cartesian.velocity + randn(3,1) .* uncertainty_cartesian.velocity;
end


function [r_ijk, v_ijk] = keplerian2ijk(a, ecc, incl, raan, argp, nu, mu)
    % Convert Keplerian orbital elements to Cartesian coordinates
    % Inputs:
    %   a    - semi-major axis [m]
    %   ecc  - eccentricity
    %   incl - inclination [rad]
    %   raan - right ascension of ascending node [rad]
    %   argp - argument of periapsis [rad]
    %   nu   - true anomaly [rad]
    %   mu   - gravitational parameter [m³/s²]
    % Outputs:
    %   r_ijk - position vector [m]
    %   v_ijk - velocity vector [m/s]
    
    % Calculate position and velocity in perifocal coordinate system
    p = a * (1 - ecc^2);  % semi-parameter
    r_pqw = [p * cos(nu) / (1 + ecc * cos(nu));
             p * sin(nu) / (1 + ecc * cos(nu));
             0];
    
    v_pqw = sqrt(mu / p) * [-sin(nu);
                            ecc + cos(nu);
                            0];
    
    % Rotation matrices
    R3_raan = [cos(raan), sin(raan), 0;
               -sin(raan), cos(raan), 0;
               0, 0, 1];
    
    R1_incl = [1, 0, 0;
               0, cos(incl), sin(incl);
               0, -sin(incl), cos(incl)];
    
    R3_argp = [cos(argp), sin(argp), 0;
               -sin(argp), cos(argp), 0;
               0, 0, 1];
    
    % Combined transformation matrix
    Q_pqw2ijk = R3_raan' * R1_incl' * R3_argp';
    
    % Transform to IJK coordinate system
    r_ijk = Q_pqw2ijk * r_pqw;
    v_ijk = Q_pqw2ijk * v_pqw;
end