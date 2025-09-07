%% Time Setup
startTime = datetime(2020,5,11,17,35,38);
stopTime = startTime + days(2);
sampleTime = 30; % seconds
timeVector = startTime:seconds(sampleTime):stopTime;
numTimeSteps = length(timeVector);

%% Collision Criteria
collisionThreshold = 10.0; % km
closeApproachThreshold = 20.0; % km

%% Earth Parameters
earthRadius = 6.378137e6; % meters
mu_earth = 3.986004418e14;

xmlFolder = '/home/roboticsadmin/MATLABProjects/Debris_Folder';
ommFiles = dir(fullfile(xmlFolder, '*.xml'));

fprintf('Reading space debris data from folder: %s\n', xmlFolder);
fprintf('Found %d XML files.\n', length(ommFiles));

numSamples = 1000;
lat = 10;
lon = -50;

sc_mean = satelliteScenario(startTime, stopTime, sampleTime);

%% Load ISS Data
ISS_omm_file = 'ISS.xml';
ISS_data = readOMMFile(ISS_omm_file);
fprintf('  Loaded: %s (ID: %s)\n', ISS_data.objectName, ISS_data.objectId);

[ISS_nominal, ISS_uncertainty, P_ISS] = convertOMMToEquinoctialElements(ISS_data);

%% ISS Gaussian Sum Filter Parameters
mu = [ISS_nominal.a;
    ISS_nominal.h;
    ISS_nominal.k;
    ISS_nominal.p;
    ISS_nominal.q;
    ISS_nominal.l];

P = P_ISS;

fprintf('\n=== Processing Multiple N Values ===\n');

[Weights, Particle_Positions, Particle_Velocities, Keplerian_Particles, Equinoctial_Particles] = GenerateGSF(mu, P, numSamples, startTime, stopTime, sc_mean);
%Weights = Weights';
Particle_ECI = [Particle_Positions; Particle_Velocities];

numParticles = numSamples;

% Pre-allocate the array for parallel processing
next_Particle_ECI = zeros(6, numTimeSteps, numParticles);

% Parallel loop for particle propagation
fprintf('Starting parallel particle propagation...\n');
parfor i = 1:numParticles
    % Use the i-th particle's orbital elements
    Keplerian_Particles_i = equinoctialToKeplerian(Equinoctial_Particles(:, i));
    
    % Create a new scenario for this particle (each worker needs its own)
    sc_local = satelliteScenario(startTime, stopTime, sampleTime);
    
    sat = satellite(sc_local, Keplerian_Particles_i(1), Keplerian_Particles_i(2), ...
        Keplerian_Particles_i(3), Keplerian_Particles_i(4), ...
        Keplerian_Particles_i(5), Keplerian_Particles_i(6), "OrbitPropagator", "two-body-keplerian");
    
    % Get states for the new scenario (this will start from current_time)
    [positions, velocities, times] = states(sat);
    
    % Store the propagated states
    next_Particle_ECI(:, :, i) = [positions; velocities];
end
fprintf('Parallel particle propagation completed.\n');

Particle_weighted_mean = zeros(6, numTimeSteps);
Particle_weighted_covariance = zeros(6, 6, numTimeSteps);
Particle_weighted_skewness = zeros(6, numTimeSteps);
Particle_weighted_kurtosis = zeros(6, numTimeSteps);

% The statistics calculation loop - keeping as regular for loop since it's sequential
fprintf('Computing weighted statistics...\n');
for t = 1:numTimeSteps
    % Extract particles at time step t: 6 x numParticles
    particles_t = squeeze(next_Particle_ECI(:, t, :)); % 6 x numParticles
    
    % Ensure weights are normalized and in correct shape
    w = ones(1, numSamples) / numSamples;
    
    % Weighted mean: 6 x 1 vector
    Particle_weighted_mean(:, t) = particles_t * w';
    
    % Center the particles
    particles_centered = particles_t - Particle_weighted_mean(:, t);
    
    % Weighted covariance: 6 x 6 matrix
    Particle_weighted_covariance(:, :, t) = (particles_centered .* w) * particles_centered' / (1 - sum(w.^2));
    
    % Weighted skewness and kurtosis for each component
    for j = 1:6
        x = particles_centered(j, :); % 1 x numParticles - already centered!
        
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
            
            % Correct excess kurtosis calculation (subtract 3 for excess kurtosis)
            % or keep as-is for raw kurtosis
            Particle_weighted_kurtosis(j, t) = mu4 / (mu2^2); % Raw kurtosis
            % For excess kurtosis, use: mu4 / (mu2^2) - 3;
        end
    end
end

% Extract mean positions and covariances
Particle_ISS_positions = Particle_weighted_mean(1:3, :);
Particle_ISS_position_covs = Particle_weighted_covariance(1:3, 1:3, :);
variances = squeeze([Particle_ISS_position_covs(1,1,end); ...
                     Particle_ISS_position_covs(2,2,end); ...
                     Particle_ISS_position_covs(3,3,end)]);

% Convert to uncertainties (standard deviations)
uncertainties = sqrt(variances);

% ---------------------------------------------------------
% Uncertainty in the MEAN (standard error of the mean)
% If you used N particles:
N = size(Particle_weighted_mean, 2);  % number of MC samples
mean_cov = Particle_ISS_position_covs(:,:,end) / N;   % covariance of the mean
mean_std = sqrt(diag(mean_cov));      % SE per coordinate
% ---------------------------------------------------------

disp('Mean position at final time:');
disp(Particle_ISS_positions(:, end));

disp('1-sigma ensemble spread (m):');
disp(uncertainties);

disp('1-sigma uncertainty in the mean (standard error, m):');
disp(mean_std);

ISS_positions_km = Particle_ISS_positions / 1000;
earthRadius_km = earthRadius / 1000;
time_hours = hours(timeVector - startTime);

% Debug: Check dimensions
fprintf('Updated_Particles dimensions: [%d, %d, %d]\n', size(next_Particle_ECI));
fprintf('Expected: [6, %d, %d]\n', numTimeSteps, numSamples);

% Position magnitude plot
figure('Position', [100, 100, 1400, 1000]);
subplot(2, 3, 1);
position_magnitude = sqrt(sum(ISS_positions_km.^2, 1));
plot(time_hours, position_magnitude, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Position Magnitude');
hold on;

% Calculate uncertainty bounds - this inner loop could also be parallelized if needed
std_dev_magnitude = zeros(1, numTimeSteps);
for i = 1:numTimeSteps
    combined_variance = sum(diag(Particle_ISS_position_covs(:, :, i))) / (1000^2);
    std_dev_magnitude(i) = sqrt(combined_variance);
end

% Add uncertainty bounds
upper_bound = position_magnitude + std_dev_magnitude;
lower_bound = position_magnitude - std_dev_magnitude;
fill([time_hours, fliplr(time_hours)], [upper_bound, fliplr(lower_bound)], ...
    'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

xlabel('Time (hours)');
ylabel('Position Magnitude (km)');
title('ISS Position Magnitude vs Time (±1σ)');
legend('Location', 'best');
grid on;

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
function [Weights, positions, velocities, Keplerian_Particles, Equinoctial_Particles] = GenerateGSF(mu, P, numSamples, startTime, stopTime, sc_mean)
    sampleTime = 30; % seconds
    timeVector = startTime:seconds(sampleTime):stopTime;
    numTimeSteps = length(timeVector);
    data = mvnrnd(mu', P, numSamples);  % Note: mu' to make it row vector
    Equinoctial_Particles = data';
    %Weights = mvnpdf(data, mu', P);
    %Weights = weights/sum(weights);
    Weights = ones(1, numSamples) / numSamples;
    % Convert each sample to Keplerian
    Keplerian_Particles = zeros(6, numSamples);
    positions = zeros(3, numSamples);
    velocities = zeros(3, numSamples);
    for sample = 1:numSamples
        data_keplerian = equinoctialToKeplerian(data(sample, :)');
        Keplerian_Particles(:, sample) = data_keplerian;
        [positions(:, sample), velocities(:, sample)] = keplerian2ijk(data_keplerian(1), data_keplerian(2), data_keplerian(3), data_keplerian(4), data_keplerian(5), data_keplerian(6));
    end
   
    
end














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


function [nominal, uncertainty, P] = convertOMMToEquinoctialElements(ommData)
    nominal = struct();
    uncertainty = struct();
    
    % Inputs from OMM
    a = ommData.semiMajorAxis * 1000;
    ecc = ommData.eccentricity;
    incl = deg2rad(ommData.inclination);
    raan = deg2rad(ommData.raOfAscNode);
    argp = deg2rad(ommData.argOfPericenter);
    
    % Mean anomaly to true anomaly
    M = deg2rad(ommData.meanAnomaly);
    E = M;
    for i = 1:10
        E = E - (E - ecc*sin(E) - M) / (1 - ecc*cos(E));  % Newton-Raphson
    end
    nu = 2 * atan(sqrt((1 + ecc)/(1 - ecc)) * tan(E/2));

    % Equinoctial elements
    nominal.a = a;
    nominal.h = ecc * sin(raan + argp);
    nominal.k = ecc * cos(raan + argp);
    nominal.p = tan(incl / 2) * sin(raan);
    nominal.q = tan(incl / 2) * cos(raan);
    nominal.l = mod(nu + raan + argp, 2*pi);
    %{
    perc.a = 0.0002;   % 5% of a
    perc.h = 0.006;   % 1% of h
    perc.k = 0.006;   % 1% of k
    perc.p = 0.006;   % 1% of p
    perc.q = 0.006;   % 1% of q
    perc.l = 0.00006;  % 0.5% of l
    %}
    
    perc.a = 0.00005;   % 5% of a
    perc.h = 0.00005;   % 1% of h
    perc.k = 0.00005;   % 1% of k
    perc.p = 0.00005;   % 1% of p
    perc.q = 0.00005;   % 1% of q
    perc.l = 0.00005;  % 0.5% of l
    
    %{
    perc.a = 0.0001;   % 5% of a
    perc.h = 0.005;   % 1% of h
    perc.k = 0.005;   % 1% of k
    perc.p = 0.005;   % 1% of p
    perc.q = 0.005;   % 1% of q
    perc.l = 0.00001;  % 0.5% of l
    %}
    % Calculate uncertainties as percentages of nominal values
    uncertainty.a = perc.a * nominal.a;
    uncertainty.h = perc.h * nominal.h;
    uncertainty.k = perc.k * nominal.k;
    uncertainty.p = perc.p * nominal.p;
    uncertainty.q = perc.q * nominal.q;
    uncertainty.l = perc.l * nominal.l;
    

    %{
    uncertainty.a = 1000;
    uncertainty.h = 1e-3;
    uncertainty.k = 1e-3;
    uncertainty.p = 1e-3;
    uncertainty.q = 1e-3;
    uncertainty.l = 1e-2*pi/180;
    %}

    P = diag([uncertainty.a^2, ...
              uncertainty.h^2, ...
              uncertainty.k^2, ...
              uncertainty.p^2, ...
              uncertainty.q^2, ...
              uncertainty.l^2]);
   
end


function keplerian = equinoctialToKeplerian(equinoctial)

% Extract equinoctial elements
a      = equinoctial(1);
h      = equinoctial(2);
k      = equinoctial(3);
p      = equinoctial(4);
q      = equinoctial(5);
lambda = equinoctial(6);

e = sqrt(h^2 + k^2);

i = 2 * atan( sqrt(p^2 + q^2) );

Omega = atan2(p, q);
Omega = mod(Omega, 2*pi);

omega = atan2(h, k) - Omega;
omega = mod(omega, 2*pi);

nu = lambda - Omega - omega;
nu = mod(nu, 2*pi);

E = 2 * atan( tan(nu / 2) / sqrt((1 + e) / (1 - e)) );
E = mod(E, 2*pi);

M = E - e * sin(E);
M = mod(M, 2*pi);

i     = rad2deg(i);
Omega = rad2deg(Omega);
omega = rad2deg(omega);
M     = rad2deg(M);

keplerian = [a, e, i, Omega, omega, M];

end

