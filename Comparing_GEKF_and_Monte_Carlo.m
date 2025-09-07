%% Time Setup
startTime = datetime(2020,5,11,17,35,38);
stopTime = startTime + days(2);
sampleTime = 120; % seconds
timeVector = startTime:seconds(sampleTime):stopTime;
numTimeSteps = length(timeVector);

%% Earth Parameters
earthRadius = 6.378137e6; % meters
mu_earth = 3.986004418e14;

xmlFolder = '/home/roboticsadmin/MATLABProjects/Debris_Folder';
ommFiles = dir(fullfile(xmlFolder, '*.xml'));

fprintf('Reading space debris data from folder: %s\n', xmlFolder);
fprintf('Found %d XML files.\n', length(ommFiles));

numSamples = 100;
lat = 10;
lon = -50;

sc_mean = satelliteScenario(startTime, stopTime, sampleTime);

%% Load ISS Data
ISS_omm_file = 'CZ-4C 2.xml';
ISS_data = readOMMFile(ISS_omm_file);
fprintf('  Loaded: %s (ID: %s)\n', ISS_data.objectName, ISS_data.objectId);

[ISS_nominal, ISS_uncertainty, P_ISS] = convertOMMToEquinoctialElements(ISS_data);

%% ISS EKF Parameters
mu = [ISS_nominal.a;
    ISS_nominal.h;
    ISS_nominal.k;
    ISS_nominal.p;
    ISS_nominal.q;
    ISS_nominal.l];

P = P_ISS;

%% Define different numbers of Gaussian components to test
N_values = [1, 20, 40, 60, 80, 100];

fprintf('\n=== Processing Multiple N Values with EKF ===\n');

% Store results for each N value
EKF_results = cell(length(N_values), 1);
bhattacharyya_distances = cell(length(N_values), 1);

% Generate reference using the highest N value (N=30)
ref_positions = Particle_ISS_positions(1:3, :);
ref_position_covs = Particle_ISS_position_covs(1:3, 1:3, :);

% Process each N value
for idx = 1:length(N_values)
    N = N_values(idx);
    fprintf('Processing N = %d Gaussians with EKF...\n', N);
    
    [weighted_covariance, weighted_mean, EKF_covariances, EKF_means] = GenerateEKF(mu, P, N, startTime, stopTime, mu_earth);
    
    % Store results
    EKF_results{idx}.weighted_covariance = weighted_covariance;
    EKF_results{idx}.weighted_mean = weighted_mean;
    EKF_results{idx}.N = N;
    
    % Extract position data
    ISS_positions = weighted_mean(1:3, 1:24:end);
    ISS_position_covs = weighted_covariance(1:3, 1:3, 1:24:end);
    
    % Calculate Bhattacharyya distance against reference
    d = zeros(1, numTimeSteps);
    for k = 1:numTimeSteps
        mu1 = ISS_positions(:, k);
        Sigma1 = ISS_position_covs(:, :, k);
        mu2 = ref_positions(:, k);
        Sigma2 = ref_position_covs(:, :, k);
        
        d(k) = bhattacharyya_distance(mu1, Sigma1, mu2, Sigma2);
    end
    
    bhattacharyya_distances{idx} = d;
    
    fprintf('  Mean Bhattacharyya distance: %.6f\n', mean(d));
    fprintf('  Max Bhattacharyya distance: %.6f\n', max(d));
end

%% Plotting
timeVector = datetime(timeVector, 'TimeZone', '');
time_hours = hours(timeVector - startTime);

% Convert reference positions to km for plotting
ref_positions_km = ref_positions / 1000;
earthRadius_km = earthRadius / 1000;

%% Enhanced Plot 1: Reference ISS Trajectory 
figure('Position', [100, 100, 1400, 1000]);
subplot(2, 1, 1);
position_magnitude = sqrt(sum(ref_positions_km.^2, 1)); % Euclidean distance
plot(time_hours, position_magnitude, 'k-', 'LineWidth', 2, 'DisplayName', 'Reference (N=30)');
hold on;

% Plot trajectories for different N values
for idx = 1:length(N_values)
    positions = EKF_results{idx}.weighted_mean(1:3, 1:24:end) / 1000; % Convert to km
    pos_magnitude = sqrt(sum(positions.^2, 1));
    disp(size(time_hours))
    disp(size(pos_magnitude))
    plot(time_hours, pos_magnitude, ...
         'LineWidth', 1.5, 'DisplayName', sprintf('N=%d', N_values(idx)));
end

% Add uncertainty bounds for reference
std_dev_magnitude = zeros(1, numTimeSteps);
for i = 1:numTimeSteps
    combined_variance = sum(diag(ref_position_covs(:, :, i))) / (1000^2); % Convert to km^2
    std_dev_magnitude(i) = sqrt(combined_variance);
end
upper_bound = position_magnitude + std_dev_magnitude;
lower_bound = position_magnitude - std_dev_magnitude;
fill([time_hours, fliplr(time_hours)], [upper_bound, fliplr(lower_bound)], ...
     'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

xlabel('Time (hours)');
ylabel('Position Magnitude (km)');
title('ISS Trajectory: Comparison of Different EKF Component Numbers');
legend('Location', 'best');
grid on;

%% Plot 2: Bhattacharyya Distances for Different N Values
subplot(2, 1, 2);
hold on;

for idx = 1:length(N_values)
    d = bhattacharyya_distances{idx};
    plot(time_hours, d, ...
         'LineWidth', 2, 'DisplayName', sprintf('N=%d', N_values(idx)));
end

xlabel('Time (hours)');
ylabel('Bhattacharyya Distance');
title('Bhattacharyya Distance vs Reference (N=30) Over Time - EKF');
legend('Location', 'best');
grid on;

%% Additional Analysis Plot: Convergence Analysis
figure('Position', [200, 200, 1200, 800]);

% Plot 1: Mean Bhattacharyya Distance vs N
subplot(2, 2, 1);
mean_distances = cellfun(@mean, bhattacharyya_distances);
semilogx(N_values, mean_distances, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Mean Bhattacharyya Distance');
title('EKF Convergence: Mean Distance vs N');
grid on;

% Plot 2: Max Bhattacharyya Distance vs N  
subplot(2, 2, 2);
max_distances = cellfun(@max, bhattacharyya_distances);
semilogx(N_values, max_distances, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Maximum Bhattacharyya Distance');
title('EKF Convergence: Maximum Distance vs N');
grid on;

% Plot 3: Standard deviation of distances vs N
subplot(2, 2, 3);
std_distances = cellfun(@std, bhattacharyya_distances);
semilogx(N_values, std_distances, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Std Dev of Bhattacharyya Distance');
title('EKF Variability: Std Dev vs N');
grid on;

% Plot 4: Distance evolution over time for all N values
subplot(2, 2, 4);
hold on;
for idx = 1:length(N_values)
    d = bhattacharyya_distances{idx};
    plot(time_hours, d, ...
         'LineWidth', 1.5, 'DisplayName', sprintf('N=%d', N_values(idx)));
end
xlabel('Time (hours)');
ylabel('Bhattacharyya Distance');
title('EKF Distance Evolution: All N Values');
legend('Location', 'best');
grid on;

%% Time-specific analysis
% Define time points of interest (1, 2, 3, ... hours)
target_hours = 1:1:floor(max(time_hours));
time_indices = zeros(size(target_hours));

% Find closest time indices for each target hour
for i = 1:length(target_hours)
    [~, time_indices(i)] = min(abs(time_hours - target_hours(i)));
end

% Extract distances at specific time points
distances_at_times = zeros(length(N_values), length(target_hours));
for idx = 1:length(N_values)
    d = bhattacharyya_distances{idx};
    distances_at_times(idx, :) = d(time_indices);
end

% Create new figure for time-specific analysis
figure('Position', [300, 300, 1400, 900]);

% Plot 1: Bar chart showing distances at different hours
subplot(2, 2, 1);
bar_width = 0.8;
x_positions = 1:length(target_hours);
bar_handle = bar(x_positions, distances_at_times', bar_width, 'grouped');

% Color the bars
for idx = 1:length(N_values)
    bar_handle(idx).DisplayName = sprintf('N=%d', N_values(idx));
end

xlabel('Time (hours)');
ylabel('Bhattacharyya Distance');
title('EKF: Bhattacharyya Distance at Hourly Intervals');
legend('Location', 'best');
grid on;
set(gca, 'XTick', x_positions, 'XTickLabel', target_hours);

%% Display comprehensive statistics
fprintf('\n=== Comprehensive EKF Bhattacharyya Distance Statistics ===\n');
fprintf('%-8s %-12s %-12s %-12s %-12s\n', 'N', 'Mean', 'Std Dev', 'Min', 'Max');
fprintf('%-8s %-12s %-12s %-12s %-12s\n', '---', '----', '-------', '---', '---');

for idx = 1:length(N_values)
    d = bhattacharyya_distances{idx};
    fprintf('%-8d %-12.6f %-12.6f %-12.6f %-12.6f\n', ...
        N_values(idx), mean(d), std(d), min(d), max(d));
end

% Compute improvement percentages
fprintf('\n=== EKF Improvement Analysis ===\n');
base_mean = mean(bhattacharyya_distances{1}); % N=1 as baseline
for idx = 2:length(N_values)
    current_mean = mean(bhattacharyya_distances{idx});
    improvement = ((base_mean - current_mean) / base_mean) * 100;
    fprintf('N=%d vs N=1: %.2f%% improvement in mean distance\n', N_values(idx), improvement);
end

function d = bhattacharyya_distance(mu1, Sigma1, mu2, Sigma2)
% Computes the Bhattacharyya distance between two multivariate Gaussians.
% 
% Inputs:
%   mu1, mu2     - Mean vectors (nx1)
%   Sigma1, Sigma2 - Covariance matrices (nxn)
%
% Output:
%   d - Bhattacharyya distance (scalar)

    % Check dimension consistency
    assert(all(size(mu1) == size(mu2)), 'Mean vectors must be the same size.');
    assert(all(size(Sigma1) == size(Sigma2)), 'Covariance matrices must be the same size.');
    
    % Average covariance
    Sigma = 0.5 * (Sigma1 + Sigma2);
    
    % Regularize in case of near-singular matrix
    epsilon = 1e-10;
    Sigma = Sigma + epsilon * eye(size(Sigma));
    
    % First term (Mahalanobis-like term)
    diff = mu1 - mu2;
    term1 = 0.125 * (diff') * (Sigma \ diff);
    
    % Second term (determinant term)
    det_Sigma1 = det(Sigma1 + epsilon * eye(size(Sigma1)));
    det_Sigma2 = det(Sigma2 + epsilon * eye(size(Sigma2)));
    det_Sigma  = det(Sigma);
    
    % Handle numerical issues with determinants
    if det_Sigma <= 0 || det_Sigma1 <= 0 || det_Sigma2 <= 0
        term2 = 0;
    else
        term2 = 0.5 * log(det_Sigma / sqrt(det_Sigma1 * det_Sigma2));
    end
    
    % Bhattacharyya distance
    d = term1 + term2;
end
function [ISS_mixture_cov, ISS_mixture_mean, ISS_EKF_cov, ISS_EKF_mean, timeVector, numTimeSteps] = GenerateEKF(mu, P, N, startTime, stopTime, mu_earth)
    sampleTime = 5; % seconds (matching original GSF code)
    timeVector = startTime:seconds(sampleTime):stopTime;
    numTimeSteps = length(timeVector);
    
    numSamples = 10000;

    % Generate samples from multivariate normal distribution
    data = mvnrnd(mu, P, numSamples);

    % Fit Gaussian Mixture Model
    options = statset('MaxIter', 1000, 'Display', 'off');
    fprintf('Number of samples: %d\n', size(data, 1));
    fprintf('Number of dimensions: %d\n', size(data, 2));
    fprintf('Requested GMM components: %d\n', N);
    gm = fitgmdist(data, N, 'Options', options, 'RegularizationValue', 1e-20);
    
    ISS_gaussian_means = gm.mu';            
    ISS_gaussian_covs = gm.Sigma;       
    Weights = gm.ComponentProportion;
    
    % Initialize output arrays
    ISS_EKF_mean = zeros(6, numTimeSteps, N);
    ISS_EKF_cov = zeros(6, 6, numTimeSteps, N);
    
    % Process each Gaussian component
    for comp = 1:N
        fprintf('Processing EKF component %d of %d\n', comp, N);
        
        
        % Convert equinoctial to Keplerian elements
        Kep = equinoctialToKeplerian(ISS_gaussian_means(:, comp)');
       
        % Get Jacobian for covariance transformation
        J = jacobian_equinoctial_to_keplerian(ISS_gaussian_means(:, comp));
        P_kep = J * ISS_gaussian_covs(:, :, comp) * J';
        

        %[Kep, P_kep] = convertEquinoctialtoKep(ISS_gaussian_means(:, comp), ISS_gaussian_covs(:, :, comp));
        % Create satellite scenario
        sc_mean = satelliteScenario(startTime, stopTime, sampleTime);
        sat = satellite(sc_mean, Kep(1), Kep(2), ...
                       Kep(3), Kep(4), ...
                       Kep(5), Kep(6), "OrbitPropagator", "two-body-keplerian");
        
        % Get states for the new scenario
        [positions, velocities, ~] = states(sat);
        ISS_EKF_mean(:, :, comp) = [positions; velocities];
        
        % Initialize arrays for Jacobians and state transition matrices
        jacob_kep2eci = zeros(6, 6, numTimeSteps);
        Phi = zeros(6, 6, numTimeSteps);
        
        % Calculate Jacobians and state transition matrices for all time steps
        for t = 1:numTimeSteps
            % Get Keplerian elements and Jacobian from ECI
            [~, jacob_eci2kep] = car2kep_ell(positions(:, t), velocities(:, t), mu_earth, true);
            unit_conversion = diag([1, 1, 180/pi, 180/pi, 180/pi, 180/pi]);
            jacob_eci2kep_deg = unit_conversion * jacob_eci2kep;
            jacob_kep2eci(:, :, t) = inv(jacob_eci2kep_deg);
            
            % Calculate state transition matrix for two-body dynamics
            F = twoBodyJacobian(positions(:, t), mu_earth);
            Phi(:, :, t) = expm(F * sampleTime);
        end
        
        % Initialize covariance propagation
        P_eci = zeros(6, 6, numTimeSteps);
        P_eci(:, :, 1) = jacob_kep2eci(:, :, 1) * P_kep * jacob_kep2eci(:, :, 1)';
        
        % Process noise (tune as needed)
        Q = eye(6) * 1e-6;  % Small process noise for two-body dynamics
        
        % Propagate covariance through time
        for t = 1:numTimeSteps-1
            % Use the state transition matrix to propagate covariance
            P_eci(:, :, t+1) = Phi(:, :, t) * P_eci(:, :, t) * Phi(:, :, t)';
        end
        
        % Store results for this component
        ISS_EKF_cov(:, :, :, comp) = P_eci;
    end

    % Calculate mixture statistics
    ISS_mixture_mean = zeros(6, numTimeSteps);
    ISS_mixture_cov = zeros(6, 6, numTimeSteps);
    
    for t = 1:numTimeSteps
        % Calculate weighted mean
        mean_mix_t = zeros(6, 1);
        for j = 1:N
            mean_j = ISS_EKF_mean(:, t, j);
            w_j = Weights(j);
            mean_mix_t = mean_mix_t + w_j * mean_j;
        end
        ISS_mixture_mean(:, t) = mean_mix_t;
        
        % Calculate mixture covariance
        cov_mix_t = zeros(6, 6);
        for j = 1:N
            mean_j = ISS_EKF_mean(:, t, j); 
            cov_j = ISS_EKF_cov(:, :, t, j);
            w_j = Weights(j);
            
            diff = mean_j - mean_mix_t;
            cov_mix_t = cov_mix_t + w_j * (cov_j + diff * diff'); 
        end
        ISS_mixture_cov(:, :, t) = cov_mix_t;
    end
end


function [mu_kep, P_kep] = convertEquinoctialtoKep(mu_est, P_est)
    numSamples = 10000;
    stateDim = 6;
    samples = mvnrnd(mu_est, P_est, numSamples)';
    kep_samples = zeros(6, numSamples);
    
    for k = 1:numSamples
        kep_samples(:, k) = equinoctialToKeplerian(samples(:, k));
    end
    
    % Initialize mean vector
    mu_kep = zeros(6, 1);
    
    % For non-angular components (a, e, i) - use regular mean
    % Components 1-3: semi-major axis, eccentricity, inclination
    mu_kep(1:3) = mean(kep_samples(1:3, :), 2);
    
    % For angular components (Ω, ω, M) - use circular mean
    % Components 4-6: RAAN, argument of periapsis, mean anomaly
    for i = 4:6
        % Convert degrees to radians for circular statistics
        angles_rad = deg2rad(kep_samples(i, :));
        
        % Compute circular mean
        mean_cos = mean(cos(angles_rad));
        mean_sin = mean(sin(angles_rad));
        circular_mean_rad = atan2(mean_sin, mean_cos);
        
        % Convert back to degrees and ensure [0, 360) range
        mu_kep(i) = mod(rad2deg(circular_mean_rad), 360);
    end
    
    % Compute covariance matrix with special handling for angular components
    P_kep = zeros(6, 6);
    
    % Create a matrix to store transformed samples for covariance calculation
    transformed_samples = kep_samples;
    
    % For angular components, unwrap them relative to the circular mean
    for i = 4:6
        ref_angle = mu_kep(i);
        for k = 1:numSamples
            % Unwrap angle relative to reference
            angle_diff = kep_samples(i, k) - ref_angle;
            angle_diff = mod(angle_diff + 180, 360) - 180; % Wrap to [-180, 180]
            transformed_samples(i, k) = ref_angle + angle_diff;
        end
    end
    
    % Now compute covariance using the unwrapped samples
    P_kep = cov(transformed_samples');

end

% Helper function for two-body Jacobian (you may need to implement this)
function F = twoBodyJacobian(r, mu)
    % Two-body dynamics Jacobian matrix
    % r: position vector (3x1)
    % mu: gravitational parameter
    
    r_norm = norm(r);
    r_norm3 = r_norm^3;
    r_norm5 = r_norm^5;
    
    % Initialize F matrix
    F = zeros(6, 6);
    
    % Upper right block (velocity coupling)
    F(1:3, 4:6) = eye(3);
    
    % Lower left block (position coupling for acceleration)
    F(4:6, 1:3) = -mu/r_norm3 * eye(3) + 3*mu/r_norm5 * (r * r');
end


function [kep, jacob] = car2kep_ell(pos, vel, mu, cjac)
% CAR2KEP_ELL Cartesian to classical Keplerian orbital elements
%
% Syntax:
%   [kep, jacob] = car2kep_ell(pos, vel)
%   [kep, jacob] = car2kep_ell(pos, vel, mu)
%   [kep, jacob] = car2kep_ell(pos, vel, mu, cjac)
%
% Description:
%   Converts cartesian orbital elements to classical Keplerian orbital elements.
%   The transformation jacobian is optionally computed.
%   This function only handles elliptical orbits.
%
% Inputs:
%   pos  - Position [X;Y;Z] [m] (3xN matrix)
%   vel  - Velocity [Vx;Vy;Vz] [m/s] (3xN matrix)
%   mu   - (optional) Gravitational constant [m^3/s^2] (default: 3.986004418e14 for Earth)
%   cjac - (optional) Set to false to avoid computing jacobian (default: true)
%
% Outputs:
%   kep   - Classical Keplerian orbital elements [sma;e;inc;raan;pom;nu] [m,rad] (6xN)
%           where nu is the true anomaly (changed from mean anomaly)
%           Note: RAAN and argument of perigee are swapped from typical order
%   jacob - (optional) Transformation jacobian (6x6xN)
%
% Example:
%   pos = [7000e3; 1000e3; -500e3];
%   vel = [1e3; 2e3; 7e3];
%   [kep, jacob] = car2kep_ell(pos, vel);

% Default values
if nargin < 3 || isempty(mu)
    mu = 3.986004418e14; % Earth's gravitational parameter [m^3/s^2]
end

if nargin < 4
    cjac = true;
end

% Orbital thresholds
EPS_CIR = 1e-10;   % Circular orbit threshold
EPS_EQUA = 1e-10;  % Equatorial orbit threshold

N = size(pos, 2);

% Initialize jacobian
if nargout > 1
    jacob = zeros(6, 6, N);
else
    cjac = false;
end

% Helper functions
function result = norm_vec(v)
    result = sqrt(sum(v.^2, 1));
end

function result = cross_vec(a, b)
    result = [a(2,:).*b(3,:) - a(3,:).*b(2,:);
              a(3,:).*b(1,:) - a(1,:).*b(3,:);
              a(1,:).*b(2,:) - a(2,:).*b(1,:)];
end

function result = dot_vec(a, b)
    result = sum(a.*b, 1);
end

function result = reduce_angle(angles)
    result = mod(angles, 2*pi);
end

% Main calculations
r = norm_vec(pos);
V = norm_vec(vel);
W = cross_vec(pos, vel);
pos_vel = dot_vec(pos, vel);

% Semi-major axis
a = r ./ (2 - r .* V.^2 / mu);

% Inclination
inc = atan2(sqrt(W(1,:).^2 + W(2,:).^2), W(3,:));

% Right ascension of ascending node (RAAN)
node = cross_vec([0; 0; 1], W);
gom = atan2(node(2,:), node(1,:));

% Handle equatorial orbits
Ieq = find(sin(inc) < EPS_EQUA);
gom(Ieq) = 0;
node(1, Ieq) = 1;
node(2, Ieq) = 0;
node(3, Ieq) = 0;

% Eccentricity
esinE = pos_vel ./ sqrt(mu .* a);
ecosE = r .* V.^2 / mu - 1;
e = sqrt(ecosE.^2 + esinE.^2);
Icir = find(e < EPS_CIR);

% Eccentric anomaly and true anomaly
E = atan2(esinE, ecosE);

% Argument of perigee
cos_alpha_v = dot_vec(pos, node);
sin_alpha_v = dot_vec(pos, cross_vec(W, node)) ./ norm_vec(W);
alpha_v = atan2(sin_alpha_v, cos_alpha_v);

% True anomaly from eccentric anomaly
nu = 2 * atan2(sqrt(1+e) .* sin(E/2), sqrt(1-e) .* cos(E/2));
pom = alpha_v - nu;

% Handle circular orbits - for circular orbits, use argument of latitude
% In circular orbits, true anomaly becomes the argument of latitude
pom_circular = pom;
pom(Icir) = 0;
nu(Icir) = nu(Icir) + pom_circular(Icir);

% Reduce angles to [0, 2π)
kep = [a; e; inc; reduce_angle(gom); reduce_angle(pom); reduce_angle(nu)];

% Jacobian computation
if cjac && nargout > 1
    % Set eccentricity to NaN for circular orbits to handle divisions
    ecir = e;
    ecir(Icir) = NaN;
    
    sinE = sin(E);
    cosE = cos(E);
    n = sqrt(mu ./ a.^3);
    
    % P and Q vectors (rotation matrix columns)
    P1 = cosE ./ r .* pos(1,:) - sinE ./ (n .* a) .* vel(1,:);
    P2 = cosE ./ r .* pos(2,:) - sinE ./ (n .* a) .* vel(2,:);
    P3 = cosE ./ r .* pos(3,:) - sinE ./ (n .* a) .* vel(3,:);
    
    Q1 = 1 ./ sqrt(1-e.^2) .* (sinE ./ r .* pos(1,:) + (cosE - e) ./ (n .* a) .* vel(1,:));
    Q2 = 1 ./ sqrt(1-e.^2) .* (sinE ./ r .* pos(2,:) + (cosE - e) ./ (n .* a) .* vel(2,:));
    Q3 = 1 ./ sqrt(1-e.^2) .* (sinE ./ r .* pos(3,:) + (cosE - e) ./ (n .* a) .* vel(3,:));
    
    % Trigonometric expressions
    cosi = P1.*Q2 - P2.*Q1;
    sini = sqrt(P3.^2 + Q3.^2);
    sini(Ieq) = NaN; % Handle equatorial orbits
    
    cosp = Q3 ./ sini;
    sinp = P3 ./ sini;
    cosg = P1 .* cosp - Q1 .* sinp;
    sing = P2 .* cosp - Q2 .* sinp;
    
    % True anomaly derivatives in terms of eccentric anomaly
    % nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2))
    % dnu/dE = sqrt(1-e^2)/(1-e*cos(E))
    % dnu/de = sin(E)/(1-e*cos(E))
    
    dnu_dE = sqrt(1-e.^2) ./ (1 - e.*cosE);
    dnu_de = sinE ./ (1 - e.*cosE);
    
    % Compute partial derivatives
    for j = 1:3
        % Semi-major axis derivatives
        dadr = 2 * a.^2 ./ r.^3 .* pos(j,:);
        dadv = 2 ./ (n.^2 .* a) .* vel(j,:);
        
        % Eccentricity components derivatives
        decosEdr = V.^2 ./ (r * mu) .* pos(j,:);
        desinEdr = vel(j,:) ./ (n .* a.^2) - pos_vel ./ (n .* a .* r.^3) .* pos(j,:);
        decosEdv = 2 * r / mu .* vel(j,:);
        desinEdv = pos(j,:) ./ (n .* a.^2) - pos_vel ./ (n .* a * mu) .* vel(j,:);
        
        % Eccentricity derivatives
        dedr = sinE .* desinEdr + cosE .* decosEdr;
        dedv = sinE .* desinEdv + cosE .* decosEdv;
        
        % Eccentric anomaly derivatives
        dEdr = 1 ./ ecir .* (cosE .* desinEdr - sinE .* decosEdr);
        dEdv = 1 ./ ecir .* (cosE .* desinEdv - sinE .* decosEdv);
        
        % True anomaly derivatives (replacing mean anomaly derivatives)
        dnudr = dnu_dE .* dEdr + dnu_de .* dedr;
        dnudv = dnu_dE .* dEdv + dnu_de .* dedv;
        
        % P vector derivatives
        dPkdr = zeros(3, N);
        for k = 1:3
            dPkdr(k,:) = (j==k) * cosE ./ r - pos(k,:) .* cosE ./ r.^3 .* pos(j,:) - ...
                         vel(k,:) .* sinE ./ (2*n.*a.^2) .* dadr - ...
                         (sinE ./ r .* pos(k,:) + cosE ./ (n .* a) .* vel(k,:)) .* dEdr;
        end
        
        dPkdv = zeros(3, N);
        for k = 1:3
            dPkdv(k,:) = (j==k) * (-sinE ./ (n.*a)) - vel(k,:) .* sinE ./ (2*n.*a.^2) .* dadv - ...
                         (sinE ./ r .* pos(k,:) + cosE ./ (n .* a) .* vel(k,:)) .* dEdv;
        end
        
        % Q vector derivatives
        dQkdr = zeros(3, N);
        for k = 1:3
            dQkdr(k,:) = 1 ./ sqrt(1-e.^2) .* ...
                         ((j==k) * sinE ./ r - pos(k,:) .* sinE ./ r.^3 .* pos(j,:) + ...
                          vel(k,:) .* (cosE-e) ./ (2*n.*a.^2) .* dadr + ...
                          (cosE ./ r .* pos(k,:) - sinE ./ (n .* a) .* vel(k,:)) .* dEdr + ...
                          (esinE ./ (r .* (1-e.^2)) .* pos(k,:) + (ecosE-1) ./ (n.*a.*(1-e.^2)) .* vel(k,:)) .* dedr);
        end
        
        dQkdv = zeros(3, N);
        for k = 1:3
            dQkdv(k,:) = 1 ./ sqrt(1-e.^2) .* ...
                         ((j==k) * (cosE-e) ./ (n.*a) + vel(k,:) .* (cosE-e) ./ (2*n.*a.^2) .* dadv + ...
                          (cosE ./ r .* pos(k,:) - sinE ./ (n .* a) .* vel(k,:)) .* dEdv + ...
                          (esinE ./ (r .* (1-e.^2)) .* pos(k,:) + (ecosE-1) ./ (n.*a.*(1-e.^2)) .* vel(k,:)) .* dedv);
        end
        
        % Inclination derivatives
        didr = cosi ./ sini .* (dPkdr(3,:) .* P3 + dQkdr(3,:) .* Q3) - ...
               (P1 .* dQkdr(2,:) + Q2 .* dPkdr(1,:) - P2 .* dQkdr(1,:) - Q1 .* dPkdr(2,:)) .* sini;
        didv = cosi ./ sini .* (dPkdv(3,:) .* P3 + dQkdv(3,:) .* Q3) - ...
               (P1 .* dQkdv(2,:) + Q2 .* dPkdv(1,:) - P2 .* dQkdv(1,:) - Q1 .* dPkdv(2,:)) .* sini;
        
        % Argument of perigee derivatives
        dpomdr = 1 ./ sini .* (dPkdr(3,:) .* cosp - dQkdr(3,:) .* sinp);
        dpomdv = 1 ./ sini .* (dPkdv(3,:) .* cosp - dQkdv(3,:) .* sinp);
        
        % RAAN derivatives
        dgomdr = (dPkdr(2,:) .* cosp - dQkdr(2,:) .* sinp) .* cosg - ...
                 (dPkdr(1,:) .* cosp - dQkdr(1,:) .* sinp) .* sing - dpomdr .* cosi;
        dgomdv = (dPkdv(2,:) .* cosp - dQkdv(2,:) .* sinp) .* cosg - ...
                 (dPkdv(1,:) .* cosp - dQkdv(1,:) .* sinp) .* sing - dpomdv .* cosi;
        
        % Assign to jacobian matrix
        jacob(1, j, :) = dadr;
        jacob(2, j, :) = dedr;
        jacob(3, j, :) = didr;
        jacob(4, j, :) = dgomdr;  % Changed from dpomdr to dgomdr (RAAN now 4th)
        jacob(5, j, :) = dpomdr;  % Changed from dgomdr to dpomdr (pom now 5th)
        jacob(6, j, :) = dnudr;  % Changed from dMdr to dnudr
        
        jacob(1, j+3, :) = dadv;
        jacob(2, j+3, :) = dedv;
        jacob(3, j+3, :) = didv;
        jacob(4, j+3, :) = dgomdv;  % Changed from dpomdv to dgomdv (RAAN now 4th)
        jacob(5, j+3, :) = dpomdv;  % Changed from dgomdv to dpomdv (pom now 5th)
        jacob(6, j+3, :) = dnudv;  % Changed from dMdv to dnudv
    end
end

end

%{
function F = twoBodyJacobian(r, mu)
    r_norm = norm(r);
    F = [zeros(3,3), eye(3);
         mu/r_norm^5 * (3*r*r' - r_norm^2*eye(3)), zeros(3,3)];
end
%}

function J = jacobian_equinoctial_to_keplerian(E)
% JACOBIAN_EQUINOCTIAL_TO_KEPLERIAN Compute the Jacobian matrix converting
% Equinoctial elements to Keplerian elements.
%
% Inputs:
% E - 6x1 vector of equinoctial elements: [a; h; k; p; q; lambda]
%
% Outputs:
% J - 6x6 Jacobian matrix d(Keplerian)/d(Equinoctial)
% Keplerian output order: [a; e; i; Omega; omega; nu] 
% where i, Omega, omega, nu are in DEGREES
%
% Unpack equinoctial elements
a = E(1);
h = E(2);
k = E(3);
p = E(4);
q = E(5);
lambda = E(6);

% Compute intermediates
e = sqrt(h^2 + k^2);
omega_plus_Omega = atan2(h, k);
tan_i_2 = sqrt(p^2 + q^2);
i = 2 * atan(tan_i_2);
Omega = atan2(p, q);
omega = omega_plus_Omega - Omega;

% Conversion factor from radians to degrees
rad2deg = 180/pi;

% Handle singularities
if e < 1e-12
    % For circular orbits, set h and k derivatives carefully
    de_dh = 0;
    de_dk = 0;
else
    de_dh = h/e;
    de_dk = k/e;
end

if tan_i_2 < 1e-12
    % For equatorial orbits
    di_dp = 0;
    di_dq = 0;
    dOmega_dp = 0;
    dOmega_dq = 0;
else
    % ∂i/∂p and ∂i/∂q - CORRECTED
    % d/dp[2*atan(sqrt(p^2+q^2))] = 2 * 1/(1+(p^2+q^2)) * p/sqrt(p^2+q^2)
    di_dp = 2 * p / ((1 + p^2 + q^2) * tan_i_2);
    di_dq = 2 * q / ((1 + p^2 + q^2) * tan_i_2);
    
    % ∂Ω/∂p and ∂Ω/∂q - CORRECTED signs
    % d/dp[atan2(p,q)] = q/(p^2+q^2), d/dq[atan2(p,q)] = -p/(p^2+q^2)
    dOmega_dp = q / (p^2 + q^2);
    dOmega_dq = -p / (p^2 + q^2);
end

if h^2 + k^2 < 1e-12
    % For circular orbits
    datan2hk_dh = 0;
    datan2hk_dk = 0;
else
    % ∂/∂h[atan2(h,k)] = k/(h^2+k^2), ∂/∂k[atan2(h,k)] = -h/(h^2+k^2)
    datan2hk_dh = k / (h^2 + k^2);
    datan2hk_dk = -h / (h^2 + k^2);
end

% Partial derivatives
% ∂a/∂E = [1 0 0 0 0 0] (no change - semi-major axis stays in same units)
da_dE = [1, 0, 0, 0, 0, 0];

% ∂e/∂E (no change - eccentricity is dimensionless)
de_dE = [0, de_dh, de_dk, 0, 0, 0];

% ∂i/∂E - MULTIPLY BY rad2deg for degrees output
di_dE = rad2deg * [0, 0, 0, di_dp, di_dq, 0];

% ∂Ω/∂E - MULTIPLY BY rad2deg for degrees output  
dOmega_dE = rad2deg * [0, 0, 0, dOmega_dp, dOmega_dq, 0];

% ∂ω/∂E = ∂(ω+Ω)/∂E - ∂Ω/∂E - MULTIPLY BY rad2deg for degrees output
domega_dE = rad2deg * [0, datan2hk_dh, datan2hk_dk, -dOmega_dp, -dOmega_dq, 0];

% For True Anomaly (nu): - MULTIPLY BY rad2deg for degrees output
% nu = lambda - (omega + Omega) = lambda - atan2(h,k)
dnu_dE = rad2deg * [0, -datan2hk_dh, -datan2hk_dk, 0, 0, 1/rad2deg];
% Note: ∂nu/∂lambda = 1, but since lambda is in radians and nu output is in degrees,
% we need 1 * rad2deg / rad2deg = 1 in the final term

% Assemble Jacobian
J = [da_dE;
     de_dE;
     di_dE;
     dOmega_dE;
     domega_dE;
     dnu_dE];
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
    a = ommData.semiMajorAxis * 1000; % Convert km to m
    ecc = ommData.eccentricity;
    incl = deg2rad(ommData.inclination);
    raan = deg2rad(ommData.raOfAscNode);
    argp = deg2rad(ommData.argOfPericenter);
    
    % Mean anomaly to true anomaly
    M = deg2rad(ommData.meanAnomaly);
    E = M; % Initial guess
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
    nominal.l = mod(nu + raan + argp, 2*pi);  % Ensure l is in [0, 2π)
    
    perc.a = 0.0001;   % 5% of a
    perc.h = 0.005;   % 1% of h
    perc.k = 0.005;   % 1% of k
    perc.p = 0.005;   % 1% of p
    perc.q = 0.005;   % 1% of q
    perc.l = 0.00001;  % 0.5% of l
    
    % Calculate uncertainties as percentages of nominal values
    uncertainty.a = perc.a * nominal.a;
    uncertainty.h = perc.h * nominal.h;
    uncertainty.k = perc.k * nominal.k;
    uncertainty.p = perc.p * nominal.p;
    uncertainty.q = perc.q * nominal.q;
    uncertainty.l = perc.l * nominal.l;


    % Diagonal covariance matrix
    P = diag([uncertainty.a^2, ...
              uncertainty.h^2, ...
              uncertainty.k^2, ...
              uncertainty.p^2, ...
              uncertainty.q^2, ...
              uncertainty.l^2]);
end


function keplerian = equinoctialToKeplerian(equinoctial)
% Converts equinoctial orbital elements to classical Keplerian elements
%
% Input:
%   equinoctial = [a, h, k, p, q, lambda]
%     a      : semi-major axis
%     h, k   : e*sin(omega+Omega), e*cos(omega+Omega)
%     p, q   : tan(i/2)*sin(Omega), tan(i/2)*cos(Omega)
%     lambda : mean longitude = Omega + omega + M
%
% Output:
%   keplerian = [a, e, i, Omega, omega, M]
%     All angles in degrees

% Extract equinoctial elements
a      = equinoctial(1);
h      = equinoctial(2);
k      = equinoctial(3);
p      = equinoctial(4);
q      = equinoctial(5);
lambda = equinoctial(6);

% Eccentricity
e = sqrt(h^2 + k^2);

% Inclination
i = 2 * atan( sqrt(p^2 + q^2) );

% RAAN (Omega)
Omega = atan2(p, q);
Omega = mod(Omega, 2*pi);

% Argument of periapsis (omega)
omega = atan2(h, k) - Omega;
omega = mod(omega, 2*pi);

% Compute true anomaly (nu)
nu = lambda - Omega - omega;
nu = mod(nu, 2*pi);

% Compute eccentric anomaly (E) from true anomaly (nu)
E = 2 * atan( tan(nu / 2) / sqrt((1 + e) / (1 - e)) );
E = mod(E, 2*pi);

% Compute mean anomaly (M)
M = E - e * sin(E);
M = mod(M, 2*pi);

% Convert angles to degrees
i     = rad2deg(i);
Omega = rad2deg(Omega);
omega = rad2deg(omega);
M     = rad2deg(M);

% Output Keplerian elements
keplerian = [a, e, i, Omega, omega, M];

end





