%% Time Setup
startTime = datetime(2020,5,11,17,35,38);
stopTime = startTime + days(1);
sampleTime = 30; % seconds
timeVector = startTime:seconds(sampleTime):stopTime;
numTimeSteps = length(timeVector);

%% Earth Parameters
earthRadius = 6.378137e6; % meters
mu_earth = 3.986004418e14;

xmlFolder = '/home/roboticsadmin/MATLABProjects/Debris_Folder';
ommFiles = dir(fullfile(xmlFolder, '*.xml'));

fprintf('Reading space debris data from folder: %s\n', xmlFolder);
fprintf('Found %d XML files.\n', length(ommFiles));

N = 1;
numSamples = 10000;
lat = 10;
lon = -50;

sc_mean = satelliteScenario(startTime, stopTime, sampleTime);

%% Load ISS Data
ISS_omm_file = 'SL24.xml';
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

[weighted_covariance, weighted_mean, gaussian_ECI_means, gaussian_ECI_covs, ISS_gaussian_means, ISS_gaussian_covs, Weights] = GenerateGSF(mu, P, N, numSamples, startTime, stopTime);

measurement_data = mvnrnd(mu, P, numSamples);
options = statset('MaxIter', 1000, 'Display', 'off');
gm = fitgmdist(measurement_data, 1, 'Options', options, 'RegularizationValue', 1e-20);
measurement_means = gm.mu';   

Initial_kep = equinoctialToKeplerian(measurement_means);
Initial_sat = satellite(sc_mean, Initial_kep(1), Initial_kep(2), Initial_kep(3), Initial_kep(4), Initial_kep(5), Initial_kep(6), "OrbitPropagator", "two-body-keplerian");

Initial_gs = groundStation(sc_mean, lat, lon);
ac = access(Initial_sat, Initial_gs);
intvls = accessIntervals(ac);

Nmin = 0.2 * N;
fprintf('\n=== Processing All Collision Times ===\n');
fprintf('Found %d access intervals\n', height(intvls));

timeVector = datetime(timeVector, 'TimeZone', '');
intvls.StartTime = datetime(intvls.StartTime, 'TimeZone', '');
intvls.EndTime = datetime(intvls.EndTime, 'TimeZone', '');

for interval_idx = 1:height(intvls)
    interval_start = intvls.StartTime(interval_idx);
    interval_end = intvls.EndTime(interval_idx);
    
    interval_mask = (timeVector >= interval_start) & (timeVector <= interval_end);
    interval_time_indices = find(interval_mask);
end

ISS_positions = weighted_mean(1:3, :);
ISS_position_covs = weighted_covariance(1:3, 1:3, :);

ISS_positions_km = ISS_positions / 1000;
earthRadius_km = earthRadius / 1000;

time_hours = hours(timeVector - startTime);

figure('Position', [100, 100, 1400, 1000]);
subplot(2, 3, 1);
position_magnitude = sqrt(sum(ISS_positions_km.^2, 1));
plot(time_hours, position_magnitude, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Position Magnitude');
hold on;

std_dev_magnitude = zeros(1, numTimeSteps);
for i = 1:numTimeSteps
    combined_variance = sum(diag(ISS_position_covs(:, :, i))) / (1000^2);
    std_dev_magnitude(i) = sqrt(combined_variance);
end

upper_bound = position_magnitude + std_dev_magnitude;
lower_bound = position_magnitude - std_dev_magnitude;
fill([time_hours, fliplr(time_hours)], [upper_bound, fliplr(lower_bound)], ...
     'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

xlabel('Time (hours)');
ylabel('Position Magnitude (km)');
legend('Location', 'best');
grid on;

predicted_measurements = [];
measurement_residuals = [];
update_indices = [];
update_counter = 0;
all_determinants = [];
mahalanobis_distances = [];
Neff_vals = [];
sc_measurement = satelliteScenario(startTime, stopTime, sampleTime);
gs_measurement = groundStation(sc_measurement, lat, lon);

for interval_idx = 1:height(intvls)
    fprintf('\nProcessing interval %d of %d...\n', interval_idx, height(intvls));
    
    interval_start = intvls.StartTime(interval_idx);
    interval_end = intvls.EndTime(interval_idx);
    
    interval_mask = (timeVector >= interval_start) & (timeVector <= interval_end);
    interval_time_indices = find(interval_mask);
    interval_time_indices = interval_time_indices(1:6:end);

    if isempty(interval_time_indices)
        fprintf('  Warning: No time points found in interval, skipping...\n');
        continue;
    end
    
    fprintf('  Interval from %s to %s\n', datetime(interval_start), datetime(interval_end));
    fprintf('  Processing %d observation times in this interval\n', length(interval_time_indices));
    
    for obs_idx = 1:length(interval_time_indices)
        timeIndex = interval_time_indices(obs_idx);
        measurement_time = timeVector(timeIndex);
        
        update_counter = update_counter + 1;
        update_indices(update_counter) = timeIndex;
    
        fprintf('    Observation %d/%d at time: %s (index %d)\n', ...
                obs_idx, length(interval_time_indices), datetime(measurement_time), timeIndex);

        [measured_position, measured_velocity] = states(Initial_sat, measurement_time);

        z_true = [measured_position; measured_velocity];
        
        R = diag([10^2,10^2,10^2, 2e-2^2,2e-2^2,2e-2^2]);
     
        noise = sqrt(diag(R)) .* randn(size(R,1),1);
        
        z = z_true + noise;
        
        mu_est = zeros(6, N);
        P_est = zeros(6, 6, N);
        y_all = cell(N,1);
        z_pred_all = zeros(6, N);
        S_all = cell(N,1);
        
        for i = 1:N
            ECI_cov = gaussian_ECI_covs(:,:, timeIndex, i);
            ECI_mean = gaussian_ECI_means(:, timeIndex, i);
            z_pred = ECI_mean;
            z_pred_all(:, i) = z_pred;
            H_augmented = eye(6);
            [mu_est(:, i), P_est(:, :, i), S, y] = MakeEstimates(ECI_mean, ECI_cov, z, R, z_pred, H_augmented);
            y_all{i} = y;
            S_all{i} = S;
        end

        z_pred_weighted = zeros(6, 1);
        for i = 1:N
            z_pred_weighted = z_pred_weighted + Weights(i) * z_pred_all(:, i);
        end
        fprintf('Predicted: Pos_x = %.6f, Vel_x = %.6f', z_pred_weighted(1), z_pred_weighted(4));
        fprintf('Measured : Pos_x = %.6f, Vel_x = %.6f', z(1), z(4));

        y_weighted = zeros(6, 1);
        y_weighted(1) = z(1) - z_pred_weighted(1);
        y_weighted(2) = z(2) - z_pred_weighted(2);
        y_weighted(3) = z(3) - z_pred_weighted(3);
        y_weighted(4) = z(4) - z_pred_weighted(4);
        y_weighted(5) = z(5) - z_pred_weighted(5);
        y_weighted(6) = z(6) - z_pred_weighted(6);

        predicted_measurements = [predicted_measurements, z_pred_weighted];
        measurement_residuals = [measurement_residuals, y_weighted];
        all_determinants = [all_determinants, det(weighted_covariance(:, :, timeIndex))];

        S_mixture = zeros(6, 6);
        for i = 1:N
            S_mixture = S_mixture + Weights(i) * S_all{i};
            y_diff = z_pred_all(:, i) - z_pred_weighted;
            S_mixture = S_mixture + Weights(i) * (y_diff * y_diff');
        end
        mahalanobis_dist_overall = sqrt(y_weighted' * (S_mixture \ y_weighted));
        mahalanobis_distances(update_counter) = mahalanobis_dist_overall;
        disp('Mahalanobis distance of measurement residual')
        disp(mahalanobis_dist_overall);
        
        l = compute_weight(z_pred_all, z, inv(S_mixture), N);
        Weights = Weights .* l;
        Weights = Weights / sum(Weights);
      
        [keep, Neff] = stratified_resample_unique(Weights);
        Neff_vals = [Neff_vals, Neff];

        disp(Neff)

        if Neff < Nmin
            fprintf('    Neff (%.2f) < Nmin (%.2f), performing resampling...\n', Neff, Nmin);
            [keep_indices, ~] = stratified_resample_unique(Weights);
            kept_means = mu_est(:, keep_indices);
            kept_covs = P_est(:, :, keep_indices);
            kept_weights = ones(length(keep_indices), 1) / length(keep_indices);
            fprintf('    Selected %d components from original %d\n', length(keep_indices), N);
            x_mean = kept_means * kept_weights;
            P_mixture = zeros(size(kept_covs, 1), size(kept_covs, 2));
            for i = 1:length(keep_indices)
                diff = kept_means(:, i) - x_mean;
                P_mixture = P_mixture + kept_weights(i) * (kept_covs(:, :, i) + diff * diff');
            end
            P_mixture = P_mixture + 1e-8 * eye(size(P_mixture));
            try
                sampled_Particles = mvnrnd(x_mean, P_mixture, numSamples);
                options = statset('MaxIter', 1000, 'Display', 'off');
                gm = fitgmdist(sampled_Particles, N, 'Options', options, 'RegularizationValue', 1e-6);
                sample_gaussian_means = gm.mu';
                sample_gaussian_covs = gm.Sigma;
                Weights = gm.ComponentProportion';
                Weights = Weights / sum(Weights);
                fprintf('    Resampling successful. New weights sum: %.6f\n', sum(Weights));
            catch ME
                fprintf('    Warning: Resampling failed (%s). Using original estimates.\n', ME.message);
                sample_gaussian_means = mu_est;
                sample_gaussian_covs = P_est;
                Weights = Weights / sum(Weights);
            end
        else
            fprintf('    Neff (%.2f) >= Nmin (%.2f), no resampling needed\n', Neff, Nmin);
            sample_gaussian_means = mu_est;
            sample_gaussian_covs = P_est;
            Weights = Weights / sum(Weights);
        end

        alpha = 1e-3;
        beta = 2;
        kappa = -3;
        n = 6; 
        n_sigma = 2*n + 1; 
        Q_reduced = diag([1e-2, 1e-2, 1e-2, 1e-4, 1e-4, 1e-4]);
        ISS_all_sigma_points = zeros(n, n_sigma, N);
        ISS_all_Wm = zeros(n_sigma, N);
        ISS_all_Wc = zeros(n_sigma, N);
        
        for i = 1:N
            P_prop = sample_gaussian_covs(:, :, i)+ Q_reduced;
            [sigma_points, Wm, Wc] = generate_sigma_points(sample_gaussian_means(:, i), P_prop, alpha, beta, kappa);
            ISS_all_sigma_points(:, :, i) = sigma_points;
            ISS_all_Wm(:, i) = Wm;
            ISS_all_Wc(:, i) = Wc;
        end

        keplerian_sigma_points = zeros(n, n_sigma, N);
        for i = 1:N
            for j = 1:n_sigma
                [kep_j, jacob_eci2kep] = car2kep_ell(ISS_all_sigma_points(1:3, j, i), ISS_all_sigma_points(4:6, j, i), mu_earth, true);
                kep_j(3:6) = rad2deg(kep_j(3:6));
                keplerian_sigma_points(:, j, i) = kep_j';
            end
        end
        
        ISS_positions = cell(n_sigma, N);
        ISS_velocities = cell(n_sigma, N);
        remaining_times = timeVector(timeIndex:end);
        num_remaining_steps = length(remaining_times);
        ISS_predicted_mean_cartesian = zeros(6, num_remaining_steps, N);
        ISS_predicted_cov_cartesian = zeros(6, 6, num_remaining_steps, N);
        
        for j = 1:N
            for i = 1:n_sigma
                a = keplerian_sigma_points(1, i, j);
                e = keplerian_sigma_points(2, i, j);
                inc = keplerian_sigma_points(3, i, j);
                raan = keplerian_sigma_points(4, i, j);
                argP = keplerian_sigma_points(5, i, j);
                theta = keplerian_sigma_points(6, i, j);
                
                sc_ISS = satelliteScenario(measurement_time, stopTime, sampleTime);
                sat = satellite(sc_ISS, a, e, inc, raan, argP, theta, "OrbitPropagator", "two-body-keplerian");
                
                [pos, vel, ~] = states(sat);
                ISS_positions{i, j} = pos;
                ISS_velocities{i, j} = vel;
            end
           
            for t = 1:num_remaining_steps
                sigma_points_cartesian = zeros(6, n_sigma);
                for i = 1:n_sigma
                    r_t = ISS_positions{i, j}(:, t);
                    v_t = ISS_velocities{i, j}(:, t);
                    sigma_points_cartesian(:, i) = [r_t; v_t];
                end
                
                Wm_j = ISS_all_Wm(:, j);
                Wc_j = ISS_all_Wc(:, j);
                
                mean_jt = sigma_points_cartesian * Wm_j; 
                cov_jt = zeros(6, 6);
                for i = 1:n_sigma
                    diff = sigma_points_cartesian(:, i) - mean_jt;
                    cov_jt = cov_jt + Wc_j(i) * (diff * diff');
                end
                
                ISS_predicted_mean_cartesian(:, t, j) = mean_jt;
                ISS_predicted_cov_cartesian(:, :, t, j) = cov_jt;
            end
            
            gaussian_ECI_means(:, timeIndex:end, j) = ISS_predicted_mean_cartesian(:, :, j);
            gaussian_ECI_covs(:, :, timeIndex:end, j) = ISS_predicted_cov_cartesian(:, :, :, j);
            gaussian_ECI_means(:, timeIndex, j) = sample_gaussian_means(:, j);
            gaussian_ECI_covs(:, :, timeIndex, j) = sample_gaussian_covs(:, :, j);
        end

        for t = timeIndex:size(gaussian_ECI_means, 2)
            mean_mix_t = zeros(6, 1);
            for j = 1:N
                mean_j = gaussian_ECI_means(:, t, j);
                w_j = Weights(j);
                mean_mix_t = mean_mix_t + w_j * mean_j;
            end
            weighted_mean(:, t) = mean_mix_t;
            
            cov_mix_t = zeros(6, 6);
            for j = 1:N
                mean_j = gaussian_ECI_means(:, t, j); 
                cov_j = gaussian_ECI_covs(:, :, t, j);
                w_j = Weights(j);
                diff = mean_j - mean_mix_t;
                cov_mix_t = cov_mix_t + w_j * (cov_j + diff * diff'); 
            end
            weighted_covariance(:, :, t) = cov_mix_t;
        end
    end
end

function l = compute_weight(z, zObs, R, N)
    RSqrtm = sqrtm(R);
    nu = RSqrtm * (z - zObs);
    log_l = -0.5 * sum(nu.^2, 1);
    max_log_l = max(log_l);  
    l = exp(log_l - max_log_l);
    l = l / sum(l);
end


%% Plot First Particle Trajectory Over Time
fprintf('\n=== Plotting First Particle Trajectory ===\n');
ISS_positions = weighted_mean(1:3, :);
ISS_position_covs = weighted_covariance(1:3, 1:3, :);
ISS_positions_km = ISS_positions / 1000;
earthRadius_km = earthRadius / 1000;
time_hours = hours(timeVector - startTime);


% Position magnitude plot
figure('Position', [100, 100, 1400, 1000]);
subplot(2, 2, 1);
position_magnitude = sqrt(sum(ISS_positions_km.^2, 1));
disp(size(position_magnitude))
disp(size(time_hours))
plot(time_hours, position_magnitude, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Position Magnitude');
hold on;
% Calculate uncertainty bounds
std_dev_magnitude = zeros(1, numTimeSteps);
for i = 1:numTimeSteps
    combined_variance = sum(diag(ISS_position_covs(:, :, i))) / (1000^2);
    std_dev_magnitude(i) = sqrt(combined_variance);
end

% Add uncertainty bounds
upper_bound = position_magnitude + std_dev_magnitude;
lower_bound = position_magnitude - std_dev_magnitude;
fill([time_hours, fliplr(time_hours)], [upper_bound, fliplr(lower_bound)], ...
    'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

xlabel('Time (hours)');
ylabel('Position Magnitude (km)');
title('Position Magnitude vs Time (±1σ)');
legend('Location', 'best');
grid on;


subplot(2, 2, 2);
if ~isempty(measurement_residuals)
    % Get the number of updates
    num_updates = size(measurement_residuals, 2);
    update_numbers = 1:num_updates;
    
    % Plot all components with different colors
    hold on;
    bar_width = 0.12; % Width of each bar group
    colors = {'r', 'g', 'b', 'c', 'm', 'y'}; % Colors for each component
    labels = {'Px', 'Py', 'Pz', 'Vx', 'Vy', 'Vz'};
    
    for i = 1:6
        residual_norms = measurement_residuals(i, :);
        % Offset each bar slightly for grouped bar effect
        x_pos = update_numbers + (i-3.5) * bar_width;
        bar(x_pos, residual_norms, bar_width, 'FaceColor', colors{i}, 'DisplayName', labels{i});
    end
    
    xlabel('Update Number');
    ylabel('Residual Norm');
    title('Innovation Residuals (All Components)');
    legend('Location', 'best');
    grid on;
    hold off;
end


subplot(2, 2, 3);
uncertainty_magnitude = zeros(1, numTimeSteps);
for i = 1:numTimeSteps
    uncertainty_magnitude(i) = sqrt(trace(ISS_position_covs(:, :, i))) / 1000;
end

plot(time_hours, uncertainty_magnitude, 'b-', 'LineWidth', 2);
hold on;

% Only plot update markers if we have updates and they're within bounds
if ~isempty(update_indices)
    for i = 1:length(update_indices)
        idx = update_indices(i);
        if idx <= length(time_hours) && idx <= length(uncertainty_magnitude)
            plot(time_hours(idx), uncertainty_magnitude(idx), 'ro', 'MarkerSize', 8, ...
                 'MarkerFaceColor', 'r');
        end
    end
end

xlabel('Time (hours)');
ylabel('Position Uncertainty (km)');
title('Position Uncertainty Trace with Updates');
legend('Uncertainty', 'Measurement Updates', 'Location', 'best');
grid on;

% Log Determinants plot section
subplot(2, 2, 4);
if ~isempty(update_indices)
    log_determinants_at_updates = zeros(1, length(update_indices));
    for i = 1:length(update_indices)
        idx = update_indices(i);
        if idx <= numTimeSteps
            det_val = det(ISS_position_covs(:, :, idx)); % Keep in m^6
            log_determinants_at_updates(i) = log(det_val); % Natural logarithm
        end
    end
    bar(1:length(log_determinants_at_updates), log_determinants_at_updates, 'g');
    xlabel('Update Number');
    ylabel('ln(Covariance Determinant) (m^6)');
    title('Natural Log Covariance Determinant at Updates');
    grid on;
end

%{
%% Plot 11: Log Determinants
subplot(2, 5, 9);
%residual_norms = sqrt(sum(predicted_measurements.^2, 1));
residual_norms = log(abs(all_determinants));
bar(1:length(residual_norms), residual_norms);
xlabel('Update Number');
ylabel('Residual Norm');
title('log Determinants');
grid on;
%}
%{
subplot(2, 5, 10);
bar(1:length(residual_norms), mahalanobis_distances);
xlabel('Update Number');
ylabel('Mahalanobis Distance');
title('Overall Innovation (Mahalanobis Norm)');
grid on;
%}
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


function [mu_est, P_est, S, y] = MakeEstimates(mu_prior, P_prior, z, R, z_hat, H)
    y = z - z_hat;
    S = H * P_prior * H' + R;
    K = P_prior * H' / S;
    mu_est = mu_prior + K * y;
    
    I = eye(size(P_prior));
    A = I - K * H;
    P_est = A * P_prior;

end

function [ISS_mixture_cov, ISS_mixture_mean, ISS_gaussian_ECI_means, ISS_gaussian_ECI_covs, ISS_gaussian_means, ISS_gaussian_covs, Weights] = GenerateGSF(mu, P, N, numSamples, startTime, stopTime)
sampleTime = 30; % seconds
timeVector = startTime:seconds(sampleTime):stopTime;
numTimeSteps = length(timeVector);

alpha = 1e-3;
beta = 2;
kappa = -3;
n = 6; 
n_sigma = 2*n + 1; 

data = mvnrnd(mu, P, numSamples);


options = statset('MaxIter', 1000, 'Display', 'off');
fprintf('Number of samples: %d\n', size(data', 1));
fprintf('Number of dimensions: %d\n', size(data', 2));
fprintf('Requested GMM components: %d\n', N);
gm = fitgmdist(data, N, 'Options', options, 'RegularizationValue', 1e-20);

N_components = gm.NumComponents; 


ISS_gaussian_means = gm.mu';            
ISS_gaussian_covs = gm.Sigma;       
Weights = gm.ComponentProportion;
%Weights = ones(1, N_components) / N;
ISS_all_sigma_points = zeros(n, n_sigma, N_components);
ISS_all_Wm = zeros(n_sigma, N_components);
ISS_all_Wc = zeros(n_sigma, N_components);

for i = 1:N_components
    mu_i = ISS_gaussian_means(:, i);   
    P_i = ISS_gaussian_covs(:, :, i);  
    
    [sigma_points_i, Wm_i, Wc_i] = generate_sigma_points(mu_i, P_i, alpha, beta, kappa);
    

    keplerian_sigma_points_i = zeros(n, n_sigma);
    for j = 1:n_sigma
        kep_j = equinoctialToKeplerian(sigma_points_i(:, j));
        keplerian_sigma_points_i(:, j) = kep_j';
    end
    
    ISS_all_sigma_points(:, :, i) = keplerian_sigma_points_i;
    ISS_all_Wm(:, i) = Wm_i;
    ISS_all_Wc(:, i) = Wc_i;
end


ISS_propagated_points = cell(1, n_sigma, N_components);
sc_ISS = satelliteScenario(startTime, stopTime, sampleTime);

for j = 1:N_components
    for i = 1:n_sigma
        a = ISS_all_sigma_points(1, i, j);
        e = ISS_all_sigma_points(2, i, j);
        inc = ISS_all_sigma_points(3, i, j);
        raan = ISS_all_sigma_points(4, i, j);
        argP = ISS_all_sigma_points(5, i, j);
        theta = ISS_all_sigma_points(6, i, j);
        
        sat = satellite(sc_ISS, ...
            a, ...
            e, ...
            inc, ...
            raan, ...
            argP, ...
            theta, ...
            "OrbitPropagator", "two-body-keplerian");
        
        ISS_propagated_points{1, i, j} = sat;
    end
end

ISS_positions = cell(1, n_sigma, N_components);
ISS_velocities = cell(1, n_sigma, N_components);

for j = 1:N_components
    for i = 1:n_sigma
        sat = ISS_propagated_points{1, i, j};
        [pos, vel, ~] = states(sat);
        ISS_positions{1, i, j} = pos;
        ISS_velocities{1, i, j} = vel;
    end
end

ISS_predicted_mean_cartesian = zeros(6, numTimeSteps, N_components);
ISS_predicted_cov_cartesian = zeros(6, 6, numTimeSteps, N_components);
ISS_sigma_points_cartesian = zeros(6, n_sigma, numTimeSteps, N_components);

for j = 1:N_components
    for t = 1:numTimeSteps
        for i = 1:n_sigma
            r_t = ISS_positions{1,i,j}(:,t);
            v_t = ISS_velocities{1,i,j}(:,t);
            ISS_sigma_points_cartesian(:,i, t, j) = [r_t; v_t];
        end
        
        Wm_j = ISS_all_Wm(:,j);
        Wc_j = ISS_all_Wc(:,j);
        
 
        mean_jt = ISS_sigma_points_cartesian(:, :, t, j) * Wm_j; 
        

        cov_jt = zeros(6, 6);
        for i = 1:n_sigma
            diff = ISS_sigma_points_cartesian(:,i, t, j) - mean_jt;
            cov_jt = cov_jt + Wc_j(i) * (diff * diff');
        end
        
        ISS_predicted_mean_cartesian(:,t,j) = mean_jt;
        ISS_predicted_cov_cartesian(:,:,t,j) = cov_jt;
    end
end

ISS_gaussian_ECI_means = ISS_predicted_mean_cartesian; 
ISS_gaussian_ECI_covs = ISS_predicted_cov_cartesian;  
gaussian_weights = Weights;               

ISS_mixture_mean = zeros(6, numTimeSteps);
ISS_mixture_cov = zeros(6, 6, numTimeSteps);

for t = 1:numTimeSteps
    mean_mix_t = zeros(6,1);
    for j = 1:N_components
        mean_j = ISS_gaussian_ECI_means(:,t,j);
        w_j = gaussian_weights(j);
        mean_mix_t = mean_mix_t + w_j * mean_j;
    end
    ISS_mixture_mean(:,t) = mean_mix_t;
    
    cov_mix_t = zeros(6, 6);
    for j = 1:N_components
        mean_j = ISS_gaussian_ECI_means(:,t,j); 
        cov_j = ISS_gaussian_ECI_covs(:,:,t,j);
        w_j = gaussian_weights(j);
        
        diff = mean_j - mean_mix_t;

        cov_mix_t = cov_mix_t + w_j * (cov_j + diff * diff'); 
    end
    ISS_mixture_cov(:,:,t) = cov_mix_t;
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
    uncertainty.a = 1000;
    uncertainty.h = 1e-3;
    uncertainty.k = 1e-3;
    uncertainty.p = 1e-3;
    uncertainty.q = 1e-3;
    uncertainty.l = 1e-2*pi/180;
    %}
    perc.a = 0.0001;   % 5% of a
    perc.h = 0.005;   % 1% of h
    perc.k = 0.005;   % 1% of k
    perc.p = 0.005;   % 1% of p
    perc.q = 0.005;   % 1% of q
    perc.l = 0.00001;  % 0.5% of l
    %{
    perc.a = 0.0003;   % 5% of a
    perc.h = 0.01;   % 1% of h
    perc.k = 0.01;   % 1% of k
    perc.p = 0.005;   % 1% of p
    perc.q = 0.005;   % 1% of q
    perc.l = 0.0001;  % 0.5% of l
    %}
    % Calculate uncertainties as percentages of nominal values
    uncertainty.a = perc.a * nominal.a;
    uncertainty.h = perc.h * nominal.h;
    uncertainty.k = perc.k * nominal.k;
    uncertainty.p = perc.p * nominal.p;
    uncertainty.q = perc.q * nominal.q;
    uncertainty.l = perc.l * nominal.l;
    
    P = diag([uncertainty.a^2, ...
              uncertainty.h^2, ...
              uncertainty.k^2, ...
              uncertainty.p^2, ...
              uncertainty.q^2, ...
              uncertainty.l^2]);
end


function [sigma_points, Wm, Wc] = generate_sigma_points(mu, P, alpha, beta, kappa)
    mu = mu(:);
    n = length(mu);
    lambda = alpha^2 * (n + kappa) - n;
    gamma = sqrt(n + lambda);
    
    % Initialize outputs
    sigma_points = zeros(n, 2 * n + 1);
    Wm = ones(1, 2 * n + 1) * 1 / (2 * (n + lambda));
    Wc = Wm;
    Wm(1) = lambda / (n + lambda);
    Wc(1) = Wm(1) + (1 - alpha^2 + beta);
    sigma_points(:, 1) = mu;
    
    % Method: Force positive definiteness using eigendecomposition
    [V, D] = eig(P);
    D = diag(D);
    
    % Replace negative/zero eigenvalues with small positive values
    min_eigenvalue = 1e-6;
    D(D <= 0) = min_eigenvalue;
    
    % Reconstruct positive definite matrix
    P_pd = V * diag(D) * V';
    
    % Ensure symmetry (numerical errors can cause asymmetry)
    P_pd = (P_pd + P_pd') / 2;
    
    % Compute square root
    S = gamma * V * diag(sqrt(D));
    
    % Generate sigma points
    for i = 1:n
        sigma_points(:, i + 1) = mu + S(:, i);     % Positive direction
        sigma_points(:, i + 1 + n) = mu - S(:, i); % Negative direction
    end
end
%{
function [sigma_points, Wm, Wc] = generate_sigma_points(mu, P, alpha, beta, kappa)
    mu = mu(:);
    n = length(mu);
    
    lambda = alpha^2 * (n + kappa) - n;
    gamma = sqrt(n + lambda);
    
    sigma_points = zeros(n, 2 * n + 1);
    
    Wm = ones(1, 2 * n + 1) * 1 / (2 * (n + lambda));
    Wc = Wm;
    
    Wm(1) = lambda / (n + lambda);
    Wc(1) = Wm(1) + (1 - alpha^2 + beta);
    
    sigma_points(:, 1) = mu;
    
    S = gamma * chol(P, 'lower');
    for i = 1:n
        sigma_points(:, i + 1) = mu + S(:, i);        % Positive direction
        sigma_points(:, i + 1 + n) = mu - S(:, i);    % Negative direction
    end
end
%}
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
nu  = rad2deg(nu);
keplerian = [a, e, i, Omega, omega, nu];

end



