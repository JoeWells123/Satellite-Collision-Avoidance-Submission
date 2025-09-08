global N lat lon numSamples startTime stopTime sampleTime timeVector numTimeSteps sc_mean mu_earth

%% Time Setup
startTime = datetime(2025, 9, 7, 10, 0, 0);
stopTime = startTime + hours(24);
sampleTime = 10; % seconds
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
lon = -40;

sc_mean = satelliteScenario(startTime, stopTime, sampleTime);

numDebris = length(ommFiles);
debris_weighted_mean = zeros(6, numTimeSteps, numDebris);
debris_weighted_covariance = zeros(6, 6, numTimeSteps, numDebris);

fprintf('\n=== Processing %d Debris Objects ===\n', numDebris);

for debrisIdx = 1:numDebris
    fprintf('\n--- Processing Debris Object %d/%d ---\n', debrisIdx, numDebris);
    
    thisFile = ommFiles(debrisIdx).name;
    debrisFilePath = fullfile(xmlFolder, thisFile);

    debris_data = readOMMFile(debrisFilePath);
    fprintf('  Loaded: %s (ID: %s)\n', debris_data.objectName, debris_data.objectId);

    sat = satellite(sc_mean, 'SL24.xml', 'OrbitPropagator','two-body-keplerian');

    elements1 = orbitalElements(sat);

    
    [nominal, uncertainty, debris_P] = convertOMMToSGP4Elements(elements1);
    
    
    mu_debris = [nominal.SemiMajorAxis;
        nominal.eccentricity;
        nominal.inclination;
        nominal.raan;
        nominal.argp;
        nominal.trueAnomaly];
    
    [debris_weighted_covariance_preupdate, debris_weighted_mean_preupdate, debris_gaussian_ECI_means, debris_gaussian_ECI_covs, ~, ~, debris_Weights, ~, ~, ~] = GenerateGSF(mu_debris, debris_P, N, numSamples, startTime, stopTime, 'SL24.xml');
    %[debris_weighted_covariance_preupdate, debris_weighted_mean_preupdate, debris_gaussian_ECI_means, debris_gaussian_ECI_covs, debris_ISS_gaussian_means, debris_ISS_gaussian_covs, debris_Weights] = GenerateGSF(mu_debris, debris_P);

    [debris_weighted_covariance(:, :, :, debrisIdx), debris_weighted_mean(:, :, debrisIdx), debris_measurement_residuals, debris_update_indices] = UpdateGSF(debris_weighted_covariance_preupdate, debris_weighted_mean_preupdate, debris_gaussian_ECI_means, debris_gaussian_ECI_covs, debris_Weights, mu_debris, debris_P, 'SL24.xml');

    %% Plot Debris Results
    fprintf('\n=== Plotting Debris Object %d ===\n', debrisIdx);
    plotResults(debris_weighted_covariance(:, :, :, debrisIdx), debris_weighted_mean(:, :, debrisIdx), debris_measurement_residuals, debris_update_indices, sprintf('Debris Object %d', debrisIdx), debrisIdx);
end

fprintf('\n=== Completed Processing All %d Debris Objects ===\n', numDebris);

%% Load ISS Data
ISS_omm_file = 'Trial.xml';
ISS_data = readOMMFile(ISS_omm_file);
fprintf('  Loaded: %s (ID: %s)\n', ISS_data.objectName, ISS_data.objectId);

sat = satellite(sc_mean, 'COS.xml', 'OrbitPropagator','two-body-keplerian');

elements1 = orbitalElements(sat);


[nominal, uncertainty, P_ISS] = convertOMMToSGP4Elements(elements1);


mu = [nominal.SemiMajorAxis;
    nominal.eccentricity;
    nominal.inclination;
    nominal.raan;
    nominal.argp;
    nominal.trueAnomaly];

P = P_ISS;

fprintf('\n=== Processing ISS ===\n');

[weighted_covariance, weighted_mean, gaussian_ECI_means, gaussian_ECI_covs, ~, ~, Weights, Satellite_sigma_point_files, Wm_all, Wc_all] = GenerateGSF(mu, P, N, numSamples, startTime, stopTime, 'COS.xml');
%[weighted_covariance, weighted_mean, gaussian_ECI_means, gaussian_ECI_covs, ISS_gaussian_means, ISS_gaussian_covs, Weights] = GenerateGSF(mu, P);

[weighted_covariance, weighted_mean, ISS_measurement_residuals, ISS_update_indices] = UpdateGSF(weighted_covariance, weighted_mean, gaussian_ECI_means, gaussian_ECI_covs, Weights, mu, P, 'COS.xml');

%% Plot ISS Results
fprintf('\n=== Plotting ISS Results ===\n');
plotResults(weighted_covariance, weighted_mean, ISS_measurement_residuals, ISS_update_indices, 'ISS', 0);

%% Plotting function
function plotResults(weighted_covariance, weighted_mean, measurement_residuals, update_indices, objectName, figureOffset)
    global timeVector startTime numTimeSteps
    
    ISS_positions = weighted_mean(1:3, :);
    ISS_position_covs = weighted_covariance(1:3, 1:3, :);
    ISS_positions_km = ISS_positions / 1000;
    time_hours = hours(timeVector - startTime);
    
    % Create figure with offset to avoid overlapping
    figure('Position', [100 + figureOffset*50, 100 + figureOffset*50, 1400, 1000]);
    
    % Position magnitude plot
    subplot(2, 2, 1);
    position_magnitude = sqrt(sum(ISS_positions_km.^2, 1));
    fprintf('Position magnitude size: %dx%d\n', size(position_magnitude));
    fprintf('Time hours size: %dx%d\n', size(time_hours));
    
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
    title(sprintf('%s: Position Magnitude vs Time (±1σ)', objectName));
    legend('Location', 'best');
    grid on;
    
    % Innovation residuals plot
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
        title(sprintf('%s: Innovation Residuals (All Components)', objectName));
        legend('Location', 'best');
        grid on;
        hold off;
    else
        text(0.5, 0.5, 'No measurement residuals available', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
        title(sprintf('%s: No Measurement Residuals', objectName));
    end
    
    % Position uncertainty trace plot
    subplot(2, 2, 3);
    uncertainty_magnitude = zeros(1, numTimeSteps);
    for i = 1:numTimeSteps
        uncertainty_magnitude(i) = sqrt(trace(ISS_position_covs(:, :, i))) / 1000;
    end
    
    plot(time_hours, uncertainty_magnitude, 'b-', 'LineWidth', 2);
    hold on;
    
    % Only plot update markers if we have updates and they're within bounds
    if ~isempty(update_indices)
        valid_updates = 0;
        for i = 1:length(update_indices)
            idx = update_indices(i);
            if idx <= length(time_hours) && idx <= length(uncertainty_magnitude)
                plot(time_hours(idx), uncertainty_magnitude(idx), 'ro', 'MarkerSize', 8, ...
                     'MarkerFaceColor', 'r');
                valid_updates = valid_updates + 1;
            end
        end
        if valid_updates > 0
            legend('Uncertainty', 'Measurement Updates', 'Location', 'best');
        else
            legend('Uncertainty', 'Location', 'best');
        end
    else
        legend('Uncertainty', 'Location', 'best');
    end
    
    xlabel('Time (hours)');
    ylabel('Position Uncertainty (km)');
    title(sprintf('%s: Position Uncertainty Trace with Updates', objectName));
    grid on;
    
    % Log Determinants plot
    subplot(2, 2, 4);
    if ~isempty(update_indices)
        valid_determinants = [];
        valid_update_nums = [];
        
        for i = 1:length(update_indices)
            idx = update_indices(i);
            if idx <= numTimeSteps
                det_val = det(ISS_position_covs(:, :, idx)); % Keep in m^6
                if det_val > 0  % Only plot positive determinants
                    valid_determinants(end+1) = log(det_val); % Natural logarithm
                    valid_update_nums(end+1) = i;
                end
            end
        end
        
        if ~isempty(valid_determinants)
            bar(valid_update_nums, valid_determinants, 'g');
            xlabel('Update Number');
            ylabel('ln(Covariance Determinant) (m^6)');
            title(sprintf('%s: Natural Log Covariance Determinant at Updates', objectName));
            grid on;
        else
            text(0.5, 0.5, 'No valid determinants to plot', ...
                'HorizontalAlignment', 'center', 'Units', 'normalized');
            title(sprintf('%s: No Valid Determinants', objectName));
        end
    else
        text(0.5, 0.5, 'No updates available', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
        title(sprintf('%s: No Updates Available', objectName));
    end
    
    % Add overall title to the figure
    sgtitle(sprintf('%s Analysis Results', objectName));
end


relative_positions = debris_weighted_mean - weighted_mean;
distance = sqrt(sum(relative_positions(1:3,:).^2, 1));
distance_km = distance / 1000;
[min_distance, min_idx] = min(distance_km);
min_time = timeVector(min_idx);

relative_positions_sgp4 = debris_weighted_mean - weighted_mean;
distance = sqrt(sum(relative_positions_sgp4.^2, 1)); 
disp(min(distance))
figure('Position', [100, 100, 1200, 800]);
hold on;
time_hours = hours(timeVector - startTime);
plot(time_hours, distance, 'b-', 'LineWidth', 1.8, 'DisplayName', 'SGP4 vs Two-Body (LEO)');

% Labels & aesthetics
xlabel('Time (hours)', 'FontSize', 12);
ylabel('Distance (m)', 'FontSize', 12);
title('Distance between SGP4/SDP4 and Two-Body Keplerian Models', 'FontSize', 14);

legend('Location', 'northwest');
grid on;
box on;

dr = (dot(relative_positions(1:3, min_idx)/1000, relative_positions(4:6, min_idx)/1000)/min_distance)*sampleTime;
disp(dr)
format long
disp(min_distance)
disp(min_time)

function [weighted_covariance, weighted_mean, measurement_residuals, update_indices] = UpdateGSF(weighted_covariance, weighted_mean, gaussian_ECI_means, gaussian_ECI_covs, Weights, mu, P, File)
global N lat lon numSamples startTime stopTime sampleTime timeVector numTimeSteps sc_mean mu_earth
measurement_data = mvnrnd(mu, P, numSamples);
options = statset('MaxIter', 1000, 'Display', 'off');
gm = fitgmdist(measurement_data, 1, 'Options', options, 'RegularizationValue', 1e-20);
measurement_means = gm.mu';   

elems = struct( ...
    'semiMajorAxis', measurement_means(1)/1000, ...
    'eccentricity', measurement_means(2), ...
    'inclination', measurement_means(3), ...
    'raan', measurement_means(4), ...
    'argPeriapsis', measurement_means(5), ...
    'trueAnomaly', measurement_means(6) ...
);
updateOMMFile(File, 'measurement.xml', elems);

Initial_sat = satellite(sc_mean, 'measurement.xml', "OrbitPropagator", "two-body-keplerian");
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

predicted_measurements = [];
measurement_residuals = [];
update_indices = [];
update_counter = 0;
all_determinants = [];
mahalanobis_distances = [];
Neff_vals = [];
sc_measurement = satelliteScenario(startTime, stopTime, sampleTime);
gs_measurement = groundStation(sc_measurement, lat, lon);

all_sigma_point_files = {};

for interval_idx = 1:height(intvls)
    fprintf('\nProcessing interval %d of %d...\n', interval_idx, height(intvls));
    
    interval_start = intvls.StartTime(interval_idx);
    interval_end = intvls.EndTime(interval_idx);
    
    interval_mask = (timeVector >= interval_start) & (timeVector <= interval_end);
    interval_time_indices = find(interval_mask);
    interval_time_indices = interval_time_indices(1:30:end);

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
        
        R = diag([5^2,5^2,5^2, 1e-1^2,1e-1^2,1e-1^2]);
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
                [a, ecc, incl, RAAN, argp, nu, truelon, arglat, lonper] = ijk2keplerian(ISS_all_sigma_points(1:3, j, i), ISS_all_sigma_points(4:6, j, i));
                keplerian_sigma_points(:, j, i) = [a, ecc, incl, RAAN, argp, nu]';
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

end
function l = compute_weight(z, zObs, R, N)

    % Compute weights in log-space for numerical stability
    RSqrtm = sqrtm(R);
    nu = RSqrtm * (z - zObs);

    % Compute log-likelihood
    log_l = -0.5 * sum(nu.^2, 1);

    % To avoid underflow, normalize in log-space
    max_log_l = max(log_l);  
    l = exp(log_l - max_log_l);
    
    % Optionally normalize to sum to 1
    l = l / sum(l);
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

function [ISS_mixture_cov, ISS_mixture_mean, ISS_gaussian_means, ISS_gaussian_covs, ISS_mu_components, ISS_P_components, gaussian_weights, updated_files, Wm_all, Wc_all] = GenerateGSF(mu, P, N, numSamples, startTime, stopTime, File)
    sampleTime = 10; % seconds
    timeVector = startTime:seconds(sampleTime):stopTime;
    numTimeSteps = length(timeVector);
    
    % Unscented transform parameters
    alpha = 1e-3;
    beta = 2;
    kappa = -3;
    n = 6;
    n_sigma = 2*n + 1;
    
    % Generate samples and fit Gaussian Mixture Model
    data = mvnrnd(mu, P, numSamples);
    options = statset('MaxIter', 10000, 'Display', 'off');
    fprintf('Number of samples: %d\n', size(data', 1));
    fprintf('Number of dimensions: %d\n', size(data', 2));
    fprintf('Requested GMM components: %d\n', N);
    gm = fitgmdist(data, N, 'Options', options, 'RegularizationValue', 1e-30);
    
    N_components = gm.NumComponents;
    
    % Extract GMM parameters
    ISS_mu_components = gm.mu';
    ISS_P_components = gm.Sigma;
    
    % Set equal weights for all components
    gaussian_weights = ones(1, N_components) / N_components;
    
    % Initialize storage for sigma points and weights for each component
    ISS_all_sigma_points = zeros(n, n_sigma, N_components);
    Wm_all = zeros(n_sigma, N_components);
    Wc_all = zeros(n_sigma, N_components);
    
    % Generate sigma points for each Gaussian component
    for i = 1:N_components
        mu_i = ISS_mu_components(:, i);
        P_i = ISS_P_components(:, :, i);
        
        [sigma_points_i, Wm_i, Wc_i] = generate_sigma_points(mu_i, P_i, alpha, beta, kappa);
        
        ISS_all_sigma_points(:, :, i) = sigma_points_i;
        Wm_all(:, i) = Wm_i;
        Wc_all(:, i) = Wc_i;
    end
    
    % Initialize cell array to store updated filenames for each component and sigma point
    updated_files = cell(N_components, n_sigma);
    
    % Propagate sigma points for each component
    ISS_propagated_points = cell(N_components, n_sigma);
    sc_ISS = satelliteScenario(startTime, stopTime, sampleTime);
    
    for j = 1:N_components
        for i = 1:n_sigma
            elems = struct( ...
                'semiMajorAxis', ISS_all_sigma_points(1, i, j)/1000, ...
                'eccentricity', ISS_all_sigma_points(2, i, j), ...
                'inclination', ISS_all_sigma_points(3, i, j), ...
                'raan', ISS_all_sigma_points(4, i, j), ...
                'argPeriapsis', ISS_all_sigma_points(5, i, j), ...
                'trueAnomaly', ISS_all_sigma_points(6, i, j) ...   
            );
            
            % Create unique filename for each sigma point and component
            updated_filename = sprintf('ISS_updated_component_%d_sigma_%d.xml', j, i);
            updated_files{j, i} = updated_filename; % Store filename in the list
            updateOMMFile(File, updated_filename, elems);
            sat = satellite(sc_ISS, updated_filename, 'OrbitPropagator','two-body-keplerian');
            ISS_propagated_points{j, i} = sat;
        end
    end
    
    % Collect propagated positions and velocities for each component
    ISS_positions = cell(N_components, n_sigma);
    ISS_velocities = cell(N_components, n_sigma);
    
    for j = 1:N_components
        for i = 1:n_sigma
            sat = ISS_propagated_points{j, i};
            [pos, vel, ~] = states(sat);
            ISS_positions{j, i} = pos;
            ISS_velocities{j, i} = vel;
        end
    end
    
    % Compute mean and covariance in Cartesian coordinates for each component
    ISS_gaussian_means = zeros(6, numTimeSteps, N_components);
    ISS_gaussian_covs = zeros(6, 6, numTimeSteps, N_components);
    
    for j = 1:N_components
        for t = 1:numTimeSteps
            sigma_points_t = zeros(6, n_sigma);
            for i = 1:n_sigma
                sigma_points_t(:, i) = [ISS_positions{j, i}(:, t); ISS_velocities{j, i}(:, t)];
            end
            
            % Mean for component j at time t
            Wm_j = Wm_all(:, j);
            mean_t = sigma_points_t * Wm_j;
            ISS_gaussian_means(:, t, j) = mean_t;
            
            % Covariance for component j at time t
            Wc_j = Wc_all(:, j);
            cov_t = zeros(6, 6);
            for i = 1:n_sigma
                diff = sigma_points_t(:, i) - mean_t;
                cov_t = cov_t + Wc_j(i) * (diff * diff');
            end
            ISS_gaussian_covs(:, :, t, j) = cov_t;
        end
    end
    
    % Compute mixture mean and covariance
    ISS_mixture_mean = zeros(6, numTimeSteps);
    ISS_mixture_cov = zeros(6, 6, numTimeSteps);
    
    for t = 1:numTimeSteps
        % Mixture mean at time t
        mean_mix_t = zeros(6, 1);
        for j = 1:N_components
            mean_j = ISS_gaussian_means(:, t, j);
            w_j = gaussian_weights(j);
            mean_mix_t = mean_mix_t + w_j * mean_j;
        end
        ISS_mixture_mean(:, t) = mean_mix_t;
        
        % Mixture covariance at time t
        cov_mix_t = zeros(6, 6);
        for j = 1:N_components
            mean_j = ISS_gaussian_means(:, t, j);
            cov_j = ISS_gaussian_covs(:, :, t, j);
            w_j = gaussian_weights(j);
            
            diff = mean_j - mean_mix_t;
            cov_mix_t = cov_mix_t + w_j * (cov_j + diff * diff');
        end
        ISS_mixture_cov(:, :, t) = cov_mix_t;
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



function [nominal, uncertainty, P] = convertOMMToSGP4Elements(elements1)
    nominal = struct();
    uncertainty = struct();
   

    nominal.SemiMajorAxis = elements1.SemiMajorAxis; % Convert km to m
    nominal.eccentricity = elements1.Eccentricity;
    nominal.inclination = elements1.Inclination; % degrees
    nominal.raan = elements1.RightAscensionOfAscendingNode; % degrees
    nominal.argp = elements1.ArgumentOfPeriapsis; % degrees
    nominal.trueAnomaly = elements1.TrueAnomaly;

    %{
    perc.semiMajorAxis = 0.00001;   % 5% of a
    perc.eccentricity = 0.0001;   % 1% of h
    perc.inclination = 0.001;   % 1% of k
    perc.raan = 0.00001;   % 1% of p
    perc.argp = 0.001;   % 1% of q
    perc.trueAnomaly = 0.0001;  % 0.5% of l
    %}
    perc.semiMajorAxis = 0.0000001;   % 5% of a
    perc.eccentricity = 0.00001;   % 1% of h
    perc.inclination = 0.00001;   % 1% of k
    perc.raan = 0.00001;   % 1% of p
    perc.argp = 0.00001;   % 1% of q
    perc.trueAnomaly = 0.00001;  % 0.5% of l

    
    uncertainty.SemiMajorAxis = perc.semiMajorAxis * nominal.SemiMajorAxis;
    uncertainty.eccentricity = perc.eccentricity * nominal.eccentricity;
    uncertainty.inclination = perc.inclination * nominal.inclination;
    uncertainty.raan = perc.raan * nominal.raan;
    uncertainty.argPeriapsis = perc.argp * nominal.argp;
    uncertainty.trueAnomaly = perc.trueAnomaly * nominal.trueAnomaly;
    


    P = diag([uncertainty.SemiMajorAxis^2, ...
        uncertainty.eccentricity^2, ...
        uncertainty.inclination^2, ...
        uncertainty.raan^2, ...
        uncertainty.argPeriapsis^2, ...
        uncertainty.trueAnomaly^2]);
end

function updateOMMFile(inputFile, outputFile, elements)
% updateOMMFile  Update Keplerian orbital elements in an OMM XML file
%
%   updateOMMFile(inputFile, outputFile, elements)
%
%   INPUTS:
%     inputFile   - path to original OMM XML file
%     outputFile  - path to save updated file
%     elements    - struct with fields:
%                      semiMajorAxis [km] (primary)
%                      eccentricity
%                      inclination [deg]
%                      raan [deg]
%                      argPeriapsis [deg]
%                      trueAnomaly [deg]

    if ~exist(inputFile, 'file')
        error('Input file not found: %s', inputFile);
    end

    mu_earth = 3.986004418e14; % m^3/s^2
    
    % Check if semi-major axis is provided
    if ~isfield(elements, 'semiMajorAxis') || isempty(elements.semiMajorAxis)
        error('semiMajorAxis must be provided in elements struct.');
    end
    
    % Check if true anomaly is provided
    if ~isfield(elements, 'trueAnomaly') || isempty(elements.trueAnomaly)
        error('trueAnomaly must be provided in elements struct.');
    end
    
    % Calculate mean motion from semi-major axis
    a_m = elements.semiMajorAxis * 1000; % km -> m
    n_rad_s = sqrt(mu_earth / a_m^3);    % rad/s
    elements.meanMotion = n_rad_s * 86400 / (2*pi); % rad/s -> rev/day
    
    % Calculate period from semi-major axis
    elements.period = 2*pi / n_rad_s / 60; % minutes
    
    % Convert true anomaly to mean anomaly
    elements.meanAnomaly = trueToMeanAnomaly(elements.trueAnomaly, elements.eccentricity);
    
    try
        xmlDoc = xmlread(inputFile);

        % Update scalar tags
        setXMLNodeValue(xmlDoc, 'MEAN_MOTION', elements.meanMotion);
        setXMLNodeValue(xmlDoc, 'ECCENTRICITY', elements.eccentricity);
        setXMLNodeValue(xmlDoc, 'INCLINATION', elements.inclination);
        setXMLNodeValue(xmlDoc, 'RA_OF_ASC_NODE', elements.raan);
        setXMLNodeValue(xmlDoc, 'ARG_OF_PERICENTER', elements.argPeriapsis);
        setXMLNodeValue(xmlDoc, 'MEAN_ANOMALY', elements.meanAnomaly);

        % Update USER_DEFINED SEMIMAJOR_AXIS and PERIOD
        userNodes = xmlDoc.getElementsByTagName('USER_DEFINED');
        smaUpdated = false;
        periodUpdated = false;
        for k = 0:userNodes.getLength-1
            node = userNodes.item(k);
            paramName = char(node.getAttribute('parameter'));
            if strcmp(paramName, 'SEMIMAJOR_AXIS')
                node.getFirstChild.setNodeValue(num2str(elements.semiMajorAxis,'%.6f'));
                smaUpdated = true;
            elseif strcmp(paramName, 'PERIOD')
                node.getFirstChild.setNodeValue(num2str(elements.period,'%.6f'));
                periodUpdated = true;
            end
        end

        % Create SEMIMAJOR_AXIS if missing
        if ~smaUpdated
            root = xmlDoc.getDocumentElement;
            userParamsNodes = xmlDoc.getElementsByTagName('USER_DEFINED_PARAMETERS');
            if userParamsNodes.getLength() == 0
                userParamsSection = xmlDoc.createElement('USER_DEFINED_PARAMETERS');
                root.appendChild(userParamsSection);
            else
                userParamsSection = userParamsNodes.item(0);
            end
            newSmaNode = xmlDoc.createElement('USER_DEFINED');
            newSmaNode.setAttribute('parameter','SEMIMAJOR_AXIS');
            newSmaNode.setAttribute('units','km');
            newSmaNode.appendChild(xmlDoc.createTextNode(num2str(elements.semiMajorAxis,'%.6f')));
            userParamsSection.appendChild(newSmaNode);
        end

        % Create PERIOD if missing
        if ~periodUpdated
            root = xmlDoc.getDocumentElement;
            userParamsNodes = xmlDoc.getElementsByTagName('USER_DEFINED_PARAMETERS');
            if userParamsNodes.getLength() == 0
                userParamsSection = xmlDoc.createElement('USER_DEFINED_PARAMETERS');
                root.appendChild(userParamsSection);
            else
                userParamsSection = userParamsNodes.item(0);
            end
            newPeriodNode = xmlDoc.createElement('USER_DEFINED');
            newPeriodNode.setAttribute('parameter','PERIOD');
            newPeriodNode.setAttribute('units','min');
            newPeriodNode.appendChild(xmlDoc.createTextNode(num2str(elements.period,'%.6f')));
            userParamsSection.appendChild(newPeriodNode);
        end

        % Fix namespace issues
        root = xmlDoc.getDocumentElement;
        if ~root.hasAttribute('xmlns:xsi')
            root.setAttribute('xmlns:xsi','http://www.w3.org/2001/XMLSchema-instance');
        end
        if root.hasAttribute('xsi:noNamespaceSchemaLocation')
            root.removeAttribute('xsi:noNamespaceSchemaLocation');
        end

        saveXML(xmlDoc, outputFile);

    catch ME
        error('Failed to update OMM file: %s', ME.message);
    end
end

%% --- Helper Functions ---
function success = setXMLNodeValue(xmlDoc, tagName, newValue)
    success = false;
    nodes = xmlDoc.getElementsByTagName(tagName);
    if nodes.getLength() > 0
        node = nodes.item(0);
        if ~isempty(node.getFirstChild)
            node.getFirstChild.setNodeValue(num2str(newValue,'%.6f'));
            success = true;
        end
    end
end

function saveXML(xmlDoc, filename)
    import javax.xml.transform.*
    import javax.xml.transform.dom.*
    import javax.xml.transform.stream.*
    transformer = TransformerFactory.newInstance.newTransformer;
    transformer.setOutputProperty(OutputKeys.INDENT,'yes');
    transformer.setOutputProperty('{http://xml.apache.org/xslt}indent-amount','2');
    source = DOMSource(xmlDoc);
    result = StreamResult(java.io.File(filename));
    transformer.transform(source,result);
end

function M = trueToMeanAnomaly(nu, e)
    % Convert true anomaly to mean anomaly
    % nu: true anomaly in degrees
    % e: eccentricity
    % M: mean anomaly in degrees
    
    nu_rad = deg2rad(nu);
    
    % Calculate eccentric anomaly from true anomaly
    E = 2 * atan(sqrt((1-e)/(1+e)) * tan(nu_rad/2));
    
    % Calculate mean anomaly from eccentric anomaly
    M_rad = E - e * sin(E);
    
    % Convert back to degrees and ensure positive
    M = rad2deg(M_rad);
    M = mod(M, 360);
end
