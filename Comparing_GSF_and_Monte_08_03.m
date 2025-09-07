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

%% ISS Gaussian Sum Filter Parameters
mu = [ISS_nominal.a;
    ISS_nominal.h;
    ISS_nominal.k;
    ISS_nominal.p;
    ISS_nominal.q;
    ISS_nominal.l];

P = P_ISS;

%% Define different numbers of Gaussian components to test
N_values = [10];
%N_values = [1, 5, 10, 15];
fprintf('\n=== Processing Multiple N Values ===\n');

% Store results for each N value
GSF_results = cell(length(N_values), 1);
bhattacharyya_distances = cell(length(N_values), 1);
euclidean_distances = cell(length(N_values), 1); % NEW: Store Euclidean distances
computation_times = zeros(length(N_values), 1); % Store computation times

% Generate particle filter reference (or use the highest N as reference)
ref_positions = Particle_ISS_positions(1:3, :);
ref_position_covs = Particle_ISS_position_covs(1:3, 1:3, :);
ref_position_skewness = Particle_ISS_positions_skewness(1:3, :);
ref_position_kurtosis = Particle_ISS_positions_kurtosis(1:3, :);
% Process each N value
for idx = 1:length(N_values)
    N = N_values(idx);
    fprintf('Processing N = %d Gaussians...\n', N);
    
    % Start timing
    tic;
    
    [weighted_covariance, weighted_mean, gaussian_ECI_means, gaussian_ECI_covs, ISS_gaussian_means, ISS_gaussian_covs, Weights] = GenerateGSF(mu, P, N, numSamples, startTime, stopTime);
   
    % End timing
    computation_times(idx) = toc;
    
    fprintf('  Computation time for N=%d: %.3f seconds\n', N, computation_times(idx));
    
    timeVector = datetime(timeVector, 'TimeZone', '');

    % Store results
    GSF_results{idx}.weighted_covariance = weighted_covariance;
    GSF_results{idx}.weighted_mean = weighted_mean;
    GSF_results{idx}.N = N;
    
    % Extract position data
    ISS_positions = weighted_mean(1:3, :);
    ISS_position_covs = weighted_covariance(1:3, 1:3, :);
    
    % Calculate Bhattacharyya distance against reference
    d_bhatt = zeros(1, numTimeSteps);
    % NEW: Calculate Euclidean distance against reference
    d_eucl = zeros(1, numTimeSteps);
    
    for k = 1:numTimeSteps
        % Bhattacharyya distance calculation (existing)
        mu1 = ISS_positions(:, k);
        Sigma1 = ISS_position_covs(:, :, k);
        mu2 = ref_positions(:, k);
        Sigma2 = ref_position_covs(:, :, k);
        
        d_bhatt(k) = bhattacharyya_distance(mu1, Sigma1, mu2, Sigma2);
        
        % NEW: Euclidean distance calculation (in meters)
        d_eucl(k) = norm(ISS_positions(:, k) - ref_positions(:, k));
    end
    
    bhattacharyya_distances{idx} = d_bhatt;
    euclidean_distances{idx} = d_eucl; % NEW: Store Euclidean distances
    
    fprintf('  Mean Bhattacharyya distance: %.6f\n', mean(d_bhatt));
    fprintf('  Max Bhattacharyya distance: %.6f\n', max(d_bhatt));
    fprintf('  Mean Euclidean distance: %.3f m\n', mean(d_eucl)); % NEW
    fprintf('  Max Euclidean distance: %.3f m\n', max(d_eucl));   % NEW
    fprintf('  RMS Euclidean error: %.3f m\n', sqrt(mean(d_eucl.^2))); % NEW
end

%% Display computation time statistics
fprintf('\n=== Computation Time Analysis ===\n');
fprintf('%-8s %-15s\n', 'N', 'Time (seconds)');
fprintf('%-8s %-15s\n', '---', '-------------');
for idx = 1:length(N_values)
    fprintf('%-8d %-15.3f\n', N_values(idx), computation_times(idx));
end

%% Plotting
timeVector = datetime(timeVector, 'TimeZone', '');
time_hours = hours(timeVector - startTime);

% Convert reference positions to km for plotting
ref_positions_km = ref_positions / 1000;
earthRadius_km = earthRadius / 1000;

%% Enhanced Plot 1: Reference ISS Trajectory 
figure('Position', [100, 100, 1400, 1000]);
% Plot 1: Mean Bhattacharyya Distance vs N
subplot(2, 2, 1);
mean_distances = cellfun(@mean, bhattacharyya_distances);
semilogx(N_values, mean_distances, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Mean Bhattacharyya Distance');
title('GS-UKF Convergence: Mean Distance vs N');
grid on;

% Plot 2: Max Bhattacharyya Distance vs N  
subplot(2, 2, 2);
max_distances = cellfun(@max, bhattacharyya_distances);
semilogx(N_values, max_distances, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Maximum Bhattacharyya Distance');
title('GS-UKF Convergence: Maximum Distance vs N');
grid on;

% Plot 3: Standard deviation of distances vs N
subplot(2, 2, 3);
std_distances = cellfun(@std, bhattacharyya_distances);
semilogx(N_values, std_distances, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Std Dev of Bhattacharyya Distance');
title('GS-UKF Variability: Std Dev vs N');
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
title('GS-UKF Distance Evolution: All N Values');
legend('Location', 'best');
grid on;


figure('Position', [150, 150, 1400, 1000]);

% Plot 1: Mean Euclidean Distance vs N
subplot(2, 2, 1);
mean_distances = cellfun(@mean, euclidean_distances);
semilogx(N_values, mean_distances, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Mean Euclidean Distance (m)');
title('GS-UKF Convergence: Mean Euclidean Distance vs N');
grid on;

% Plot 2: Max Euclidean Distance vs N
subplot(2, 2, 2);
max_distances = cellfun(@max, euclidean_distances);
semilogx(N_values, max_distances, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Max Euclidean Distance (m)');
title('GS-UKF Convergence: Max Euclidean Distance vs N');
grid on;

% Plot 3: Std Dev of Euclidean Distance vs N
subplot(2, 2, 3);
std_distances = cellfun(@std, euclidean_distances);
semilogx(N_values, std_distances, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Std Dev of Euclidean Distance (m)');
title('GS-UKF Variability: Std Dev vs N');
grid on;

% Plot 4: Euclidean Distance Evolution Over Time
subplot(2, 2, 4);
hold on;
for idx = 1:length(N_values)
    d = euclidean_distances{idx};
    plot(time_hours, d, 'LineWidth', 1.5, 'DisplayName', sprintf('N=%d', N_values(idx)));
end
xlabel('Time (hours)');
ylabel('Euclidean Distance (m)');
title('GS-UKF Distance Evolution: All N Values');
legend('Location', 'best');
grid on;


%% NEW: Computation Time Analysis Plot
figure('Position', [50, 50, 1400, 900]);

% Plot 1: Computation Time vs Number of Gaussians
subplot(2, 2, 1);
plot(N_values, computation_times, 'bo-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
xlabel('Number of Gaussian Components (N)');
ylabel('Computation Time (seconds)');
title('Computation Time vs Number of Gaussians');
grid on;

% Add annotations with time values
for idx = 1:length(N_values)
    text(N_values(idx), computation_times(idx) + max(computation_times)*0.02, ...
         sprintf('%.2fs', computation_times(idx)), ...
         'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end

% Plot 2: Semi-log plot of computation time
subplot(2, 2, 2);
semilogy(N_values, computation_times, 'ro-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('Number of Gaussian Components (N)');
ylabel('Computation Time (seconds) - Log Scale');
title('Computation Time vs N (Log Scale)');
grid on;

% Plot 3: Time complexity analysis (if more than 2 points)
subplot(2, 2, 3);
if length(N_values) > 2
    % Fit polynomial to estimate complexity
    p = polyfit(log(N_values), log(computation_times), 1);
    complexity_order = p(1);
    
    % Plot actual times and fitted line
    loglog(N_values, computation_times, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    hold on;
    N_fit = linspace(min(N_values), max(N_values), 100);
    time_fit = exp(p(2)) * N_fit.^complexity_order;
    loglog(N_fit, time_fit, 'g--', 'LineWidth', 2);
    
    xlabel('Number of Gaussian Components (N) - Log Scale');
    ylabel('Computation Time (seconds) - Log Scale');
    title(sprintf('Time Complexity Analysis (O(N^{%.2f}))', complexity_order));
    legend('Actual Times', sprintf('Fitted: O(N^{%.2f})', complexity_order), 'Location', 'best');
    grid on;
else
    % Simple linear plot if only 2 points
    plot(N_values, computation_times, 'go-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    xlabel('Number of Gaussian Components (N)');
    ylabel('Computation Time (seconds)');
    title('Time Complexity Analysis');
    grid on;
    
    % Calculate growth rate
    if length(N_values) == 2
        growth_rate = (computation_times(2) - computation_times(1)) / (N_values(2) - N_values(1));
        text(mean(N_values), mean(computation_times), ...
             sprintf('Growth: %.3f s/component', growth_rate), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', ...
             'BackgroundColor', 'white', 'EdgeColor', 'black');
    end
end

% Plot 4: Performance efficiency (Time per Gaussian component)
subplot(2, 2, 4);
time_per_component = computation_times ./ N_values;
bar(N_values, time_per_component, 'FaceColor', [0.7, 0.3, 0.9], 'EdgeColor', 'black', 'LineWidth', 1.5);
xlabel('Number of Gaussian Components (N)');
ylabel('Time per Component (seconds)');
title('Computational Efficiency (Time per Gaussian)');
grid on;

% Add value labels on bars
for idx = 1:length(N_values)
    text(N_values(idx), time_per_component(idx) + max(time_per_component)*0.02, ...
         sprintf('%.3f', time_per_component(idx)), ...
         'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
end

%% Additional Analysis Plot: Convergence Analysis with Timing
figure('Position', [200, 200, 1200, 800]);

% Plot 1: Mean Bhattacharyya Distance vs N
subplot(2, 2, 1);
mean_bhatt_distances = cellfun(@mean, bhattacharyya_distances);
semilogx(N_values, mean_bhatt_distances, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Mean Bhattacharyya Distance');
title('Convergence: Mean Bhattacharyya Distance vs N');
grid on;

% Plot 2: NEW - Mean Euclidean Distance vs N
subplot(2, 2, 2);
mean_eucl_distances = cellfun(@mean, euclidean_distances);
semilogx(N_values, mean_eucl_distances, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Mean Euclidean Distance Error (m)');
title('Convergence: Mean Euclidean Error vs N');
grid on;

% Plot 3: Performance Trade-off (Euclidean Accuracy vs Time)
subplot(2, 2, 3);
scatter(computation_times, mean_eucl_distances, 100, N_values, 'filled', 'MarkerEdgeColor', 'black');
colorbar;
xlabel('Computation Time (seconds)');
ylabel('Mean Euclidean Distance Error (m)');
title('Performance Trade-off: Euclidean Accuracy vs Time');
grid on;

% Add labels for each point
for idx = 1:length(N_values)
    text(computation_times(idx) + max(computation_times)*0.02, mean_eucl_distances(idx), ...
         sprintf('N=%d', N_values(idx)), 'FontSize', 10, 'FontWeight', 'bold');
end

% Plot 4: Distance evolution over time for all N values (Euclidean)
subplot(2, 2, 4);
hold on;
for idx = 1:length(N_values)
    d_eucl_km = euclidean_distances{idx} / 1000; % Convert to km for plotting
    plot(time_hours, d_eucl_km, ...
         'LineWidth', 1.5, 'DisplayName', sprintf('N=%d (%.2fs)', N_values(idx), computation_times(idx)));
end
xlabel('Time (hours)');
ylabel('Euclidean Distance Error (km)');
title('Euclidean Error Evolution with Computation Times');
legend('Location', 'best');
grid on;

%% NEW: Bhattacharyya Distance at Specific Time Intervals
% Define time points of interest (1, 2, 3, ... hours)
target_hours = 1:1:floor(max(time_hours));
time_indices = zeros(size(target_hours));

% Find closest time indices for each target hour
for i = 1:length(target_hours)
    [~, time_indices(i)] = min(abs(time_hours - target_hours(i)));
end

% Extract distances at specific time points
bhatt_distances_at_times = zeros(length(N_values), length(target_hours));
eucl_distances_at_times = zeros(length(N_values), length(target_hours)); % NEW
for idx = 1:length(N_values)
    d_bhatt = bhattacharyya_distances{idx};
    d_eucl = euclidean_distances{idx};
    bhatt_distances_at_times(idx, :) = d_bhatt(time_indices);
    eucl_distances_at_times(idx, :) = d_eucl(time_indices); % NEW
end

% Create new figure for time-specific analysis
figure('Position', [300, 300, 1400, 900]);

% Plot 1: Bar chart showing Bhattacharyya distances at different hours
subplot(2, 2, 1);
bar_width = 0.8;
x_positions = 1:length(target_hours);
bar_handle = bar(x_positions, bhatt_distances_at_times', bar_width, 'grouped');

% Color the bars
for idx = 1:length(N_values)
    bar_handle(idx).DisplayName = sprintf('N=%d', N_values(idx));
end

xlabel('Time (hours)');
ylabel('Bhattacharyya Distance');
title('Bhattacharyya Distance at Hourly Intervals');
legend('Location', 'best');
grid on;
set(gca, 'XTick', x_positions, 'XTickLabel', target_hours);

% Plot 2: NEW - Bar chart showing Euclidean distances at different hours
subplot(2, 2, 2);
bar_handle2 = bar(x_positions, eucl_distances_at_times', bar_width, 'grouped');

% Color the bars
for idx = 1:length(N_values)
    bar_handle2(idx).DisplayName = sprintf('N=%d', N_values(idx));
end

xlabel('Time (hours)');
ylabel('Euclidean Distance Error (m)');
title('Euclidean Distance Error at Hourly Intervals');
legend('Location', 'best');
grid on;
set(gca, 'XTick', x_positions, 'XTickLabel', target_hours);

% Plot 3: NEW - Correlation between Bhattacharyya and Euclidean distances
subplot(2, 2, 3);
colors = lines(length(N_values));
for idx = 1:length(N_values)
    scatter(bhattacharyya_distances{idx}, euclidean_distances{idx}, 50, colors(idx,:), 'filled', ...
           'DisplayName', sprintf('N=%d', N_values(idx)));
    hold on;
end
xlabel('Bhattacharyya Distance');
ylabel('Euclidean Distance Error (m)');
title('Correlation: Bhattacharyya vs Euclidean Distance');
legend('Location', 'best');
grid on;

% Plot 4: NEW - Error reduction over time
subplot(2, 2, 4);
if length(N_values) > 1
    % Show improvement of higher N compared to N=1
    base_eucl = euclidean_distances{1}; % N=1 baseline
    for idx = 2:length(N_values)
        current_eucl = euclidean_distances{idx};
        error_reduction = (base_eucl - current_eucl) ./ base_eucl * 100; % Percentage improvement
        plot(time_hours, error_reduction, 'LineWidth', 2, ...
             'DisplayName', sprintf('N=%d vs N=1', N_values(idx)));
        hold on;
    end
    xlabel('Time (hours)');
    ylabel('Error Reduction (%)');
    title('Euclidean Error Reduction Compared to N=1');
    legend('Location', 'best');
    grid on;
end

%% Display comprehensive statistics
fprintf('\n=== Comprehensive Distance Statistics ===\n');
fprintf('%-8s %-12s %-12s %-12s %-12s %-12s %-12s %-15s\n', ...
        'N', 'Bhatt Mean', 'Bhatt Max', 'Eucl Mean', 'Eucl Max', 'Eucl RMS', 'Eucl Std', 'Comp. Time (s)');
fprintf('%-8s %-12s %-12s %-12s %-12s %-12s %-12s %-15s\n', ...
        '---', '----------', '----------', '---------', '---------', '---------', '--------', '--------------');

for idx = 1:length(N_values)
    d_bhatt = bhattacharyya_distances{idx};
    d_eucl = euclidean_distances{idx};
    fprintf('%-8d %-12.6f %-12.6f %-12.3f %-12.3f %-12.3f %-12.3f %-15.3f\n', ...
        N_values(idx), mean(d_bhatt), max(d_bhatt), ...
        mean(d_eucl), max(d_eucl), sqrt(mean(d_eucl.^2)), std(d_eucl), ...
        computation_times(idx));
end

% Compute improvement percentages
fprintf('\n=== Improvement Analysis ===\n');
base_bhatt_mean = mean(bhattacharyya_distances{1}); % N=1 as baseline
base_eucl_mean = mean(euclidean_distances{1});
base_time = computation_times(1);

for idx = 2:length(N_values)
    current_bhatt_mean = mean(bhattacharyya_distances{idx});
    current_eucl_mean = mean(euclidean_distances{idx});
    current_time = computation_times(idx);
    
    bhatt_improvement = ((base_bhatt_mean - current_bhatt_mean) / base_bhatt_mean) * 100;
    eucl_improvement = ((base_eucl_mean - current_eucl_mean) / base_eucl_mean) * 100;
    time_increase = ((current_time - base_time) / base_time) * 100;
    
    bhatt_efficiency = bhatt_improvement / time_increase;
    eucl_efficiency = eucl_improvement / time_increase;
    
    fprintf('N=%d vs N=1:\n', N_values(idx));
    fprintf('  Bhattacharyya improvement: %.2f%%\n', bhatt_improvement);
    fprintf('  Euclidean improvement: %.2f%%\n', eucl_improvement);
    fprintf('  Time increase: %.2f%%\n', time_increase);
    fprintf('  Bhattacharyya efficiency ratio: %.4f\n', bhatt_efficiency);
    fprintf('  Euclidean efficiency ratio: %.4f\n', eucl_efficiency);
    fprintf('\n');
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
    epsilon = 1e-6;
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

function [ISS_mixture_cov, ISS_mixture_mean, ISS_gaussian_means, ISS_gaussian_covs, ISS_mu_components, ISS_P_components, gaussian_weights] = GenerateGSF(mu, P, N, numSamples, startTime, stopTime)
sampleTime = 30; % seconds
timeVector = startTime:seconds(sampleTime):stopTime;
numTimeSteps = length(timeVector);

alpha = 1e-2;
beta = 2;
kappa = -3;
n = 6; 
n_sigma = 2*n + 1; 
numSamples = 100000;
data = mvnrnd(mu, P, numSamples);
options = statset('MaxIter', 10000, 'Display', 'off');
fprintf('Number of samples: %d\n', size(data', 1));
fprintf('Number of dimensions: %d\n', size(data', 2));
fprintf('Requested GMM components: %d\n', N);
gm = fitgmdist(data, N, 'Options', options, 'RegularizationValue', 1e-30);

N_components = gm.NumComponents; 


ISS_mu_components = gm.mu';            
ISS_P_components = gm.Sigma;       

Weights = gm.ComponentProportion;
%Weights = ones(1, N_components) / N;
ISS_all_sigma_points = zeros(n, n_sigma, N_components);
ISS_all_Wm = zeros(n_sigma, N_components);
ISS_all_Wc = zeros(n_sigma, N_components);
%disp(P)
for i = 1:N_components
    mu_i = ISS_mu_components(:, i);   
    P_i = ISS_P_components(:, :, i);  
    
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
%disp(ISS_all_sigma_points)
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


ISS_predicted_cov_keplerian = zeros(6, 6, N_components);

for j = 1:N_components
        
        Wm_j = ISS_all_Wm(:,j);
        Wc_j = ISS_all_Wc(:,j);
        
 
        mean_jt = ISS_all_sigma_points(:, :, j) * Wm_j; 
        %disp(mean_jt)

        cov_jt = zeros(6, 6);
        for i = 1:n_sigma
            diff = ISS_all_sigma_points(:,i, j) - mean_jt;
            cov_jt = cov_jt + Wc_j(i) * (diff * diff');
        end
        
        ISS_predicted_cov_keplerian(:,:,j) = cov_jt;
end
disp(mean_jt')
disp(ISS_predicted_cov_keplerian)




ISS_gaussian_means = ISS_predicted_mean_cartesian; 
ISS_gaussian_covs = ISS_predicted_cov_cartesian;  
gaussian_weights = Weights;               

ISS_mixture_mean = zeros(6, numTimeSteps);
ISS_mixture_cov = zeros(6, 6, numTimeSteps);

for t = 1:numTimeSteps
    mean_mix_t = zeros(6,1);
    for j = 1:N_components
        mean_j = ISS_gaussian_means(:,t,j);
        w_j = gaussian_weights(j);
        mean_mix_t = mean_mix_t + w_j * mean_j;
    end
    ISS_mixture_mean(:,t) = mean_mix_t;
    
    cov_mix_t = zeros(6, 6);
    for j = 1:N_components
        mean_j = ISS_gaussian_means(:,t,j); 
        cov_j = ISS_gaussian_covs(:,:,t,j);
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
    perc.a = 0.00005;   % 5% of a
    perc.h = 0.005;   % 1% of h
    perc.k = 0.005;   % 1% of k
    perc.p = 0.005;   % 1% of p
    perc.q = 0.005;   % 1% of q
    perc.l = 0.00001;  % 0.5% of l
    %}
    perc.a = 0.0001;   % 5% of a
    perc.h = 0.005;   % 1% of h
    perc.k = 0.005;   % 1% of k
    perc.p = 0.005;   % 1% of p
    perc.q = 0.005;   % 1% of q
    perc.l = 0.00001;  % 0.5% of l
    %{
    perc.a = 0.00001;   % 5% of a
    perc.h = 0.001;   % 1% of h
    perc.k = 0.001;   % 1% of k
    perc.p = 0.001;   % 1% of p
    perc.q = 0.001;   % 1% of q
    perc.l = 0.00001;  % 0.5% of l
    %}
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