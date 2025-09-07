global time_hours timeVector
startTime = datetime(2025, 9, 7, 10, 0, 0);
stopTime = startTime + days(1);
sampleTime = 10; % seconds
timeVector = startTime:seconds(sampleTime):stopTime;
numTimeSteps = length(timeVector);
numDebris = 1; % Number of debris objects
time_hours = hours(timeVector - startTime);
N = 1;
high_prob_threshold = 1e-4;
n_sigma = 13;
[collisionprobs_total, lognoprob, noprob] = propagate(weighted_mean, debris_weighted_mean, weighted_covariance, debris_weighted_covariance, numTimeSteps, numDebris);





fprintf('Calculating total collision probability across all debris...\n');
for j = 1:numTimeSteps
    lognoprob(j) = 0; 
    for i = 1:numDebris
      
        lognoprob(j) = lognoprob(j) + log1p(-collisionprobs_total(1, j, i));
    end
    noprob(j) = exp(lognoprob(j));
end
total_collision_prob_per_timestep = 1 - noprob; 
max_total_prob = max(total_collision_prob_per_timestep);


fprintf('\n=== TOTAL COLLISION PROBABILITY ANALYSIS ===\n');
fprintf('Maximum total collision probability: %.6e\n', max_total_prob);


[~, max_idx] = max(total_collision_prob_per_timestep);
fprintf('Maximum occurs at time step: %d\n', max_idx);
fprintf('Time: %s\n', char(timeVector(max_idx)));


high_risk_total = sum(total_collision_prob_per_timestep > 1e-4);
if high_risk_total > 0
    fprintf('HIGH RISK: %d time steps with total Pc > 10^-4\n', high_risk_total);
end

%{
figure('Position', [100, 100, 1200, 600]);
plot(time_hours, total_collision_prob_per_timestep, 'b-', 'LineWidth', 2);
hold on;
yline(1e-4, 'm--', 'LineWidth', 1.5, 'Label', 'Pc = 10^{-4}');
xlabel('Time Step');
ylabel('Total Collision Probability');
title('Total Collision Probability over Time - Linear Scale');
grid on;
hold off;
%}
find_minimum_constant_binary_search(total_collision_prob_per_timestep, weighted_mean, debris_weighted_mean, weighted_covariance, debris_weighted_covariance)


function optimal_constant = find_minimum_constant_binary_search(total_collision_prob_per_timestep, weighted_mean, debris_weighted_mean, weighted_covariance, debris_weighted_covariance)
    % Binary search to find minimum constant that achieves Pc < 1e-4
    
    % Define search bounds
    constant_min = 0;        % Lower bound
    constant_max = 0.1;      % Upper bound (adjust based on your system)
    tolerance = 1e-7;        % Search precision
    target_prob = 1e-4;      % Target maximum probability
    
    % Verify that constant_max achieves the target
    max_prob_upper = evaluate_collision_probability(constant_max, total_collision_prob_per_timestep, weighted_mean, debris_weighted_mean, weighted_covariance, debris_weighted_covariance);
    if max_prob_upper >= target_prob
        error('Upper bound constant_max = %.6f does not achieve target probability', constant_max);
    end
    
    % Binary search
    iteration = 0;
    max_iterations = 50;
    
    fprintf('Starting binary search...\n');
    fprintf('Target: max_total_prob < %.2e\n', target_prob);
    fprintf('Initial bounds: [%.6f, %.6f]\n', constant_min, constant_max);
    
    while (constant_max - constant_min) > tolerance && iteration < max_iterations
        iteration = iteration + 1;
        constant_mid = (constant_min + constant_max) / 2;
        
        % Evaluate collision probability at midpoint
        max_prob_mid = evaluate_collision_probability(constant_mid, total_collision_prob_per_timestep, weighted_mean, debris_weighted_mean, weighted_covariance, debris_weighted_covariance);
        format long
        fprintf('Iteration %d: constant = %.8f, max_prob = %.8e\n', ...
                iteration, constant_mid, max_prob_mid);
        
        if max_prob_mid < target_prob
            % Mid value achieves target, try smaller constant
            constant_max = constant_mid;
        else
            % Mid value doesn't achieve target, need larger constant
            constant_min = constant_mid;
        end
    end
    
    optimal_constant = constant_max;
    final_prob = evaluate_collision_probability(optimal_constant, total_collision_prob_per_timestep, weighted_mean, debris_weighted_mean, weighted_covariance, debris_weighted_covariance);
    
    fprintf('\nOptimization complete!\n');
    fprintf('Optimal constant: %.6f\n', optimal_constant);
    fprintf('Achieved max_total_prob: %.6e\n', final_prob);
end


function max_total_prob = evaluate_collision_probability(constant, total_collision_prob_per_timestep, weighted_mean, debris_weighted_mean, weighted_covariance, debris_weighted_covariance)
global timeVector
startTime =datetime(2025, 9, 7, 10, 0, 0);
stopTime = startTime + days(1);
sampleTime = 10; % seconds
timeVector = startTime:seconds(sampleTime):stopTime;
numTimeSteps = length(timeVector);
numDebris = 1; % Number of debris objects
N = 1;
n_sigma = 13;


high_prob_threshold = 1e-4;

first_exceedance_idx = find(total_collision_prob_per_timestep > high_prob_threshold, 1, 'first');

hour_index = 3600;

transfer_index = first_exceedance_idx - hour_index;

sc_transfer = satelliteScenario(timeVector(transfer_index), stopTime, sampleTime);

transfer_mean = weighted_mean(:, transfer_index);

[a, ecc, incl, RAAN, argp, nu, truelon, arglat, lonper] = ijk2keplerian(transfer_mean(1:3), transfer_mean(4:6));;

sat_transfer = satellite(sc_transfer, a, ecc, incl, RAAN, argp, nu, 'OrbitPropagator','two-body-keplerian');

[pre_positions, pre_velocities] = states(sat_transfer, timeVector(transfer_index));

delta_v = constant * pre_velocities;

mu_earth = 3.986e14;

delta_a = (2 * a^2 / mu_earth) * dot(pre_velocities, delta_v);

sat_post_transfer = satellite(sc_transfer, a + delta_a, ecc, incl, RAAN, argp, nu, 'OrbitPropagator','two-body-keplerian');

[post_positions, post_velocities] = states(sat_post_transfer);

new_weighted_mean = [weighted_mean(:, 1:transfer_index - 1), [post_positions;post_velocities]];
new_weighted_covariance = weighted_covariance;
%{
ISS_propagated_points = cell(1, n_sigma);
for i = 1:n_sigma
    sat_sigma = satellite(sc_transfer, Satellite_sigma_point_files{i});

    elements1 = orbitalElements(sat_sigma);

    [nominal, ~, ~] = convertOMMToSGP4Elements(elements1);


    mu_sigma = [nominal.MeanMotion;
        nominal.eccentricity;
        nominal.inclination;
        nominal.raan;
        nominal.argp;
        nominal.meanAnomaly];

    [pre_positions, pre_velocities] = states(sat_sigma, timeVector(transfer_index));

    delta_v = constant * pre_velocities;
    
    mu_earth = 3.986e14;
    
    a = (mu_earth / nominal.MeanMotion^2)^(1/3); 
    delta_mM = - (3 * a * mu_sigma(1) / mu_earth) * dot(pre_velocities, delta_v);


    elems = struct( ...
        'meanMotion', (mu_sigma(1)+delta_mM) * 240, ...   
        'eccentricity', mu_sigma(2), ...
        'inclination', mu_sigma(3), ...
        'raan', mu_sigma(4), ...
        'argPeriapsis', mu_sigma(5), ...
        'meanAnomaly', mu_sigma(6) ...
    );
    updateOMMFile(Satellite_sigma_point_files{i}, 'Transfer_updated.xml', elems);
    sat = satellite(sc_transfer, 'Transfer_updated.xml');
    ISS_propagated_points{i} = sat;
end
% Collect propagated positions and velocities
ISS_positions = cell(1, n_sigma);
ISS_velocities = cell(1, n_sigma);
for i = 1:n_sigma
    sat = ISS_propagated_points{i};
    [pos, vel, ~] = states(sat);
    ISS_positions{i} = pos;
    ISS_velocities{i} = vel;
end

remaining_times = timeVector(transfer_index:end);
num_remaining_steps = length(remaining_times);
% Compute mean and covariance in Cartesian coordinates
ISS_gaussian_means = zeros(6, num_remaining_steps);
ISS_gaussian_covs = zeros(6, 6, num_remaining_steps);

for t = 1:num_remaining_steps
    sigma_points_t = zeros(6, n_sigma);
    for i = 1:n_sigma
        sigma_points_t(:, i) = [ISS_positions{i}(:, t); ISS_velocities{i}(:, t)];
    end
    
    % Mean
    mean_t = sigma_points_t * Wm';
    ISS_gaussian_means(:, t) = mean_t;
    
    % Covariance
    cov_t = zeros(6, 6);
    for i = 1:n_sigma
        diff = sigma_points_t(:, i) - mean_t;
        cov_t = cov_t + Wc(i) * (diff * diff');
    end
    ISS_gaussian_covs(:, :, t) = cov_t;
end

new_weighted_mean = [weighted_mean(:, 1:transfer_index - 1), ISS_gaussian_means];
new_weighted_covariance = cat(3, weighted_covariance(:, :, 1:transfer_index - 1), ISS_gaussian_covs);
%}


ISS_positions = new_weighted_mean(1:3, :);

ISS_position_covs = new_weighted_covariance(1:3, 1:3, :);
ISS_positions_km = ISS_positions / 1000;
time_hours = hours(timeVector - startTime);

%{
% Create figure with offset to avoid overlapping
figure('Position', [100, 100, 1200, 800]);

% Position magnitude plot
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
title('Position Magnitude vs Time');
legend('Location', 'best');
grid on;
%}
disp(size(new_weighted_mean))
disp(size(debris_weighted_mean))

[collisionprobs_total, lognoprob, noprob] = propagate(new_weighted_mean, debris_weighted_mean, new_weighted_covariance, debris_weighted_covariance, numTimeSteps, numDebris);

fprintf('Calculating total collision probability across all debris...\n');
for j = 1:numTimeSteps
    lognoprob(j) = 0; 
    for i = 1:numDebris
      
        lognoprob(j) = lognoprob(j) + log1p(-collisionprobs_total(1, j, i));
    end
    noprob(j) = exp(lognoprob(j));
end
total_collision_prob_per_timestep = 1 - noprob; 
max_total_prob = max(total_collision_prob_per_timestep);


figure('Position', [100, 100, 1200, 600]);
plot(timeVector, total_collision_prob_per_timestep, 'b-', 'LineWidth', 2);
hold on;

% Add horizontal reference line
yline(1e-4, 'm--', 'LineWidth', 1.5, 'Label', 'Pc = 10^{-4}');

% Mark the maximum point
[max_total_prob, max_idx] = max(total_collision_prob_per_timestep);
plot(timeVector(max_idx), max_total_prob, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red');

% Add text annotation attached to the maximum point
text(timeVector(max_idx), max_total_prob, ...
     sprintf('  Max: %.2e', max_total_prob), ...
     'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
     'FontSize', 12, 'FontWeight', 'bold', 'Color', 'red');


xlabel('Time (hours)');
ylabel('Total Collision Probability');
title('Total Collision Probability over Time - Linear Scale');
grid on;
legend('Total Pc', 'Pc = 10^{-4}', 'Maximum Point', 'Location', 'best');
hold off;

end

function [collisionprobs_total, lognoprob, noprob] = propagate(weighted_mean, debris_weighted_mean, weighted_covariance, debris_weighted_covariance, numTimeSteps, numDebris)
global time_hours timeVector
% Initialize arrays to store collision probabilities for all debris
collisionprobs_total = zeros(1, numTimeSteps, numDebris);
log_collisionprobs_total = zeros(1, numTimeSteps, numDebris);
ratios_total = zeros(1, numTimeSteps, numDebris);

% Initialize arrays to store sigma ratios for all debris
sigma_x_ratios_total = zeros(1, numTimeSteps, numDebris);
sigma_z_ratios_total = zeros(1, numTimeSteps, numDebris);

% Initialize total collision probability arrays
lognoprob = zeros(1, numTimeSteps);
noprob = zeros(1, numTimeSteps);


for i = 1:numDebris
    % Use weighted means instead of individual Gaussian components
    ISS_predicted_mean_cartesian_positions = weighted_mean(1:3, :);     % [3 x numTimeSteps]
    ISS_predicted_mean_cartesian_velocities = weighted_mean(4:6, :);    % [3 x numTimeSteps]
    debris_predicted_mean_cartesian_positions = debris_weighted_mean(1:3, :, i);  % [3 x numTimeSteps]
    debris_predicted_mean_cartesian_velocities = debris_weighted_mean(4:6, :, i); % [3 x numTimeSteps]

    % Use weighted covariances instead of individual Gaussian components
    ISS_predicted_cov_cartesian_positions = weighted_covariance(1:3, 1:3, :);     % [3 x 3 x numTimeSteps]
    ISS_predicted_cov_cartesian_velocities = weighted_covariance(4:6, 4:6, :);    % [3 x 3 x numTimeSteps]
    debris_predicted_cov_cartesian_positions = debris_weighted_covariance(1:3, 1:3, :, i);  % [3 x 3 x numTimeSteps]
    debris_predicted_cov_cartesian_velocities = debris_weighted_covariance(4:6, 4:6, :, i); % [3 x 3 x numTimeSteps]

    % Compute relative positions and velocities
    relative_positions = ISS_predicted_mean_cartesian_positions - debris_predicted_mean_cartesian_positions;
    relative_velocity = ISS_predicted_mean_cartesian_velocities - debris_predicted_mean_cartesian_velocities;

    distance = sqrt(sum(relative_positions.^2, 1)); % size: [1, numTimeSteps]
    %{
    figure('Position', [100, 100, 1200, 800]);
    plot(1:numTimeSteps, distance);

    xlabel('Time Step');
    ylabel('Distance (m)');
    title(sprintf('Distance between ISS and Debris %d', i));
    grid on;

    fprintf('ISS mean position at t=1: [%.2f, %.2f, %.2f]\n', weighted_mean(1:3, 1));
    fprintf('Debris mean position at t=1: [%.2f, %.2f, %.2f]\n', debris_weighted_mean(1:3, 1, i));
    distance_test = norm(weighted_mean(1:3, 1) - debris_weighted_mean(1:3, 1, i));
    fprintf('Distance between mixture means: %.2f km\n', distance_test/1000);
    %}
    R = 5;
    % Compute collision probabilities using weighted means and covariances
    [log_collisionprobs, collisionprobs, ratios, sigma_x_ratios, sigma_z_ratios] = compute_collision_probability_over_time_weighted(...
        ISS_predicted_mean_cartesian_positions, ISS_predicted_cov_cartesian_positions, ISS_predicted_mean_cartesian_velocities, ...
        debris_predicted_mean_cartesian_positions, debris_predicted_cov_cartesian_positions, debris_predicted_mean_cartesian_velocities, R);

    % Store results for total probability calculation
    collisionprobs_total(1, :, i) = collisionprobs;
    log_collisionprobs_total(1, :, i) = log_collisionprobs;
    ratios_total(1, :, i) = ratios;
    sigma_x_ratios_total(1, :, i) = sigma_x_ratios;
    sigma_z_ratios_total(1, :, i) = sigma_z_ratios;

    % Define validity threshold (8.5 sigma)
    validity_threshold = 8.5;
    
    % Create validity mask - only show where approximation is valid
    valid_indices = ratios <= validity_threshold;
    
    % Create arrays for plotting (NaN where invalid)
    valid_log_probs = log_collisionprobs;
    valid_log_probs(~valid_indices) = NaN;
    
    valid_collision_probs = collisionprobs;
    valid_collision_probs(~valid_indices) = NaN;
    
    valid_ratios = ratios;
    valid_ratios(~valid_indices) = NaN;
    
    figure('Position', [100, 100, 1200, 800]);
    % Plot log probabilities (only valid regions)
    subplot(2,1,1);
    semilogy(timeVector, valid_collision_probs, 'b-', 'LineWidth', 1.5);
    hold on;
    semilogy(timeVector, collisionprobs, 'r', 'LineWidth', 0.5, 'Color', [0.8 0.8 0.8]);
    
    % Add markers for high collision probabilities (>10^-4) in valid regions
    high_prob_threshold = 1e-4;
    high_prob_indices = (valid_collision_probs > high_prob_threshold) & ~isnan(valid_collision_probs);
    
    if any(high_prob_indices)
semilogy(timeVector(high_prob_indices), valid_collision_probs(high_prob_indices), ...
         'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
    end
    
    % Add horizontal reference line at 10^-4
    yline(high_prob_threshold, 'm--', 'LineWidth', 1.5, 'Label', 'Pc = 10^{-4}');
    
    xlabel('Time (hours)');
    ylabel('Collision Probability');
    title('Collision Probability over Time (Valid regions only)');
    if any(high_prob_indices)
        legend('Valid (≤8.5σ)', 'All data', 'High Risk (>10^{-4})', 'Pc = 10^{-4}', 'Location', 'best');
    else
        legend('Valid (≤8.5σ)', 'All data', 'Pc = 10^{-5}', 'Location', 'best');
    end
    grid on;
    hold off;

    subplot(2,1,2);
    plot(timeVector, ratios, 'k-', 'LineWidth', 1);
    hold on;
    yline(validity_threshold, 'r--', 'LineWidth', 2, 'Label', '8.5σ threshold');
    % Highlight valid regions
    valid_y = ratios;
    valid_y(~valid_indices) = NaN;
    plot(timeVector, valid_y, 'g-', 'LineWidth', 2);
    xlabel('Time Step');
    ylabel('Miss Distance / Uncertainty Radius (σ)');
    title('Approximation Validity Check');
    legend('All ratios', '8.5σ threshold', 'Valid regions', 'Location', 'best');
    grid on;
    hold off;
    
end
end

function [log_collision_probs, collision_probs, ratios, sigma_x_ratios, sigma_z_ratios] = compute_collision_probability_over_time_weighted(...
    means1_pos, covs1_pos, vels1, ...
    means2_pos, covs2_pos, vels2, R)

T = size(means1_pos, 2); % Number of time steps
log_collision_probs = zeros(1, T);
collision_probs = zeros(1, T);
ratios = zeros(1, T);
sigma_x_ratios = zeros(1, T);
sigma_z_ratios = zeros(1, T);

for i = 1:T
    % Extract states and covariances at time step i
    mu1 = means1_pos(:, i); % [3x1] ISS position
    mu2 = means2_pos(:, i); % [3x1] Debris position
    P1 = covs1_pos(:, :, i); % [3x3] ISS position covariance
    P2 = covs2_pos(:, :, i); % [3x3] Debris position covariance
    v1 = vels1(:, i); % [3x1] ISS velocity
    v2 = vels2(:, i); % [3x1] Debris velocity
   


    % Relative quantities
    rel_mu = mu1 - mu2;
    rel_vel = v1 - v2;
    
    if norm(rel_vel) < 1e-8
        % Objects have nearly identical velocities - high collision risk
        log_collision_probs(i) = 0; % log(1) = 0, meaning Pc = 1
        collision_probs(i) = 1;
        ratios(i) = 0;
        sigma_x_ratios(i) = NaN;
        sigma_z_ratios(i) = NaN;
        continue;
    end
    
    % Construct encounter frame
    ye = rel_vel / norm(rel_vel);
    x_temp = rel_mu - dot(rel_mu, ye) * ye;
    
    if norm(x_temp) < 1e-12
        % Objects are on collision course
        log_collision_probs(i) = 0; % log(1) = 0, meaning Pc = 1
        collision_probs(i) = 1;
        ratios(i) = 0;
        sigma_x_ratios(i) = NaN;
        sigma_z_ratios(i) = NaN;
        continue;
    end
    
    xe = x_temp / norm(x_temp);
    ze = cross(xe, ye);
    Me = [xe, ye, ze]; % 3x3 encounter frame transformation matrix
    
    % Transform covariances to encounter frame
    Pe1 = Me' * P1 * Me;
    Pe2 = Me' * P2 * Me;
    Pe = Pe1 + Pe2;
    
    % Drop along-track (ye) component → 2D in encounter plane
    Pe_2D = Pe([1 3], [1 3]); % xe-ze plane
    mu_rel_2D = Me' * rel_mu;
    mu_2D = mu_rel_2D([1 3]);
    
    % Check if covariance matrix is positive definite
    if det(Pe_2D) <= 0 || any(eig(Pe_2D) <= 0)
        log_collision_probs(i) = -Inf;
        collision_probs(i) = 0;
        ratios(i) = Inf;
        sigma_x_ratios(i) = NaN;
        sigma_z_ratios(i) = NaN;
        continue;
    end
    
    % Rotate covariance to principal axes
    [V, D] = eig(Pe_2D);
    P_diag = D;
    mu_rot = V' * mu_2D;
    sigma_x = sqrt(abs(P_diag(1,1)));
    sigma_z = sqrt(abs(P_diag(2,2)));
    dx = mu_rot(1);
    dz = mu_rot(2);
    
    % Check for numerical issues
    if sigma_x < 1e-12 || sigma_z < 1e-12
        log_collision_probs(i) = -Inf;
        collision_probs(i) = 0;
        ratios(i) = Inf;
        sigma_x_ratios(i) = NaN;
        sigma_z_ratios(i) = NaN;
        continue;
    end
    
    % Compute collision probability using double integral
    % Integrate bivariate normal PDF over circular region of radius R
    
    % Define the bivariate normal PDF in rotated coordinates
    pdf_func = @(x, z) (1 / (2 * pi * sigma_x * sigma_z)) * ...
               exp(-0.5 * ((x - dx).^2 / sigma_x^2 + (z - dz).^2 / sigma_z^2));
    
    % Integration limits for circular region
    % For each x, z ranges from -sqrt(R^2 - x^2) to +sqrt(R^2 - x^2)
    z_lower = @(x) -sqrt(max(0, R^2 - x.^2));
    z_upper = @(x) sqrt(max(0, R^2 - x.^2));
    
    try
        % Compute the double integral
        Pc = integral2(pdf_func, -R, R, z_lower, z_upper, ...
                      'AbsTol', 1e-12, 'RelTol', 1e-9);
        
        % Handle numerical precision issues
        Pc = max(0, min(1, Pc)); % Clamp to [0, 1]
        
        if Pc > 0
            log_Pc = log(Pc);
        else
            log_Pc = -Inf;
        end
        
    catch ME
        % If integration fails, fall back to very small probability
        warning('Integration failed at time step %d: %s', i, ME.message);
        Pc = 0;
        log_Pc = -Inf;
    end
    
    % Compute ratio to check approximation validity
    miss_distance = norm(mu_2D);
    uncertainty_radius = sqrt(sigma_x^2 + sigma_z^2);
    ratio = miss_distance / uncertainty_radius;
    
    % Compute sigma/R ratios
    sigma_x_ratios(i) = sigma_x / R;
    sigma_z_ratios(i) = sigma_z / R;
    
    % Store results
    log_collision_probs(i) = log_Pc;
    collision_probs(i) = Pc;
    ratios(i) = ratio;
end
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


function updateOMMFile(inputFile, outputFile, elements)
% updateOMMFile  Update Keplerian orbital elements in an OMM XML file
%
%   updateOMMFile(inputFile, outputFile, elements)
%
%   INPUTS:
%     inputFile   - path to original OMM XML file
%     outputFile  - path to save updated file
%     elements    - struct with fields:
%                      meanMotion [rev/day] (primary)
%                      eccentricity
%                      inclination [deg]
%                      raan [deg]
%                      argPeriapsis [deg]
%                      meanAnomaly [deg]

    if ~exist(inputFile, 'file')
        error('Input file not found: %s', inputFile);
    end

    mu_earth = 3.986004418e14; % m^3/s^2
    
    % If mean motion is provided, calculate semi-major axis and period from it
    if isfield(elements, 'meanMotion') && ~isempty(elements.meanMotion)
        n_rad_s = elements.meanMotion * 2*pi / 86400; % rev/day -> rad/s
        a_m = (mu_earth / n_rad_s^2)^(1/3);          % semi-major axis in meters
        elements.semiMajorAxis = a_m / 1000;         % km
        elements.period = 2*pi / n_rad_s / 60;       % minutes

    else
        error('meanMotion must be provided in elements struct.');
    end
    
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


function [nominal, uncertainty, P] = convertOMMToSGP4Elements(elements1)
    nominal = struct();
    uncertainty = struct();
   

    nominal.MeanMotion = elements1.MeanMotion; % Convert km to m
    nominal.eccentricity = elements1.Eccentricity;
    nominal.inclination = elements1.Inclination; % degrees
    nominal.raan = elements1.RightAscensionOfAscendingNode; % degrees
    nominal.argp = elements1.ArgumentOfPeriapsis; % degrees
    nominal.meanAnomaly = elements1.MeanAnomaly;

    %{
    perc.semiMajorAxis = 0.00001;   % 5% of a
    perc.eccentricity = 0.0001;   % 1% of h
    perc.inclination = 0.001;   % 1% of k
    perc.raan = 0.00001;   % 1% of p
    perc.argp = 0.001;   % 1% of q
    perc.trueAnomaly = 0.0001;  % 0.5% of l
    %}
    perc.semiMajorAxis = 0.000001;   % 5% of a
    perc.eccentricity = 0.00001;   % 1% of h
    perc.inclination = 0.00001;   % 1% of k
    perc.raan = 0.00001;   % 1% of p
    perc.argp = 0.00001;   % 1% of q
    perc.trueAnomaly = 0.00001;  % 0.5% of l

    uncertainty.MeanMotion = perc.semiMajorAxis * nominal.MeanMotion;
    uncertainty.eccentricity = perc.eccentricity * nominal.eccentricity;
    uncertainty.inclination = perc.inclination * nominal.inclination;
    uncertainty.raan = perc.raan * nominal.raan;
    uncertainty.argPeriapsis = perc.argp * nominal.argp;
    uncertainty.meanAnomaly = perc.trueAnomaly * nominal.meanAnomaly;

    P = diag([uncertainty.MeanMotion^2, ...
        uncertainty.eccentricity^2, ...
        uncertainty.inclination^2, ...
        uncertainty.raan^2, ...
        uncertainty.argPeriapsis^2, ...
        uncertainty.meanAnomaly^2]);
end