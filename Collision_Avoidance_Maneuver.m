global time_hours timeVector
startTime = datetime(2025, 9, 7, 10, 0, 0);
stopTime = startTime + days(1);
sampleTime = 10;
timeVector = startTime:seconds(sampleTime):stopTime;
numTimeSteps = length(timeVector);
numDebris = 1;
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

find_minimum_constant_binary_search(total_collision_prob_per_timestep, weighted_mean, debris_weighted_mean, weighted_covariance, debris_weighted_covariance)

function optimal_constant = find_minimum_constant_binary_search(total_collision_prob_per_timestep, weighted_mean, debris_weighted_mean, weighted_covariance, debris_weighted_covariance)
    constant_min = 0;
    constant_max = 0.1;
    tolerance = 1e-7;
    target_prob = 1e-4;
    
    max_prob_upper = evaluate_collision_probability(constant_max, total_collision_prob_per_timestep, weighted_mean, debris_weighted_mean, weighted_covariance, debris_weighted_covariance);
    if max_prob_upper >= target_prob
        error('Upper bound constant_max = %.6f does not achieve target probability', constant_max);
    end
    
    iteration = 0;
    max_iterations = 50;
    
    fprintf('Starting binary search...\n');
    fprintf('Target: max_total_prob < %.2e\n', target_prob);
    fprintf('Initial bounds: [%.6f, %.6f]\n', constant_min, constant_max);
    
    while (constant_max - constant_min) > tolerance && iteration < max_iterations
        iteration = iteration + 1;
        constant_mid = (constant_min + constant_max) / 2;
        max_prob_mid = evaluate_collision_probability(constant_mid, total_collision_prob_per_timestep, weighted_mean, debris_weighted_mean, weighted_covariance, debris_weighted_covariance);
        format long
        fprintf('Iteration %d: constant = %.8f, max_prob = %.8e\n', iteration, constant_mid, max_prob_mid);
        
        if max_prob_mid < target_prob
            constant_max = constant_mid;
        else
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
sampleTime = 10;
timeVector = startTime:seconds(sampleTime):stopTime;
numTimeSteps = length(timeVector);
numDebris = 1;
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

ISS_positions = new_weighted_mean(1:3, :);
ISS_position_covs = new_weighted_covariance(1:3, 1:3, :);
ISS_positions_km = ISS_positions / 1000;
time_hours = hours(timeVector - startTime);

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
yline(1e-4, 'm--', 'LineWidth', 1.5, 'Label', 'Pc = 10^{-4}');
[max_total_prob, max_idx] = max(total_collision_prob_per_timestep);
plot(timeVector(max_idx), max_total_prob, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red');
text(timeVector(max_idx), max_total_prob, sprintf('  Max: %.2e', max_total_prob), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'red');
xlabel('Time (hours)');
ylabel('Total Collision Probability');
title('Total Collision Probability over Time - Linear Scale');
grid on;
legend('Total Pc', 'Pc = 10^{-4}', 'Maximum Point', 'Location', 'best');
hold off;
end

function [collisionprobs_total, lognoprob, noprob] = propagate(weighted_mean, debris_weighted_mean, weighted_covariance, debris_weighted_covariance, numTimeSteps, numDebris)
global time_hours timeVector
collisionprobs_total = zeros(1, numTimeSteps, numDebris);
log_collisionprobs_total = zeros(1, numTimeSteps, numDebris);
ratios_total = zeros(1, numTimeSteps, numDebris);
sigma_x_ratios_total = zeros(1, numTimeSteps, numDebris);
sigma_z_ratios_total = zeros(1, numTimeSteps, numDebris);
lognoprob = zeros(1, numTimeSteps);
noprob = zeros(1, numTimeSteps);

for i = 1:numDebris
    ISS_predicted_mean_cartesian_positions = weighted_mean(1:3, :);
    ISS_predicted_mean_cartesian_velocities = weighted_mean(4:6, :);
    debris_predicted_mean_cartesian_positions = debris_weighted_mean(1:3, :, i);
    debris_predicted_mean_cartesian_velocities = debris_weighted_mean(4:6, :, i);
    ISS_predicted_cov_cartesian_positions = weighted_covariance(1:3, 1:3, :);
    ISS_predicted_cov_cartesian_velocities = weighted_covariance(4:6, 4:6, :);
    debris_predicted_cov_cartesian_positions = debris_weighted_covariance(1:3, 1:3, :, i);
    debris_predicted_cov_cartesian_velocities = debris_weighted_covariance(4:6, 4:6, :, i);

    relative_positions = ISS_predicted_mean_cartesian_positions - debris_predicted_mean_cartesian_positions;
    relative_velocity = ISS_predicted_mean_cartesian_velocities - debris_predicted_mean_cartesian_velocities;
    distance = sqrt(sum(relative_positions.^2, 1));
    R = 5;

    [log_collisionprobs, collisionprobs, ratios, sigma_x_ratios, sigma_z_ratios] = compute_collision_probability_over_time_weighted(...
        ISS_predicted_mean_cartesian_positions, ISS_predicted_cov_cartesian_positions, ISS_predicted_mean_cartesian_velocities, ...
        debris_predicted_mean_cartesian_positions, debris_predicted_cov_cartesian_positions, debris_predicted_mean_cartesian_velocities, R);

    collisionprobs_total(1, :, i) = collisionprobs;
    log_collisionprobs_total(1, :, i) = log_collisionprobs;
    ratios_total(1, :, i) = ratios;
    sigma_x_ratios_total(1, :, i) = sigma_x_ratios;
    sigma_z_ratios_total(1, :, i) = sigma_z_ratios;

    validity_threshold = 8.5;
    valid_indices = ratios <= validity_threshold;
    valid_log_probs = log_collisionprobs;
    valid_log_probs(~valid_indices) = NaN;
    valid_collision_probs = collisionprobs;
    valid_collision_probs(~valid_indices) = NaN;
    valid_ratios = ratios;
    valid_ratios(~valid_indices) = NaN;

    figure('Position', [100, 100, 1200, 800]);
    subplot(2,1,1);
    semilogy(timeVector, valid_collision_probs, 'b-', 'LineWidth', 1.5);
    hold on;
    semilogy(timeVector, collisionprobs, 'r', 'LineWidth', 0.5, 'Color', [0.8 0.8 0.8]);
    high_prob_threshold = 1e-4;
    high_prob_indices = (valid_collision_probs > high_prob_threshold) & ~isnan(valid_collision_probs);
    if any(high_prob_indices)
        semilogy(timeVector(high_prob_indices), valid_collision_probs(high_prob_indices), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
    end
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
    means1_pos, covs1_pos, vels1, means2_pos, covs2_pos, vels2, R)

T = size(means1_pos, 2);
log_collision_probs = zeros(1, T);
collision_probs = zeros(1, T);
ratios = zeros(1, T);
sigma_x_ratios = zeros(1, T);
sigma_z_ratios = zeros(1, T);

for i = 1:T
    mu1 = means1_pos(:, i);
    mu2 = means2_pos(:, i);
    P1 = covs1_pos(:, :, i);
    P2 = covs2_pos(:, :, i);
    v1 = vels1(:, i);
    v2 = vels2(:, i);

    rel_mu = mu1 - mu2;
    rel_vel = v1 - v2;
    
    if norm(rel_vel) < 1e-8
        log_collision_probs(i) = 0;
        collision_probs(i) = 1;
        ratios(i) = 0;
        sigma_x_ratios(i) = NaN;
        sigma_z_ratios(i) = NaN;
        continue;
    end
    
    ye = rel_vel / norm(rel_vel);
    x_temp = rel_mu - dot(rel_mu, ye) * ye;
    
    if norm(x_temp) < 1e-12
        log_collision_probs(i) = 0;
        collision_probs(i) = 1;
        ratios(i) = 0;
        sigma_x_ratios(i) = NaN;
        sigma_z_ratios(i) = NaN;
        continue;
    end
    
    xe = x_temp / norm(x_temp);
    ze = cross(xe, ye);
    Me = [xe, ye, ze];
    Pe1 = Me' * P1 * Me;
    Pe2 = Me' * P2 * Me;
    Pe = Pe1 + Pe2;
    Pe_2D = Pe([1 3], [1 3]);
    mu_rel_2D = Me' * rel_mu;
    mu_2D = mu_rel_2D([1 3]);
    
    if det(Pe_2D) <= 0 || any(eig(Pe_2D) <= 0)
        log_collision_probs(i) = -Inf;
        collision_probs(i) = 0;
        ratios(i) = Inf;
        sigma_x_ratios(i) = NaN;
        sigma_z_ratios(i) = NaN;
        continue;
    end
    
    [V, D] = eig(Pe_2D);
    P_diag = D;
    mu_rot = V' * mu_2D;
    sigma_x = sqrt(abs(P_diag(1,1)));
    sigma_z = sqrt(abs(P_diag(2,2)));
    dx = mu_rot(1);
    dz = mu_rot(2);
    
    if sigma_x < 1e-12 || sigma_z < 1e-12
        log_collision_probs(i) = -Inf;
        collision_probs(i) = 0;
        ratios(i) = Inf;
        sigma_x_ratios(i) = NaN;
        sigma_z_ratios(i) = NaN;
        continue;
    end
    
    pdf_func = @(x, z) (1 / (2 * pi * sigma_x * sigma_z)) * exp(-0.5 * ((x - dx).^2 / sigma_x^2 + (z - dz).^2 / sigma_z^2));
    z_lower = @(x) -sqrt(max(0, R^2 - x.^2));
    z_upper = @(x) sqrt(max(0, R^2 - x.^2));
    
    try
        Pc = integral2(pdf_func, -R, R, z_lower, z_upper, 'AbsTol', 1e-12, 'RelTol', 1e-9);
        Pc = max(0, min(1, Pc));
        if Pc > 0
            log_Pc = log(Pc);
        else
            log_Pc = -Inf;
        end
    catch ME
        warning('Integration failed at time step %d: %s', i, ME.message);
        Pc = 0;
        log_Pc = -Inf;
    end
    
    miss_distance = norm(mu_2D);
    uncertainty_radius = sqrt(sigma_x^2 + sigma_z^2);
    ratio = miss_distance / uncertainty_radius;
    sigma_x_ratios(i) = sigma_x / R;
    sigma_z_ratios(i) = sigma_z / R;
    log_collision_probs(i) = log_Pc;
    collision_probs(i) = Pc;
    ratios(i) = ratio;
end
end
