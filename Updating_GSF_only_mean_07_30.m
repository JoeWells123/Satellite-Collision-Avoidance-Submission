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

N = 1;
numSamples = 1000;
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

fprintf('\n=== Processing Multiple Alpha Values ===\n');

% Define alpha values to test
alpha_values = [1, 1e-1, 1e-2, 1e-3];
num_alphas = length(alpha_values);

% Storage for results
results = struct();
for i = 1:num_alphas
    alpha_str = sprintf('alpha_%g', alpha_values(i));
    alpha_str = strrep(alpha_str, '.', 'p'); % Replace . with p for valid field names
    alpha_str = strrep(alpha_str, '-', 'neg'); % Replace - with neg
    
    results.(alpha_str) = struct();
end

% Process each alpha value
for alpha_idx = 1:num_alphas
    current_alpha = alpha_values(alpha_idx);
    alpha_str = sprintf('alpha_%g', current_alpha);
    alpha_str = strrep(alpha_str, '.', 'p');
    alpha_str = strrep(alpha_str, '-', 'neg');
    
    fprintf('\n--- Processing Alpha = %g ---\n', current_alpha);
    
    [weighted_covariance, weighted_mean, ISS_gaussian_means, ISS_gaussian_covs, initial_ISS_gaussian_means, initial_ISS_gaussian_covs, Weights] = GenerateGSF(mu, P, N, numSamples, startTime, stopTime, current_alpha);
    
    % Store results
    results.(alpha_str).weighted_covariance = weighted_covariance;
    results.(alpha_str).weighted_mean = weighted_mean;
    results.(alpha_str).alpha_value = current_alpha;
    
    % Calculate determinant and trace for position covariance (3x3)
    % Convert to km units (divide by 1000^2 for covariance)
    position_covs = weighted_covariance(1:3, 1:3, :) / (1000^2); % Convert m^2 to km^2
    det_values = zeros(1, numTimeSteps);
    trace_values = zeros(1, numTimeSteps);
    
    for t = 1:numTimeSteps
        det_values(t) = det(position_covs(:, :, t)); % Now in km^6
        trace_values(t) = trace(position_covs(:, :, t)); % Now in km^2
    end
    
    results.(alpha_str).determinant = det_values;
    results.(alpha_str).trace = trace_values;
end

timeVector = datetime(timeVector, 'TimeZone', '');
time_hours = hours(timeVector - startTime);

%% Plot Determinant and Trace for Different Alpha Values
figure('Position', [100, 100, 1400, 800]);

% Subplot 1: Determinant
subplot(2, 1, 1);
colors = {'b-', 'r-', 'g-', 'm-'};
hold on;

field_names = fieldnames(results);
for i = 1:length(field_names)
    alpha_str = field_names{i};
    current_alpha = results.(alpha_str).alpha_value;
    det_values = results.(alpha_str).determinant;
    
    plot(time_hours, det_values, colors{i}, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('\\alpha = %g', current_alpha));
end

xlabel('Time (hours)');
ylabel('Determinant of Position Covariance Matrix (km^6)'); % Updated units
title('Determinant of Position Covariance Matrix vs Time for Different Alpha Values');
legend('Location', 'best');
grid on;
set(gca, 'YScale', 'log'); % Log scale for better visualization

% Subplot 2: Trace
subplot(2, 1, 2);
hold on;

for i = 1:length(field_names)
    alpha_str = field_names{i};
    current_alpha = results.(alpha_str).alpha_value;
    trace_values = results.(alpha_str).trace;
    
    plot(time_hours, trace_values, colors{i}, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('\\alpha = %g', current_alpha));
end

xlabel('Time (hours)');
ylabel('Trace of Position Covariance Matrix (km^2)'); % Updated units
title('Trace of Position Covariance Matrix vs Time for Different Alpha Values');
legend('Location', 'best');
grid on;

sgtitle('Covariance Matrix Analysis for Different Unscented Transform Alpha Parameters');

%% Additional Plot: Position Magnitude for Default Alpha
% Use the last processed alpha (1e-3) for position magnitude plot
default_result = results.(field_names{end});
ISS_positions = default_result.weighted_mean(1:3, :);
ISS_position_covs = default_result.weighted_covariance(1:3, 1:3, :);
ISS_positions_km = ISS_positions / 1000;

figure('Position', [200, 200, 1400, 600]);
position_magnitude = sqrt(sum(ISS_positions_km.^2, 1));

plot(time_hours, position_magnitude, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Position Magnitude');
hold on;

% Calculate uncertainty in magnitude (1σ bounds)
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
title(sprintf('Position Magnitude against Time (\\alpha = %g)', default_result.alpha_value));
legend('Location', 'best');
grid on;

function [ISS_mixture_cov, ISS_mixture_mean, ISS_gaussian_means, ISS_gaussian_covs, ISS_mu_components, ISS_P_components, gaussian_weights] = GenerateGSF(mu, P, N, numSamples, startTime, stopTime, alpha)
% Added alpha as an input parameter with default value
if nargin < 7
    alpha = 1e-3; % Default alpha value
end

sampleTime = 120; % seconds
timeVector = startTime:seconds(sampleTime):stopTime;
numTimeSteps = length(timeVector);

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
fprintf('Using alpha = %g\n', alpha);
gm = fitgmdist(data, N, 'Options', options, 'RegularizationValue', 1e-30);

N_components = gm.NumComponents; 

ISS_mu_components = gm.mu';            
ISS_P_components = gm.Sigma;       

Weights = ones(1, N_components) / N;
ISS_all_sigma_points = zeros(n, n_sigma, N_components);
ISS_all_Wm = zeros(n_sigma, N_components);
ISS_all_Wc = zeros(n_sigma, N_components);

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
    
    cov_jt = zeros(6, 6);
    for i = 1:n_sigma
        diff = ISS_all_sigma_points(:,i, j) - mean_jt;
        cov_jt = cov_jt + Wc_j(i) * (diff * diff');
    end
    
    ISS_predicted_cov_keplerian(:,:,j) = cov_jt;
end

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
    uncertainty.a = 100;
    uncertainty.h = 1e-4;
    uncertainty.k = 1e-4;
    uncertainty.p = 1e-4;
    uncertainty.q = 1e-4;
    uncertainty.l = 1e-3*pi/180;
    %}
    
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



