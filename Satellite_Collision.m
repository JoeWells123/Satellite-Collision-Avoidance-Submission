%% Time Setup
%startTime = datetime(2025, 8, 29, 19, 59, 48.154);
startTime = datetime(2025, 9, 7, 0, 0, 0);
stopTime = startTime + hours(64);
sampleTime = 10; % second
timeVector = startTime:seconds(sampleTime):stopTime;
numTimeSteps = length(timeVector);
time_hours = hours(timeVector - startTime);

sc_mean = satelliteScenario(startTime, stopTime, sampleTime);
sat1 = satellite(sc_mean, 'SL24.xml', 'OrbitPropagator','two-body-keplerian');
sat2 = satellite(sc_mean, 'COS.xml', 'OrbitPropagator','two-body-keplerian');
sat1.ShowLabel = false;
sat2.ShowLabel = false;
sat1.MarkerColor = "blue";
sat2.MarkerColor = "green";
%play(sc_mean);

[positions1, ~] = states(sat1);
[positions2, ~] = states(sat2);
relative_positions_sgp4 = positions1 - positions2;
distance = sqrt(sum(relative_positions_sgp4.^2, 1)); 
disp(min(distance))


figure('Position', [100, 100, 1200, 800]);
hold on;

plot(time_hours, distance, 'b-', 'LineWidth', 1.8, 'DisplayName', 'SGP4 vs Two-Body (LEO)');

% Labels & aesthetics
xlabel('Time (hours)', 'FontSize', 12);
ylabel('Distance (m)', 'FontSize', 12);
title('Distance between SGP4/SDP4 and Two-Body Keplerian Models', 'FontSize', 14);

legend('Location', 'northwest');
grid on;
box on;

% Optional: make y-axis log scale if values vary a lot
% set(gca, 'YScale', 'log');

hold off;
%}
%{
sat_low_satellite = satellite(sc_mean, 'Cosmos 2251.xml', 'Name', 'Satellite (two-body-keplerian)-Blue', OrbitPropagator='two-body-keplerian');
sat_low_satellite.MarkerColor = "blue";
sat_low_debris = satellite(sc_mean, 'Cosmos 2251.xml', 'Name', 'Satellite (SGP4)-Green', OrbitPropagator = 'sgp4');
sat_low_debris.MarkerColor = "green";
sat_high_satellite = satellite(sc_mean, 'Delta 1.xml', 'Name', 'Satellite (two-body-keplerian)-Blue', OrbitPropagator='two-body-keplerian');
sat_high_satellite.MarkerColor = "blue";
sat_high_debris = satellite(sc_mean, 'Delta 1.xml', 'Name', 'Satellite (SDP4)-Magenta', OrbitPropagator = 'sdp4');
sat_high_debris.MarkerColor = "magenta";
play(sc_mean);


% Extract states
positions_two_body_1 = states(sat_low_satellite);
positions_two_body_2 = states(sat_high_satellite);
positions_sgp4       = states(sat_low_debris);
positions_sdp4       = states(sat_high_debris);

% Compute relative distances
relative_positions_sgp4 = positions_two_body_1 - positions_sgp4;
relative_positions_sdp4 = positions_two_body_2 - positions_sdp4;

distance_sgp4 = sqrt(sum(relative_positions_sgp4.^2, 1)); % [1, numTimeSteps]
distance_sdp4 = sqrt(sum(relative_positions_sdp4.^2, 1)); % [1, numTimeSteps]

% Plot
figure('Position', [100, 100, 1200, 800]);
hold on;

plot(time_hours, distance_sgp4, 'b-', 'LineWidth', 1.8, 'DisplayName', 'SGP4 vs Two-Body (LEO)');
plot(time_hours, distance_sdp4, 'r--', 'LineWidth', 1.8, 'DisplayName', 'SDP4 vs Two-Body (HEO)');

% Labels & aesthetics
xlabel('Time (hours)', 'FontSize', 12);
ylabel('Distance (m)', 'FontSize', 12);
title('Distance between SGP4/SDP4 and Two-Body Keplerian Models', 'FontSize', 14);

legend('Location', 'northwest');
grid on;
box on;

% Optional: make y-axis log scale if values vary a lot
% set(gca, 'YScale', 'log');

hold off;
%}
