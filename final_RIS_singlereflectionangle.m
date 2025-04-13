% Input values
N = 64;                  % Number of RIS elements
lambda = 1;              % Wavelength (normalized)
d = lambda / 2;          % Spacing between RIS elements
k = 2 * pi / lambda;     % Wave number
alpha = 45;              % Incident angle (degrees) 
beta = 60;               % Desired reflection angle (degrees) (target)

% Convert angles to radians
alpha_rad = deg2rad(alpha);
beta_rad = deg2rad(beta);

% Phase shifts from Tx to RIS and RIS to target
h1 = exp(1i * k * d * (0:N-1)' * cos(alpha_rad));
h2 = exp(1i * k * d * (0:N-1)' * cos(beta_rad));

% Optimal RIS phase shift
theta_opt = angle(h1) - angle(h2);
Theta = diag(exp(1i * theta_opt));

% Observation angles
theta = linspace(0, 180, 1000);  
theta_rad = deg2rad(theta);

% Compute array factor over all angles
array_factor = zeros(size(theta));
for idx = 1:length(theta)
    h_obs = exp(1i * k * d * (0:N-1)' * cos(theta_rad(idx)));
    array_factor(idx) = abs(h1' * Theta * h_obs)^2;
end

% Normalize and convert to dB
array_factor_norm = array_factor / max(array_factor);
array_factor_db = 10 * log10(array_factor_norm + eps);

% Plot beam pattern
figure; 
plot(theta, array_factor_db, 'LineWidth', 1.5);
hold on;

% Add red star at target angle
[~, target_idx] = min(abs(theta - beta));
plot(beta, array_factor_db(target_idx), 'r*', 'MarkerSize', 10, 'LineWidth', 2);

% Labels and title
grid on;
xlabel('Observation Angle (°)');
ylabel('Normalized Gain (dB)');
title(sprintf('RIS Beam Pattern | Incident Angle = %d°, Target Angle = %d°', alpha, beta));
ylim([-50, 0]);
legend('Beampattern', 'Target Angle', 'Location', 'best');
hold off;

% Display main lobe info
[max_gain, max_idx] = max(array_factor_norm);
main_lobe_angle = theta(max_idx);

disp(['Main Lobe Angle: ', num2str(main_lobe_angle), '°']);
disp(['Maximum Gain: ', num2str(10*log10(max_gain)), ' dB']);
