function RISExhaustiveBeamSteering(N, alpha, beta, target_angle)
    % Input values:
    lambda = 1;              % Wavelength (normalized to 1)
    d = lambda / 2;          % Spacing between RIS elements
    k = 2 * pi / lambda;     % Wave number
    
    % Convert angles to radians:
    alpha_rad = deg2rad(alpha);
    beta_rad = deg2rad(beta);

    % Tx to RIS phase shift
    h1 = exp(1i * k * d * (0:N-1)' * cos(alpha_rad));

    % Observation angles:
    theta = linspace(0, 180, 1000);  % Observation angles in degrees
    theta_rad = deg2rad(theta);      % Convert to radians

    % Find the exact index corresponding to the target angle:
    target_idx = find(theta == target_angle);

    if isempty(target_idx)
        % If no exact match for the target angle, find the closest value:
        [~, target_idx] = min(abs(theta - target_angle)); 
    end

    % Array to store powers at the target angle for each reflection angle:
    powers_at_target = zeros(size(beta));
    
    % Matrix to store the full power distribution
    power_matrix = zeros(length(beta), length(theta));

    for b_idx = 1:length(beta)
        % Desired RIS to Rx phase shift for current reflection angle:
        h2 = exp(1i * k * d * (0:N-1)' * cos(beta_rad(b_idx)));
        theta_opt = angle(h1) - angle(h2); % Optimal phase shifts
        Theta = diag(exp(1i * theta_opt)); % RIS phase configuration

        % Calculate the beam pattern for this reflection angle:
        array_factor = zeros(size(theta));
        for idx = 1:length(theta)
            h_obs = exp(1i * k * d * (0:N-1)' * cos(theta_rad(idx)));
            array_factor(idx) = abs(h1' * Theta * h_obs)^2; % Beam pattern
        end

        % Normalize the beam pattern:
        array_factor_norm = array_factor / max(array_factor);
        power_matrix(b_idx, :) = array_factor_norm; % Store full pattern

        % Extract the power at the target angle:
        powers_at_target(b_idx) = array_factor_norm(target_idx);
    end

    % Find the observation angle that resulted in the highest power at the target angle
    [~, max_power_idx] = max(powers_at_target);
    best_observation_angle = beta(max_power_idx);
    max_power_value = 10 * log10(powers_at_target(max_power_idx));
    
    % Display the results
    fprintf('Observation Angle that resulted in highest power at target angle: %.2f degrees\n', best_observation_angle);
    fprintf('Mean Absolute Error: %.2f degrees\n', abs(best_observation_angle - target_angle));
    fprintf('Maximum Power at Target Angle: %.2f dB\n', max_power_value);
    
    % Plot powers at the target angle using vertical lines:
    figure;
    hold on;
    for b_idx = 1:length(beta)
        x = b_idx; % X-position corresponds to the reflection angle index
        y = 10 * log10(powers_at_target(b_idx)); % Power in dB at the target angle

        % Plot the vertical line from power point to x-axis (solid line)
        plot([x, x], [y, 0], 'k-', 'LineWidth', 1.5); % Solid vertical line
        plot(x, y, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5); % Point at the top of the line
    end
    hold off;

    grid on;
    xticks(1:length(beta)); % Set x-axis ticks to reflection angle indices
    xticklabels(string(beta)); % Label ticks with reflection angles
    xlabel('Reflection Angle Index (β)');
    ylabel(['Normalized Gain at ', num2str(target_angle), '° (dB)']);
    title(['Gain at User Angle ', num2str(target_angle), '° for Each Reflection Angle']);
    legend('Power at target angle', 'Location', 'best');

    % Add an additional bar plot for the power at the target angle:
    figure;
    bar(beta, 10 * log10(powers_at_target), 'FaceColor', [0.1, 0.5, 0.2], 'EdgeColor', 'k', 'LineWidth', 1.5);

    % Configure the bar chart:
    grid on;
    xlabel('Reflection Angle (°)');
    ylabel(['Normalized Gain at ', num2str(target_angle), '° (dB)']);
    title(['Bar Plot of Gain at User Angle ', num2str(target_angle), '°']);
end

% Define input values
alpha = 45;               % Incident angle in degrees
beta = 10:10:170;      % Reflection angles (degrees)
N = 16;                   % Number of RIS elements
target_angle = 60;        % Desired observation angle in degrees

% Call the function to plot the RIS beam patterns and power at the target angle
RISExhaustiveBeamSteering(N, alpha, beta, target_angle);
