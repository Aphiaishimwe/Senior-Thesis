function RISBeamSteering(N, alpha, beta, target_angle, max_iterations)
    % Input values:
    lambda = 1;              % Wavelength (normalized to 1)
    d = lambda / 2;          % Spacing between RIS elements
    k = 2 * pi / lambda;     % Wave number

    % Convert angles to radians:
    alpha_rad = deg2rad(alpha);

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

    % Initialize variables for hierarchical refinement
    current_beta = beta; % Start with the given reflection angles
    max_power = -Inf;    % Initialize max power to a very low value
    iteration_powers = {}; % Store powers for each iteration
    iteration_betas = {};  % Store reflection angles for each iteration
    best_angle = 0;

    % Loop through the maximum number of iterations
    for iteration = 1:max_iterations
        % Array to store powers at the target angle for each reflection angle:
        powers_at_target = zeros(size(current_beta));

        % Calculate powers for current reflection angles
        for b_idx = 1:length(current_beta)
            % Desired RIS to Rx phase shift for current reflection angle:
            beta_rad = deg2rad(current_beta(b_idx));
            h2 = exp(1i * k * d * (0:N-1)' * cos(beta_rad));
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

            % Extract the power at the target angle:
            powers_at_target(b_idx) = array_factor_norm(target_idx);
        end

        % Store powers and reflection angles for this iteration
        iteration_powers{iteration} = 10 * log10(powers_at_target);
        iteration_betas{iteration} = current_beta;

        % Plot power vs reflection angle for the current iteration
        figure; % Create a new figure for each iteration
        plot(current_beta, iteration_powers{iteration}, 'bo-', 'LineWidth', 1.5);
        grid on;
        xlabel('Reflection Angle (°)');
        ylabel('Power (dB)');
        title(['Iteration ', num2str(iteration), ': Power vs Reflection Angle']);
        drawnow; % Force MATLAB to display the current plot

        % If it's the last iteration, add an additional bar plot
        if iteration == max_iterations
            figure; % Create a new figure for the bar plot
            bar(current_beta, iteration_powers{iteration}, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'k', 'LineWidth', 1.5);
            grid on;
            xlabel('Reflection Angle (°)');
            ylabel('Power (dB)');
            title(['Bar Plot: Iteration ', num2str(iteration), ' - Power vs Reflection Angle']);
        end

        % Find the maximum power and its corresponding reflection angle:
        [max_power, max_idx] = max(powers_at_target);
        best_angle = current_beta(max_idx); % Best reflection angle

        % Refine the reflection angles around the best angle:
        step = (current_beta(2) - current_beta(1)) / 2; % Halve the step size
        new_beta = best_angle + (-step:step/2:step); % Generate new search segment
        current_beta = new_beta(new_beta >= 0 & new_beta <= 180);% Enforce valid angular range

    end

    % Print the final best observation angle and its corresponding power
    fprintf('Final Best Observation Angle: %.2f degrees\n', best_angle);
    fprintf('The mean squared error: %.4f degrees\n', (best_angle - target_angle).^2);
    fprintf('Final Maximum Power: %.2f dB\n', 10 * log10(max_power));

    % Print the angle corresponding to max_idx and its power
    fprintf('Max Power Found at Angle: %.2f degrees with Power: %.2f dB\n', best_angle, 10 * log10(max_power));
end

% Define input values
alpha = 30;               % Incident angle in degrees
beta = linspace(10, 170, 3);         % Reflection angles (degrees)
N = 16;                   % Number of RIS elements
target_angle = 20;        % Desired observation angle in degrees
max_iterations = 3;

% Call the function to perform hierarchical beam steering
RISBeamSteering(N, alpha, beta, target_angle, max_iterations);
