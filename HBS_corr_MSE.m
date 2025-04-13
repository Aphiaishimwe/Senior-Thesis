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
    iteration_powers = {}; % Store powers for each iteration
    iteration_betas = {};  % Store reflection angles for each iteration
    snr_values = -20:10:30; % Define SNR range
    selected_iterations = 1:5:max_iterations; % Iterations to analyze
    
    % Storage for MSE and Correlation for selected iterations
    mse_results = zeros(length(selected_iterations), length(snr_values));
    correlation_results = zeros(length(selected_iterations), length(snr_values));
    iteration_labels = strings(1, length(selected_iterations));
    iter_idx = 1;

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

        % Store powers for this iteration
        iteration_powers{iteration} = powers_at_target;

        % Set estimated signature from first iteration
        if iteration == 1
            estimated_signature = powers_at_target;
        end

        % Compute MSE and Correlation for different SNR levels at selected iterations
        if ismember(iteration, selected_iterations)
            for s = 1:length(snr_values)
                noisy_signature = awgn(estimated_signature, snr_values(s), 'measured');
                mse_results(iter_idx, s) = mean((noisy_signature - estimated_signature).^2);
                correlation_results(iter_idx, s) = corr(estimated_signature', noisy_signature');
            end
            
            % Store iteration labels
            iteration_labels(iter_idx) = ['Iteration ', num2str(iteration)];
            iter_idx = iter_idx + 1;
        end
    end

    % Final combined MSE vs SNR plot
    figure;
    hold on;
    cmap = lines(length(selected_iterations));
    for i = 1:length(selected_iterations)
        plot(snr_values, mse_results(i, :), 'o-', 'Color', cmap(i, :), 'LineWidth', 1.5, 'MarkerSize', 6);
    end
    grid on;
    xlabel('SNR (dB)');
    ylabel('Mean Squared Error (MSE)');
    title('MSE vs. SNR for Selected Iterations, Hierarchical');
    legend(iteration_labels, 'Location', 'bestoutside');
    hold off;
    
    % Final combined Correlation vs SNR plot
    figure;
    hold on;
    for i = 1:length(selected_iterations)
        plot(snr_values, correlation_results(i, :), 's-', 'Color', cmap(i, :), 'LineWidth', 1.5, 'MarkerSize', 6);
    end
    grid on;
    xlabel('SNR (dB)');
    ylabel('Correlation with Estimated Signature');
    title('Correlation vs. SNR for Selected Iterations, Hierarchical');
    legend(iteration_labels, 'Location', 'bestoutside');
    hold off;
    
    % Print the total number of iterations it took:
    disp(['Total number of iterations: ', num2str(max_iterations)]);
end

% Define input values
alpha = 30;               % Incident angle in degrees
beta = 10:40:170;          % Reflection angles (degrees)
N = 64;                   % Number of RIS elements
target_angle = 30;        % Desired observation angle in degrees
max_iterations = 26;

% Call the function to perform hierarchical beam steering
RISBeamSteering(N, alpha, beta, target_angle, max_iterations);
