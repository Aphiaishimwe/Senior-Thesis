function RISRandomPhases(N, alpha, target_angle, num_iterations)
    % Inputs:
    %   N: Number of RIS elements
    %   alpha: Incident angle (degrees)
    %   target_angle: Desired observation angle (degrees)
    %   num_iterations: Number of random phase trials

    if nargin < 4
        num_iterations = 5; % Default to 5 beam patterns
    end

    lambda = 1;              % Wavelength (normalized)
    d = lambda / 2;          % RIS element spacing
    k = 2 * pi / lambda;     % Wave number

    % Convert angles to radians:
    alpha_rad = deg2rad(alpha);
    
    % Incident wave phase shifts (Tx to RIS):
    h1 = exp(1i * k * d * (0:N-1)' * cos(alpha_rad));

    % Observation angles:
    theta = linspace(0, 90, 50); % 50 beam indices
    theta_rad = deg2rad(theta);

    % Precompute target angle index:
    [~, target_idx] = min(abs(theta - target_angle));

    % Define SNR levels:
    SNR_levels = 0:5:30;

    % Initialize storage:
    clean_signature = zeros(1, length(theta));
    noisy_signatures = zeros(length(SNR_levels), length(theta));
    correlation_values = zeros(1, length(SNR_levels)); % Store correlation values

    % Generate one clean power signature:
    theta_opt = 2 * pi * rand(N, 1); 
    Theta = diag(exp(1i * theta_opt));

    for idx = 1:length(theta)
        h_obs = exp(1i * k * d * (0:N-1)' * cos(theta_rad(idx)));
        clean_signature(idx) = abs(h1' * Theta * h_obs)^2;
    end
    
    % Normalize to dB:
    clean_signature = 10 * log10(clean_signature / max(clean_signature));

    % Add noise at different SNR levels and compute correlation:
    for s = 1:length(SNR_levels)
        noisy_signatures(s, :) = awgn(clean_signature, SNR_levels(s), 'measured');
        correlation_values(s) = corr(clean_signature', noisy_signatures(s, :)'); % Compute correlation
    end

    % Plot clean vs noisy signatures:
    figure;
    hold on;
    cmap = lines(length(SNR_levels));
    
    plot(theta, clean_signature, 'k-', 'LineWidth', 2, 'DisplayName', 'Clean Signature'); % Clean version

    for s = 1:length(SNR_levels)
        plot(theta, noisy_signatures(s, :), 'Color', cmap(s, :), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('SNR = %d dB', SNR_levels(s)));
    end
    
    hold off;
    grid on;
    xlabel('Beam Index');
    ylabel('Normalized Gain (dB)');
    title(['Power Signature at ', num2str(target_angle), '° with Noise']);
    legend('show');

    % Plot Correlation vs SNR:
    figure;
    plot(SNR_levels, correlation_values, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    xlabel('SNR (dB)');
    ylabel('Correlation Coefficient');
    title('Correlation between Noisy and Clean Signatures');
    ylim([-1, 1]);

    % Save results to CSV:
    data_table = array2table([SNR_levels(:), correlation_values(:)], ...
        'VariableNames', {'SNR_dB', 'Correlation'});
    writetable(data_table, 'correlation_vs_snr.csv');
    fprintf('Correlation values saved to CSV file.\n');
end

% Example usage:
alpha = 30;          % Incident angle (degrees)
target_angle = 30;   % User angle (degrees)
N = 16;              % RIS elements
num_iterations = 50; % Number of beam patterns

% Call the function:
RISRandomPhases(N, alpha, target_angle, num_iterations);
