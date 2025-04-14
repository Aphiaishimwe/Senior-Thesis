function RISRandomPhases_MultipathVirtualAP_GeneralizedK(N, alpha, num_iterations, num_runs, K_values)
    rng(40);

    lambda = 1;
    d = lambda / 2;
    k_wave = 2 * pi / lambda;

    theta = linspace(0, 90, 100);            % Angular scan grid
    theta_rad = deg2rad(theta);

    angle_sweep = 0:1:90;
    num_angles = length(angle_sweep);
    num_K = length(K_values);

    mean_abs_errors = zeros(1, num_K);
    SNR_dB = 5;  % Additive White Gaussian Noise SNR

    % === Precompute LOS steering vectors for all angles ===
    h_obs_LOS_all = zeros(N, length(theta));
    for idx = 1:length(theta)
        h_obs_LOS_all(:, idx) = exp(1i * k_wave * d * (0:N-1)' * cos(theta_rad(idx)));
    end

    % === Load RIS phase configuration list ===
    theta_filename = sprintf('ThetaList_N%d_Iter%d.mat', N, num_iterations);
    if isfile(theta_filename)
        load(theta_filename, 'Theta_list');
        fprintf(" Loaded RIS Theta_list from %s\n", theta_filename);
    else
        error("Theta_list file not found. Please generate it first using generate_fixed_theta_list.m");
    end

    % === Rician fading config ===
    fs = 3.84e6;
    pathDelays = [0 200 800 1200 2300 3700]*1e-9;
    avgPathGains = [0 -0.9 -4.9 -8 -7.8 -23.9];
    fD = 50;

    for k_idx = 1:num_K
        K = K_values(k_idx);
        fprintf(" Simulating for K = %.2f\n", K);
        abs_errors_all_angles = zeros(1, num_angles);

        for a_idx = 1:num_angles
            angle_ref = angle_sweep(a_idx);
            alpha_rad = deg2rad(alpha);
            h1 = exp(1i * k_wave * d * (0:N-1)' * cos(alpha_rad));
            [~, ref_idx] = min(abs(theta - angle_ref));

            % === Clean Signature ===
            power_matrix_clean = zeros(num_iterations, length(theta));
            for iter = 1:num_iterations
                Theta = Theta_list{iter};
                for idx = 1:length(theta)
                    h_obs_LOS = h_obs_LOS_all(:, idx);
                    power_matrix_clean(iter, idx) = abs(h1' * Theta * h_obs_LOS)^2;
                end
            end

            power_matrix_clean = power_matrix_clean ./ vecnorm(power_matrix_clean, 2, 2);
            power_matrix_clean_db = 10 * log10(power_matrix_clean + eps);
            reference_signature = power_matrix_clean_db(:, ref_idx);

            %  Add AWGN to reference signature
            reference_signature = awgn(reference_signature, SNR_dB, 'measured');

            peak_angles = zeros(num_runs, 1);

            for run = 1:num_runs
                % === Rician channel ===
                ricianChan = comm.RicianChannel( ...
                    'SampleRate', fs, ...
                    'PathDelays', pathDelays, ...
                    'AveragePathGains', avgPathGains, ...
                    'KFactor', K, ...
                    'MaximumDopplerShift', fD, ...
                    'ChannelFiltering', true, ...
                    'DopplerSpectrum', doppler('Bell', 8), ...
                    'RandomStream', 'mt19937ar with seed', ...
                    'Seed', 100 + k_idx + run, ...
                    'PathGainsOutputPort', true ...
                );
                ricianChan(1 + 1i);  % Prime

                power_matrix = zeros(num_iterations, length(theta));
                for iter = 1:num_iterations
                    Theta = Theta_list{iter};
                    for idx = 1:length(theta)
                        h_obs_LOS = h_obs_LOS_all(:, idx);
                        [faded_vector, ~] = ricianChan(h_obs_LOS);
                        h_obs_total = faded_vector(:);
                        power_matrix(iter, idx) = abs(h1' * Theta * h_obs_total)^2;
                    end
                end

                power_matrix = power_matrix ./ vecnorm(power_matrix, 2, 2);
                power_matrix_db = 10 * log10(power_matrix + eps);

                % === Correlate each angle signature with noisy reference ===
                correlations = zeros(1, length(theta));
                for idx = 1:length(theta)
                    noisy_sig = awgn(power_matrix_db(:, idx), SNR_dB, 'measured');
                    correlations(idx) = corr(reference_signature, noisy_sig);
                end

                [~, peak_idx] = max(correlations);
                peak_angles(run) = theta(peak_idx);
            end

            abs_errors_all_angles(a_idx) = mean(abs(peak_angles - angle_ref));
        end

        mean_abs_errors(k_idx) = mean(abs_errors_all_angles);
    end

    % === Plot MAE vs K ===
    figure;
    plot(K_values, mean_abs_errors, 'mo-', 'LineWidth', 2, 'MarkerFaceColor', 'm');
    grid on;
    xlabel('Rician K-Factor');
    ylabel('Mean Absolute Error (°)');
    title(sprintf('MAE vs K-Factor (SNR = %d dB)', SNR_dB));

    % === Save MAE vs K results to CSV ===
    results_table = table(K_values(:), mean_abs_errors(:), ...
        'VariableNames', {'K_Factor', 'Mean_Absolute_Error'});
    
    csv_filename = sprintf('MAE_vs_lowK_SNR_%ddB.csv', SNR_dB);
    writetable(results_table, csv_filename);
    fprintf("📁 MAE results saved to CSV: %s\n", csv_filename);
end

% === Run Parameters ===
N = 16;
alpha = 30;
num_iterations = 20;
num_runs = 1;
K_values = 0.1:0.1:1.2;

RISRandomPhases_MultipathVirtualAP_GeneralizedK(N, alpha, num_iterations, num_runs, K_values);
