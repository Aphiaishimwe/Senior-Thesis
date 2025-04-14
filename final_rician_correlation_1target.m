function RISRandomPhases_MultipathVirtualAP_SweepK(N, alpha, num_iterations, angle_ref, num_runs, K_values)
    rng(40);

    lambda = 1;
    d = lambda / 2;
    k = 2 * pi / lambda;

    alpha_rad = deg2rad(alpha);
    h1 = exp(1i * k * d * (0:N-1)' * cos(alpha_rad));

    theta = linspace(0, 90, 100);
    theta_rad = deg2rad(theta);
    [~, ref_angle_index] = min(abs(theta - angle_ref));

    num_K = length(K_values);
    correlation_accum = zeros(num_K, length(theta));
    power_accum = zeros(num_K, num_iterations);
    clean_power_accum = zeros(1, num_iterations);
    clean_correlation = zeros(1, length(theta));
    peak_angle_accum = zeros(num_runs, num_K);

    theta_filename = sprintf('ThetaList_N%d_Iter%d.mat', N, num_iterations);
    if isfile(theta_filename)
        load(theta_filename, 'Theta_list');
        fprintf("Loaded RIS Theta_list from %s\n", theta_filename);
    else
        error("Theta_list file not found. Please generate it first using generate_fixed_theta_list.m");
    end

    for run = 1:num_runs
        % === Clean LOS Power ===
        power_matrix_clean = zeros(num_iterations, length(theta));
        for iter = 1:num_iterations
            Theta = Theta_list{iter};
            for idx = 1:length(theta)
                h_obs_LOS = exp(1i * k * d * (0:N-1)' * cos(theta_rad(idx)));
                power_matrix_clean(iter, idx) = abs(h1' * Theta * h_obs_LOS)^2;
            end
        end

        power_matrix_clean = power_matrix_clean ./ vecnorm(power_matrix_clean, 2, 2);
        power_matrix_clean_db = 10 * log10(power_matrix_clean + eps);
        reference_signature = power_matrix_clean_db(:, ref_angle_index);
        % Inject AWGN using MATLAB's built-in function
        % reference_signature = awgn(reference_signature, 5, 'measured');
        clean_power_accum = clean_power_accum + reference_signature';

        for idx = 1:length(theta)
            clean_correlation(idx) = clean_correlation(idx) + corr(reference_signature, power_matrix_clean_db(:, idx));
        end

        % === Loop over K values ===
        for k_idx = 1:num_K
            K = K_values(k_idx);

            % Rician fading config (ITU Pedestrian B)
            fs = 3.84e9;
            pathDelays = [0 200 800 1200 2300 3700]*1e-9;
            avgPathGains = [0 -0.9 -4.9 -8 -7.8 -23.9];
            fD = 50;

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
                    h_obs_LOS = exp(1i * k * d * (0:N-1)' * cos(theta_rad(idx)));

                    h_obs_total = zeros(N, 1);
                    for n = 1:N
                        [faded_sig, ~] = ricianChan(h_obs_LOS(n));
                        h_obs_total(n) = faded_sig;
                    end

                    power_matrix(iter, idx) = abs(h1' * Theta * h_obs_total)^2;
                end
            end

            power_matrix = power_matrix ./ vecnorm(power_matrix, 2, 2);
            power_matrix_db = 10 * log10(power_matrix + eps);
            power_accum(k_idx, :) = power_accum(k_idx, :) + power_matrix_db(:, ref_angle_index)';

            for idx = 1:length(theta)
                correlation_accum(k_idx, idx) = correlation_accum(k_idx, idx) + corr(reference_signature, power_matrix_db(:, idx));
            end

            [~, peak_idx] = max(correlation_accum(k_idx, :) / run);
            peak_angle_accum(run, k_idx) = theta(peak_idx);
        end
    end

    avg_corr = correlation_accum / num_runs;
    avg_clean_corr = clean_correlation / num_runs;
    avg_power = power_accum / num_runs;
    avg_clean_power = clean_power_accum / num_runs;
    mean_peak_angles = mean(peak_angle_accum, 1);
    abs_errors = abs(mean_peak_angles - angle_ref);
    mae = mean(abs_errors);

    % === Plot 1: Correlation vs Angle ===
    figure;
    hold on;
    cmap = lines(num_K);
    legend_entries = cell(num_K+1, 1);
    for k_idx = 1:num_K
        plot(theta, avg_corr(k_idx, :), 'Color', cmap(k_idx,:), 'LineWidth', 1.5);
        legend_entries{k_idx} = sprintf('K = %.1f (Peak Corr = %.3f @ %.0f°)', ...
            K_values(k_idx), max(avg_corr(k_idx, :)), mean_peak_angles(k_idx));
    end
    [max_clean_corr, clean_peak_idx] = max(avg_clean_corr);
    clean_peak_angle = theta(clean_peak_idx);
    plot(theta, avg_clean_corr, 'k--', 'LineWidth', 2);
    legend_entries{end} = sprintf('Clean LOS (Peak Corr = %.3f @ %.0f°)', max_clean_corr, clean_peak_angle);
    legend(legend_entries, 'Location', 'eastoutside');
    grid on;
    xlabel('Observation Angle (°)');
    ylabel('Correlation with Clean Signature');
    title(sprintf('Avg Corr vs Angle | Ref = %d°, MAE = %.2f°', angle_ref, mae));
    ylim([-1 1]);
    hold off;

    % === Plot 2: Power vs Beam Index ===
    figure;
    hold on;
    plot(1:num_iterations, avg_clean_power, 'k-', 'LineWidth', 2, 'DisplayName', 'Clean LOS');
    for k_idx = 1:num_K
        plot(1:num_iterations, avg_power(k_idx, :), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('K = %.1f', K_values(k_idx)));
    end
    grid on;
    xlabel('Beam Index (Iteration)');
    ylabel('Power at Reference Angle (dB)');
    title(sprintf('Power vs Beam Index @ %d° (Avg over %d Runs)', angle_ref, num_runs));
    legend('show', 'Location', 'southwest');
    hold off;

    % === Plot 3: Abs Error vs K ===
    figure;
    plot(K_values, abs_errors, 'mo-', 'LineWidth', 2);
    hold on;
    yline(mae, 'r--', 'LineWidth', 2, ...
        'Label', sprintf('Mean Abs Error = %.2f°', mae), ...
        'LabelHorizontalAlignment', 'center', ...
        'LabelVerticalAlignment', 'bottom');
    grid on;
    xlabel('Rician K-Factor');
    ylabel('Absolute Error (°)');
    title(sprintf('Angle Estimation Error vs K | Ref = %d°', angle_ref));
    hold off;
end
% Parameters
N = 16;
alpha = 30;
num_iterations = 20;
angle_ref = 30;
num_runs = 1;
K_values = 5:5:20;

% Run the simulation
RISRandomPhases_MultipathVirtualAP_SweepK(N, alpha, num_iterations, angle_ref, num_runs, K_values);
