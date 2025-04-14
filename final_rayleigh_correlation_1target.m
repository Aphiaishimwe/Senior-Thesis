function RISRandomPhases_MultipathVirtualAP_SweepG(N, alpha, num_iterations, angle_ref, num_runs)
    lambda = 1;
    d = lambda / 2;
    k = 2 * pi / lambda;
    rng(50);

    alpha_rad = deg2rad(alpha);
    h1 = exp(1i * k * d * (0:N-1)' * cos(alpha_rad));

    theta = linspace(0, 90, 100);
    theta_rad = deg2rad(theta);
    [~, ref_angle_index] = min(abs(theta - angle_ref));

    nlos_linear_values = 1:2:10;
    num_nlos = length(nlos_linear_values);

    correlation_accum = zeros(num_nlos, length(theta));
    power_accum = zeros(num_nlos, num_iterations);
    clean_power_accum = zeros(1, num_iterations);
    clean_correlation = zeros(1, length(theta));
    peak_angle_accum = zeros(num_runs, num_nlos);

    % === Load RIS Theta list ===
    theta_filename = sprintf('ThetaList_N%d_Iter%d.mat', N, num_iterations);
    if isfile(theta_filename)
        load(theta_filename, 'Theta_list');
        fprintf("Loaded RIS Theta_list from %s\n", theta_filename);
    else
        error("Theta_list file not found. Please generate it using generate_fixed_theta_list.m");
    end

    for run = 1:num_runs
        fprintf("Run %d of %d\n", run, num_runs);

        % === Clean LOS power ===
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
        clean_power_accum = clean_power_accum + reference_signature';

        for idx = 1:length(theta)
            clean_correlation(idx) = clean_correlation(idx) + corr(reference_signature, power_matrix_clean_db(:, idx));
        end

        % === Loop over NLOS gains (g values) ===
        for g_idx = 1:num_nlos
            g = nlos_linear_values(g_idx);
            power_matrix = zeros(num_iterations, length(theta));

            for iter = 1:num_iterations
                Theta = Theta_list{iter};
                for idx = 1:length(theta)
                    % h_obs_LOS = exp(1i * k * d * (0:N-1)' * cos(theta_rad(idx)));
                    rayleigh_coeff = (randn(N,1) + 1i * randn(N,1)) / sqrt(2);  % Normalized
                    h_obs_NLOS = rayleigh_coeff .* exp(1i * k * d * (0:N-1)' * cos(deg2rad(angle_ref)));
                    h_obs_total = g * h_obs_NLOS;
                    power_matrix(iter, idx) = abs(h1' * Theta * h_obs_total)^2;
                end
            end

            power_matrix = power_matrix ./ vecnorm(power_matrix, 2, 2);
            power_matrix_db = 10 * log10(power_matrix + eps);
            power_accum(g_idx, :) = power_accum(g_idx, :) + power_matrix_db(:, ref_angle_index)';

            for idx = 1:length(theta)
                correlation_accum(g_idx, idx) = correlation_accum(g_idx, idx) + corr(reference_signature, power_matrix_db(:, idx));
            end

            [~, peak_idx] = max(correlation_accum(g_idx, :) / run);
            peak_angle_accum(run, g_idx) = theta(peak_idx);
        end
    end

    % === Compute statistics ===
    avg_corr = correlation_accum / num_runs;
    avg_clean_corr = clean_correlation / num_runs;
    avg_power = power_accum / num_runs;
    avg_clean_power = clean_power_accum / num_runs;
    mean_peak_angles = mean(peak_angle_accum, 1);
    abs_errors = abs(mean_peak_angles - angle_ref);
    mae = mean(abs_errors);

    %% Plot 1: Correlation vs Angle
    figure;
    hold on;
    cmap = lines(num_nlos);
    legend_entries = cell(num_nlos+1, 1);
    for g_idx = 1:num_nlos
        plot(theta, avg_corr(g_idx, :), 'Color', cmap(g_idx,:), 'LineWidth', 1.5);
        legend_entries{g_idx} = sprintf('NLOS Gain = %.1f (Peak Corr = %.3f @ %.0f°)', ...
            nlos_linear_values(g_idx), max(avg_corr(g_idx, :)), mean_peak_angles(g_idx));
    end
    [max_clean_corr, clean_peak_idx] = max(avg_clean_corr);
    clean_peak_angle = theta(clean_peak_idx);
    plot(theta, avg_clean_corr, 'k--', 'LineWidth', 2);
    legend_entries{end} = sprintf('Clean LOS (Peak Corr = %.3f @ %.0f°)', max_clean_corr, clean_peak_angle);
    legend(legend_entries, 'Location', 'eastoutside');
    grid on;
    xlabel('Observation Angle (°)');
    ylabel('Correlation with Clean Signature');
    title(sprintf('Avg Correlation vs Angle | Ref = %d°, MAE = %.2f°', angle_ref, mae));
    ylim([-1 1]);
    hold off;

    %% Plot 2: Power vs Beam Index
    figure;
    hold on;
    plot(1:num_iterations, avg_clean_power, 'k-', 'LineWidth', 2, 'DisplayName', 'Clean LOS');
    for g_idx = 1:num_nlos
        plot(1:num_iterations, avg_power(g_idx, :), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('NLOS Gain = %.1f', nlos_linear_values(g_idx)));
    end
    grid on;
    xlabel('Beam Index (Iteration)');
    ylabel('Power at Reference Angle (dB)');
    title(sprintf('Power vs Beam Index at %d° (Avg over %d Runs)', angle_ref, num_runs));
    legend('show', 'Location', 'southwest');
    ylim([-30 -5]);
    hold off;

    %% Plot 3: Absolute Error vs NLOS Gain
    figure;
    plot(nlos_linear_values, abs_errors, 'bo-', 'LineWidth', 2);
    grid on;
    xlabel('NLOS Gain (g)');
    ylabel('Absolute Error (°)');
    title(sprintf('MAE vs NLOS Gain | Ref = %d°', angle_ref));
end

% === Example usage ===
alpha = 30;
N = 16;
num_iterations = 50;
angle_ref = 30;
num_runs = 5;
RISRandomPhases_MultipathVirtualAP_SweepG(N, alpha, num_iterations, angle_ref, num_runs);
