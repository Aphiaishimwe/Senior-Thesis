function RISRandomPhases_MultipathVirtualAP_SweepG(N, alpha, num_iterations, num_runs)
    lambda = 1;
    d = lambda / 2;
    k = 2 * pi / lambda;
    rng(50)


    theta = linspace(0, 90, 100);
    theta_rad = deg2rad(theta);
    angle_sweep = 0:1:90;
    num_angles = length(angle_sweep);

    g_values = 1:1:10;  % NLOS gain values
    num_g = length(g_values);
    mean_abs_errors = zeros(1, num_g);

    for g_idx = 1:num_g
        g = g_values(g_idx);
        fprintf("Simulating for NLOS gain g = %.1f\n", g);
        abs_errors_all_angles = zeros(1, num_angles);

        for a_idx = 1:num_angles
            angle_ref = angle_sweep(a_idx);
            alpha_rad = deg2rad(alpha);
            angle_ref_rad = deg2rad(angle_ref);
            h1 = exp(1i * k * d * (0:N-1)' * cos(alpha_rad));

            % === Generate Clean Signature ===
            Theta_list = cell(num_iterations, 1);
            for iter = 1:num_iterations
                theta_opt = 2 * pi * rand(N, 1);
                Theta_list{iter} = diag(exp(1i * theta_opt));
            end

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

            [~, ref_idx] = min(abs(theta - angle_ref));
            reference_signature = power_matrix_clean_db(:, ref_idx);

            % === Rayleigh fading + MAE estimation ===
            peak_angles = zeros(num_runs, 1);
            for run = 1:num_runs
                power_matrix = zeros(num_iterations, length(theta));
                for iter = 1:num_iterations
                    Theta = Theta_list{iter};
                    for idx = 1:length(theta)
                        % h_obs_LOS = exp(1i * k * d * (0:N-1)' * cos(theta_rad(idx)));
                        rayleigh_coeff = (randn(N,1) + 1i * randn(N,1)) / sqrt(2); % normalized Rayleigh
                        h_obs_NLOS = rayleigh_coeff .* exp(1i * k * d * (0:N-1)' * cos(angle_ref_rad));
                        h_obs_total = g * h_obs_NLOS;
                        power_matrix(iter, idx) = abs(h1' * Theta * h_obs_total)^2;
                    end
                end

                power_matrix = power_matrix ./ vecnorm(power_matrix, 2, 2);
                power_matrix_db = 10 * log10(power_matrix + eps);

                % Correlate each angular profile with reference signature
                correlations = zeros(1, length(theta));
                for idx = 1:length(theta)
                    correlations(idx) = corr(reference_signature, power_matrix_db(:, idx));
                end

                [~, peak_idx] = max(correlations);
                peak_angles(run) = theta(peak_idx);
            end

            % MAE for current reference angle
            abs_errors_all_angles(a_idx) = mean(abs(peak_angles - angle_ref));
        end

        % Average MAE over all reference angles for current g
        mean_abs_errors(g_idx) = mean(abs_errors_all_angles);
    end

    %% === Plot MAE vs g ===
    figure;
    plot(g_values, mean_abs_errors, 'mo-', 'LineWidth', 2, 'MarkerFaceColor', 'm');
    grid on;
    xlabel('NLOS Gain (g)');
    ylabel('Mean Absolute Error (°)');
    title('MAE vs NLOS Gain (Rayleigh Fading)');

    %% === Save MAE to CSV ===
    results_table = table(g_values(:), mean_abs_errors(:), ...
        'VariableNames', {'NLOS_Gain', 'Mean_Absolute_Error'});
    csv_filename = 'MAE_vs_G.csv';
    writetable(results_table, csv_filename);
    fprintf("MAE results saved to: %s\n", csv_filename);
end

% === Example Usage ===
alpha = 30;
N = 16;
num_iterations = 50;
num_runs = 1;
RISRandomPhases_MultipathVirtualAP_SweepG(N, alpha, num_iterations, num_runs);
