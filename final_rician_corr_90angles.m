function RISRandomPhases_MultipathVirtualAP_SweepG(N, alpha, num_iterations, num_runs)
    lambda = 1;
    d = lambda / 2;
    k = 2 * pi / lambda;

    theta = linspace(0, 90, 100);         % Angular scan grid
    theta_rad = deg2rad(theta);
    angle_sweep = 0:1:90;                 % Reference angles
    num_angles = length(angle_sweep);

    g_values = 1:1:10;                    % NLOS gain values
    num_g = length(g_values);
    mean_abs_errors = zeros(1, num_g);

    % === Load RIS Theta list ===
    theta_filename = sprintf('ThetaList_N%d_Iter%d.mat', N, num_iterations);
    if isfile(theta_filename)
        load(theta_filename, 'Theta_list');
        fprintf(" Loaded RIS Theta_list from %s\n", theta_filename);
    else
        error("Theta_list file not found. Please generate it using generate_fixed_theta_list.m");
    end

    for g_idx = 1:num_g
        g = g_values(g_idx);
        fprintf("Simulating for NLOS gain g = %.1f\n", g);

        abs_errors_all_angles = zeros(1, num_angles);

        for a_idx = 1:num_angles
            angle_ref = angle_sweep(a_idx);
            angle_ref_rad = deg2rad(angle_ref);
            alpha_rad = deg2rad(alpha);
            h1 = exp(1i * k * d * (0:N-1)' * cos(alpha_rad));

            % === Clean signature at reference angle ===
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

            % Use clean signature at correct ref angle
            [~, ref_idx] = min(abs(theta - angle_ref));
            reference_signature = power_matrix_clean_db(:, ref_idx);

            % === Rayleigh + NLOS fading ===
            peak_angles = zeros(num_runs, 1);
            for run = 1:num_runs
                power_matrix = zeros(num_iterations, length(theta));
                for iter = 1:num_iterations
                    Theta = Theta_list{iter};
                    for idx = 1:length(theta)
                        h_obs_LOS = exp(1i * k * d * (0:N-1)' * cos(theta_rad(idx)));
                        rayleigh_coeff = (randn(N,1) + 1i * randn(N,1)) / sqrt(2);
                        h_obs_NLOS = rayleigh_coeff .* exp(1i * k * d * (0:N-1)' * cos(angle_ref_rad));
                        h_obs_total = h_obs_LOS + g * h_obs_NLOS;
                        power_matrix(iter, idx) = abs(h1' * Theta * h_obs_total)^2;
                    end
                end

                power_matrix = power_matrix ./ vecnorm(power_matrix, 2, 2);
                power_matrix_db = 10 * log10(power_matrix + eps);

                correlations = zeros(1, length(theta));
                for idx = 1:length(theta)
                    correlations(idx) = corr(reference_signature, power_matrix_db(:, idx));
                end

                [~, peak_idx] = max(correlations);
                peak_angles(run) = theta(peak_idx);
            end

            abs_errors_all_angles(a_idx) = mean(abs(peak_angles - angle_ref));
        end

        mean_abs_errors(g_idx) = mean(abs_errors_all_angles);
    end

    %% === Plot MAE vs g ===
    figure;
    plot(g_values, mean_abs_errors, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');
    grid on;
    xlabel('NLOS Gain (g)');
    ylabel('Mean Absolute Error (°)');
    title('MAE vs NLOS Gain (Rayleigh Fading, Fixed Theta)');

    %% === Save MAE to CSV ===
    T = table(g_values(:), mean_abs_errors(:), ...
        'VariableNames', {'NLOS_Gain', 'Mean_Absolute_Error'});
    writetable(T, 'MAE_vs_G_Correlation_FixedTheta.csv');
    fprintf(" MAE results saved to MAE_vs_G_Correlation_FixedTheta.csv\n");
end

% === Example usage ===
N = 16;
alpha = 30;
num_iterations = 50;
num_runs = 5;

RISRandomPhases_MultipathVirtualAP_SweepG(N, alpha, num_iterations, num_runs);
