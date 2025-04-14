function RIS_GenerateRayleighNNData_CSV(N, alpha, num_iterations, num_runs, g_values, save_filename)
    % Constants
    lambda = 1;
    d = lambda / 2;
    k = 2 * pi / lambda;
    rng(50);

    % Sweep settings
    angle_sweep = 0:1:90;                 % True AoA labels
    num_angles = length(angle_sweep);
    num_g = length(g_values);
    total_samples = num_g * num_angles * num_runs;
    num_features = num_iterations;        % Each feature is power from a beam config

    % Preallocate data
    X_all = zeros(total_samples, num_features);
    g_all = zeros(total_samples, 1);
    angle_labels = zeros(total_samples, 1);

    sample_idx = 1;

    % === Load fixed RIS phase configurations ===
    theta_filename = sprintf('ThetaList_N%d_Iter%d.mat', N, num_iterations);
    if isfile(theta_filename)
        load(theta_filename, 'Theta_list');
        fprintf("Loaded RIS Theta_list from %s\n", theta_filename);
    else
        error("Theta_list file not found. Generate it using generate_fixed_theta_list.m");
    end

    % === Loop over NLOS gains, angles, and runs ===
    for g_idx = 1:num_g
        g = g_values(g_idx);
        fprintf("Simulating for NLOS gain g = %.1f\n", g);

        for a_idx = 1:num_angles
            angle_ref = angle_sweep(a_idx);
            angle_ref_rad = deg2rad(angle_ref);
            alpha_rad = deg2rad(alpha);
            h1 = exp(1i * k * d * (0:N-1)' * cos(alpha_rad));

            for run = 1:num_runs
                power_vector = zeros(1, num_iterations);

                for iter = 1:num_iterations
                    Theta = Theta_list{iter};

                    % Observation at the reference angle
                    h_obs_LOS = exp(1i * k * d * (0:N-1)' * cos(angle_ref_rad));
                    rayleigh_coeff = (randn(N,1) + 1i * randn(N,1)) / sqrt(2);
                    h_obs_NLOS = rayleigh_coeff .* exp(1i * k * d * (0:N-1)' * cos(angle_ref_rad));
                    h_obs_total = g * h_obs_NLOS;

                    power = abs(h1' * Theta * h_obs_total)^2;
                    power_vector(iter) = power;
                end

                % Normalize and convert to dB
                power_vector = power_vector ./ norm(power_vector, 2);
                % power_vector_db = 10 * log10(power_vector + eps);

                % Store sample
                X_all(sample_idx, :) = power_vector;
                g_all(sample_idx) = g;
                angle_labels(sample_idx) = angle_ref;

                sample_idx = sample_idx + 1;
            end
        end
    end

    % === Combine into a table and save to CSV ===
    feature_names = arrayfun(@(i) sprintf('Iter_%d', i), 1:num_iterations, 'UniformOutput', false);
    column_names = [feature_names, {'NLOS_Gain', 'True_Angle'}];
    data_matrix = [X_all, g_all, angle_labels];
    data_table = array2table(data_matrix, 'VariableNames', column_names);

    % Optional: sort for easier viewing
    data_table = sortrows(data_table, {'True_Angle', 'NLOS_Gain'});

    writetable(data_table, save_filename);
    fprintf("DNN training data saved to %s\n", save_filename);
end

% === Example Usage ===
N = 16;
alpha = 30;
num_iterations = 50;
num_runs = 5;
g_values = 1:1:10;
save_filename = 'DNN_RayleighTrainingData_FixedTheta.csv';

RIS_GenerateRayleighNNData_CSV(N, alpha, num_iterations, num_runs, g_values, save_filename);
