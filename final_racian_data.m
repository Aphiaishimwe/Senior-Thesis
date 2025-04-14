function generate_ris_dataset(N, num_iterations, num_runs, K, ref_angles, nlos_gains)
% Generates RIS power signature dataset for training a neural network.
% Each row = [beam powers (num_iterations), NLOS gain g, reference angle]
%
% Inputs:
% - N: number of RIS elements
% - num_iterations: number of random RIS phase configurations (beams)
% - num_runs: number of runs per (angle, g) pair
% - K: Rician K-factor
% - ref_angles: vector of reference angles in degrees (e.g., 0:10:90)
% - nlos_gains: vector of NLOS gain values (e.g., 0:0.2:4)

    if nargin < 6
        error('Must provide all six inputs: N, num_iterations, num_runs, K, ref_angles, nlos_gains');
    end
    rng(40);

    lambda = 1;
    d = lambda / 2;
    k = 2 * pi / lambda;
    alpha_rad = deg2rad(30);
    h1 = exp(1i * k * d * (0:N-1)' * cos(alpha_rad));

    total_rows = length(ref_angles) * length(nlos_gains) * num_runs;
    all_data = zeros(total_rows, num_iterations + 2); % +2 for g and angle
    row_idx = 1;

    for angle_idx = 1:length(ref_angles)
        angle_ref = ref_angles(angle_idx);
        angle_rad = deg2rad(angle_ref);
        

        for run = 1:num_runs
            % Generate all random phase profiles for this run
            Theta_list = cell(num_iterations, 1);
            for iter = 1:num_iterations
                theta_opt = 2 * pi * rand(N, 1);
                Theta_list{iter} = diag(exp(1i * theta_opt));
            end

            for g_idx = 1:length(nlos_gains)
                g = nlos_gains(g_idx);
                signature = zeros(1, num_iterations);

                for iter = 1:num_iterations
                    Theta = Theta_list{iter};

                    % Observation channel (LOS + g * Rician fading)
                    h_obs_LOS = exp(1i * k * d * (0:N-1)' * cos(angle_rad));
                    LOS = exp(1i * k * d * (0:N-1)' * cos(angle_rad));
                    NLOS = (randn(N,1) + 1i * randn(N,1)) / sqrt(2);
                    rician_fading = sqrt(K / (K + 1)) * LOS + sqrt(1 / (K + 1)) * NLOS;

                    h_obs_total = h_obs_LOS + g * rician_fading;

                    % Beamforming power response
                    signature(iter) = abs(h1' * Theta * h_obs_total)^2;
                end

                % Append row: [signature, g, angle]
                all_data(row_idx, :) = [signature, g, angle_ref];
                row_idx = row_idx + 1;
            end
        end
    end

    % === Assign column headers ===
    beam_headers = arrayfun(@(x) sprintf('beam_%d', x), 1:num_iterations, 'UniformOutput', false);
    column_headers = [beam_headers, {'g', 'angle'}];

    % === Convert to table and save with headers ===
    data_table = array2table(all_data, 'VariableNames', column_headers);
    filename = sprintf('ris_rician_low_50iter_dataset_K%.1f.csv', K);
    writetable(data_table, filename);
    fprintf('Dataset saved to "%s" with %d rows and %d columns (with headers).\n', ...
            filename, size(all_data,1), size(all_data,2));
end
N = 16;
num_iterations = 50;
num_runs = 5;
K = 2;
ref_angles = 0:1:90;
nlos_gains = 0:0.1:4;

generate_ris_dataset(N, num_iterations, num_runs, K, ref_angles, nlos_gains);
