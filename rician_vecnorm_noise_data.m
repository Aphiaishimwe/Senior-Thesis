function generate_ris_dataset_K_AngleSweep(N, num_iterations, num_runs, K_values, ref_angles)
% Generates RIS dataset using ITU Pedestrian B multipath Rician fading model
% Each row = [beam powers (num_iterations), K-factor, angle_ref]

    if nargin < 5
        error('Must provide N, num_iterations, num_runs, K_values, ref_angles');
    end
    rng(40);  % For reproducibility

    lambda = 1;
    d = lambda / 2;
    k = 2 * pi / lambda;
    SNR_dB = 5;  % Desired SNR for added noise

    total_rows = length(K_values) * length(ref_angles) * num_runs;
    all_data = zeros(total_rows, num_iterations + 2);  % +2 for K and angle
    row_idx = 1;

    % Load RIS beam pattern list
    theta_filename = sprintf('ThetaList_N%d_Iter%d.mat', N, num_iterations);
    if isfile(theta_filename)
        load(theta_filename, 'Theta_list');
        fprintf("✅ Loaded RIS Theta_list from %s\n", theta_filename);
    else
        error("❌ Theta_list file not found. Please generate it using generate_fixed_theta_list.m");
    end

    for angle_idx = 1:length(ref_angles)
        angle_ref = ref_angles(angle_idx);
        angle_rad = deg2rad(angle_ref);

        % h1: AP → RIS (fixed incident angle)
        alpha_rad = deg2rad(30);
        h1 = exp(1i * k * d * (0:N-1)' * cos(alpha_rad));

        % h_los: LOS channel from RIS to user at given angle
        h_los = exp(1i * k * d * (0:N-1)' * cos(angle_rad));

        for k_idx = 1:length(K_values)
            K = K_values(k_idx);

            % === ITU Pedestrian B Rician Fading ===
            fs = 3.84e6;
            pathDelays = [0 200 800 1200 2300 3700]*1e-6;
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
                'Seed', 73 + k_idx + angle_idx, ...
                'PathGainsOutputPort', true ...
            );

            ricianChan(1 + 1i);  % Prime the channel with complex input

            for run = 1:num_runs
                signature = zeros(1, num_iterations);

                for iter = 1:num_iterations
                    Theta = Theta_list{iter};

                    % Simulate fading channel per RIS element
                    h_obs_total = zeros(N, 1);
                    for n = 1:N
                        [faded_sig, ~] = ricianChan(h_los(n));
                        h_obs_total(n) = faded_sig;
                    end

                    % Compute power at user after beamforming
                    signature(iter) = abs(h1' * Theta * h_obs_total)^2;
                end

                % Normalize and add AWGN to the signature
                signature = signature ./ norm(signature, 2);  % Normalize power vector
                signature = awgn(signature, SNR_dB, 'measured');  % Add noise at SNR

                % Store in dataset: [signature, K, angle]
                all_data(row_idx, :) = [signature, K, angle_ref];
                row_idx = row_idx + 1;
            end
        end
    end

    % Save to CSV
    beam_headers = arrayfun(@(x) sprintf('beam_%d', x), 1:num_iterations, 'UniformOutput', false);
    column_headers = [beam_headers, {'K', 'angle'}];
    data_table = array2table(all_data, 'VariableNames', column_headers);

    filename = sprintf('ris_less1_test_dataset_KAngleSweep_ITUrician_SNR%ddB_N%d_iter%d.csv', SNR_dB, N, num_iterations);
    writetable(data_table, filename);
    fprintf('✅ Dataset saved to "%s" with %d rows and %d columns.\n', ...
        filename, size(data_table, 1), size(data_table, 2));
end

N = 16;
num_iterations = 20;
num_runs = 1;
K_values = 0.1:0.1:1.2;
ref_angles = 0:1:90;

generate_ris_dataset_K_AngleSweep(N, num_iterations, num_runs, K_values, ref_angles);
