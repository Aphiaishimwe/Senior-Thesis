function mae = RISBeamSteering(N, alpha, beta)
    % Beamsteering to multiple target angles, returning MAE
    lambda = 1;
    d = lambda / 2;
    k = 2 * pi / lambda;

    alpha_rad = deg2rad(alpha);
    h1 = exp(1i * k * d * (0:N-1)' * cos(alpha_rad));

    % Define full observation grid
    theta = linspace(0, 180, 1000);
    theta_rad = deg2rad(theta);

    % Target angles to evaluate (in degrees)
    target_angles = 10:1:170;
    angle_errors = zeros(size(target_angles));

    for t = 1:length(target_angles)
        target_angle = target_angles(t);
        target_idx = find(theta == target_angle);
        if isempty(target_idx)
            [~, target_idx] = min(abs(theta - target_angle));
        end

        powers_at_target = zeros(size(beta));

        for b_idx = 1:length(beta)
            beta_rad = deg2rad(beta(b_idx));
            h2 = exp(1i * k * d * (0:N-1)' * cos(beta_rad));
            theta_opt = angle(h1) - angle(h2);
            Theta = diag(exp(1i * theta_opt));

            array_factor = zeros(size(theta));
            for idx = 1:length(theta)
                h_obs = exp(1i * k * d * (0:N-1)' * cos(theta_rad(idx)));
                array_factor(idx) = abs(h1' * Theta * h_obs)^2;
            end

            array_factor_norm = array_factor / max(array_factor);
            powers_at_target(b_idx) = array_factor_norm(target_idx);
        end

        [~, max_power_idx] = max(powers_at_target);
        best_observation_angle = beta(max_power_idx);
        angle_errors(t) = abs(best_observation_angle - target_angle);

        fprintf('Target = %d°, Best Reflection = %.1f°, Error = %.2f°\n', ...
            target_angle, best_observation_angle, angle_errors(t));
    end

    mae = mean(angle_errors);
    fprintf('✅ Mean Absolute Error (MAE) across %d target angles = %.2f°\n', ...
        length(target_angles), mae);
end

% === Example Run ===
N = 16;
alpha = 45;
beta = linspace(10, 170, 30);  % Adjustable beam steering angles
RISBeamSteering(N, alpha, beta);
