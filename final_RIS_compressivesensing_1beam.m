function RISRandomPhases(N, alpha, target_angle)
    % Inputs:
    %   N: Number of RIS elements
    %   alpha: Incident angle (degrees)
    %   target_angle: Desired observation angle (degrees)
    rng(40)

    lambda = 1;              % Wavelength (normalized)
    d = lambda / 2;          % RIS element spacing
    k = 2 * pi / lambda;     % Wave number
    d_los = 2;

    % Convert angles to radians:
    alpha_rad = deg2rad(alpha);
    
    % 1. Incident wave phase shifts (Tx to RIS channel):
    h1 = exp(1i * k * d * (0:N-1)' * cos(alpha_rad));

    % 2. Generate RANDOM optimal phase shifts (0 to 2π):
    theta_opt = 2 * pi * rand(N, 1); % Uniform random phases
    RIS_config = diag(exp(1i * theta_opt)); % RIS configuration matrix

    % 3. Calculate beam pattern across all observation angles:
    observ_angles = linspace(0, 90, 1000); % Observation angles
    observ_angles_rad = deg2rad(observ_angles);
    array_factor = zeros(size(observ_angles));
    
    for j = 1:length(observ_angles)
        h_obs = exp(1i * k * d_los * (0:N-1)' * cos(observ_angles_rad(j)));
        array_factor(j) = abs(h1' * RIS_config * h_obs)^2; % Signal power
    end
    
    % Normalize the array factor
    array_factor_norm = array_factor / vecnorm(array_factor, 2);
    
     % plot the beampattern
    figure;
    plot(observ_angles, 10*log10(array_factor_norm), 'LineWidth', 1.5);
    hold on;
    % xline(target_angle, '--r', 'LineWidth', 1.5); % Mark target angle
    hold off;
    
    grid on;
    xlabel('Observation Angle (°)');
    ylabel('Normalized Gain (dB)');
    title('Beam Pattern with Random RIS Phases');
    ylim([-60, 0]);
    legend('Random Beam', 'Target Angle', 'Location', 'best');
end

% Example usage:
alpha = 30;          % Incident angle (degrees)
target_angle = 35;   % User angle (degrees)
N = 16;              % RIS elements

% Call the function:
RISRandomPhases(N, alpha, target_angle);