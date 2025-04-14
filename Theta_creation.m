% generate_fixed_theta_list.m
% Creates a fixed set of RIS random phase configurations and saves them to disk

clc;
clear;

% === Parameters ===
N = 16;
num_iterations = 100;
theta_filename = sprintf('ThetaList_N%d_Iter%d.mat', N, num_iterations);

% === Generate Theta_list ===
Theta_list = cell(num_iterations, 1);
for iter = 1:num_iterations
    theta_opt = 2 * pi * rand(N, 1);                 % Random phase per element
    Theta_list{iter} = diag(exp(1i * theta_opt));    % Diagonal phase matrix
end

% === Save to file ===
save(theta_filename, 'Theta_list');
fprintf('Saved fixed Theta_list to "%s" with N = %d, Iter = %d\n', ...
        theta_filename, N, num_iterations);
