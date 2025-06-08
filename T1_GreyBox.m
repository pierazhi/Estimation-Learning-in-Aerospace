% Piercarlo Fontana - System Identification: Longitudinal Dynamics
clc;
clearvars;
close all;

rng(42);  % For reproducibility

% ------------------------------------------------------------------------
% Define Initial Continuous-Time Model (A, B, C, D)
% State: [u; q; theta], Input: normalized pitching moment, Output: state
% -------------------------------------------------------------------------

Xu = -0.1068;     Xq = 0.1192;
Mu = -5.9755;     Mq = -2.6478;
Xd = -10.1647;    Md = 450.71;

A = [Xu, Xq, -9.81;
     Mu, Mq,  0;
      0,  1,  0];

B = [Xd;
     Md;
      0];

C = eye(3);  % Output all states
D = zeros(3,1);

% ------------------------------------------------------------------------
% Add Measurement Noise and Delays
% -------------------------------------------------------------------------

noise.Enabler = 1;

noise.pos_stand_dev      = noise.Enabler * 0.0011;            % [m]
noise.vel_stand_dev      = noise.Enabler * 0.01;              % [m/s]
noise.attitude_stand_dev = noise.Enabler * deg2rad(0.33);     % [rad]
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(1);        % [rad/s]

delay.position_filter = 1;
delay.attitude_filter = 1;
delay.mixer           = 1;

% ------------------------------------------------------------------------
% Load Controller Parameters and Input Signal (Excitation)
% -------------------------------------------------------------------------

parameters_controller;  % Loads sample_time

load ExcitationM        % Contains [time, normalized pitching moment]
SetPoint = [0, 0];       % Placeholder for reference

t = ExcitationM(:,1);
simulation_time = t(end) - t(1);

% ------------------------------------------------------------------------
% Run Simulink Simulation ('Identification.slx')
% -------------------------------------------------------------------------

load_system('Identification');
set_param('Identification', "FastRestart", "off");

simulation = sim('Identification', 'SrcWorkspace', 'current');

% Cleanup
if exist('slprj', 'dir')
    rmdir('slprj', 's');
end

time        = 0:sample_time:simulation_time;
input       = simulation.Mtot;   % Input: normalized pitching moment
true_output = simulation.q;      % Output: pitch rate [rad/s]

% ------------------------------------------------------------------------
% Plot Input and Output
% -------------------------------------------------------------------------

figure;
plot(ExcitationM(:,1), ExcitationM(:,2), 'b', 'LineWidth', 1.5);
title('Input Signal: Normalized Pitching Moment', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('$\delta_{lon}$', 'Interpreter', 'latex', 'FontSize', 14);
grid on; xlim tight; set(gca, 'FontSize', 12);

figure;
plot(time, true_output, 'r', 'LineWidth', 1.5);
title('System Output: Pitch Rate Response', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('$q(t)$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 14);
grid on; xlim tight; set(gca, 'FontSize', 12);

% ------------------------------------------------------------------------
% Grey-Box System Identification using Frequency-Domain Data
% -------------------------------------------------------------------------

true_sys = ss(A, B, C, D, sample_time);  % True system (discrete-time)
true_tf  = tf(true_sys);

sys_id = iddata(true_output, input, sample_time);  % Create ID data object
sys_ft = fft(sys_id);                              % Frequency domain

% Initial parameter guess
parameters = {
    'Xu',  -0.1;
    'Xq',   0.1;
    'Mu',  -6.0;
    'Mq',  -2.5;
    'Xd', -10.0;
    'Md', 450.0
};

% Grey-box model setup
sys_init   = idgrey(@statespacefun, parameters, 'c');
grey_id = greyest(sys_ft, sys_init);  % Estimate parameters

% Extract and reassign estimated values
param_names = {'X_u', 'X_q', 'M_u', 'M_q', 'X_d', 'M_d'}';
params = getpvec(grey_id);
Xu = params(1); Xq = params(2);
Mu = params(3); Mq = params(4);
Xd = params(5); Md = params(6);

cov_matrix = getcov(grey_id);  % Asymptotic covariance
std_devs = sqrt(diag(cov_matrix));  % Standard deviation (1-sigma)

% True parameter values
true_params = [-0.1068, 0.1192, -5.9755, -2.6478, -10.1647, 450.71];
labels = {'Xu', 'Xq', 'Mu', 'Mq', 'Xd', 'Md'};

% Compute absolute and percentage errors
abs_errors = abs(params' - true_params);
percent_errors = 100 * abs_errors ./ abs(true_params);

% Print header
fprintf('\n--- State-Space Parameter Comparison ---\n\n');
fprintf('%-10s %15s %15s %15s %15s\n', ...
    'Parameter', 'Estimated', 'True', 'Abs Error', '% Error');
fprintf('%s\n', repmat('-', 1, 80));

% Print values in aligned columns
for i = 1:length(true_params)
    fprintf('%-10s %+15.4f %+15.4f %15.4f %14.2f%%\n', ...
        labels{i}, params(i), true_params(i), abs_errors(i), percent_errors(i));
end

% Compute RMSE
rmse = sqrt(mean((params' - true_params).^2));
rmse_percent = sqrt(mean(percent_errors.^2));

fprintf('\n%-25s %.6f\n', 'Root-Mean-Square Error (RMSE):', rmse);
fprintf('%-25s %.2f%%\n\n', 'RMSE (Percentage):', rmse_percent);

g = 9.81;

A_grey_id = grey_id.A;
B_grey_id = grey_id.B;
C_grey_id = [1, 0, 0;
     0, 1, 0;
     0, 0, 1];
D_grey_id = [0;
     0;
     0];

grey_idsys = ss(A_grey_id, B_grey_id, C_grey_id, D_grey_id, sample_time);
grey_idtf  = tf(grey_idsys);

Ts = sample_time;

lambda = eig(A_grey_id);         % Eigenvalues of the system
wn_all = abs(imag(lambda));
wn_all = wn_all(wn_all > 1e-2);   % Filter out near-zero (non-oscillatory modes)

dominant_omega = max(wn_all);     % Dominant mode in rad/s
T_natural = 2*pi / dominant_omega;
 
%% ------------------------------------------------------------------------
% Task 1c - Validation Using a 3211 Sequence
% ------------------------------------------------------------------------

clc; close all;

% Load input excitation (3211 maneuver)
load Excitation3211.mat
t3211 = Excitation3211(:, 1);
u3211 = Excitation3211(:, 2);
simulation_time = t3211(end) - t3211(1);
Ts = sample_time;  % ensure sample_time is loaded beforehand

% ------------------------------------------------------------------------
% TRUE SYSTEM Simulation (ground truth with known parameters)
% ------------------------------------------------------------------------

Xu = -0.1068;  Xq =  0.1192;
Mu = -5.9755;  Mq = -2.6478;
Xd = -10.1647; Md = 450.71;

A_true = [Xu, Xq, -9.81;
          Mu, Mq,  0;
           0,  1,  0];
B_true = [Xd; Md; 0];
C_true = eye(3);
D_true = zeros(3, 1);

% Assign to Simulink base workspace
assignin('base', 'A', A_true);
assignin('base', 'B', B_true);
assignin('base', 'C', C_true);
assignin('base', 'D', D_true);

% Clean and run Simulink
bdclose('all');
load_system('Validation');
set_param('Validation', "FastRestart", "off");
true_validation = sim('Validation', 'SrcWorkspace', 'current');
if exist('slprj', 'dir'), rmdir('slprj', 's'); end

% Extract true system I/O
true_input_validation  = true_validation.Mtot;
true_output_validation = true_validation.q;

% Align time vector
time = t3211(1:length(true_output_validation));

% Create true system object for freq/time domain comparisons
true_sys = ss(A_true, B_true, C_true, D_true, Ts);

% ------------------------------------------------------------------------
% GREY-BOX SYSTEM Simulation (identified model)
% ------------------------------------------------------------------------

A_grey_val = grey_id.A;
B_grey_val = grey_id.B;
C_grey_val = eye(3);
D_grey_val = zeros(3,1);

assignin('base', 'A', A_grey_val);
assignin('base', 'B', B_grey_val);
assignin('base', 'C', C_grey_val);
assignin('base', 'D', D_grey_val);

bdclose('all');
load_system('Validation');
set_param('Validation', "FastRestart", "off");
estimated_validation = sim('Validation', 'SrcWorkspace', 'current');
if exist('slprj', 'dir'), rmdir('slprj', 's'); end

grey_valout = estimated_validation.q;

grey_valsys = ss(A_grey_val, B_grey_val, C_grey_val, D_grey_val, Ts);

% ------------------------------------------------------------------------
% Evaluation Metrics - Grey Box
% ------------------------------------------------------------------------

disp('--- Grey Box Validation Evaluation ---');
VAF_GREY = vaf(true_output_validation, grey_valout);
fit_grey = fit(true_output_validation, grey_valout);
pec_grey = pec(true_output_validation, grey_valout);
fprintf('VAF: %.2f%% | FIT: %.2f%% | PEC: %.2f%%\n', ...
        VAF_GREY, fit_grey, pec_grey);

% ------------------------------------------------------------------------
% Plots
% ------------------------------------------------------------------------

% 3211 Input
figure;
plot(t3211, u3211, 'r', 'LineWidth', 1.5);
ylabel('3211 Excitation Sequence', 'Interpreter', 'latex');
title('3211 Excitation Input');
xlabel('Time [s]');
grid on; xlim tight; set(gca, 'FontSize', 12);

% Normalized Moment Input
figure;
plot(t3211, true_input_validation, 'b', 'LineWidth', 1.5);
ylabel('$\delta_{lon}$ (input)', 'Interpreter', 'latex');
title('Normalized Pitching Moment (3211)');
xlabel('Time [s]');
grid on; xlim tight; set(gca, 'FontSize', 12);

% Output Comparison
figure;
plot(time, true_output_validation, 'k', 'LineWidth', 1.5); hold on;
plot(time, grey_valout, 'r--', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Pitch Rate $q(t)$ [rad/s]', 'Interpreter', 'latex');
title(sprintf('3211 Validation – Grey-Box vs True Output (VAF: %.3f%%)', VAF_GREY));
legend('True Output', 'Grey-Box Model');
grid on; xlim tight; set(gca, 'FontSize', 12);

% Bode Plot
figure;
bodemag(true_sys, 'k-', grey_valsys, 'r--', {0.1, 100});
legend('True System', 'Identified Model', 'Location', 'Best');
title('Bode Magnitude Comparison');
grid on;

% Pole Comparison (clean)
true_poles = pole(true_sys);
id_poles = pole(grey_valsys);

figure;
hold on;
plot(real(true_poles), imag(true_poles), 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'True Poles');
plot(real(id_poles), imag(id_poles), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Identified Poles');
legend('show');
title('Pole Comparison');
xlabel('Real Axis'); ylabel('Imaginary Axis');
grid on; axis equal; set(gca, 'FontSize', 12);

%%
function [A, B, C, D] = statespacefun(Xu, Xq, Mu, Mq, Xd, Md, ~)

g = 9.81;

A = [Xu Xq -g; Mu Mq 0; 0 1 0];

B = [Xd Md 0]';

C = [0 1 0];

D = 0;

end
