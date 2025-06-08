% T3.1 Compare the Grey
clearvars; close all; clc

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
opt = greyestOptions();
opt.SearchMethod = 'lm';
opt.SearchOptions.MaxIterations = 100;
opt.SearchOptions.Tolerance = 1e-6;
opt.Display = 'on';

sys_init   = idgrey(@statespacefun, parameters, 'c');
grey_id = greyest(sys_ft, sys_init, opt);  % Estimate parameters

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

% ------------------------------------------------------------------------
% Black-Box System Identification using Frequency-Domain Data
% -------------------------------------------------------------------------

n = 3;    % order of system
p = 750;   % past window size
f = 300;   % future window size

[Aest, Best, Cest, Dest, Kest, Sest] = estimate_ABCDK(input', true_output', f, p, n);
black_sys = ss(Aest, Best, Cest, Dest, sample_time);
black_idout = lsim(black_sys, input, time);

guess = [0 0 0 -2 -10 300]';
%guess = [Aest(1,1) Aest(1,2) Aest(2,1) Aest(2,2) Best(1) Best(2)]';
freq_vec = 2*pi*logspace(-1, 2, length(time));  

% Structure black-box into grey-box

frd_data = frd(black_sys, freq_vec);
sys_black = idgrey(@longitdynRot, guess, 'c');

opt = greyestOptions();
opt.SearchMethod = 'lm';
opt.SearchOptions.MaxIterations = 100;
opt.SearchOptions.Tolerance = 1e-6;
opt.Display = 'on';

% Estimate structured model
structuredModel = greyest(frd_data, sys_black, opt);
structured_tf = tf(structuredModel);

% Extract matrices and parameters
A_blackstruct =  structuredModel.A;
B_blackstruct  = structuredModel.B;
C_blackstruct  = structuredModel.C;
D_blackstruct  = structuredModel.D;
est_params = getpvec(structuredModel);

structured_sys = ss(A_blackstruct, B_blackstruct, C_blackstruct, D_blackstruct, Ts);

% True parameter values
true_params = [-0.1068, 0.1192, -5.9755, -2.6478, -10.1647, 450.71];
labels = {'Xu', 'Xq', 'Mu', 'Mq', 'Xd', 'Md'};

% Compute absolute and percentage errors
abs_errors = abs(est_params' - true_params);
percent_errors = 100 * abs_errors ./ abs(true_params);

% Print header
fprintf('\n--- State-Space Parameter Comparison ---\n\n');
fprintf('%-10s %15s %15s %15s %15s\n', ...
    'Parameter', 'Estimated', 'True', 'Abs Error', '% Error');
fprintf('%s\n', repmat('-', 1, 80));

% Print values in aligned columns
for i = 1:length(true_params)
    fprintf('%-10s %+15.4f %+15.4f %15.4f %14.2f%%\n', ...
        labels{i}, est_params(i), true_params(i), abs_errors(i), percent_errors(i));
end

% Compute RMSE
rmse = sqrt(mean((est_params' - true_params).^2));
rmse_percent = sqrt(mean(percent_errors.^2));

fprintf('\n%-25s %.6f\n', 'Root-Mean-Square Error (RMSE):', rmse);
fprintf('%-25s %.2f%%\n\n', 'RMSE (Percentage):', rmse_percent);

% Plots

figure;
pzplot(true_sys, 'k', grey_idsys, structuredModel);
legend('True', 'Grey-box', 'Structured Black-box');
title('Pole-Zero Map Comparison');
grid on;

w = logspace(-1, 2, 500);  % Frequency in rad/s
H_grey = squeeze(freqresp(grey_idsys(2,:), w));       % Returns complex FRF
H_black = squeeze(freqresp(structured_sys(1,:), w));
H_true = squeeze(freqresp(true_sys(2,:), w));

figure;
subplot(2,1,1);
semilogx(w, 20*log10(abs(H_true)), 'k',...
         w, 20*log10(abs(H_grey)), 'b--', ...
         w, 20*log10(abs(H_black)), 'r-');
ylabel('Magnitude [dB]');
legend('True', 'Grey-box', 'Structured Black-box');
grid on;

subplot(2,1,2);
semilogx(w, angle(H_true)*180/pi, 'k', ...
         w, angle(H_grey)*180/pi, 'b--', ...
         w, angle(H_black)*180/pi, 'r-');
ylabel('Phase [deg]');
xlabel('Frequency [rad/s]');
grid on;

%% T3.2 GREY Monte Carlo
clc;
clearvars;
close all;

N = 100;
params_grey = zeros(N, 6);  % Preallocate
wait_h = waitbar(0, 'Running Monte Carlo for Grey-Box Identification...');
    
for i = 1:N  % For reproducibility

    rng(42);
    waitbar(i/N, wait_h, sprintf('Run %d of %d', i, N));
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

    seed.x = i;
    seed.vx = i + 1;
    seed.theta = i + 2;
    seed.q = i + 3;
    
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
    
    load_system('MonteCarlo');
    set_param('MonteCarlo', "FastRestart", "off");
    
    simulation = sim('MonteCarlo', 'SrcWorkspace', 'current');
    
    % Cleanup
    if exist('slprj', 'dir')
        rmdir('slprj', 's');
    end
    
    time        = 0:sample_time:simulation_time;
    input       = simulation.Mtot;   % Input: normalized pitching moment
    true_output = simulation.q;      % Output: pitch rate [rad/s]
    
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
    opt = greyestOptions();
    opt.SearchMethod = 'lm';
    opt.SearchOptions.MaxIterations = 100;
    opt.SearchOptions.Tolerance = 1e-6;
    opt.Display = 'off';
    
    sys_init   = idgrey(@statespacefun, parameters, 'c');
    grey_id = greyest(sys_ft, sys_init, opt);  % Estimate parameters
    
    % Extract and reassign estimated values
    param_names = {'X_u', 'X_q', 'M_u', 'M_q', 'X_d', 'M_d'}';
    params = getpvec(grey_id);
    Xu = params(1); Xq = params(2);
    Mu = params(3); Mq = params(4);
    Xd = params(5); Md = params(6);
    
    cov_matrix = getcov(grey_id);  % Asymptotic covariance
    std_devs = sqrt(diag(cov_matrix));  % Standard deviation (1-sigma)

    params_grey(i, :) = getpvec(grey_id);  % size [N × 6]

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

end

%% Compute stats
mean_params = mean(params_grey);
std_params = std(params_grey);
bias = mean_params - true_params;
rmse = sqrt(mean((params_grey - true_params).^2));

% Print table
fprintf('\n--- Grey-Box Parameter Statistics (N = %d) ---\n\n', size(params_grey,1));
fprintf('%-10s %12s %12s %12s %12s\n', 'Param', 'True', 'Mean', 'Bias', 'RMSE');
fprintf('%s\n', repmat('-', 1, 62));
for j = 1:6
    fprintf('%-10s %12.4f %12.4f %12.4f %12.4f\n', ...
        labels{j}, true_params(j), mean_params(j), bias(j), rmse(j));
end

figure;
for j = 1:6
    subplot(2,3,j);
    histogram(params_grey(:,j), 30);
    xlabel(labels{j});
    ylabel('Count');
    title(sprintf('%s Distribution', labels{j}));
    grid on;
end
sgtitle('Grey-Box Parameter Histograms (Monte Carlo)');

w = logspace(-1, 2, 1000);  % rad/s
Ts = sample_time;
H_grey_all = zeros(length(w), N);  % one output channel: pitch rate

for i = 1:N
    p = params_grey(i,:);
    A = [p(1), p(2), -9.81;
         p(3), p(4),  0;
         0,    1,     0];
    B = [p(5); p(6); 0];
    Cq = [0 1 0];  % pitch rate output
    sys_i = ss(A, B, Cq, 0, Ts);
    H_grey_all(:,i) = squeeze(freqresp(sys_i, w));
end

mag = abs(H_grey_all);
mean_mag = mean(mag, 2);
std_mag = std(mag, 0, 2);

figure;
semilogx(w, 20*log10(mean_mag), 'b', 'LineWidth', 1.5); hold on;
semilogx(w, 20*log10(mean_mag + std_mag), 'b--');
semilogx(w, 20*log10(mean_mag - std_mag), 'b--');
xlabel('Frequency [rad/s]');
ylabel('Magnitude [dB]');
title('Grey-Box FRF Dispersion (Pitch Rate)');
grid on;
legend('Mean', '+1σ', '-1σ');
mag = abs(H_grey_all);
mean_mag = mean(mag, 2);
std_mag = std(mag, 0, 2);

phase = angle(H_grey_all);
mean_phase = mean(phase, 2);
std_phase = std(phase, 0, 2);

figure;
semilogx(w, mean_phase*180/pi, 'r', 'LineWidth', 1.5); hold on;
semilogx(w, (mean_phase + std_phase)*180/pi, 'r--');
semilogx(w, (mean_phase - std_phase)*180/pi, 'r--');
xlabel('Frequency [rad/s]');
ylabel('Phase [deg]');
title('Grey-Box FRF Phase Dispersion');
grid on;
legend('Mean', '+1σ', '-1σ');

for i = 1:N
    % Reconstruct system matrix A for this run
    p = params_grey(i,:);
    A = [p(1), p(2), -9.81;
         p(3), p(4),  0;
         0,    1,     0];
    poles_grey(i,:) = eig(A);  % 3 eigenvalues per system
end

figure;
subplot(1,2,1);
histogram(real(poles_grey(:, 1)), 100);
xlabel('Re(λ)');
ylabel('Count');
title('Grey-Box: First Pole Real Part Distribution');
grid on;

subplot(1,2,2);
histogram(imag(poles_grey(:, 1)), 100);
xlabel('Im(λ)');
ylabel('Count');
title('Grey-Box: Pole Imaginary Part Distribution');
grid on;

figure;
plot(real(poles_grey), imag(poles_grey), 'x');
xlabel('Re(λ)');
ylabel('Im(λ)');
title('Pole Locations (Complex Plane)');
grid on;

%% T3.2 BLACK Monte Carlo
clc;
clearvars;
close all;

N = 100;
params_black = zeros(N, 6);  % Preallocate
wait_h = waitbar(0, 'Running Monte Carlo for Grey-Box Identification...');
    
for i = 1:N  % For reproducibility

    rng(42);
    waitbar(i/N, wait_h, sprintf('Run %d of %d', i, N));
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

    seed.x = i;
    seed.vx = i + 1;
    seed.theta = i + 2;
    seed.q = i + 3;
    
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
    
    load_system('MonteCarlo');
    set_param('MonteCarlo', "FastRestart", "off");
    
    simulation = sim('MonteCarlo', 'SrcWorkspace', 'current');
    
    % Cleanup
    if exist('slprj', 'dir')
        rmdir('slprj', 's');
    end
    
    time        = 0:sample_time:simulation_time;
    input       = simulation.Mtot;   % Input: normalized pitching moment
    true_output = simulation.q;      % Output: pitch rate [rad/s]
    
    n = 3;    % order of system
    p = 750;   % past window size
    f = 300;   % future window size
    
    [Aest, Best, Cest, Dest, Kest, Sest] = estimate_ABCDK(input', true_output', f, p, n);
    black_sys = ss(Aest, Best, Cest, Dest, sample_time);
    black_idout = lsim(black_sys, input, time);
    
    guess = [0 0 0 -2 -10 300]';
    freq_vec = 2*pi*logspace(-1, 2, length(time));  
    
    % Structure black-box into grey-box
    
    frd_data = frd(black_sys, freq_vec);
    sys_black = idgrey(@longitdynRot, guess, 'c');
    
    opt = greyestOptions();
    opt.SearchMethod = 'lm';
    opt.SearchOptions.MaxIterations = 100;
    opt.SearchOptions.Tolerance = 1e-6;
    opt.Display = 'off';
    
    % Estimate structured model
    structuredModel = greyest(frd_data, sys_black, opt);
    
    % Extract matrices and parameters
    A_blackstruct =  structuredModel.A;
    B_blackstruct  = structuredModel.B;
    C_blackstruct  = structuredModel.C;
    D_blackstruct  = structuredModel.D;
    est_params = getpvec(structuredModel);

    params_black(i, :) = getpvec(structuredModel);  % size [N × 6]

    % True parameter values
    true_params = [-0.1068, 0.1192, -5.9755, -2.6478, -10.1647, 450.71];
    labels = {'Xu', 'Xq', 'Mu', 'Mq', 'Xd', 'Md'};
    
    % Compute absolute and percentage errors
    abs_errors = abs(est_params' - true_params);
    percent_errors = 100 * abs_errors ./ abs(true_params);
    
    % Print header
    fprintf('\n--- State-Space Parameter Comparison ---\n\n');
    fprintf('%-10s %15s %15s %15s %15s\n', ...
        'Parameter', 'Estimated', 'True', 'Abs Error', '% Error');
    fprintf('%s\n', repmat('-', 1, 80));
    
    % Print values in aligned columns
    for i = 1:length(true_params)
        fprintf('%-10s %+15.4f %+15.4f %15.4f %14.2f%%\n', ...
            labels{i}, est_params(i), true_params(i), abs_errors(i), percent_errors(i));
    end
    
    % Compute RMSE
    rmse = sqrt(mean((est_params' - true_params).^2));
    rmse_percent = sqrt(mean(percent_errors.^2));
    
    fprintf('\n%-25s %.6f\n', 'Root-Mean-Square Error (RMSE):', rmse);
    fprintf('%-25s %.2f%%\n\n', 'RMSE (Percentage):', rmse_percent);

    
    g = 9.81;
    
    A_black_id = structuredModel.A;
    B_black_id = structuredModel.B;
    C_black_id = [1, 0, 0;
         0, 1, 0;
         0, 0, 1];
    D_black_id = [0;
         0;
         0];

end

%% Compute stats
mean_params = mean(params_black);
std_params = std(params_black);
bias = mean_params - true_params;
rmse = sqrt(mean((params_black - true_params).^2));

% Print table
fprintf('\n--- Black-Box Parameter Statistics (N = %d) ---\n\n', size(params_black,1));
fprintf('%-10s %12s %12s %12s %12s\n', 'Param', 'True', 'Mean', 'Bias', 'RMSE');
fprintf('%s\n', repmat('-', 1, 62));
for j = 1:6
    fprintf('%-10s %12.4f %12.4f %12.4f %12.4f\n', ...
        labels{j}, true_params(j), mean_params(j), bias(j), rmse(j));
end

figure;
for j = 1:6
    subplot(2,3,j);
    histogram(params_black(:,j), 30);
    xlabel(labels{j});
    ylabel('Count');
    title(sprintf('%s Distribution', labels{j}));
    grid on;
end
sgtitle('Black-Box Parameter Histograms (Monte Carlo)');

w = logspace(-1, 2, 1000);  % rad/s
Ts = sample_time;
H_grey_all = zeros(length(w), N);  % one output channel: pitch rate

for i = 1:N
    p = params_black(i,:);
    A = [p(1), p(2), -9.81;
         p(3), p(4),  0;
         0,    1,     0];
    B = [p(5); p(6); 0];
    Cq = [0 1 0];  % pitch rate output
    sys_i = ss(A, B, Cq, 0, Ts);
    H_black_all(:,i) = squeeze(freqresp(sys_i, w));
end

mag = abs(H_black_all);
mean_mag = mean(mag, 2);
std_mag = std(mag, 0, 2);

figure;
semilogx(w, 20*log10(mean_mag), 'b', 'LineWidth', 1.5); hold on;
semilogx(w, 20*log10(mean_mag + std_mag), 'b--');
semilogx(w, 20*log10(mean_mag - std_mag), 'b--');
xlabel('Frequency [rad/s]');
ylabel('Magnitude [dB]');
title('Black-Box FRF Dispersion (Pitch Rate)');
grid on;
legend('Mean', '+1σ', '-1σ');
mag = abs(H_black_all);
mean_mag = mean(mag, 2);
std_mag = std(mag, 0, 2);

phase = angle(H_black_all);
mean_phase = mean(phase, 2);
std_phase = std(phase, 0, 2);

figure;
semilogx(w, mean_phase*180/pi, 'r', 'LineWidth', 1.5); hold on;
semilogx(w, (mean_phase + std_phase)*180/pi, 'r--');
semilogx(w, (mean_phase - std_phase)*180/pi, 'r--');
xlabel('Frequency [rad/s]');
ylabel('Phase [deg]');
title('Black-Box FRF Phase Dispersion');
grid on;
legend('Mean', '+1σ', '-1σ');

for i = 1:N
    % Reconstruct system matrix A for this run
    p = params_black(i,:);
    A = [p(1), p(2), -9.81;
         p(3), p(4),  0;
         0,    1,     0];
    poles_black(i,:) = eig(A);  % 3 eigenvalues per system
end

figure;
subplot(1,2,1);
histogram(real(poles_black(:, 1)), 100);
xlabel('Re(λ)');
ylabel('Count');
title('Black-Box: First Pole Real Part Distribution');
grid on;

subplot(1,2,2);
histogram(imag(poles_black(:, 1)), 100);
xlabel('Im(λ)');
ylabel('Count');
title('Black-Box: Pole Imaginary Part Distribution');
grid on;

figure;
plot(real(poles_black), imag(poles_black), 'x');
xlabel('Re(λ)');
ylabel('Im(λ)');
title('Pole Locations (Complex Plane)');
grid on;

%% Functions

function [A, B, C, D] = longitdynRot(THETA, Ts,  ~)
    
    Xu = THETA(1);
    Xq = THETA(2);
    Mu = THETA(3);
    Mq = THETA(4);
    Xd = THETA(5);
    Md = THETA(6);
    
    A = [ Xu, Xq, -9.81;
          Mu, Mq,   0;  
          0,   1,   0  ]; 
    
    B = [Xd; Md; 0];
      
    C = [0, 1, 0]; 
    
    D = 0;
    
    
end

function [A, B, C, D] = statespacefun(Xu, Xq, Mu, Mq, Xd, Md, ~)

g = 9.81;

A = [Xu Xq -g; Mu Mq 0; 0 1 0];

B = [Xd Md 0]';

C = [0 1 0];

D = 0;

end

function [Aest, Best, Cest, Dest, Kest, Sest] = estimate_ABCDK(u,y,f,p,n)

r = size(u,1);  % Number of inputs
l = size(y,1);  % Number of outputs
N = size(y,2);

m = r + l;
z = [u; y];
Z = zeros(p*m, N-p);

for i = 1:p
    Z((i-1)*m+1:i*m, :) = z(:, i:N+i-p-1);
end

Y = y(:, p+1:N);
U = u(:, p+1:N);

ZU = [Z; U];

Q = Y * pinv(ZU);

Qp_minus = zeros(f*l,p*m);

for i = 1:f
    Qp_minus((i-1)*l+1:i*l,p*m-(p-i+1)*m+1:p*m) = Q(:,1:(p-i+1)*m);
end

[~,S,V] = svd(Qp_minus*Z(1:p*m,:),'econ');
Sest = diag(S)';
X = diag(sqrt(Sest))*V';
X = X(1:n, :);

%S = estimate_order(u,y,f,p);

u = u(:, p+1:p+size(X,2));
y = y(:, p+1:p+size(X,2));

CD = y(:, 1:end-1)*pinv(vertcat(X(:, 1:end-1), u(:, 1:end-1)));

e = y - CD*[X;u];
FROBE_ABK = vertcat(vertcat(X(:, 1:end-1), u(:, 1:end-1)), e(:, 1:end-1));
ABK = X(:, 2:end)*pinv(FROBE_ABK);

Cest = CD(:, 1:n);
Dest = CD(:, n+1:n+r);
Best = ABK(:, n+1:n+r);
Kest = ABK(1:n, n+r+1:n+r+1);
Aest = ABK(:, 1:n);

end

function S = estimate_order(u,y,f,p)

r = size(u,1);  % Number of inputs
l = size(y,1);  % Number of outputs
N = size(y,2);

m = r + l;
z = [u; y];
Z = zeros(p*m, N-p);

for i = 1:p
    Z((i-1)*m+1:i*m, :) = z(:, i:N+i-p-1);
end

Y = y(:, p+1:N);
U = u(:, p+1:N);

ZU = [Z; U];

Q = Y * pinv(ZU);

Qp_minus = zeros(f*l,p*m);

for i = 1:f
    Qp_minus((i-1)*l+1:i*l,p*m-(p-i+1)*m+1:p*m) = Q(:,1:(p-i+1)*m);
end

[~,S,~] = svd(Qp_minus*Z(1:p*m,:),'econ');
S = diag(S)';

%Plot singular values on a logarithmic scale
figure;
hold on
scatter(1:length(S), S, 'filled');
set(gca, 'YScale', 'log'); % Use log scale for better visibility of gaps
xlabel('Index');
ylabel('Singular Value (log scale)');
title('Singular Value Spectrum');

end