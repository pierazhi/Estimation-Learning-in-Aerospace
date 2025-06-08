% Piercarlo Fontana - T2
close all; clear; clc;

%% First-Order Thermal System Identification using PBSID + Structuring
close all; clear; clc;

rng('default')

% TRUE PARAMETERS
tau_true = 1;     % [s] thermal time constant
K_true   = 3.0;     % [°C/unit input]

A_th = -1 / tau_true;
B_th = K_true / tau_true;
C_th = 1;
D_th = 0;

Ts = 0.01;                     % sampling time
N  = 1000;                     % number of samples
t  = (0:N-1)' * Ts;

sys_th = ss(A_th, B_th, C_th, D_th);     % continuous system
u = randn(N,1);                          % random input
y0 = lsim(sys_th, u, t);                 % ideal output

noise_std = 0.001;                        % small measurement noise
y = y0 + noise_std * randn(N,1);

% PBSID PARAMETERS
f = 10; p = 20; n = 1;

% PBSID (Black-Box Identification)
[Ai, Bi, Ci, Di] = estimate_ABCDK(u', y', f, p, n);
sys_id = ss(Ai, Bi, Ci, Di, Ts);

% Simulate identified model
y_id = lsim(sys_id, u, t);

% Plot
% figure; plot(t, y, 'b', t, y_id, 'r--');
% xlabel('Time [s]'); ylabel('Output');
% legend('Measured', 'Identified'); grid on;
% title('PBSID Model vs Measured Output');

% STRUCTURING METHOD
param0 = [Ai, Bi];  % initial guess for [tau, K]
freq_vec = 2*pi*logspace(-1, 1, length(u));  % up to ~100 rad/s

[structuredModel, data_BB_freq] = structuring(sys_id, param0, @thermal_dynamics, freq_vec);

% Estimated parameters
est_params = getpvec(structuredModel);

% --- TRUE MATRICES ---
A_true = -1 / tau_true;
B_true =  K_true / tau_true;
C_true = 1;
D_true = 0;

% --- ESTIMATED MATRICES ---
A_est = structuredModel.A;
B_est = structuredModel.B;
C_est = structuredModel.C;
D_est = structuredModel.D;

fprintf('\n--- State-Space Matrices Comparison ---\n\n');
fprintf('A_est = %.4f \t A_true = %.4f\n', A_est, A_true);
fprintf('B_est = %.4f \t B_true = %.4f\n', B_est, B_true);
fprintf('C_est = %.4f \t C_true = %.4f\n', C_est, C_true);
fprintf('D_est = %.4f \t D_true = %.4f\n', D_est, D_true);

fprintf('\nEstimated tau = %.4f (true: %.2f)\n', est_params(1), tau_true);
fprintf('Estimated  K  = %.4f (true: %.2f)\n\n', est_params(2), K_true);

% Compare Bode plots
% figure;
% bode(data_BB_freq, structuredModel);
% legend('Black-box (FRD)', 'Structured Model');
% title('Frequency-Domain Fit');

% Generate PRBS input

u_chirp = chirp(t, 0.1, t(end), 10);

% Simulate true system
y_true_prbs = lsim(sys_th, u_chirp, t);

% Reconstruct structured model from estimated parameters
tau_hat = est_params(1);
K_hat   = est_params(2);
A_struct = -1 / tau_hat;
B_struct = K_hat / tau_hat;
C_struct = 1;
D_struct = 0;
sys_structured = ss(A_struct, B_struct, C_struct, D_struct);

% Simulate structured model on PRBS input
y_val_prbs = lsim(sys_structured, u_chirp, t);

figure;
plot(t, u_chirp, 'b');
xlabel('Time [s]');
ylabel('Output');
title('Input');
grid on;

% Plot comparison
figure;
plot(t, y_true_prbs, 'k', t, y_val_prbs, 'r--');
xlabel('Time [s]');
ylabel('Output');
legend('True System', 'Structured Model');
title('Validation: Structured Model vs True System');
grid on;

% === Frequency-Domain Comparison (Bode Plot) ===
figure;
bode(sys_th, 'k', sys_structured, 'r--');
legend('True System', 'Identified Model');
title('Validation: Bode Plot - True vs Structured System');
grid on;

% === Pole-Zero Map Comparison ===
figure;
pzmap(sys_th);
hold on;
pzmap(sys_structured);
legend('True System', 'Identified Model');
title('Validation: Pole-Zero Map - True vs Structured System');
grid on;

% Error metrics
disp('--- Black Box Validation Evaluation ---');
VAF_BLACK = vaf(y_true_prbs, y_val_prbs);
fit_black = fit(y_true_prbs, y_val_prbs);
pec_black = pec(y_true_prbs, y_val_prbs);


%% Identify a black box model for the longitudinal dynamics of the UAV & estimate its parameters

close all; clearvars; 

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

% Black-Box Identification

% p_values = 10:5:100;
% f_values = 5:5:50;
% AIC_matrix = NaN(length(p_values), length(f_values)); % Initialize with NaNs
% 
% h = waitbar(0, 'Searching over (p,f)...');
% counter = 1;
% total = length(p_values) * length(f_values);
% 
% for i = 1:length(p_values)
%     p = p_values(i);
% 
%     for j = 1:length(f_values)
%         f = f_values(j);
% 
%         if f >= p
%             continue; % Enforce p > f
%         end
% 
%         try
%             [Aest, Best, Cest, Dest, Kest, Sest] = estimate_ABCDK(input', true_output', f, p, 3);
%             black_sys = ss(Aest, Best, Cest, Dest, sample_time);
% 
%             time = (0:size(input,1)-1)' * sample_time;
%             black_valout = lsim(black_sys, input, time);
% 
%             residuals = black_valout - true_output;
%             V = var(residuals);
% 
%             n = size(Aest, 1); % number of states
%             m = size(Best, 2); % number of inputs
%             k = n * p + n * m;
% 
%             AIC_matrix(i,j) = 2 * k + length(black_valout) * log(max(V, eps)); % avoid log(0)
%         catch
%             AIC_matrix(i,j) = NaN; % Handle estimation failures gracefully
%         end
% 
%         waitbar(counter / total, h, sprintf('Progress: %d/%d', counter, total));
%         counter = counter + 1;
%     end
% end
% 
% close(h);
% 
% % Find the optimal (p, f) with minimum AIC
% [minAIC, idx] = min(AIC_matrix(:));
% [i_opt, j_opt] = ind2sub(size(AIC_matrix), idx);
% optimal_p = p_values(i_opt);
% optimal_f = f_values(j_opt);
% optimal_AIC = minAIC;
% 
% disp(['Optimal past window size p: ', num2str(optimal_p), ', corresponding f: ', num2str(optimal_f), ', AIC = ', num2str(optimal_AIC)]);
% 
% figure;
% contourf(f_values, p_values, AIC_matrix, 100); 
% xlabel('Future Horizon (f)');
% ylabel('Past Window Size (p)');
% title('AIC Porkchop Plot');
% colorbar;
% grid on;
% set(gca, 'YDir', 'normal'); % so p increases upward
% hold on;
% plot(f_values(j_opt), p_values(i_opt), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
% text(f_values(j_opt), p_values(i_opt), sprintf('  min AIC = %.2f', minAIC), 'Color', 'red');
% hold off;

n = 3;    % order of system

p = 70;   % past window size
f = 35;   % future window size

[Aest, Best, Cest, Dest] = PBSID(input', true_output', f, p, n)
black_sys = ss(Aest, Best, Cest, Dest, sample_time);
black_idout = lsim(black_sys, input, time);

guess = [0 0 0 -2 -10 300]';
freq_vec = 2*pi*logspace(-1, 2, length(time)); 

%% Structure black-box into grey-box

frd_data = frd(black_sys, freq_vec);
sys_black = idgrey(@longitdynRot, guess, 'c');

opt = greyestOptions();
opt.SearchMethod = 'lm';
opt.SearchOptions.MaxIterations = 100;
opt.SearchOptions.Tolerance = 1e-6;
opt.Display = 'on';

% Estimate structured model
structuredModel = greyest(frd_data, sys_black, opt);

% Extract matrices and parameters
A_blackstruct =  structuredModel.A;
B_blackstruct  = structuredModel.B;
C_blackstruct  = structuredModel.C;
D_blackstruct  = structuredModel.D;
est_params = getpvec(structuredModel);

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


%% Black-Box Validation

% Load input excitation (3211 maneuver)
load Excitation3211.mat
t3211 = Excitation3211(:, 1);
u3211 = Excitation3211(:, 2);
simulation_time = t3211(end) - t3211(1);
Ts = sample_time;  % ensure sample_time is loaded beforehand

% ------------------------------------------------------------------------
% TRUE SYSTEM Simulation 
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

assignin('base', 'A', structuredModel.A);
assignin('base', 'B', structuredModel.B);
assignin('base', 'C', eye(3));
assignin('base', 'D', zeros(3, 1));

bdclose('all');
load_system('Validation');
set_param('Validation', "FastRestart", "off");
estimated_validation = sim('Validation', 'SrcWorkspace', 'current');
if exist('slprj', 'dir'), rmdir('slprj', 's'); end

black_valout = estimated_validation.q;
black_valsys = ss(A_blackstruct, B_blackstruct, C_blackstruct, D_blackstruct, Ts);

% ------------------------------------------------------------------------
% Evaluation Metrics - Black Box
% ------------------------------------------------------------------------

disp('--- Black Box Validation Evaluation ---');
VAF_BLACK = vaf(true_output_validation, black_valout);
fit_black = fit(true_output_validation, black_valout);
pec_black = pec(true_output_validation, black_valout);
fprintf('VAF: %.2f%% | FIT: %.2f%% | PEC: %.2f%%\n', ...
        VAF_BLACK, fit_black, pec_black);

% Output Comparison
figure;
plot(time, true_output_validation, 'k', 'LineWidth', 1.5); hold on;
plot(time, black_valout, 'r--', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Pitch Rate $q(t)$ [rad/s]', 'Interpreter', 'latex');
title(sprintf('3211 Validation – Black-Box vs True Output (VAF: %.3f%%)', VAF_BLACK));
legend('True Output', 'Grey-Box Model');
grid on; xlim tight; set(gca, 'FontSize', 12);

% Bode Plot
figure;
bodemag(true_sys, 'k-', black_valsys, 'r--', {0.1, 100});
legend('True System', 'Identified Model', 'Location', 'Best');
title('Bode Magnitude Comparison');
grid on;

% Pole Comparison (clean)
true_poles = pole(true_sys);
id_poles = pole(black_valsys);

figure;
hold on;
plot(real(true_poles), imag(true_poles), 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'True Poles');
plot(real(id_poles), imag(id_poles), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Identified Poles');
legend('show');
title('Pole Comparison');
xlabel('Real Axis'); ylabel('Imaginary Axis');
grid on; axis equal; set(gca, 'FontSize', 12);


%% Functions

function [A, B, C, D] = longitdynRot(THETA, Ts)
    
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

function [A_BB,B_BB,C_BB,D_BB,K_BB] = PBSID(u,y,f,p,n)
    
    % PBSID-varx
    [S,X] = dordvarx(u,y,f,p);
    figure;
    hold on
    scatter(1:length(S), S, 'filled');
    set(gca, 'YScale', 'log'); % Use log scale for better visibility of gaps
    xlabel('Index');
    ylabel('Singular Value (log scale)');
    title('Singular Value Spectrum');
    x = dmodx(X,n);
    [A_BB,B_BB,C_BB,D_BB,K_BB] = dx2abcdk(x,u,y,f,p);
    
end
    
function X = dmodx(X,n)
    
    % check number of arguments
    if nargin < 2
        error('DMODX requires at least two input arguments.');
    end
    
    % check for batches
    if iscell(X)
        batch = length(X);
    else
        batch = 1;
    end
    
    % do for all batches
    for k = 1:batch
        if batch == 1
            x = X;
        else
            x = X{k};
    end
    
    % check the dimensions of the inputs
    mx = size(x,1);
    if (n < 1) || isempty(n)
        error('System order of zero or lower does not make sense!')
    end
    if mx < n
        error('The number of rows of matrix X must be equal or higher then the order n.')
    end
    
    % use only n rows
    x = x(1:n,:);
    
        if batch == 1
            X = x;
        else
            X{k} = x;
        end
    end
    
    
end    
    
function [S,X,VARX,U,Zps] = dordvarx(u,y,f,p,reg,opt,weight,noD)
    
    
    % check number of input arguments
    if nargin < 4
    error('DORDVARX requires four or five input arguments.')
    end
    
    % assign default values to unspecified parameters
    if (nargin < 8) || isempty(noD)
    noD = 0;
    end
    if (nargin < 7) || isempty(weight)
    weight = 0;
    end
    if (nargin < 6) || isempty(opt)
    opt = 'gcv';
    end
    if (nargin < 5) || isempty(reg)
    reg = 'none';
    end
    
    % check the size of the windows
    if f > p
        error('Future window size f must equal or smaller then past window p. (f <= p)')
    end
    
    % % check for batches
    % if iscell(y)
    %     batch = length(y);
    %     ZZ = cell(1,batch);
    %     X = cell(1,batch);
    %     yb = y;
    %     ub = u;
    % else
    %     batch = 1;
    % end
    % 
    % % do for all batches
    % for k = 1:batch
    %     if batch > 1
    %         y = yb{k};
    %         u = ub{k};
    %     end
    % 
    % check dimensions of inputs
    if size(y,2) < size(y,1)
        y = y';
    end
    N = size(y,2);
    l = size(y,1);
    
    if isempty(u)
        r = 0;
        u = zeros(0,N);
    else

    if size(u,2) < size(u,1)
        u = u';
    end
    r = size(u,1);
        if ~isequal(N,length(u))
            error('The number of rows of vectors/matrices u and y must be the same.')
        end
    end
    if l == 0
        error('DORDVARX requires an output vector y.')
    end
    
    % store the past and future vectors
    m = r+l;      % since I have a SISO r = 1 (number of inputs) and l = 1 (numer of outputs)
    z = [u; y];
    
    % creating past data matrix
    Z = zeros(p*m,N-p); 
    for i = 1:p
        Z((i-1)*m+1:i*m,:) = z(:,i:N+i-p-1); 
    end
    
    % solve VARX problem
    Y = y(:,p+1:N);
    U = u(:,p+1:N);
    
    PHI = [Z; U];
    VARX = Y * pinv(PHI);  %function to solve the pseudo inverse
    
    
    % construct LambdaKappa (GAMMA * DELTA)
    LK = zeros(f*l,p*m);

    if weight == 0
        for i = 1:f
            LK((i-1)*l+1:i*l,p*m-(p-i+1)*m+1:p*m) = VARX(:,1:(p-i+1)*m);
        end
        elseif weight == 1
        for i = 0:f-1
            LK(i*l+1:(i+1)*l,i*m+1:p*m) = VARX(:,1:(p-i)*m);
            if i ~= 0
                for j = 0:i-1
                    LK(i*l+1:(i+1)*l,:) = LK(i*l+1:(i+1)*l,:) + VARX(:,(p-i+j)*m+r+(1:l))*LK(j*l+1:(j+1)*l,:);
                end
            end
        end
    end
    
    % singular value decomposition
    
    [U,S,V] = svd(LK*Z(1:p*m,:),'econ');
    S = diag(S)';
    X = diag(sqrt(S))*V';
    
    if nargout > 3
        if weight == 0
            U = diag(1./sqrt(S))*U';
        elseif weight == 1
            LK0 = zeros(f*l,p*m);
        for i = 1:f
            LK0((i-1)*l+1:i*l,p*m-(p-i+1)*m+1:p*m) = VARX(:,1:(p-i+1)*m);
        end
            U = diag(1./sqrt(S))*U'*(LK*pinv(LK0));
        end
    end
    
end
    
function [A, B, C, D, K] = dx2abcdk(x, u, y, f, p, c)

    if nargin < 5
        error('DX2ABCDK requires at least five input arguments.');
    end
    if nargin == 5 || isempty(c)
        c = 'none';
    end

    % Ensure proper orientation
    if size(y, 2) < size(y, 1), y = y'; end
    if size(x, 2) < size(x, 1), x = x'; end
    if size(u, 2) < size(u, 1), u = u'; end

    % Sizes
    [l, N] = size(y);
    n = size(x, 1);
    r = size(u, 1);

    % Sanity checks
    if l == 0, error('Output vector y cannot be empty.'); end
    if n == 0, error('State vector x cannot be empty.'); end

    if N ~= size(u, 2)
        error('Input and output time series must have the same length.');
    end
    if rank(x) < n
        error('State matrix x is not full rank (rank(x) < n).');
    end
    if f > p
        error('Future window size f must be less than or equal to past window size p.');
    end

    % Adjust for windowed data
    u = u(:, p+1:p+size(x, 2));
    y = y(:, p+1:p+size(x, 2));

    % Estimate C and D
    CD = y(:, 1:end-1) * pinv([x(:, 1:end-1); u(:, 1:end-1)]);
    e = y - CD * [x; u];

    % Estimate A, B, K
    z = [x(:, 1:end-1); u(:, 1:end-1); e(:, 1:end-1)];
    ABK = x(:, 2:end) * pinv(z);

    A = ABK(:, 1:n);
    B = ABK(:, n+1:n+r);
    K = ABK(:, n+r+1:end);

    C = CD(:, 1:n);
    D = CD(:, n+1:n+r); 
end

function [sys_est_grey_from_black,data_BB_freq] = structuring(sysBB, theta0, greyfun,freq_vec)
    
    % Structuring with Transfer Function
    sys_est_cont = sysBB;
    frd_data = frd(sys_est_cont, freq_vec);
    data_BB_freq = frd_data;
    sys_grey = idgrey(greyfun, theta0, 'c');

    opt = greyestOptions();
    opt.SearchMethod = 'lm';
    opt.SearchOptions.MaxIterations = 100;
    opt.SearchOptions.Tolerance = 1e-6;
    opt.EnforceStability = true;
    opt.Display = 'on';

    sys_est_grey_from_black = greyest(frd_data, sys_grey, opt); 
    
    figure;
    bode(sys_est_cont, sys_est_grey_from_black)
    
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

S = estimate_order(u,y,f,p);

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

function [A, Acl, B, K, C, ABK, s_singular, X_p_p_l] = estimateModel(U, Y, Markov, Z_0_pm1_l, past, future, order_estimate)
    % Function to estimate state-space model parameters using PBSID
    % Inputs:
    %   U: Input matrix (m x timeSteps)
    %   Y: Output matrix (r x timeSteps)
    %   Markov: Markov parameter matrix
    %   Z_0_pm1_l: Data matrix from past steps
    %   past: Number of past time steps
    %   future: Number of future time steps
    %   order_estimate: Estimated system order
    % Outputs:
    %   A: State-space A matrix (including feedback)
    %   Acl: Closed-loop A matrix (without feedback)
    %   B: State-space B matrix
    %   K: Feedback matrix
    %   C: Output matrix
    %   s_singular: Singular values from SVD
    %   X_p_p_l: Estimated state sequence

    [m, timeSteps] = size(U);
    r = size(Y, 1);
    l = timeSteps - past - 1;
    n = order_estimate;

    % Initialize Qpm1
    Qpm1 = zeros(future * r, past * (m + r));
    for i = 1:future
        Qpm1((i-1)*r+1:i*r, (i-1)*(m+r)+1:end) = Markov(:, 1:(m+r)*(past-i+1));
    end

    % Estimate state sequence
    Qpm1_times_Z_0_pm1_l = Qpm1 * Z_0_pm1_l;
    [Usvd, Ssvd, Vsvd] = svd(Qpm1_times_Z_0_pm1_l, 'econ');
    s_singular = diag(Ssvd);
    % Estimated state sequence
    X_p_p_l = diag(sqrt(s_singular(1:n))) * Vsvd(:, 1:n)';

    % Partition state sequence
    X_pp1_pp1_lm1 = X_p_p_l(:, 2:end);
    X_p_p_lm1 = X_p_p_l(:, 1:end-1);

    % Form Z_p_p_lm1 and Y_p_p_l
    Z_p_p_lm1 = zeros(m + r, l);
    Z_p_p_lm1(1:m, :) = U(:, past+1:past+l);
    Z_p_p_lm1(m+1:end, :) = Y(:, past+1:past+l);

    Y_p_p_l = Y(:, past+1:end);

    % Concatenate S matrix and compute ABK
    S = [X_p_p_lm1; Z_p_p_lm1];
    ABK = X_pp1_pp1_lm1 * pinv(S);

    % Estimate C, Acl, B, K, and A matrices
    C = Y_p_p_l * pinv(X_p_p_l);

    Acl = ABK(1:n, 1:n);
    B = ABK(1:n, n+1:n+m);
    K = ABK(1:n, n+m+1:n+m+r);
    A = Acl + K * C;
end

function [A, B, C, D] = thermal_dynamics(param, Ts)
    tau = param(1);
    K   = param(2);

    A = -1 / tau;
    B = K / tau;
    C = 1;
    D = 0;
end