clc; close all; clearvars

rng(42);  % For reproducibility

% ------------------------------------------------------------------------
% True Model Parameters
% ------------------------------------------------------------------------
Xu = -0.1068;     Xq = 0.1192;
Mu = -5.9755;     Mq = -2.6478;
Xd = -10.1647;    Md = 450.71;
A = [Xu, Xq, -9.81; Mu, Mq,  0; 0,  1,  0];
B = [Xd; Md; 0];
C = eye(3);  
D = zeros(3,1);

% ------------------------------------------------------------------------
% Noise + Delay
% ------------------------------------------------------------------------
noise.Enabler = 1;
noise.pos_stand_dev      = noise.Enabler * 0.0011;
noise.vel_stand_dev      = noise.Enabler * 0.01;
noise.attitude_stand_dev = noise.Enabler * deg2rad(0.33);
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(1);
delay.position_filter = 1;
delay.attitude_filter = 1;
delay.mixer = 1;

% ------------------------------------------------------------------------
% Input and Simulink Simulation
% ------------------------------------------------------------------------
parameters_controller;
load ExcitationM
SetPoint = [0, 0];
t = ExcitationM(:,1);
simulation_time = t(end) - t(1);

load_system('Identification');
set_param('Identification', "FastRestart", "off");
simulation = sim('Identification', 'SrcWorkspace', 'current');
if exist('slprj', 'dir'); rmdir('slprj', 's'); end

time  = (0:length(simulation.Mtot)-1)' * sample_time;
input = simulation.Mtot;
output = simulation.q;

% ------------------------------------------------------------------------
% Cross-validation on Full Data for Structuring
% ------------------------------------------------------------------------
n_vals = 3;
p_vals = 0:25:400;
f_vals = 0:10:300;
freq_vec = 2*pi*logspace(-1, 2, 1000);
guess = [0 0 0 -2 -10 300]';

best_rmse = Inf;
best_combo = [NaN NaN NaN];
best_est_params = [];
RMSE_matrix = NaN(length(p_vals), length(f_vals));  % Rows = p, Columns = f

fprintf('Running parameter-based cross-validation over (n, p, f)...\n');
total = length(p_vals) * length(f_vals);
count = 1;
h = waitbar(0, 'Searching for best (n, p, f)...');

tic
for pi = 1:length(p_vals)
    p = p_vals(pi);
    for fi = 1:length(f_vals)
        f = f_vals(fi);
        waitbar(count / total, h, sprintf('Trying p=%d, f=%d...', p, f));
        count = count + 1;

        if f >= p
            continue;
        end

        try
            [Aest, Best, Cest, Dest] = estimate_ABCDK(input', output', f, p, n_vals);
            sys_id = ss(Aest, Best, Cest, Dest, sample_time);
            
            frd_data = frd(sys_id, freq_vec);
            sys_black = idgrey(@longitdynRot, guess, 'c');
            opt = greyestOptions();
            opt.SearchMethod = 'lm';
            opt.SearchOptions.MaxIterations = 100;
            opt.SearchOptions.Tolerance = 1e-6;
            opt.Display = 'off';
            structuredModel = greyest(frd_data, sys_black, opt);
            est_params = getpvec(structuredModel);

            % Relative RMSE
            true_params = [-0.1068, 0.1192, -5.9755, -2.6478, -10.1647, 450.71];
            rel_errors = (est_params' - true_params) ./ true_params;
            rel_rmse = sqrt(mean(rel_errors.^2));
            if rel_rmse * 100 > 500  % Arbitrary threshold, e.g., 500%
                fprintf('⚠️ Discarding run with unstable parameters: RMSE = %.2f%%\n', rel_rmse * 100);
                continue;
            end

            RMSE_matrix(pi, fi) = rel_rmse * 100;  % store as percentage

            % Print summary
            labels = {'Xu', 'Xq', 'Mu', 'Mq', 'Xd', 'Md'};
            abs_errors = abs(est_params' - true_params);
            percent_errors = 100 * abs_errors ./ abs(true_params);
            fprintf('\n--- Parameter Comparison (p = %d, f = %d) ---\n\n', p, f);
            fprintf('%-10s %15s %15s %15s %15s\n', 'Parameter', 'Estimated', 'True', 'Abs Error', '% Error');
            fprintf('%s\n', repmat('-', 1, 80));
            for j = 1:length(true_params)
                fprintf('%-10s %+15.4f %+15.4f %15.4f %14.2f%%\n', ...
                    labels{j}, est_params(j), true_params(j), abs_errors(j), percent_errors(j));
            end
            fprintf('\nRelative RMSE (this run): %.2f%%\n', rel_rmse * 100);

            if rel_rmse < best_rmse
                best_rmse = rel_rmse;
                best_combo = [n_vals, p, f];
                best_est_params = est_params;
            end

        catch
            fprintf('❌ Failed for p=%d, f=%d\n', p, f);
            continue;
        end
    end
end
elapsed_time = toc;
close(h);

%% Final report
fprintf('\n✅ Best (n, p, f) combination: (%d, %d, %d)\n', best_combo);
fprintf('🔍 Best relative RMSE: %.2f%%\n', best_rmse * 100);
fprintf('\n🧾 Best Estimated Parameters:\n');
fprintf('  %-10s %12s\n', 'Parameter', 'Estimate');
fprintf('  %-10s %12s\n', '---------', '--------');
param_labels = {'Xu', 'Xq', 'Mu', 'Mq', 'Xd', 'Md'};
for k = 1:length(best_est_params)
    fprintf('  %-10s %12.5f\n', param_labels{k}, best_est_params(k));
end
fprintf('Elapsed Time: %.2f min\n', elapsed_time/60);

% Porkchop plot of Relative RMSE
figure;
contourf(f_vals, p_vals, RMSE_matrix, 100, 'LineColor', 'none');
colorbar;
xlabel('Future window size f');
ylabel('Past window size p');
title(sprintf('Relative RMSE Porkchop Plot (n = %d)', n_vals));
set(gca, 'YDir', 'normal');  % ensures p increases upward
grid on;

% Highlight best combo
hold on;
plot(best_combo(3), best_combo(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
text(best_combo(3), best_combo(2), sprintf('  Min RMSE = %.2f%%', best_rmse * 100), ...
    'Color', 'r', 'FontWeight', 'bold');
hold off;


%% Use all available data now

[Aest, Best, Cest, Dest] = estimate_ABCDK(input', output', best_combo(3), best_combo(2), best_combo(1));
black_sys = ss(Aest, Best, Cest, Dest, sample_time);
black_idout = lsim(black_sys, input, time);

%%

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