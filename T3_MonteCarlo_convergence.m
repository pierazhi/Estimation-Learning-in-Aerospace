clc;
clear;
close all;

% --- Load saved Monte Carlo parameter estimates ---
load("MonteCarlo_greyparams.mat");  % params_grey: [N × 6]
labels = {'Xu', 'Xq', 'Mu', 'Mq', 'Xd', 'Md'};

% === Settings ===
epsilon = 1e-4;
min_stable_iters = 50;

% === Initialization ===
[N, num_params] = size(params_grey);
theta_mean = NaN(N, num_params);
delta_matrix = NaN(N, num_params);
delta_overall = NaN(N, 1);
iter_converged = NaN(1, num_params);
iter_converged_overall = NaN;

% === Compute running mean and deltas ===
for i = 1:N
    theta_mean(i, :) = mean(params_grey(1:i, :), 1);
    if i > 1
        % Per-parameter relative delta
        delta_matrix(i, :) = abs(theta_mean(i,:) - theta_mean(i-1,:)) ./ abs(theta_mean(i,:));
        % Overall delta (full vector)
        delta_overall(i) = norm(theta_mean(i,:) - theta_mean(i-1,:)) / norm(theta_mean(i,:));
    end
end

% === Check per-parameter convergence ===
for j = 1:num_params
    stable_count = 0;
    for i = 2:N
        if delta_matrix(i,j) < epsilon
            stable_count = stable_count + 1;
        else
            stable_count = 0;
        end
        if stable_count >= min_stable_iters
            iter_converged(j) = i;
            break;
        end
    end
end

% === Check overall vector convergence ===
stable_count = 0;
for i = 2:N
    if delta_overall(i) < epsilon
        stable_count = stable_count + 1;
    else
        stable_count = 0;
    end
    if stable_count >= min_stable_iters
        iter_converged_overall = i;
        break;
    end
end

% === Print summary ===
fprintf('\n📋 Convergence Summary (ε = %.1e, stable for %d steps):\n', epsilon, min_stable_iters);
for j = 1:num_params
    if ~isnan(iter_converged(j))
        fprintf('✅ %s converged at iteration %d\n', labels{j}, iter_converged(j));
    else
        fprintf('❌ %s did not converge\n', labels{j});
    end
end
if ~isnan(iter_converged_overall)
    fprintf('✅ Overall parameter vector converged at iteration %d\n', iter_converged_overall);
else
    fprintf('❌ Overall parameter vector did not converge\n');
end

% === Plot: Relative delta for each parameter ===
% figure('Name','Relative Change per Parameter');
% for j = 1:num_params
%     subplot(3,2,j);
%     semilogy(delta_matrix(:,j), 'b', 'LineWidth', 1.5); hold on;
%     if ~isnan(iter_converged(j))
%         xline(iter_converged(j), '--r', 'LineWidth', 1.2);
%     end
%     title(sprintf('Convergence Δ for %s', labels{j}));
%     xlabel('Iteration'); ylabel('Rel. Δ (log)'); grid on;
% end

% === Plot: Overall relative change ===
figure('Name','Overall Parameter Convergence');
semilogy(delta_overall, 'Color', [.5 .5 .5], 'LineWidth', 1); hold on;
if ~isnan(iter_converged_overall)
    xline(iter_converged_overall, '--r', 'LineWidth', 1.5);
end
xlabel('Iteration'); ylabel('Rel. Δ');
title('Overall Convergence of Parameter Mean Vector');
grid on;
legend('Rel. Δ', 'Convergence')

% === Plot: Running means ===
figure('Name','Parameter Mean Convergence');
for j = 1:num_params
    subplot(3,2,j);
    plot(theta_mean(:,j), 'LineWidth', 1, 'Color', [.5, .5, .5]); hold on;
    % if ~isnan(iter_converged(j))
    %     xline(iter_converged(j), '--r', 'LineWidth', 1.2);
    % end
    title(sprintf('Mean Convergence: %s', labels{j}));
    xlabel('Iteration'); ylabel('Mean Value'); grid on;
end

%% BLACK CONVERGENCE
clc;
clear;
close all;

% === Load your black-box Monte Carlo results ===
load("MonteCarlo_blackparams.mat");  % should contain params_black
labels = {'Xu', 'Xq', 'Mu', 'Mq', 'Xd', 'Md'};

% === Settings ===
epsilon = 1e-4;
min_stable_iters = 50;

% === Step 1: Remove skipped runs (all-zero rows) ===
zero_rows = all(params_black == 0, 2);
params_black_clean = params_black(~zero_rows, :);

% === Initialization ===
[N, num_params] = size(params_black_clean);
theta_mean = NaN(N, num_params);
delta_matrix = NaN(N, num_params);
delta_overall = NaN(N, 1);
iter_converged = NaN(1, num_params);
iter_converged_overall = NaN;

% === Step 2: Running mean and relative deltas ===
for i = 1:N
    theta_mean(i, :) = mean(params_black_clean(1:i, :), 1);
    if i > 1
        delta_matrix(i, :) = abs(theta_mean(i,:) - theta_mean(i-1,:)) ./ abs(theta_mean(i,:));
        delta_overall(i) = norm(theta_mean(i,:) - theta_mean(i-1,:)) / norm(theta_mean(i,:));
    end
end

% === Step 3: Per-parameter convergence ===
for j = 1:num_params
    stable_count = 0;
    for i = 2:N
        if delta_matrix(i,j) < epsilon
            stable_count = stable_count + 1;
        else
            stable_count = 0;
        end
        if stable_count >= min_stable_iters
            iter_converged(j) = i;
            break;
        end
    end
end

% === Step 4: Overall convergence ===
stable_count = 0;
for i = 2:N
    if delta_overall(i) < epsilon
        stable_count = stable_count + 1;
    else
        stable_count = 0;
    end
    if stable_count >= min_stable_iters
        iter_converged_overall = i;
        break;
    end
end

% === Print Summary ===
fprintf('\n📋 Black-Box Convergence Summary (ε = %.1e, stable for %d steps):\n', epsilon, min_stable_iters);
for j = 1:num_params
    if ~isnan(iter_converged(j))
        fprintf('✅ %s converged at iteration %d\n', labels{j}, iter_converged(j));
    else
        fprintf('❌ %s did not converge\n', labels{j});
    end
end
if ~isnan(iter_converged_overall)
    fprintf('✅ Overall parameter vector converged at iteration %d\n', iter_converged_overall);
else
    fprintf('❌ Overall parameter vector did not converge\n');
end

% === Plot: Relative change per parameter ===
% figure('Name','Black-Box: Relative Change per Parameter');
% for j = 1:num_params
%     subplot(3,2,j);
%     semilogy(delta_matrix(:,j), 'b', 'LineWidth', 1.5); hold on;
%     if ~isnan(iter_converged(j))
%         xline(iter_converged(j), '--r', 'LineWidth', 1.2);
%     end
%     title(sprintf('Convergence Δ for %s', labels{j}));
%     xlabel('Iteration'); ylabel('Rel. Δ (log)'); grid on;
% end

% === Plot: Overall convergence ===
figure('Name','Black-Box: Overall Parameter Convergence');
semilogy(delta_overall, 'k', 'LineWidth', 1); hold on;
if ~isnan(iter_converged_overall)
    xline(iter_converged_overall, '--r', 'LineWidth', 1.5);
end
xlabel('Iteration'); ylabel('Rel. Δ');
title('Overall Convergence of Parameter Mean Vector');
grid on;
legend('Rel. Δ', 'Convergence')

% === Plot: Running mean per parameter ===
figure('Name','Black-Box: Parameter Mean Convergence');
for j = 1:num_params
    subplot(3,2,j);
    plot(theta_mean(:,j), 'LineWidth', 1, 'Color', 'k'); hold on;
    % if ~isnan(iter_converged(j))
    %     xline(iter_converged(j), '--r', 'LineWidth', 1.2);
    % end
    title(sprintf('Mean Convergence: %s', labels{j}));
    xlabel('Iteration'); ylabel('Mean Value'); grid on;
end


%% T3.2 GREY Monte Carlo
clc;
clearvars;
close all;


N = 500;
params_grey = zeros(N, 6);  % Preallocate

wmin = -2; wmax = 2; Npoints = 1000;
w = logspace(wmin, wmax, Npoints);
mag_all = NaN(N, length(w));
phase_all = NaN(N, length(w));

theta_mean = zeros(N, 6);
epsilon = 1e-5;
stability_counter = 0;
min_stable_iters = 25;
        
wait_h = waitbar(0, 'Starting Monte Carlo...');
start_time = tic;
for i = 1:N
    elapsed = toc(start_time);
    avg_iter_time = elapsed / i;
    est_total_time = avg_iter_time * N;
    time_left = est_total_time - elapsed;

    % Format ETA as mm:ss
    hrs  = floor(time_left / 3600);
    mins = floor(mod(time_left, 3600) / 60);
    secs = mod(round(time_left), 60);
    eta_str = sprintf('%02d:%02d:%02d', hrs, mins, secs);

    % Update waitbar with full-run ETA
    waitbar(i/N, wait_h, ...
        sprintf('Run %d of %d — ETA: %s', i, N, eta_str));
    rng(42);

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
     guess = {
        'Xu',  0;
        'Xq',   0;
        'Mu',  0;
        'Mq',  -2.5;
        'Xd', -10.0;
        'Md', 300
    };
    
    % Grey-box model setup
    opt = greyestOptions();
    opt.SearchMethod = 'lm';
    opt.SearchOptions.MaxIterations = 100;
    opt.SearchOptions.Tolerance = 1e-6;
    opt.Display = 'off';
    
    sys_init   = idgrey(@statespacefun, guess, 'c');
    grey_id = greyest(sys_ft, sys_init, opt);  % Estimate parameters
   
    [mag, phase] = bode(grey_id, w);
    mag_grey = squeeze(20*log10(mag(1,1,:)));    % pitch rate magnitude [dB]
    phase_grey = squeeze(phase(1,1,:));          % pitch rate phase [deg]

    mag_all(i, :) = mag_grey;
    phase_all(i, :) = phase_grey;

    % Extract and reassign estimated values
    param_names = {'X_u', 'X_q', 'M_u', 'M_q', 'X_d', 'M_d'}';
    params = getpvec(grey_id);
    Xu = params(1); Xq = params(2);
    Mu = params(3); Mq = params(4);
    Xd = params(5); Md = params(6);
    
    cov_matrix = getcov(grey_id);  % Asymptotic covariance
    std_devs = sqrt(diag(cov_matrix));  % Standard deviation (1-sigma)

    params_grey(i, :) = getpvec(grey_id);  % size [N × 6]
    % Compute running mean

    if i == 1
        theta_mean(i, :) = params_grey(i, :);
    else
        theta_mean(i, :) = mean(params_grey(1:i, :), 1);
        % Compute relative change in the mean vector
        delta_mean = norm(theta_mean(i, :) - theta_mean(i-1, :)) / norm(theta_mean(i, :))
    
        % Check for convergence
        if delta_mean < epsilon
            stability_counter = stability_counter + 1
        else
            stability_counter = 0;
        end
    
        if stability_counter >= min_stable_iters
            fprintf('\n✅ Parameter mean converged after %d iterations (%.2e rel change for %d steps)\n', ...
                    i, delta_mean, stability_counter);
            params_grey = params_grey(1:i, :);
            mag_all = mag_all(1:i, :);
            phase_all = phase_all(1:i, :);
            theta_mean = theta_mean(1:i, :);
            break;
        end
    end

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
    for j = 1:length(true_params)
        fprintf('%-10s %+15.4f %+15.4f %15.4f %14.2f%%\n', ...
            labels{j}, params(j), true_params(j), abs_errors(j), percent_errors(j));
    end
    
    % Compute RMSE
    rel_errors = (params' - true_params) ./ true_params;
    rel_rmse = sqrt(mean(rel_errors.^2));

    if rel_rmse * 100 > 100  % Arbitrary threshold, e.g., 500%
                fprintf('⚠️ Discarding run with unstable parameters: RMSE = %.2f%%\n', rel_rmse * 100);
                continue;
    end
    
    fprintf('\nRelative RMSE: %.2f%%\n', rel_rmse * 100);    
    g = 9.81;
    
    A_grey_id = grey_id.A;
    B_grey_id = grey_id.B;
    C_grey_id = [1, 0, 0;
         0, 1, 0;
         0, 0, 1];
    D_grey_id = [0;
         0;
         0];

    grey_sys = ss(A_grey_id, B_grey_id, C_grey_id, D_grey_id, sample_time);  % Convert to standard state-space object
    
    % Compute poles and zeros
    grey_poles(i, :) = pole(grey_sys(2)).';
    grey_zeros(i, :) = tzero(grey_sys(2)).';
end

% === End of loop: display total time taken
total_time = toc(start_time);
fprintf('\n✅ Monte Carlo complete — Total time: %.1f seconds (%.1f minutes)\n', ...
        total_time, total_time / 60);
close(wait_h);

figure;
for j = 1:6
    subplot(3,2,j);
    plot(1:size(theta_mean,1), theta_mean(:,j), 'LineWidth', 1.5);
    title(sprintf('Mean Convergence: %s', labels{j}));
    xlabel('Iteration'); ylabel('Value'); grid on;
end

%%

clc;
clearvars;
close all;

rng(42);  % Fixed seed for reproducibility

accepted = 0;
   
guess = [0 0 0 -2 -10 300]';
N = 500;
params_black = zeros(N, 6);  % Preallocate
skipped = 0;                 % Count rejected runs
wmin = -2; wmax = 2; Npoints = 1000;
w = logspace(wmin, wmax, Npoints);

mag_all = NaN(N, length(w));
phase_all = NaN(N, length(w));

theta_mean = zeros(N, 6);
epsilon = 1e-5;              % Relative convergence threshold
stability_counter = 0;       % Consecutive stability tracker
min_stable_iters = 25;       % Required consecutive convergence steps

wait_h = waitbar(0, 'Running Monte Carlo for Structured Black-Box ID...');
start_time = tic;
for i = 1:N
    elapsed = toc(start_time);
    avg_iter_time = elapsed / i;
    est_total_time = avg_iter_time * N;
    time_left = est_total_time - elapsed;

    % Format ETA as mm:ss
    hrs  = floor(time_left / 3600);
    mins = floor(mod(time_left, 3600) / 60);
    secs = mod(round(time_left), 60);
    eta_str = sprintf('%02d:%02d:%02d', hrs, mins, secs);

    % Update waitbar with full-run ETA
    waitbar(i/N, wait_h, ...
        sprintf('Run %d of %d — ETA: %s', i, N, eta_str));
    
    % --- True Parameters ---
    Xu = -0.1068;     Xq = 0.1192;
    Mu = -5.9755;     Mq = -2.6478;
    Xd = -10.1647;    Md = 450.71;

    A = [Xu, Xq, -9.81;
         Mu, Mq,  0;
          0,  1,  0];
    B = [Xd; Md; 0];
    C = eye(3);
    D = zeros(3,1);

    % --- Noise + Delay Settings ---
    noise.Enabler = 1;
    noise.pos_stand_dev      = noise.Enabler * 0.0011;
    noise.vel_stand_dev      = noise.Enabler * 0.01;
    noise.attitude_stand_dev = noise.Enabler * deg2rad(0.33);
    noise.ang_rate_stand_dev = noise.Enabler * deg2rad(1);

    delay.position_filter = 1;
    delay.attitude_filter = 1;
    delay.mixer = 1;

    seed.x = i;
    seed.vx = i + 1;
    seed.theta = i + 2;
    seed.q = i + 3;

    % --- Simulation Setup ---
    parameters_controller;
    load ExcitationM
    SetPoint = [0, 0];
    t = ExcitationM(:,1);
    simulation_time = t(end) - t(1);

    load_system('MonteCarlo');
    set_param('MonteCarlo', "FastRestart", "off");
    simulation = sim('MonteCarlo', 'SrcWorkspace', 'current');
    if exist('slprj', 'dir'); rmdir('slprj', 's'); end

    time = 0:sample_time:simulation_time;
    input = simulation.Mtot;
    true_output = simulation.q;

    % --- PBSID Identification ---
    n = 3; p = 350; f = 300;
    [Aest, Best, Cest, Dest, Kest, Sest] = estimate_ABCDK(input', true_output', f, p, n);
    black_sys = ss(Aest, Best, Cest, Dest, sample_time);

    % --- Structuring Step ---
    freq_vec = 2*pi*logspace(-1, 2, 1000);
    frd_data = frd(black_sys, freq_vec);
    sys_black = idgrey(@longitdynRot, guess, 'c');

    opt = greyestOptions();
    opt.SearchMethod = 'lm';
    opt.SearchOptions.MaxIterations = 100;
    opt.SearchOptions.Tolerance = 1e-6;
    opt.Display = 'off';

    structuredModel = greyest(frd_data, sys_black, opt);
    est_params = getpvec(structuredModel);
    [mag, phase] = bode(structuredModel, w);

    % Extract pitch rate (2nd output)
    mag_pitch = squeeze(20*log10(mag(1,1,:)));  % dB
    phase_pitch = squeeze(phase(1,1,:));        % deg

    mag_all(i, :) = mag_pitch;
    phase_all(i, :) = phase_pitch;

    % --- True and error analysis ---
    true_params = [-0.1068, 0.1192, -5.9755, -2.6478, -10.1647, 450.71];
    labels = {'Xu', 'Xq', 'Mu', 'Mq', 'Xd', 'Md'};
    abs_errors = abs(est_params' - true_params);
    percent_errors = 100 * abs_errors ./ abs(true_params);
    rel_errors = (est_params' - true_params) ./ true_params;
    rel_rmse = sqrt(mean(rel_errors.^2));
    
    % --- Outlier Rejection ---
    param_threshold = 500;        % Max absolute value of any parameter
    rel_rmse_threshold = 100;     % Relative RMSE in percent
    
    if any(abs(est_params) > param_threshold) || rel_rmse * 100 > rel_rmse_threshold
        fprintf('⚠️  Outlier detected — skipping this run (Rel. RMSE = %.2f%%)\n\n', rel_rmse * 100);
        skipped = skipped + 1;
        continue;
    end


    % Save valid run
    accepted = accepted + 1;
    params_black(accepted, :) = est_params';
    
    % Compute running mean (only on accepted samples)
    if accepted == 1
        theta_mean(accepted, :) = est_params';
    else
        theta_mean(accepted, :) = mean(params_black(1:accepted, :), 1);
        delta_mean = norm(theta_mean(accepted, :) - theta_mean(accepted - 1, :)) / norm(theta_mean(accepted, :));
    
        if delta_mean < epsilon
            stability_counter = stability_counter + 1;
        else
            stability_counter = 0;
        end
    
        if stability_counter >= min_stable_iters
            fprintf('\n✅ Parameter mean converged after %d accepted runs (%.2e rel. change for %d steps)\n', ...
                    accepted, delta_mean, stability_counter);
            
            % Trim outputs
            params_black = params_black(1:accepted, :);
            mag_all = mag_all(1:accepted, :);
            phase_all = phase_all(1:accepted, :);
            theta_mean = theta_mean(1:accepted, :);
            black_poles = black_poles(1:accepted, :);
            black_zeros = black_zeros(1:accepted, :);
            break;
        end
    end


    % --- Print summary ---
    fprintf('\n--- State-Space Parameter Comparison (Run %d) ---\n\n', i);
    fprintf('%-10s %15s %15s %15s %15s\n', 'Parameter', 'Estimated', 'True', 'Abs Error', '% Error');
    fprintf('%s\n', repmat('-', 1, 80));
    for j = 1:length(true_params)
        fprintf('%-10s %+15.4f %+15.4f %15.4f %14.2f%%\n', ...
            labels{j}, est_params(j), true_params(j), abs_errors(j), percent_errors(j));
    end

    fprintf('\n%-25s %.2f%%\n', 'Relative Root-Mean-Square Error (RMSE):', rel_rmse*100);

    % --- Save matrices if needed ---
    A_black_id = structuredModel.A;
    B_black_id = structuredModel.B;
    C_black_id = eye(3);
    D_black_id = zeros(3,1);

    black_sys = ss(A_black_id, B_black_id, C_black_id, D_black_id, sample_time);
    black_poles(i, :) = pole(black_sys(2)).';
    black_zeros(i, :) = tzero(black_sys(2)).';
end
elapsed_time = toc;

close(wait_h);
fprintf('\n✅ Monte Carlo complete. Total skipped runs: %d out of %d\n', skipped, N);
fprintf('\n✅ Elapsed Time: %.2f min\n', elapsed_time/60);

figure;
for j = 1:6
    subplot(3,2,j);
    plot(1:size(theta_mean,1), theta_mean(:,j), 'LineWidth', 1.5);
    title(sprintf('Mean Convergence: %s', labels{j}));
    xlabel('Iteration'); ylabel('Value'); grid on;
end

fprintf('\n✅ Monte Carlo complete. Accepted: %d | Skipped: %d | Total: %d\n', ...
        accepted, skipped, accepted + skipped);

%%
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

function saveAllFiguresToPath(userDefinedPath)
% SAVEALLFIGURESTOPATH Saves all open figures to a user-defined path.
%
%   saveAllFiguresToPath(userDefinedPath) saves all open MATLAB figures to the
%   specified userDefinedPath. The function generates filenames as Figure_1.png,
%   Figure_2.png, etc.

% Input Validation
if nargin == 0
    error('User-defined path is required.');
end

if ~ischar(userDefinedPath) || isempty(userDefinedPath)
    error('User-defined path must be a non-empty string.');
end

if ~exist(userDefinedPath, 'dir')
    mkdir(userDefinedPath);
    disp(['Created directory: ' userDefinedPath]);
end

% List all open figures
openFigures = findall(0, 'Type', 'figure');

if isempty(openFigures)
    disp('No open figures found.');
    return
end

% Save figures
counter = 1;
for i = 1:numel(openFigures)
    currentFigure = openFigures(i);
    
    % Generate filename and full file path
    fileName = sprintf('Figure_%d', counter);
    fullFilePath = fullfile(userDefinedPath, [fileName '.png']);
    
    % Save the figure
    try
        saveas(currentFigure, fullFilePath);
        disp(['Figure ' num2str(i) ' saved to: ' fullFilePath]);
    catch ME
        disp(['Error saving figure ' num2str(i) ': ' ME.message]);
    end
    
    counter = counter + 1;
end
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