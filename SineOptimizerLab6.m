%% waveshaper_final_THD_labstyle_compare_constraint.m
% Two-threshold waveshaper: initial (piecewise slopes) vs optimized design
% EE380 (IITK) — Triangle → Sine waveshaping experiment
% Added constraint flag: if constraint=1, t1 lies along y=x line

clear; close all; clc;

%% PARAMETERS
Fs = 200e3;         % sampling frequency (Hz)
f0 = 250;           % fundamental frequency (Hz)
Vpeak = 4;          % ±4 V triangle
nPeriods = 100;     % integer number of cycles
constraint = 1;     % 0 = free optimization, 1 = enforce t1 along y=x

samples_per_period = round(Fs/f0);
N = nPeriods * samples_per_period;
t = (0:N-1)/Fs;

% sweep ranges for optimization
t1_range = linspace(1.2,1.4,21);
t2_range = linspace(2.6,2.7,21);
Aout_range = linspace(2.5,2.5,9);

%% Input signal
tri = Vpeak * sawtooth(2*pi*f0*t, 0.5);

%% --- THD of input (normalized) ---
THD_in = compute_THD_labstyle(tri, Fs, f0);
fprintf('Input THD (normalized) = %.3f %%\n', THD_in*100);

%% --- INITIAL DESIGN (manual VTC with explicit coordinates) ---
% Define coordinate points directly (from lab manual)
t1x = 1.3;  t1y = 1.3;      % transition point 1 (x,y)
t2x = 2.7;  t2y = 2.3;      % transition point 2 (x,y)
endx = 4.0; endy = 2.6;     % endpoint saturation (x,y)
Aout_init = endy;

% Construct symmetric piecewise VTC using these coordinates
xb_init = [-endx, -t2x, -t1x, 0, t1x, t2x, endx];
yb_init = [-endy, -t2y, -t1y, 0, t1y, t2y, endy];

% Generate output using interpolation through defined points
out_init = interp1(xb_init, yb_init, tri, 'linear', 'extrap');

% Compute THD for this initial design
THD_init = compute_THD_labstyle(out_init, Fs, f0);

% Display results
fprintf('\nInitial Design (piecewise slopes using explicit coordinates):\n');
fprintf('  (t1x,t1y) = (%.3f, %.3f)\n', t1x, t1y);
fprintf('  (t2x,t2y) = (%.3f, %.3f)\n', t2x, t2y);
fprintf('  (endx,endy) = (%.3f, %.3f)\n', endx, endy);
fprintf('  Output THD = %.3f %%\n', THD_init*100);


%% --- OPTIMIZED DESIGN (search) ---
bestTHD = Inf; best_pair=[0 0]; best_Aout=0; best_out=[]; best_map=[];

for Aout = Aout_range
    for t1 = t1_range
        for t2 = t2_range
            if t1 >= t2, continue; end
            [xb,yb] = build_map(t1,t2,Vpeak,Aout,constraint);
            out = interp1(xb,yb,tri,'linear','extrap');
            THD_val = compute_THD_labstyle(out, Fs, f0);
            if THD_val < bestTHD
                bestTHD = THD_val;
                best_pair = [t1 t2];
                best_Aout = Aout;
                best_out = out;
                best_map = [xb; yb];
            end
        end
    end
end

fprintf('\nOptimized Design:\n');
fprintf('  t1 = %.3f V\n  t2 = %.3f V\n', best_pair(1), best_pair(2));
fprintf('  Aout(opt) = %.3f V\n', best_Aout);
fprintf('  Output THD = %.3f %%\n', bestTHD*100);
if constraint
    fprintf('  (Constraint active: y = x at ±t1)\n');
else
    fprintf('  (No constraint applied)\n');
end

%% --- PLOTTING ---
Lshow = 3 * samples_per_period;
[Xmag, f] = compute_fft(tri, Fs);
[Ymag_init, ~] = compute_fft(out_init, Fs);
[Ymag_best, ~] = compute_fft(best_out, Fs);
fmax = 10 * f0;

%% WINDOW 1: Initial Design
figure('Name','Initial Design (Manual VTC with slopes)','Position',[100 100 1200 850]);
subplot(3,2,1);
plot(t(1:Lshow), tri(1:Lshow), 'b', 'LineWidth', 1.3); grid on;
xlabel('Time (s)'); ylabel('V_{in} (V)'); title('Input Triangular Wave');

subplot(3,2,2);
plot(xb_init, yb_init, 'ro-', 'LineWidth', 1.3); grid on; axis equal;
xlabel('V_{in} (V)'); ylabel('V_{out} (V)');
title('Initial Design VTC (piecewise slopes)');

subplot(3,2,3);
plot(t(1:Lshow), out_init(1:Lshow), 'm', 'LineWidth', 1.3); grid on;
xlabel('Time (s)'); ylabel('V_{out} (V)');
title(sprintf('Output Approximated Sine (Initial) (THD=%.3f%%)', THD_init*100));

subplot(3,2,4); axis off;
text(0,0.8,'\bfInitial Design Parameters','FontSize',12);
text(0,0.65,sprintf('t1=(%.2f, %.2f)', t1x, t1y),'FontSize',11);
text(0,0.55,sprintf('t2=(%.2f, %.2f)', t2x, t2y),'FontSize',11);
text(0,0.45,sprintf('End=(%.2f, %.2f)', endx, endy),'FontSize',11);
text(0,0.35,sprintf('Input THD = %.3f %%', THD_in*100),'FontSize',11);
text(0,0.25,sprintf('Output THD = %.3f %%', THD_init*100),'FontSize',11);
text(0,0.10,'(THD computed up to 10th harmonic)','FontAngle','italic','FontSize',9);

subplot(3,2,5);
stem(f/f0, Xmag, 'b', 'Marker', 'none'); grid on; xlim([0 10]);
xlabel('f / f_0'); ylabel('|X(f)|'); title('FFT of Input (Normalized Frequency)');

subplot(3,2,6);
stem(f/f0, Ymag_init, 'm', 'Marker', 'none'); grid on; xlim([0 10]);
xlabel('f / f_0'); ylabel('|Y(f)|'); title('FFT of Output (Initial Design)');
sgtitle('Two-Threshold Waveshaper — Initial Design (piecewise slopes)', 'FontWeight','bold','FontSize',14);

%% WINDOW 2: Optimized Design
figure('Name','Optimized Design (Search)','Position',[100 100 1200 850]);
subplot(3,2,1);
plot(t(1:Lshow), tri(1:Lshow), 'b', 'LineWidth', 1.3); grid on;
xlabel('Time (s)'); ylabel('V_{in} (V)'); title('Input Triangular Wave');

subplot(3,2,2);
plot(best_map(1,:), best_map(2,:), 'ro-', 'LineWidth', 1.3); grid on; axis equal;
xlabel('V_{in} (V)'); ylabel('V_{out} (V)');
title(sprintf('Optimized VTC (t1=%.2fV, t2=%.2fV)', best_pair(1), best_pair(2)));

subplot(3,2,3);
plot(t(1:Lshow), best_out(1:Lshow), 'm', 'LineWidth', 1.3); grid on;
xlabel('Time (s)'); ylabel('V_{out} (V)');
title(sprintf('Output Approximated Sine (Optimized) (THD=%.3f%%)', bestTHD*100));

subplot(3,2,4); axis off;
text(0,0.8,'\bfOptimized Parameters','FontSize',12);
text(0,0.65,sprintf('t1=%.3f V', best_pair(1)),'FontSize',11);
text(0,0.55,sprintf('t2=%.3f V', best_pair(2)),'FontSize',11);
text(0,0.45,sprintf('Aout(opt)=%.3f V', best_Aout),'FontSize',11);
text(0,0.35,sprintf('Input THD = %.3f %%', THD_in*100),'FontSize',11);
text(0,0.25,sprintf('Output THD = %.3f %%', bestTHD*100),'FontSize',11);
if constraint
    text(0,0.15,'Constraint: y=x at ±t1','FontAngle','italic','FontSize',10);
else
    text(0,0.15,'No constraint applied','FontAngle','italic','FontSize',10);
end
text(0,0.08,'(THD computed up to 10th harmonic)','FontAngle','italic','FontSize',9);

subplot(3,2,5);
stem(f/f0, Xmag, 'b', 'Marker', 'none'); grid on; xlim([0 10]);
xlabel('f / f_0'); ylabel('|X(f)|'); title('FFT of Input (Normalized Frequency)');

subplot(3,2,6);
stem(f/f0, Ymag_best, 'm', 'Marker', 'none'); grid on; xlim([0 10]);
xlabel('f / f_0'); ylabel('|Y(f)|'); title('FFT of Output (Optimized Design)');
sgtitle('Two-Threshold Waveshaper — Optimized Design Results','FontWeight','bold','FontSize',14);

%% --- FUNCTIONS ---
function [xb, yb] = build_map(t1, t2, Vp, Aout, constraint)
    xb = [-Vp -t2 -t1 0 t1 t2 Vp];
    phase = (xb / Vp) * (pi/2);
    yb = Aout * sin(phase);
    if constraint == 1
        for k = 1:length(xb)
            if abs(abs(xb(k)) - t1) < 1e-6
                yb(k) = sign(xb(k)) * t1;
            end
        end
    end
end

function THD = compute_THD_labstyle(x, Fs, f0)
    x = x / max(abs(x));
    N = length(x);
    X = fft(x, 4*N);
    Xmag = abs(X(1:2*N)) / N * 2;
    f = (0:2*N-1)*(Fs/(4*N));
    [~, fundBin] = min(abs(f - f0));
    fund_amp = Xmag(fundBin);
    harm_sq = 0;
    for k = 2:10
        fk = k * f0;
        [~, idx] = min(abs(f - fk));
        harm_sq = harm_sq + Xmag(idx)^2;
    end
    THD = sqrt(harm_sq) / fund_amp;
end

function [Xmag, f] = compute_fft(x, Fs)
    N = length(x);
    X = fft(x, 4*N);
    Xmag = abs(X(1:2*N)) / N * 2;
    f = (0:2*N-1)*(Fs/(4*N));
end
