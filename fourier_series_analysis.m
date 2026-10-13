clear; close all; clc;

%% PART 1
T  = 8.0;      % period 
T1 = 2.0;      % half-width from PDF
D  = 2*T1/T;   % duty ratio = 0.5
fs = 100;     % sampling frequency
w0 = 2*pi/T;
% Multi-period timeline
t  = 0:1/fs:4*T-1/fs;
X1 = square(w0*t, D*100);

% One period for FS coefficient computation
t_pere  = 0:1/fs:T-1/fs;
X1_pere = square(w0*t_pere, D*100);

% for N=25 and N=75
[k1_25, Ck1_25] = fs_coeffs(t_pere, X1_pere, T, 25);
[k1_75, Ck1_75] = fs_coeffs(t_pere, X1_pere, T, 75);

% Reconstructions
x1_rec_25 = fs_reconstruct_Ck(t, T, k1_25, Ck1_25);
x1_rec_75 = fs_reconstruct_Ck(t, T, k1_75, Ck1_75);
% Plots
figure;
subplot(2,1,1);
plot(t, X1, 'LineWidth',1.2); hold on;
plot(t, x1_rec_25, '--'); grid on;
title('Part 1: x1(t) vs Reconstruction (N=25)');
legend('Original','FS N=25');

subplot(2,1,2);
plot(t, X1, 'LineWidth',1.2); hold on;
plot(t, x1_rec_75, '--'); grid on;
title('Part 1: x1(t) vs Reconstruction (N=75)');
legend('Original','FS N=75');


%% PART 3:
x2 = X1;           
kvec = k1_25;     
Ck2 = Ck1_25;      
% Derivative in time
x3 = [diff(x2), 0]*fs;
% FS property: multiply Ck by (j*k*w0)
Ck3 = 1j*kvec*w0 .* Ck2;

% Plots (3x2)
figure;
subplot(3,2,1); plot(t,x2); grid on; title('x2(t)');
subplot(3,2,2); plot(t,x3); grid on; title('x3(t) = d/dt x2(t)');
subplot(3,2,3); stem_fs_spectrum(kvec,Ck2,'mag','x2 (N=25)');
subplot(3,2,4); stem_fs_spectrum(kvec,Ck3,'mag','x3 (N=25)');
subplot(3,2,5); stem_fs_spectrum(kvec,Ck2,'phase','x2 (N=25)');
subplot(3,2,6); stem_fs_spectrum(kvec,Ck3,'phase','x3 (N=25)');

%% PART 4
f1 = 100; f2 = 50;
fs_cases = [1000, 200, 50];

figure;
for i=1:3
    fs_i = fs_cases(i);
    t_i = 0:1/fs_i:1-1/fs_i;
    x4 = cos(2*pi*f1*t_i) + cos(2*pi*f2*t_i);

    subplot(3,2,2*i-1);
    plot(t_i,x4); grid on; title(sprintf('Time (fs=%d Hz)',fs_i));

    X = fft(x4); Nfft = length(X);
    f = (0:Nfft-1)*(fs_i/Nfft);
    P2 = abs(X)/Nfft;
    P1 = P2(1:floor(Nfft/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);

    subplot(3,2,2*i);
    stem(f(1:numel(P1)),P1,'filled'); grid on;
    title(sprintf('Spectrum (fs=%d Hz)',fs_i));
end

%% helper functions
function [kvec, Ck] = fs_coeffs(t, x, T, N)
    w0 = 2*pi/T; dt = mean(diff(t));
    kvec = -N:N;
    Ck = zeros(size(kvec));
    for ii = 1:length(kvec)
        k = kvec(ii);
        Ck(ii) = (1/T) * sum(x .* exp(-1j*k*w0*t)) * dt;
    end
end

function xrec = fs_reconstruct_Ck(t, T, kvec, Ck)
    w0 = 2*pi/T;
    xrec = zeros(size(t));
    for ii = 1:length(kvec)
        xrec = xrec + Ck(ii) .* exp(1j*kvec(ii)*w0*t);
    end
    xrec = real(xrec);
end

function stem_fs_spectrum(kvec, Ck, kind, ttl)
    if nargin < 4, ttl = ''; end
    switch lower(kind)
        case 'mag'
            stem(kvec, abs(Ck), 'filled'); grid on;
            xlabel('k'); ylabel('|C_k|'); title(['Magnitude ', ttl]);
        case 'phase'
            stem(kvec, angle(Ck), 'filled'); grid on;
            xlabel('k'); ylabel('\angle C_k (rad)'); title(['Phase ', ttl]);
    end
end
