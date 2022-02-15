%% Estimate the pitch of a voiced signal
clc, clear all, close all


%% Load the voice signal
% Load the signal voiced_a.wav, then try with unvoiced_sh.wav
[s, Fs] = audioread('voiced_a.wav');
%[s, Fs] = audioread('unvoiced_sh.wav');
% s = randn(10000,1); Fs = 22050;

% Plot the signal
figure(1);
plot(s)
title('Signal (time domain)');


%% Define the LPC parameters
% Set the window length to 25ms and the prediction order to 25
p = 25;                 % prediction order 
fl = 25;                % frame length (ms)

% Set the number of frames to 10
N = 10;                  % number of frames
M = floor(fl*Fs / 1000); % frame length(samples)


%% Loop over all the windows
for n=0:N-1
    
    % Select window and computer auto-correlation
    sn = s(n*M + 1: n*M + M);          % Windowing of the signal 
    r = xcorr(sn);                     % compute the autocorrelation
    
    % Plot autocorrelation
    plot(r);                           % Plot of the autocorrelation function 
    
    % Compute the filter parameters using the function levinson  (see doc
    % levinson)
    a = levinson(r(M:end), p);      %Computation of the parameters of the filter 
    
    % Compute the prediction error using filter or conv
    %e = conv(sn, a);
    e = filter(a, [1], sn);
    En = sum(e.^2);
    G = sqrt(En);
    
    % Find pitch-related peaks
    [pks, locs] = findpeaks(e, 'MinPeakDistance', 100);
    
    % Compute the optimal error Dp as a * R
    %Dp = sum(a.*r(M:M+p)');             % Prediction error
    %G = sqrt(Dp);

    % Compute the shaping filter H using freqz function
    [H, w] = freqz(1, a, floor(M/2)+1); %Interpolation of the FIR filter
    f = w.*Fs/(2*pi);                   %frequency axis
    
    % Compute the DFT of the original signal
    Sn = fft(sn);                       %FFT of the signal
    
    % Plot the Magnitude of the signal spectrum, and on the same graph, the
    % the LPC spectrum estimation (remember the definition of |E(omega)|) 
    figure(3), clf
    subplot(3,1,1)
    plot(f,abs(Sn(1:ceil(length(Sn)/2)))); %Plot of the FFT of the original signal
    hold on;
    title('Input signal (blue line) and envelope (red line)');
    xlabel('Frequency [Hz]');
    plot(f, G*abs(H),'r');        %Spectral matching
    hold off
    
    % Plot the predition error e in time
    subplot(3,1,2), hold on
    t = (0:1:length(e)-1) / Fs;
    plot(t, e);               % plot prediction error (time domain)
    plot(t(locs), pks, 'v');  % plot detected peaks
    
    title('Prediction error (time)');   
    xlabel('Time[s]');
    E = fft(e, 2*length(f));
    
    % Plot the prediction error magnitude spectrum
    subplot(3,1,3)
    plot(f, abs(E(1:end/2)));           %Plot prediction error (frequency domain)
    title ('Prediction error (frequency)');
    xlabel('Frequency [Hz]');
        
    % Use the function pause to stop the for loop and check the plot
    pause();
    
end
