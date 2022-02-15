% Voiced/unvoiced binary classification
clc, clear all, close all hidden


%% 1) Load the file 'voiced.wav'
[x_v, Fs] = audioread('voiced.wav');


%% 2) Define all windowing parameters for window extraction
%%    length = 40 ms
%%    spacing = 10 ms
% parameters in samples
frame_length = floor(40e-3 * Fs);
frame_spacing = floor(10e-3 * Fs);

% generate Hamming window
win = hamming(frame_length);


%% 3) Extract ZCR and audio power for each window of the voiced sound
% number of windows
N = floor((length(x_v) - frame_length)/frame_spacing) + 1; % number of frames

% container for fetures
f_v = zeros(N, 2);

% loop over windows of x_v
for n=1:N
    
    % extract one audio window
    frame = x_v((n-1)*frame_spacing+1 : (n-1)*frame_spacing+frame_length); % select the frame
    frame = frame .* win; % apply Hamming window
    
    % compute features
    f_v(n,1) = sum(abs(diff(frame>0)))/frame_length;
    f_v(n,2) = sum(frame.^2);
end


%% 4) Plot in time on three subplots: the waveform; the ZCR; the AP
t_waveform = (0:1:length(x_v)-1) / Fs;
t_windows = (0:1:N-1) / Fs;

figure(1), clf, hold on
subplot(311), plot(t_waveform, x_v), grid on, xlabel('time [s]'), ylabel('x_v(t)')
subplot(312), plot(t_windows, f_v(:, 1)), grid on, xlabel('time [s]'), ylabel('ZCR(t)')
subplot(313), plot(t_windows, f_v(:, 2)), grid on, xlabel('time [s]'), ylabel('AP(t)')
