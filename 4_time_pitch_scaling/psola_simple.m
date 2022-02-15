%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSOLA                           
% Audio Signals course
% 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear all, close all


%% Read audio track
% Read signal
[y, fs] = audioread('A string.wav');
N = length(y);  % Signal length in samples
Ts = 1/fs;  % Sampling period

% Convert to mono if stereo
y = y(:, 1);

% Plot waveform
t = 0:1/fs:(N-1)/fs;
figure(1)
plot(t, y)
xlabel('time [s]')


%% Detect fundamental frequency from Fourier analysis
% Compute FFT magnitude
Y = abs(fft(y));

% Plot the spectrum
f = linspace(0, fs-N/fs, N);  % Define frequency axis 

figure(2), hold on
plot(f, Y);
xlabel('freq [Hz]')
grid('on')

% Detect fundamental frequency of the harmonic sound
% hint: check 'max' and 'findpeaks' functions 
[Y_max, Y_loc] = findpeaks(Y, 'MinPeakHeight', 0.3 * max(abs(Y)));  % Only peaks higher than 30% of the maximum value
Y_max = Y_max(1);  % Select first peak
Y_loc = Y_loc(1);  % Select first peak location
f_fund = Y_loc * fs / length(Y);  % Compute fundamental frequency

plot(f_fund, Y_max, '*');  % Plot peak position


%% Generate analysis train of pulses
% Compute the period associated to the fundamental frequency
T_a = 1/f_fund;
T_a_samples = round(T_a*fs); 

% Generate a train of pulses according to the computed period
train_a = zeros(N, 1);  % Start from a vector of zeros
train_a(1:T_a_samples:N-T_a_samples) = 1;  % Write 1 every T_a_samples to the end

% Compute pulses locations
idx_a = find(train_a == 1);  % Find position of elements equal to 1


%% Generate synthesis train of pulses
% Choose a pitch scaling factor
alpha = .9;

% Compute the new period associated to the scaling factor
T_b_samples = round(T_a_samples * alpha);

% Generate a train of pulses according to the computed period
train_b = zeros(N, 1);  % Start from a vector of zeros
train_b(1:T_b_samples:N-T_b_samples) = 1;  % Write 1 every T_b_samples to the end

% Compute pulses locations
idx_b = find(train_b == 1);  % Find position of elements equal to 1


%% Plot both trains of pulses
figure(3)
subplot(2,1,1), stem(train_a(1:3000)), title('Analysis peaks')
subplot(2,1,2), stem(train_b(1:3000)), title('Synthesis peaks')


%% Window
% Generate a Hann window with the correct lenght for synthesis
L = T_b_samples * 2 + 1;  % In order to honor overlap and add reconstruction condition
win = hann(L);

% Check windows ovrelapping condition
figure(4)
plot(conv(win, train_b), 'LineWidth', 3)  % All overlapped windows add to 1


%% PSOLA (simplified)
% Allocate memory for synthesis signal
y_out = zeros(N + 2*L, 1);  % Big enough buffer (can be optimized)

% Loop over all synthesis windows
for n = 1:length(idx_b)
    
    % Select synthesis window position
    start_out = idx_b(n) - floor(L/2);  % Start synthesis sample index
    stop_out = idx_b(n) + floor(L/2);   % End synthesis sample index
    
    % Select analysis window
    [~, n_a] = min(abs(idx_a - idx_b(n)));  % Find position of the closest analysis peak given the chosen synthesis one
    
    start_in = idx_a(n_a) - floor(L/2);  % Start analysis sample index
    stop_in = idx_a(n_a) + floor(L/2);   % End analysis sample index
    
    % Copy analysis window to synthesis position
    if start_in >= 0 && start_out >= 0 && stop_in <= length(y) && stop_out <= length(y_out)  % If we are not writing outside of the buffer or reading outside of the analysis signal
        x = y(start_in:stop_in) .* win;  % Select the analysis signal window
        y_out(start_out:stop_out) = y_out(start_out:stop_out) + x;  % Copy it over the correct synthesis part of the buffer using overlapp and add
    end
    
end


%% Play original and pitch-shifted signals
soundsc(y/max(y), fs)  % Original track
pause(6)
soundsc(y_out/max(y_out), fs)  % Pitch shifted track

