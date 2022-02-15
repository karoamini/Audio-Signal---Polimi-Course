%% Time-variant fractional delay line
clc, clear all, close all


%% Load the imput signal
[x, Fs] = audioread('audio.wav');
N = length(x);  % number of samples

% Convert to mono
x = mean(x, 2);

% Plot
figure(1)
t = (1:1:N)/Fs;  % define time axis
plot(t, x), xlabel('time [s]'), ylabel('x(s)'), grid('on')

% Play
% soundsc(x, Fs)


%% Aplly time-variant fractional delay (linear interpolation)
% Define parameters
delay_nominal_sec = 0.1;  % delay in seconds
delay_modulant_sec = 0.0005;  % delay in seconds

% Define modulated delay function
freq_delay = 10;  % frequency Hz
delay_fun_sec = @(n) delay_nominal_sec + delay_modulant_sec * cos(2*pi*freq_delay*n/Fs);


%% Initialize buffer and output
% Compute maximum delay
max_delay_sec = delay_nominal_sec + delay_modulant_sec;
max_delay = ceil(max_delay_sec * Fs);

% Define buffer
buff_len = 2*max_delay + 1;  % buffer size
buff = zeros(buff_len, 1);  % buffer

% Initialize output signal
y = zeros(N, 1);  % initialize output signal


%% Loop over all samples applying delay
for n = 1:N
    
    % Update buffer
    buff = [x(n); buff(1:buff_len-1)];  % update the buffer

    % Compute delay in samples
    delay_sec = delay_fun_sec(n);
    delay = delay_sec * Fs;  % delay in samples
    delay_int = floor(delay);  % integer part of the delay in samples
    delay_frac = delay - delay_int;  % fractional part of the delay

    % Apply delay
    y(n) = delay_frac * buff(delay_int) + (1-delay_frac) * buff(delay_int + 1);  % interpolate samples from buffer 
end


%% Plot
figure()
subplot(211), plot(t, x), title('Input'), xlabel('time [s]'), ylabel('x(s)'), grid('on')  % original signal
subplot(212), plot(t, y), title('Delayed'), xlabel('time [s]'), ylabel('y(s)'), grid('on')  % delayed signal


%% Play
%soundsc(x, Fs);
%pause();
%soundsc(y, Fs);

