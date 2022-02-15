%% Example for cross synthesis of two sounds
clc, clear all, close all


%% Load data
[exc, FS] = audioread('moore_guitar.wav'); % sound 1: excitation
env = audioread('Toms_diner.wav'); % sound 2: spectral env.


%% Define Hanning windowing parameters
long     = 400;    % window length for calculation of coefficients
hopsize  = long/2; % hop size
w        = hanning(long); % window


%% Define LPC parameters
order_env   = 20; % order of the LPC
order_exc   = 9;  % order for the excitation


%% Trim audio samples to the same lenght and normalize amplitude
ly = min(length(exc), length(env));

N_frames = floor( (ly - long)/hopsize ) + 1;
ly = (N_frames-1)*hopsize + long;

exc = exc(1:ly);
exc = exc/max(abs(exc));
env = env(1:ly);
env = env/max(abs(env));


%% Initialize output vector
out      = zeros(ly,1); % result sound (note: use of time difference equation
% which returns a vector of the same size of the input signal


%% Loop over each windoe and perform cross-synthesis
for j=1:N_frames
    
    % Select analysis windows and perform auto correlation
    k = hopsize*(j-1); % offset of the buffer
    r_env = xcorr(env(k+1:k+long).*w);
    r_exc = xcorr(exc(k+1:k+long).*w);
    
    % Estiamte LPC filter coefficients
    [A_env, Dp_env] = levinson(r_env(long:end), order_env);
    [A_exc, Dp_exc] = levinson(r_exc(long:end), order_exc);
    
    % Compute LPC residuals
    e_env = filter(A_env, 1, env(k+1:k+long));
    e_exc = filter(A_exc, 1, exc(k+1:k+long));
    
    % Synthesize output signal and write into buffer
    out(k+1:k+long) = out(k+1:k+long) + w.*filter(1,A_env,e_exc*sqrt(Dp_env/Dp_exc));
end


%% Normalize output signal amplitude
out_norm = .99* out/max(abs(out)); % scale for wav output


%% Listen to the result
display('press a key to listen the excitation signal');
pause
soundsc(exc, FS);

display('press a key to listen the vocal (envelope) signal');
pause
soundsc(env, FS);

display('press a key to listen the cross-synthesis resulting signal');
pause
soundsc(out_norm, FS);
