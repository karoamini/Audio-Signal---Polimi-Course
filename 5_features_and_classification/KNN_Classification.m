% Voiced/unvoiced binary classification
clc, clear all, close all hidden


%% 1) Load the files 'voiced.wav', 'unvoiced.wav' and 'test'.wav
[x_v, Fs] = audioread('voiced.wav');
[x_u] = audioread('unvoiced.wav');
%[x_t] = audioread('test.wav');
[x_t] = audioread('test_long.wav');


%% 2) Define all windowing parameters for window extraction
%%    length = 40 ms
%%    spacing = 10 ms
% parameters in samples
frame_length = floor(40e-3 * Fs);
frame_spacing = floor(10e-3 * Fs);

% generate Hamming window
win = hamming(frame_length);


%% 3) Extract ZCR and audio power for each window of the training voiced sound
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
    f_v(n,1) = sum(abs(diff(frame>0))) / frame_length;
    f_v(n,2) = sum(frame.^2);
end


%% 4) Repeat for the training unvoiced sound
% number of windows
N = floor((length(x_u) - frame_length)/frame_spacing) + 1; % number of frames

% container for fetures
f_u = zeros(N, 2);

% loop over windows of x_v
for n=1:N
    
    % extract one audio window
    frame = x_u((n-1)*frame_spacing+1 : (n-1)*frame_spacing+frame_length); % select the frame
    frame = frame .* win; % apply Hamming window
    
    % compute features
    f_u(n,1) = sum(abs(diff(frame>0)))/frame_length;
    f_u(n,2) = sum(frame.^2);
end


%% 5) Plot training features in the feature space using two different colors
figure(1), hold on
plot(f_v(:, 1), f_v(:, 2), 'vg')
plot(f_u(:, 1), f_u(:, 2), '*r')
legend('voiced', 'unvoiced')
xlabel('ZCR')
ylabel('AP')
grid('on')


%% 6) Build training (feature, label) database for KNN classification
% concatenate features
f_db = [f_v; f_u];

% concatenate labels
y_db = [zeros(length(f_v), 1); ones(length(f_u), 1)];

% plot in three graphs: ZCR, AP, labels
figure(2)
subplot(311), plot(f_db(:, 1), 'o'), grid on, xlabel('sample number'), ylabel('ZCR'), axis tight
subplot(312), plot(f_db(:, 2), 'o'), grid on, xlabel('sample number'), ylabel('AP'), axis tight
subplot(313), plot(y_db, 'o'), grid on, xlabel('sample number'), ylabel('label'), axis tight


%% 7) Extract features for the test sound
% number of windows
N = floor((length(x_t) - frame_length)/frame_spacing) + 1; % number of frames

% container for fetures
f_t = zeros(N, 2);

% loop over windows of x_v
for n=1:N
    
    % extract one audio window
    frame = x_t((n-1)*frame_spacing+1 : (n-1)*frame_spacing+frame_length); % select the frame
    frame = frame .* win; % apply Hamming window
    
    % compute features
    f_t(n,1) = sum(abs(diff(frame>0)))/frame_length;
    f_t(n,2) = sum(frame.^2);
end


%% 8) Add test features to the plot
figure(3), hold on
plot(f_v(:, 1), f_v(:, 2), 'vg')
plot(f_u(:, 1), f_u(:, 2), '*r')
plot(f_t(:, 1), f_t(:, 2), '+b')
legend('voiced', 'unvoiced', 'test')
xlabel('ZCR')
ylabel('AP')
grid('on')



%% 9) Classify each window of the test audio using KNN (K = 3)
K = 3;

% initialize contained for estimated labels
y_hat = zeros(length(f_t), 1);

% loop over all windows
for f_idx = 1:length(f_t)
    % select feature vector for this window
    f = f_t(f_idx, :);
    
    % compute distance between test feature and all training features
    distance = sqrt(sum((f_db - f).^2, 2));
    
    % sort distances and select classes from K closest features
    [val, idx] = sort(distance);
    classes = y_db(idx(1:3));
    
    % apply majority voting to estimate the class of the test window
    y_hat(f_idx) = mode(classes);
end


%% 10) Plot test waveform, features, and estimated class in time
figure(4)
subplot(411), plot(x_t), grid on, xlabel('sample number'), ylabel('x_t(n)'), axis tight
subplot(412), plot(f_t(:, 1), 'o'), grid on, xlabel('sample number'), ylabel('ZCR'), axis tight
subplot(413), plot(f_t(:, 2), 'o'), grid on, xlabel('sample number'), ylabel('AP'), axis tight
subplot(414), plot(y_hat, 'o'), grid on, xlabel('sample number'), ylabel('estimated label'), axis tight





