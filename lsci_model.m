timeRange = 60; % sec
rate = 50; % sample per sec
heartRate = [100, 420]; % beat per minutes
respiratoryRate = [10, 150]; % beat per minutes
respiratoryEffort = [50, 10]; % modulation of resp. peaks

[signal,  hBeats, rBeats] = lsci_model_fnc(timeRange, rate, heartRate, respiratoryRate, respiratoryEffort);

plot([signal; hBeats; 2*rBeats]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LSCI based respiratory and heart beat model      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signal, hBeats, rBeats] = lsci_model_fnc(timeRange, rate, heartRate, respiratoryRate, respiratoryEffort)

N = timeRange * rate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% heart pulses
t = linspace(0, timeRange-1/rate, N) / timeRange;
hBeatsClean = vco(t*2-1, heartRate/60, rate);
hBeats = hBeatsClean > 0;
hBeats = [0, diff(hBeats, 1)];
hBeats = hBeats == -1;
hBeats = imdilate(hBeats, ones(1, round(rate/12))); % beats are fixed to ~0.08 sec

% respiration pulses
rBeatsClean = vco(t*2-1, respiratoryRate/60, rate);
rBeatsClean = rBeatsClean > 0;
rBeatsClean = [diff(rBeatsClean, 1), 0];
rBeatsClean = rBeatsClean == 1;
rBeatsClean = imdilate(rBeatsClean, ones(1, round(rate/10))); % beats are fixed to ~0.2 sec

% create baseline, which corresponds to heartrate
baseLine = linspace(heartRate(1)/350*10, heartRate(2)/350*10, N);

% heartbeat attenuated and frequency filterred
Fs = rate;  % Sampling Frequency
Ord  = 8;     % Order
Fc = 250/rate;  % Cutoff Frequency based on rSPG SNR measurements
h  = fdesign.lowpass('N,F3dB', Ord, Fc, Fs);
Hd = design(h, 'butter');
hBeatsReal = filter(Hd, hBeatsClean * 0.5);

% respiration filtering
Fs = rate;  % Sampling Frequency
Ord  = 4;     % Order
Fc = 100/rate;  % Cutoff Frequency based on speckle measurements
h  = fdesign.lowpass('N,F3dB', Ord, Fc, Fs);
Hd = design(h, 'butter');

rBeatsReal = filter(Hd, double(rBeatsClean));
rBeatsReal = rBeatsReal .* linspace(respiratoryEffort(1),respiratoryEffort(2), N);
rBeatsReal = abs(rBeatsReal); % mimicing the fact, that LSCI cannot produce negative values

% measurement noise
mNoise = 0.1 * randn(1,N);

% motion artifacts
% at lower breathrate less motion is expected
K = round(timeRange); % ~1s periods of movements
mArt1 = rand(1,K) < (sqrt(linspace(respiratoryRate(1) / 100, respiratoryRate(2) / 100, K))-1);
mArt1 = interp1(1:K, 5*double(mArt1), linspace(1,K,N), 'pcip');
mArt2 = rand(1,K) < (sqrt(linspace(respiratoryRate(1) / 100, respiratoryRate(2) / 100, K))-1);
mArt2 = interp1(1:K, 3*double(mArt2), linspace(1,K,N), 'pcip');
mArt = mArt1 + mArt2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sum up
signal = rBeatsReal + hBeatsReal;
signal = signal + mArt;
signal = baseLine + signal + mNoise;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output ref peaks
hBeats = hBeats > 0;

% correction for length and filter delay
rBeatsClean = imdilate(rBeatsClean, ones(1, round(rate/5))); 
rBeats = [zeros(1,11) rBeatsClean(1:end-11) > 0];

end