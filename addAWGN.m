% Add White Gaussian Noice function
% Signal-to-noise ratio (often abbreviated SNR or S/N) is a measure to quantify how 
% much a signal has been corrupted by noise. 
% It is defined as the ratio of signal power to the noise power corrupting the signal. 
% A ratio higher than 1:1 indicates more signal than noise.


function out_signal = addAWGN(signal, targetSNR)
sigLength = length(signal); % length
awgnNoise = randn(size(signal)); % orignal noise
pwrSig = sqrt(sum(signal.^2))/sigLength; % signal power
pwrNoise = sqrt(sum(awgnNoise.^2))/sigLength; % noise power
if targetSNR ~= 0
   scaleFactor = (pwrSig/pwrNoise)/targetSNR; %find scale factor
   awgnNoise = scaleFactor*awgnNoise; 
   out_signal = signal + awgnNoise; % add noise
else
   out_signal = awgnNoise; % noise only
end