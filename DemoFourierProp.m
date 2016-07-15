% Show that the Fourier transform of the autocorrelation
% of a sequence is the magnitude squared of the Fourier transform of the
% sequence
%% Autocorrelation section

seq = [-1 -1]; % Original sequences

autoCorr = xcorr(seq)

% Tarek likes to trim the autocorrelation (good for US processing?)
% a = autoCorr(length(aseqPad)/2:length(seqPad)/2+length(seqPad)-1);
magFTACF = abs(fft(autoCorr))

%%
% Spectral power density
% Pad with zeroes (so circular convolution doesn't introduce extra terms)
% (do this earlier?)
seqPad = [zeros(1,length(seq)-1) seq ]; 
spd = abs(fft(seqPad)).^2 %fft(seqPad).*conj(fft(seqPad)) % Spectral power density
iftToGetACF = ifft(spd)

% So why do we have this implied shift in time?
% Also, why does padding zeroes in time work?


