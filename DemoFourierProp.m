% Show that the Fourier transform of the autocorrelation
% of a sequence is the magnitude squared of the Fourier transform of the
% sequence

seq = [1 -1 3];
autoCorr = xcorr(seq);
magFTACF = abs(fft(autoCorr));

% Pad with zeroes
% Theorem for transfering correlation to multiplication in frequency
% has circular convolution in it. To make the circular convolution 
% not have any extra terms this is how many zeroes you add.
seqPad = [zeros(1,length(seq)-1) seq ];
spd = abs(fft(seqPad)).^2; %fft(seqPad).*conj(fft(seqPad)) % Spectral power density
iftToGetACF = ifft(spd);

% So why do we have this implied shift in time?
% Also, why does padding zeroes in time work?


% Demo that absolute value of sum of complmentary sequences
% comes out to a constant.
fft1 = fft(xcorr([1 1 1 -1]))
fft2 = fft(xcorr([1 1 -1 1]))
abs(fft1 + fft2)
