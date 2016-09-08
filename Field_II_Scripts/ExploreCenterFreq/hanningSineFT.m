% Discover meaning of "center frequency"
% Guess: Hanning-ed sine wave has a FT with a "center frequency" = f0
% Result: it does appear to be the case!

f0 = 25; %6.67e6;            % Central frequency                        [Hz]

t = 0:1/500:10-1/500;  
impulse_response = sin(2*pi*f0*t);
impulse_response = impulse_response.*hanning(length(impulse_response))';

figure
y = fft(impulse_response);            % Compute DFT of x
m = abs(y);                               % Magnitude
f = (0:length(y)-1)*100/length(y);        % Frequency vector
plot(f,m,'o')
title('Frequency Spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude')