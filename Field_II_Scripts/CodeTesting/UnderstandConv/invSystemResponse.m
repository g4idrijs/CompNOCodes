close all
clear

% Defining ultrasound system impulse response

% Define sequence length
L = 34;

% Define number of cycles
m = 2;

% Define the corresponding m-cycle sinusoid
n = 0:L-1;
h = sin(2*pi*m/(L-1)*n);

% Impulse response of (normalized for reflection) ultrasound system
H = conv(h,h);
invH = ifft(1./fft(H));
plot(invH);

figure
plot(conv(H,invH))
