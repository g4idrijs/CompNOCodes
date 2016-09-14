% Investigate the form of the self convolution of a sequence
% of the form a[n]*sin(2*pi/(L-1)*n) for 0 <= n <= L-1
% and 0 for all other n

% I suspect there will be some kind of odd symmetry

% Define sequence length
L = 50;
n = 0:L-1;

% Define number of cycles
m = 3;

% Define the corresponding m-cycle sinusoid
h = sin(2*pi*m/(L-1)*n);
plot(n,h,'o')

% Convolve this sequence with itself, to get the ultrasound system impulse
% response
figure
plot(conv(h,h))
figure
plot(xcorr(h,h))