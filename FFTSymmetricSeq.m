% Show that the DFT of a symmetric sequence is all real

% Specifically, the sequence must be of the form:
% x_0, x_1, x_2, ..., x_2, x_1
% http://dsp.stackexchange.com/questions/2560/discrete-fourier-transform-symmetry
seq = [0 -1 2 2 -1];
fft(seq)