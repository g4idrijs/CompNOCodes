% Given two focal regions deltaF apart at depth L
% and codes of length N, sampling frequency fs, and speed of sound c
% calculates:
% time lag between return of two codes to a single transducer
% divided by the time it takes for a whole code to return
% (what fraction of a code length are the codes offset in time)
L = 40e-3;      % focal zone depth (m)
deltaF = 12e-3;  % focal zone x spacing (m)
c = 1540;       % speed of sound (m/s)
fs = 100e6;     % sampling frequency (Hz)
N = 100;        % code length (elements)

Tdelta =  (sqrt(L^2 + deltaF^2) - L)/c
Tcode = N/fs
TfracDelay = Tdelta/Tcode