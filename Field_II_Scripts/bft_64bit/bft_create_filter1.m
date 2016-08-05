%BFT_CREATE_FILTER1 Create linear phase low pass filter. Method #1
%    This function creates a low pass filter, with a cutting frequency
%    at half the sampling frequency. The design method is based on 
%    the Matlab function "FIRLS". The result of "FIRLS" is weighted
%    with a Kaiser window.
%
%USAGE : h = bft_create_filter1(Nf, Ntaps,[beta])
%
%INPUTS: Nf - Upsampling ratio                                [Integer]
%        Ntaps - Number of filter taps at the original fs     [Integer]
%        beta - Coefficient used in the creation of the Keiser window.
%               A bigger beta means a higher suppression in stop band.
%               The default value is 5.
%
%
%OUTPUT : h - The impulse response of the filter. 
%
%NOTE  :
%    Beta trades the rejection of the lowpass filter against the transition
%    width from passband to stopband.  Larger Beta means a slower
%    transition and greater stopband rejection.  See Rabiner and Gold
%    (Theory and Application of DSP) under Kaiser windows for more about
%    Beta.  The following table from Rabiner and Gold gives some feel
%    for the effect of Beta:
%
%All ripples in dB, width of transition band = D*N where N = window length
%
%               BETA    D       PB RIP   SB RIP
%               --------------------------------
%               2.120   1.50  +-0.27      -30
%               3.384   2.23    0.0864    -40
%               4.538   2.93    0.0274    -50
%               5.658   3.62    0.00868   -60
%               6.764   4.32    0.00275   -70
%               7.865   5.0     0.000868  -80
%               8.960   5.7     0.000275  -90
%               10.056  6.4     0.000087  -100
%
% The length of the window N is N = Nf*Ntaps-1
%

function h = bft_create_filter1(Nf, Ntaps, beta)

if (nargin<3)
  beta = 5;
end

L = Nf*Ntaps;

fc = 1/2/Nf*0.9;
h = Nf*firls( L-1, [0 2*fc 2*fc 1], [1 1 0 0]).*kaiser(L,beta)' ;
