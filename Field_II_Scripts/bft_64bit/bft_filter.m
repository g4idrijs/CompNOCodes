%BFT_FILTER Set a low pass filter, used for the delays in the beamforming.
%   The BFT uses the filter by upsampling the signal to a new frequency
%   fs1 = fs * Nf, and then picking the necessary sample. The length
%   of the filter must be Nf*Ntaps - 1, where Ntaps is the number of
%   samples in the original sample used to create a new sample in the 
%   upsampled one.
%
%USAGE  : bft_filter(Nf, Ntaps, h)
%
%INPUTS : Nf - Ratio between the new and  old sampling frequencies [Integer]
%         Ntaps - Number of samples from the original signal, used to 
%                 create one new sample                            [Integer]
%         h - The impulse response of the filter.
%             length(h) == Nf*Ntaps
%
%
%OUTPUT : None
%
%VERSION : 1.0 05 Sep 2000, Svetoslav Nikolov

function bft_filter(Nf, Ntaps, h)

bft(18, Nf, Ntaps, h)
