%BFT_BEAMFORM_DSTA Beamform a number of scan-lines.
%   This function must be used only when the image is dynamically focused !!
%
%USAGE : bf_lines = bft_beamform(time, rf_data, element_no)
%        bf_lines = bft_beamform(time, rf_data, origin)
%
%INPUTS : time  - The time of the first sample
%         rf_data -  A matrix with channel data. One column per channel
%         element_no - Number of element, if the receive aperture is
%                      the same as the transmitting aperture
%         origin     - Origin of transmission
%
%OUTPUT : bf_lines - A matrix with the beamformed data. 
%                    The number of samples is equal to the number of
%                    rows of the input signal RF_DATA
%
%CREATED : 12 May 2003, Svetoslav Nikolov

function bf_lines = bft_beamform(time, rf_data, element_no)

bf_lines = bft(11, time, rf_data, element_no);
