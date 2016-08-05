%BFT_SUM_APODIZATION Create a summation apodization time line.
%   This function is used in the case that the individual low
%   resolution images must be weighted during the summation
%
% USAGE : bft_sum_apodization(xdc, times, values, line_no)
%
% INPUT : xdc - Pointer to a transducer aperture
%         times - Timea after which the associated apodization is valid
%         values - Apodization values. Matrix with one row for each
%                  time value and a number of columns equal to the 
%                  number of physical elements in the aperture. 
%         line_no - Number of line. If skipped, 'line_no' is assumed
%                   to be equal to '1'
%
% OUTPUT: None
%
%VERSION: 1.0, Feb 11, 2000 Svetoslav Nikolov

function bft_sum_apodization(xdc, times, values, line_no)
if (nargin < 4) line_no = 1; end;
bft(14, xdc, times, values', line_no);

