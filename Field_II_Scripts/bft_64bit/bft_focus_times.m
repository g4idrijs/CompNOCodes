%BFT_FOCUS_TIMES Create a focus time line defined by focus delays.
%   The user supplies the delay times for each element.
%
% USAGE  : bft_focus_times(xdc, times, delays, line_no)
%
% INPUT  : xdc - Pointer to a transducer aperture.
%          times - Time after which the associated delay is valid
%          delays - Delay values. Matrix with one row for each
%                   time value and a number of columns equal to the
%                   number of physical elements in the aperture.
%          line_no - Number of line. If skipped, 'line_no' is 
%                    assumed to be equal to 1.
% OUTPUT : None
%
% VERSION: 1.0, Feb 11, 2000 Svetoslav Nikolov

function bft_focus_times(xdc, times, delays, line_no)
if (nargin < 4) line_no = 1; end;
bft(8, xdc, times, delays', line_no);
