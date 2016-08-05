%BFT_FOCUS Create a focus time line defined by focal points.
%
% USAGE  : bft_focus(xdc, times, points, line_no)
%
% INPUT  : xdc - Pointer to aperture.
%          times -  Time after which the associated focus is valid
%          points - Focus points. Vector with three columns (x,y,z) 
%                   and one row for each field point.
%          line_no - Number of line for which we set the focus. If 
%                    skipped, 'line_no' is assumed equal to '1'.
%
% OUTPUT : none
%
% VERSION: 1.0, Feb 2000, Svetoslav Nikolov

function bft_focus(xdc, times, points, line_no) 

points = points';
if(nargin < 4) line_no = 1; end; 

bft(7, xdc, times, points, line_no);
