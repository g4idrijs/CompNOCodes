%BFT_FOCUS_2WAY Create a 2way focus time line defined by focal points.
%  These focus settings are relevant only for synthetic aperture imaging.
%  This is the classical monostatic synthetic aperture focusing
%
% USAGE  : bft_focus_2wy(xdc, times, points, line_no)
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
% VERSION: 1.0, Apr 2000, Svetoslav Nikolov

function bft_focus_2way(xdc, times, points, line_no) 

points = points';
if(nargin < 4) line_no = 1; end; 

bft(16, xdc, times, points, line_no);
