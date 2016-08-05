%BFT_FOCUS_PIXEL Set the coordinates of the focal pixels
%    This type of focusing is meant to be based on pixels not
%    on lines. Therefore the term "focus-time-line" is non
%    valid. The user will get back as many focused samples as 
%    the number of points he/she has passed to this function.
%
%USAGE   : bft_focus_pixel(xdc, points, line_no)
%
%UNPIT   : xdc - Pointer to aperture.
%          points - Focus points. Vector with three columns (x,y,z) 
%                   and one row for each field point.
%          line_no - Number of line for which we set the focus. If 
%                    skipped, 'line_no' is assumed equal to '1'.
%
% OUTPUT : none
%
% VERSION: 1.0, May 2000, Svetoslav Nikolov

function bft_focus_pixel(xdc, points, line_no)

points = points';
if (nargin < 3) line_no = 1; end;
bft(17, xdc, points, line_no);
