%BFT_CENTER_FOCUS Set the center focus point for the focusing
%   This point is used as a reference for calculating the
%   focusing delay times and as a starting point for dynamic
%   focusing.
%
% USAGE : bft_center_focus(point, line_no)
%
% INPUT : point - The center point [x,y,z]              [ m ]
%         line_no - Number of line. If omitted in the parameter
%                   list 'line_no' is assumed equal to 1
%
% OUTPUT: None
%
% VERSION: 1.0, Feb 11, 2000 Svetoslav Nikolov

function bft_center_focus(point, line_no)

[m n] = size(point);

if (nargin < 2) line_no = 1; end;
if (n>m) point = point'; end;

bft(6, point,line_no);
