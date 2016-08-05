%BFT_DYNAMIC_FOCUS Set dynamic focusing for a line
%
% USAGE : bft_dynamic_focus(xdc, dir_xz, dir_zy, line_no)
% 
% INPUT : xdc     - Pointer to the transducer aperture
%         dir_zx  - Direction (angle) in radians for the dynamic
%                   focus. The direction is taken from the center for
%                   the focus of the transducer in the z-x plane.
%         dir_zy  - Direction (angle) in radians for the dynamic
%                   focus. The direction is taken from the center for
%                   the focus of the transducer in the z-y plane.
%         line_no - Number of line. If skipped, 'line_no' is assumed
%                   to be equal to '1'
%
% OUTPUT : None
% 
% VERSION: 1.0, Feb 11, 2000 Svetoslav Nikolov 
function bft_dynamic_focus(xdc, dir_xz, dir_yz, line_no)
if (nargin < 4) line_no = 1; end;
bft(10, xdc, dir_xz, dir_yz, line_no);

