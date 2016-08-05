%BFT_SCAN_PHASED - Define a phased-array sector scan.
%
%USAGE  : bft_scan_phased(xdc, sector, no_lines, depth, no_focal_zones, ...
%                          apodization, t, c)
%
%INPUTS : xdc    - Handle to transducer definition.
%         sector - Size of a sector. If sector > pi, the dimension
%                  of the sector is assumed to be in degrees, otherwise
%                  the sector is assumed to be in radians.
%         no_lines - Number of lines in the sector.
%         depth  - Depth of scan (depth = Tprf * c/2)                    [m]
%         no_focal_zones - Number of focal zones.
%         apodiztion - Matrix with apodization values.                 [0-1]
%         t - Time after which the current apodization is valid          [s]
%         c - Speed of sound
%
%NOTE  : The apodization is will take effect only for normal bemforming.
%        The SUM apodization (See BFT_SUM_APODIZATION) is set to '1'. 
%        This might cause side effects.
%
%OUTPUT : None
%
%VERSION: 1.0, 05 Jan 2001, Svetoslav Nikolov

function bft_scan_phased(xdc, sector, no_lines, depth, no_focal_zones,  apo, t, c)

if nargin ~=8
   error(nargchk(8,8,nargin));
end

if sector > pi,
   sector = sector * pi / 180;
end

no_elements = size(apo,2);

bft_no_lines(no_lines);

theta = (-sector/2:sector/(no_lines-1):sector/2);

dr = depth/ (no_focal_zones-1);
r = [dr:(depth-dr)/(no_focal_zones-1):depth];
tf = 2*r/c - dr/c;


for ii =1:no_lines
  f = r'*[sin(theta(ii)) 0 cos(theta(ii))];
  bft_focus(xdc, tf, f, ii);
  bft_center_focus([0 0 0],ii);
  bft_apodization(xdc,t, apo, ii);
  bft_sum_apodization(xdc, 0, ones(1,no_elements), ii);
end

