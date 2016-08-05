%BFT_BEAMFORM Beamform a number of scan-lines.
%   The number of the simultaneously beamformed scan-lines is set
%   by BFT_NO_LINES. If BFT_NO_LINES is not called, only one scan
%   line will be beamformed. 
%    This function can be used for synthetic transmit aperture 
%   beamforming. This is enabled by the optional argument ELEMENT_NO
%   If this argument is present, the data in the beamformed image 
%   will be offset at a number of samples, corresponding to ELEMENT_NO.
%   Then the low -resolution images can be summed in MATLAB using '+'.
%     A second approach is to omit ELEMENT_NO, and to sum the low 
%   resolution images using the command BFT_SUM_IMAGES.
%   
%   For normal beamforming ELEMENT_NO must be skipped.
%
%
%USAGE  : bf_lines = bft_beamform(time, rf_data, [element_no])
%
%INPUT  : time    - The time of the first sampled value
%         rf_data - The recorded RF data. The number of columns 
%                   is equal to the number of elements.
%         element_no - Number of element, used in transmit.
%       
%OUTPUT :bf_lines - Matrix with the beamformed data. The number 
%                   of rows of 'bf_lines' is equal to the number 
%                   of rows of 'rf_data'. The number of columns 
%                   is equal to the number of lines
%
%VERSION: 2.0, 17 Apr 2000, Svetoslav Nikolov

%VERSION: 1.0, 11 Feb 2000, Svetoslav Nikolov

function bf_lines = bft_beamform(time, rf_data, element_no) 

if (~isa(rf_data,'double')) rf_data = double(rf_data);end;

if nargin == 2,
  bf_lines = bft(11, time, rf_data);
else 
  bf_lines = bft(11, time, rf_data);
  dim = size(bf_lines);
  dummy = zeros(dim(1),dim(2));
  bf_lines = bft(13, dummy, bf_lines, element_no, time);
end  
