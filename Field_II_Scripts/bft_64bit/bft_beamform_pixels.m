%BFT_BEAMFORM_PIXELS Beamform image (line) based on pixels definitions
%    This function can be used for synthetic transmit aperture 
%   beamforming. This is enabled by the optional argument ELEMENT_NO
%   If this argument is present, the data in the beamformed image 
%   will be offset at a number of samples, corresponding to ELEMENT_NO.
%   Then the low -resolution images can be summed in MATLAB using '+'.
%
%   The output samples are ordered in the way the user has specified 
%   the focal  points.
%
%USAGE  : pixels = bft_beamform_pixels(time, rf_data, [element_no])
%
%INPUT  : time    - The time of the first sampled value
%         rf_data - The recorded RF data. The number of columns 
%                   is equal to the number of elements.
%         element_no - Number of element, used in transmit.
%       
%OUTPUT : pixels  - Matrix with the beamformed data. The 
%                   number of samples are equal to the number of
%                   desired pixels.
%
%VERSION: 1.0, 04 May 2000, Svetoslav Nikolov

function pixels = bft_beamform_pixels(time, rf_data, element_no)

if (~isa(rf_data,'double')) rf_data = double(rf_data);end;

rf_data = hilbert(rf_data);

if (nargin == 2) 
   rpixels = bft(11, time, real(rf_data));
   ipixels = bft(11, time, imag(rf_data));
else
   rpixels = bft(11, time, real(rf_data), element_no);
   ipixels = bft(11, time, imag(rf_data), element_no);
end

pixels = rpixels + i*ipixels;
