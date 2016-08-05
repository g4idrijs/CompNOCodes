%BFT_LINEAR_ARRAY - Create a linear array.
%
%USAGE  : xdc = bft_linear_array(no_elements, width,  kerf)
%   OR  : xdc = bft_linear_array(no_elements, pitch)
% 
%INPUT  : no_elements - Number of elelements in the array
%         pitch - Distance between the centers of two elements [m]
%         width - Width in x-direction                         [m]
%         kerf  - Distance between two elements                [m]
%
%   The function assumes that 'kerf' + 'width' = 'pitch'
%
%OUTPUT : xdc - Pointer to the allocated aperture
%
%VERSION: 1.0, Feb 14, 2000 Svetoslav Nikolov

function xdc = bft_linear_array(no_elements, width, kerf)

if (nargin ==3)
   pitch = width + kerf;
else
   pitch = width;
end

x = [-(no_elements-1)/2:(no_elements-1)/2]' * pitch;
y = zeros(no_elements,1);
z = zeros(no_elements,1);
xdc = bft_transducer([x, y, z]);

