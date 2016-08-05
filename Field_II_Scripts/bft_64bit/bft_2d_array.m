%BFT_2D_ARRAY - Create a 2D array
%
%USAGE  : xdc = bft_2d_array(no_ele_x, no_ele_y, pitch_x, pitch_y)
%    OR : xdc = bft_2d_array(no_ele_x, no_ele_y, width, height, kerf_x, kerf_y)
%
%INPUT  : no_ele_x - Number of rows in x direction
%         no_ele_y - Number of rows in y direction
%         pitch_x  - Center-to-center distance in x direction
%         pitch_y  - Center-to-center distance in y direction
%         width    - Size of the element in the x direction
%         height   - Size of the element in the y direction
%         kerf_x   - Distance between two elements in the x direction
%         kerf_y   - Distance between two elements in the y direction
%
% The function assumes that:
%                'kerf_x' + 'width' = 'pitch_x'
%                'kerf_y' + 'height' = 'pitch_y'
%
%OUTPUT : xdc - Pointer to the allocated aperture
%
%VERSION : 1.0, Aug 16 2000 Svetoslav Nikolov

function xdc = bft_linear_array(no_ele_x, no_ele_y, width, height, kerf_x, kerf_y)

if ((nargin~=4) & (nargin ~=6))
   error('Wrong number of input arguments')
end

if (nargin == 6)
   pitch_x = width + kerf_x;
   pitch_y = height + kerf_y;
else
   pitch_x = width;
   pitch_y = height;
end

x = [-(no_ele_x-1)/2:(no_ele_x-1)/2] * pitch_x;
y = [-(no_ele_y-1)/2:(no_ele_y-1)/2] * pitch_y;

x = ones(no_ele_y,1) * x;
y = y' * ones(1,no_ele_x);
z = zeros(no_ele_y, no_ele_x);
x = x(:);
y = y(:);
z = z(:);
xdc = bft_transducer([x,y,z]);

