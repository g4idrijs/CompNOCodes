%BFT_CONVEX_ARRAY  Create a convex array transducer
%
%USAGE  : xdc = bft_convex_array(no_elements, width, kerf, Rconvex)
%    OR : xdc = bft_convex_array(no_elements, pitch, Rconvex)
%
%INPUTS : no_elements - Number of elements in the convex array
%         width       - Width of one element
%         kerf        - Distance between 2 elements
%         Rconvex     - Convex radius
%         pitch       - Distance between the centers of two elements
% 
%OUTPUT : xdc - pointer to an array structure
%

function xdc = bft_convex_array(no_elements, width, kerf, Rconvex)

if nargin == 4
   pitch = width + kerf;
else
   pitch = width;
   Rconvex = kerf;
end

d_theta = pitch / Rconvex;
theta = [-(no_elements-1)/2:(no_elements-1)/2]*d_theta;
x = Rconvex * sin(theta);
z = Rconvex * cos(theta);
z = z-max(z);
y = zeros(1,no_elements);
xdc = bft_transducer([x', y', z']);



