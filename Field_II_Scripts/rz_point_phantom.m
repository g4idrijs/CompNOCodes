function [positions, amp] = rz_point_phantom(dz, z_start, Npoints)
% Create a computer model for a point phantom
%
% Calling [positions, amp] = pts_pha;
%
% Output:   positions - positions of the scatterers.
%           amp       - amplitude of the scatterers.
% Zemp 9 Mar 2011

% Create the point scatterer positions

positions = [zeros(1, Npoints); zeros(1, Npoints); ([1:Npoints]-1)*dz+z_start]';
amp = ones(Npoints, 1);
