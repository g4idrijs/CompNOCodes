%BFT_TRANSDUCER_SET Set new coordinates for the transducer
%
%USAGE  : bft_transducer_set(xdc, centers)
%
%INPUTS : xdc - Handle to existing transducer definition
%         centers -  Matrix with the coordinates of the centers of 
%                  the elements. It has 3 columns (x,y,z) and a
%                  number of rows equal to the number of elements.
%                  The coordinates are specified in [m]
%
%OUTPUTS: None
%
%VERSION: 1.0 May 12, 2003, Svetoslav Nikolov
%
function bft_transducer_set(xdc, centers)
bft(21, xdc, centers')
