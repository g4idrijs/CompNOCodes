%BFT_TRANSDUCER Create a new transducer definition. 
%   The transducer definition is necessary for the calculation of 
%   the delays.
%
%USAGE : xdc = bft_transducer(centers)
%
%INPUT : centers - Matrix with the coordinates of the centers of 
%                  the elements. It has 3 columns (x,y,z) and a
%                  number of rows equal to the number of elements.
%                  The coordinates are specified in [m]
%
%OUTPUT: xdc - Pointer to the memory location with the transducer 
%              definition. Do not alter this value !!!
%
%VERSION: 1.0, Feb 10, 2000 Svetoslav Nikolov
function xdc = bft_transducer(centers)
xdc = bft(4,centers');
