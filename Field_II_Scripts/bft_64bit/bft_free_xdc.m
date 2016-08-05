%BFT_FREE_XDC Free the memory allocated for a transducer definition
%
% USAGE  : bft_free_xdc(xdc)
%
% INPUT  : xdc - Pointer to the memory location returned by the
%                function BFT_TRANSDUCER
%
% OUTPUT : Nothing
%
%VERSION : 1.0, Feb 11, 2000 Svetoslav Nikolov

function bft_free_xdc(xdc)
bft(5,xdc)
