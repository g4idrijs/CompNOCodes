%BFT_NO_LINES Set the number of lines that will be beamformed in parallel.
%   After calling BFT_INIT, the number of lines that are beamformed in 
%   parallel is 1. If the user wants to beamform a whole image in one 
%   command, he/she must set the number of lines, and then specify the
%   focal zones for each of the lines.
%
%USAGE  : bft_no_lines(no_lines)
%
%INPUT  : no_lines  - Number of lines beamformed in parallel
%
%OUTPUT : None
%
%VERSION: 1.0, Feb 10, 2000 Svetoslav Nikolov

function bft_no_lines(no_lines)
bft(3,no_lines)
