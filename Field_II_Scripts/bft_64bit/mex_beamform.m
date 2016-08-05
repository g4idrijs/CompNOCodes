%MEX_BEAMFORM Compile the beamforming library for matlab
% USAGE  : mex_beamform(debug)
% 
% INPUT  : debug - if 1 "DEBUG" is defined. More messages are printed
%                     0  release version. Most of the messages are omitted
% OUTPUT : Nothing
%
function cmd = mex_beamform(debug)
if (nargin > 0)
  if (debug > 0)
    debug = '-DDEBUG';
  else
    debug = '';
  end
else
  debug = '';  
end
file_names = ['c/mex_beamform.c c/focus.c c/beamform.c c/geometry.c c/transducer.c c/motion.c'];
host = computer;
if (strcmp(host,'PCWIN') || strcmp(host,'PCWIN64'))
   cmd = ['mex ' debug ' -D__MSCVC_' ' -O -output bft ' file_names];
else
   cmd = ['mex ' debug ' -O -output bft ' file_names];
end

disp('Compiling ... ')
cmd
eval(cmd)
