%BFT_PARAM Set a paramater of the BeamForming Toolbox
%
%USAGE  : bft_param(...)
%
%INPUT  : ... Variable number of arguments, given in pairs 'name', value.
%         name - Name of the parameter (string). Currently supported:
%                -----+-----------------------+--------------+------
%                name |        Meaning        | Default value|  Unit 
%                -----+-----------------------+--------------+------
%                 'c' | Speed of sound.       | 1540         |  m/s 
%                 'fs'| Sampling frequency    | 40,000,000   |  Hz
%                -----+-----------------------+--------------+------
%         value - New value for the parameter. Must be scalar. 
%
%OUTPUT : None
%
%VERSION: 1.1, Jan 5, 2001 Svetoslav Nikolov

%VERSION: 1.0, Feb 10, 2000 Svetoslav Nikolov

function bft_param(varargin)

no_vars = length(varargin);

if rem(no_vars,2) == 1,
   error('The arguments must be given in pairs')
end

for ii = 1:2:no_vars
   name = varargin{ii};
   value = varargin{ii+1};
   bft(2,name,value);
end

