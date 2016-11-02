% custTimes == number of times to repeat code elements
% len == number of elements transducer impulse response

% highWidth = width above threshold level / length impulse response
% dBThresh = dB threshold for calculating width above threshold 
% max strength = maximum value in decoded signal (not normalized)

% 1D simulation:
% -Simulate 1D slice of firing beam at a phantom
% (simulates data received in slice containing beam)
% -Additively incorporate clutter (from other beams)
% -Set # of transmits and how to process received data

function [highWidth, maxStrength] = oneDimCustMain(varargin) % custTimes, len, plotResult, dBThresh, maxValNoise, codes

% Default noise and code values
defNoise = 50;
addpath('../../../Complementary Pairs')
defCodes = load('len10_8436113.mat');
defCodes = defCodes.('x');

% Default argument values
if(nargin == 0)
   custTimes = 1;
   len = 600;
   plotResult = 1;
   dBThresh = -45; % dB
   maxValNoise = defNoise;
   codes = defCodes;
elseif(nargin == 1)
    custTimes = varargin{1};
    len = 600;
    plotResult = 1;
    dBThresh = -45; % dB
    maxValNoise = defNoise;
    codes = defCodes;
elseif(nargin == 2)
    custTimes = varargin{1};
    len = varargin{2};
    plotResult = 1;
    dBThresh = -45; % dB
    maxValNoise = defNoise;
    codes = defCodes;
elseif(nargin == 3);
    custTimes = varargin{1};
    len = varargin{2};
    plotResult = varargin{3};
    dBThresh = -45; % dB
    maxValNoise = defNoise;
    codes = defCodes;
elseif(nargin == 4)
    custTimes = varargin{1};
    len = varargin{2};
    plotResult = varargin{3};
    dBThresh = varargin{4};
    maxValNoise = defNoise;
    codes = defCodes;
elseif(nargin == 5)
    custTimes = varargin{1};
    len = varargin{2};
    plotResult = varargin{3};
    dBThresh = varargin{4};
    maxValNoise = varargin{5};
    codes = defCodes;
elseif(nargin == 6)
    custTimes = varargin{1};
    len = varargin{2};
    plotResult = varargin{3};
    dBThresh = varargin{4};
    maxValNoise = varargin{5};
    codes = varargin{6};
end

% Start building title
titleStr = sprintf('Ultrasound System Output.  Elem rep: %.2f. \n Numel ImpResp: %.0f.',double(custTimes)/double(len), len);

%% Define system propereties
% Specify sequence values, starting at n = 0.
% Unspecified values are assumed to be zero.

% Transducer transfer function - choices
% ==============================
% Hanning'ed sine wave
hannSin = @(number_cycles, len)sin(linspace(0,number_cycles*2*pi,len)).*hanning(len)';
number_cycles = 2;

% Delta
delt = [1];
%  ==============================

% Transducer transfer function - selection
%  ==============================
hTrans = hannSin(number_cycles, len);
% hTrans = delt;
% ==============================
 
% Phantom transfer function
hPhant = [1];

% Entire system transfer function
hT = conv(conv(hTrans, hPhant),hTrans);

%% Determine raw decoded data

% Set method of transmitting and decoding data

% Uniform noise with single beam
% % % % % % % % % % % %
% Uniform noise
% a =-maxValNoise; b = maxValNoise; 
% c = @(length, transEvent) a + (b-a).*rand(length,1);

% Complementary pair - single beam
% decode1D = @(hT,c,len, custTimes,codes)twoTransComp(hT, c, len, custTimes,codes); 
% % % % % % % % % % % %

% Two beams firing - receive only on one
phant = [];
custTimes = round(len/4); % Close to ideal transmit length
decode1D = twoBeams(codes, custTimes, hTrans, phant, len)



R = decode1D(hT,c, len, custTimes,codes);

%% Determine final decoded data
envData = 1;

% Envelope
env = @(R)abs(hilbert(R));

% Identity function
id = @(R)R; 

if(envData == 1)
   F = env(R);
   titleStr = strcat(titleStr, ' Enveloped.');
else
   F = id(R); 
end

%% Determine -45 dB width (in units of impulse response length)
highWidth  = oneDWidth(F,len, dBThresh);
titleStr = strcat(titleStr, sprintf(' High width: %.2f', highWidth));

%% Get maximum value of decoded signal
maxStrength = max(abs(F));

%% Display decoded data
if(plotResult == 1)
    plotOneDim(F, titleStr, len)
end