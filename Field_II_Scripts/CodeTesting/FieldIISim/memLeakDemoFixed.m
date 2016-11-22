% This example demonstrates a workaround to the memory leak from memLeakDemo.
% The workaround is to call field_end after each single call to
% ele_waveform, and re-initialize field II before the next call to
% ele_waveform.

% Set parameters for transmit array
no_elements = 192;      
width = 0.2/1000;       
height = 5/1000;        
kerf = 0.02/1000;      

% Call ele_waveform a few times
% Note how Field II is iniatilized at the start of the loop and field_end is 
% called at the end.
for lineNo = 1:5000
    % Initialize field II each time
    % Need to call all field related things here.
    field_init(-1);  
    xmt = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);

    ele_waveform(xmt, (1:no_elements)', ones(no_elements, 256));
    
    % Release memory used by Field II
    field_end;
end

% Clears all objects in MATLAB workspace
clear all;
