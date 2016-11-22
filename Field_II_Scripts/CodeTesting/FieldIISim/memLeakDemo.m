% Demonstrate that repeatedly calling ele_waveform increases MATLAB's 
% memory usage.
% Memory usage increases by 2 GB in this example and does not return to 
% original value even after "field_end" and "clear all" are called. 
% Memory usage for MATLAB can be tracked using Window's Task Manager.

% Start Field II
field_init(-1);            

% Create the transmit array
no_elements = 192;      
width = 0.2/1000;       
height = 5/1000;        
kerf = 0.02/1000;      
xmt = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);

% Call ele_waveform a few times
% Removing this section removes the memory issue  
for lineNo = 1:5000
    ele_waveform(xmt, (1:no_elements)', ones(no_elements, 256));
end

% Release memory used by Field II
field_end;

% Clears all objects in MATLAB workspace
clear all;
