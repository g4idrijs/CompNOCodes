% Generate a triangle waveform for purpose of getting a clear
% and simple input shape that doesn't change suddenly
function triangSeq =  triangle(height, numel)

triangSeq = [linspace(0,height,round(numel/1)), linspace(height,0,round(numel/1))];