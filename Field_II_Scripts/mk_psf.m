% mk_psf.m

min_sample = min(times)*fs;
max_sample = max([max(times)*fs 2*0.121/c*fs]);

[n, m] = size(image_data);
n = n + (max_sample-min_sample);

psf = 0;
for i = 1:no_lines
    ps = [zeros(round(times(i)*fs-min_sample), 1); image_data(:, i)];
    psf(1:max(size(ps)),i) = ps;
end

clear ps;

lateral = d_x*[-(no_lines-1)/2:1:(no_lines-1)/2]*1000;

% generate axes for image
% first row in image corresponds to min(times) in seconds
timetxrx = min(times):1/fs:min(times)+(size(psf, 1)-1)*1/fs;
timetxrx = timetxrx-2*0.5*(length(conv(excitation, conv(impulse_response, impulse_response)))-1)/fs;
disttxrx = (timetxrx)*c/2;
disttxrx = (disttxrx*1000); % convert to mm
    
% imagesc(lateral, disttxrx, psf);