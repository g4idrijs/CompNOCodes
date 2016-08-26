% focusBeam.m

% Function to focus a beam with signal @excitation at focal point @focus
% This is done by calculating the time delay from each element in
% @xElements to the focal point @focus.
% This function moves the delay profile so that the element that fires last
% fires at @maxDelay. This is needed for aligning the beams when there are
% multiple focal zones in the depth direction.

function beam = focusBeam(excitation, focus, xElements, fs, c, maxDelay)
    % Calculate delays for each element
    delays = (focus(3)-sqrt((xElements-focus(1)).^2+focus(3)^2))./c * fs;
    % Adjust delays to have a maximum at @maxDelay
    delays = delays-min(delays);
    delays = delays+(maxDelay-max(delays));
    
    if min(delays) < 0
        error(...
        ['@maxDelay is not big enough to hold beam. Need to grow by ', ...
            num2str(ceil(-min(delays)))]);
    end
    
    % Calculate length of beam required to store results with delays
    maxBeamLength = ceil(max(abs(delays))+length(excitation));
    beam = zeros(maxBeamLength, length(delays));
    for i = 1:length(delays)
        % Interpolate to shift beam
        tmpBeam = interp1(1:length(excitation), excitation, (1-delays(i)):length(excitation), 'linear', 0);
        beam(1:length(tmpBeam), i) = tmpBeam;
    end
end