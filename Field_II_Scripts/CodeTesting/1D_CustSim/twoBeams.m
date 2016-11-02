% Return raw decoded data 
% codes == complementary codes array (each two rows forms a comp. pair)
    % (should have 4 rows, for two beams)
    % (first two rows = main beam, second two rows = secondary beam)
% custTimes  = number of times to repeat codes
% hTrans == transducer transfer function
% phant == phantom
    % first row = depth of phantoms in first beam
    % second row = depth of phantoms in second beam
% speedWave = number of units wave travels axially per imp. resp. time
% trLen == length of impulse response of transducer
% beamSp == spacing between beams
function R = twoBeams(codes, custTimes, hTrans, phant, speedWave, trLen, beamSp)
    %% Repeat code elements as desired (to transmit for correct time: maximize clarity)
    codes = repelem(codes,1,custTimes);
    
    %% Loop through transmit events
    % First transmit
    % Set main beam and interference beam codes
    mainCode1 = codes(1,:);
    intCode1 = codes(3,:);
    
    % Use phantom to create and add shifted copies of main code
    % (scaled down/delay by distance)
   
    % Use phantom to create and add shifted copies of interference code
    % (scaled down/delay by distance)
    
    % Second transmit
    
    % Calculate raw decoded data
    R = xcorr(htC(code1Pair1),code1Pair1) + xcorr(htC(code2Pair1),code2Pair1);
    
    % Trim zeros (very small values)
    minVal = 1e-10;
    [~,J] = find(R > minVal);
    R = R(min(J):end);
    [~,J] = find(R > minVal);
    R = R(1:max(J));
end

