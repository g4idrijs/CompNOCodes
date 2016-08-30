function joins = joins( intVect )
% Find the join lengths of a vector
% This is also the sum of continguous interval lengths
% [2,3] -> 2,3,5 join lengths
% [2,3,5] -> (2,3,5),[(5,8),(10)]

% Plan: go through each element
% At each element, start a cumulative sum that goes to the end of the
% vector
joins = zeros(1,1000); % Preallocating
i = 1;
for startEl = 1:length(intVect)
    for endEl = startEl:length(intVect)
        joins(i) = sum(intVect(startEl:endEl));
        i = i+1;
    end
end

% Remove pre-allocating zeros
joins = joins(joins~=0);

% Gives visual picture of join usage efficiency
% plot(sort(joins([2, 3, 4, 6, 8, 11, 16, 12, 24, 20, 17])))
% Efficiency of packing decreases with code length, using that output in
% the paper

end

