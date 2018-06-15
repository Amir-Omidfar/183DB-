function M = qMatrix( A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    M = vertcat([1 0 0 0], horzcat([0 0 0].', A));
end

