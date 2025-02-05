function [epsilon] = maxError(v1,v2)
% Compute the maximum Error between two arrays v1 and v2
epsilon = max(abs(v1-v2),[],'all');
end