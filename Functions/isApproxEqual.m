function [out] = isApproxEqual(v1,v2, TOL)
% Checks if two arrays v1 and v2 are apprxoximatly equal by finding 
% the max. component-based difference and compare it to a tolerance TOL.

if nargin < 3 
    TOL = 1e-9;
end

if max(abs(v1-v2),[],'all') < TOL 
    out = true;
else
    out = false;
end

end