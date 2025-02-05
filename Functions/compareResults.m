function [out] = compareResults(d1,d2)
% Compare to Results vectors d1 and d2
    
% Initalize output
out = zeros(1,4);
% Remove NaNs
idx1 = or(isnan(d1),isnan(d2));
d1(idx1==1) = [];
d2(idx1==1) = [];
% Maximal difference
out(1)           = max(abs(d1-d2),[],'all');
% Normalized Euclidian distance
out(2)           = norm(d1-d2)/norm(d1);
% Linear correlation coefficient
[out(3), out(4)] = corr(d1,d2);
disp(['The max. absolute error comapred to the reference_data is ', num2str(out(1))])
disp(['The normalized Euclidian distance to the reference_data is ', num2str(out(2))])
disp(['The linear correlation with the reference data is ', num2str(out(3)), ' (p=', num2str(out(4)),')'])
end