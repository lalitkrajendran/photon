function [bias_error, random_error, total_error] = compute_errors(x, x_ref)
% This function computes the bias, random and rms errors associated with
% the input data set. The errors are defined as follows:
% Bias/Systematic Error: E[x] - x_ref
% Random error: sqrt(E[{x - E[x]}^2])
% RMS error: sqrt(E[{x-x_ref}^2]) = sqrt(Bias^2 + Random^2) 

% INPUTS: x - actual data, 
%         x_ref - ground truth

    % convert 2d arrays to 1d arrays
%     x = x(:);
%     x_ref = x_ref(:);

    % bias error
    bias_error = nanmean(x - x_ref);    

    % random error
    random_error = std(x);
    
    % total error
    total_error = sqrt(bias_error.^2 + random_error.^2); %rms(x - x_ref);

end