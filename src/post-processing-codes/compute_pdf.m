function [x,f] = compute_pdf(u, nbins)
% This function computes the probability density function
% of a signal u by normalizing the area of histogram created using the
% number of bins specified by nbins
% Taken from Bendat and Piersol Eq. (11.85), Pg 381

% compute histogram
[counts, centers] = histogram(u, nbins);

% compute the width of each bin (data units)
bin_width = nbins / (max(u) - min(u));

% find the area of the histogram
hist_area = trapz(centers, counts);

% normalize counts by the area so the integral is equal to 1
counts_normalized = counts*bin_width / hist_area;

% convert the normalized counts to a probability estimate by dividing by
% the length of each bin interval
probability_density = counts_normalized / bin_width;

% assign counts and centers to the variables that will be returned
x = centers;
f = probability_density;
end