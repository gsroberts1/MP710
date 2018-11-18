function x_bar = SoftThresh(y,threshold)

% x_bar=SoftThresh(y,threshold)
% 
% SoftThreshold function
%   y = input array, expected in 1D vector format
%   threshold = threshold
% (c) Oliver Wieben  2017

x_bar = zeros(size(y));

% could be implemented more efficiently with find command
for n = 1:length(y)
    if y(n) < (-1)*threshold
        x_bar(n) = y(n) + threshold;
    elseif y(n) > threshold
        x_bar(n) = y(n) - threshold;
    else
        x_bar(n) = 0;
    end
end
