function x_bar = SoftThresh_C(y,threshold)

% x_bar=SoftThresh_C(y,threshold)
% 
% SoftThreshold function
%   y = input array, expected in 1D vector format
%   threshold = threshold
% (c) Oliver Wieben  2017

x_bar = zeros(size(y));

% could be implemented more efficiently with find command
for n = 1:length(y)
        if abs(y(n)) <= threshold
            x_bar(n) = 0;
        else
            x_bar(n) = (abs(y(n))-threshold)*y(n)/abs(y(n));
        end

    
    
    
    if y(n) < (-1)*threshold
        x_bar(n) = y(n) + threshold;
    elseif y(n) > threshold
        x_bar(n) = y(n) - threshold;
    else
        x_bar(n) = 0;
    end
end
