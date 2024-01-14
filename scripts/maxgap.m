function [width,lower,upper,rel] = maxgap(X)
% rows of X are samples
% columns of X are variables
n = size(X,1);
X = sort(X); % each col in ascending order
gaps = X(2:n,:) - X(1:n-1,:);
[width,idx] = max(gaps,[],1);
linIdx = idx + (0:length(idx)-1)*n;
lower = X(linIdx);
upper = X(linIdx+1);
rel = width ./ (X(n,:) - X(1,:));
end