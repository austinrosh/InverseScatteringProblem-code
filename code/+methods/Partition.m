function v = Partition(n, lgt)
% v = Partition(n)
% INPUT
%   n: non negative integer
%   lgt: optional non negative integer
% OUTPUT:
%   v: (m x lgt) non-negative integer array such as sum(v,2)==n
%       each row of v is descending sorted
%       v contains all possible combinations
%       m = P(n) in case lgt == n, where P is the partition function
%       v is (dictionnary) sorted
% Algorithm:
%    Recursive
% Example:
% >> Partition(5)
% 
% ans =
% 
%      5     0     0     0     0
%      4     1     0     0     0
%      3     2     0     0     0
%      3     1     1     0     0
%      2     2     1     0     0
%      2     1     1     1     0
%      1     1     1     1     1
if nargin < 2
    lgt = n;
end
v = PartR(lgt+1, n, Inf);
end % Partition
%% Recursive engine of integer partition
function v = PartR(n, L1, head)
rcall = isfinite(head);
if rcall
    L1 = L1-head;
end
if n <= 2
    if ~rcall
        v = L1;
    elseif head >= L1
        v = [head, L1];
    %else
       % v = zeros(0, n, class(L1));
    end
else % recursive call
    j = min(L1,head):-1:0;
    v = arrayfun(@(j) PartR(n-1, L1, j), j(:), 'UniformOutput', false);
    v = cat(1,v{:});
    if rcall
        v = [head+zeros([size(v,1),1], class(head)), v];
    end
end
end % PartR