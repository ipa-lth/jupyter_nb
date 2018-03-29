function Q_sum = calcQsum(Q_in)
% calculate the 'blown up' matrix

if ~iscell(Q_in)
    error('cell array input expected')
end

m = length(Q_in);
n = size(Q_in{1},1);
if m == 1
    Q_sum = [2*Q_in{1}, zeros(n); zeros(n), zeros(n)];
elseif m == 2
    Q_sum = [2*Q_in{1}, Q_in{2}; Q_in{2}, zeros(n)];
else
    Q_sum = [2*Q_in{1}, Q_in{2}; Q_in{2}, 2*Q_in{3}];
end
for i1 = 4:2:m
    if i1 ~= m
        Q_sum = [Q_sum, [zeros((i1/2-1)*n,n); Q_in{i1}]; ...
            zeros(n,(i1/2-1)*n), Q_in{i1}, 2*Q_in{i1+1}]; %#ok<AGROW>
    else % even number of elementary matrices --> zeros in southeasternmost position
        Q_sum = [Q_sum, [zeros((i1/2-1)*n,n); Q_in{i1}]; ...
            zeros(n,(i1/2-1)*n), Q_in{i1}, zeros(n)]; %#ok<AGROW>
    end
end
Q_sum = 1/2*Q_sum;

end % calcQsum
