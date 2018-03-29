function Q_out = transQ(Q_in,a,b)
% calculate the interval transformation (p(q) = a*q+b)

if ~iscell(Q_in)
    error('cell array input expected')
end

m = length(Q_in);
n = size(Q_in{1},1);

for i1 = 0:m-1
    Q_out{i1+1} = Q_in{i1+1}; %#ok<AGROW>
    for i2 = i1+1:m-1
        Q_out{i1+1} = Q_out{i1+1} + b^(i2-i1)*nchoosek(i2,i1)*Q_in{i2+1}; %#ok<AGROW>
    end
    Q_out{i1+1} = a^i1*Q_out{i1+1}; %#ok<AGROW>
end

Q_out;
end % transQ
