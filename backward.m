function x = backward(A,b)
[~,n] = size(A);
bn = length(b);
if n ~= bn
    error('Number of rows of A  and b must be equal');
end

x = zeros(n,1);
if A(n,n) == 0
    error('the problem does not have a solution');
    return;
end
x(n) = b(n)/A(n,n);

for i=n-1:-1:1
    sum = A(i,i+1:n)*x(i+1:n);
    x(i) = (b(i) - sum)/A(i,i);
end