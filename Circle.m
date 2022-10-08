N = 5;
H = zeros(N);
tau = -1;
for i = 1:N
    for j = 1:N
        if j == i-1 || i == j-1 || abs(j-i) == N-1
            H(i,j) = tau;
        end
    end
end
H(1,N) = tau;
H(N,1) = tau;

[V,E] = eig(H);

x = zeros(N);
for i = 1:N
    for j = 1:N
        x(i,j) = V(i,j)^2;
    end
end