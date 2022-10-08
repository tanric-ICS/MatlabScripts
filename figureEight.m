[H, fnought, alpha] = FigureEight(100);
function [H, fnought, alpha] = FigureEight(N)
H = zeros(2*N-1);
tau = 1;

for i = 1:2*N-2
    for j = 1:2*N-1
        if j == i-1 || i == j-1 || j-i == 2*N-2 || (i == N-1 && j == 2*N-1) || (i == N && j == 2*N-1)
            H(i,j) = tau;
        end
        if (i == N-1 && j == N) || (i == N && j == N-1)
            H(i,j) = 0;
        end
    end
end
H(2*N-1, [2*N-2, N-1, N, 1]) = tau;

[eigenvectors,~] = eig(H);

prob = zeros(2*N-1, 1);
for i = 1:2*N-1
    prob(i,1) = eigenvectors(i,2*N-1)^2;
end

fnought = prob(2*N-1, 1);
alpha = -log(prob(1,1)/fnought);

f = @(m) fnought*exp(-alpha*m);
domain = 0:N;
plot(domain, f(domain));
end

