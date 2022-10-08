[H,wave,prob,check] = CircleDynamics(50,10);
function [H,wave,prob,check] = CircleDynamics(N,t)
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

[V,~] = eig(H);
E = eig(H);

tzero = zeros(N,1);
tzero(ceil(N/2),1) = 1;

wave = zeros(N,1);
for j = 1:N
    wave = wave + V(:,j) * dot(V(:,j),tzero)*exp(-E(j,1)*1i*t);
end

prob = zeros(N,1);
for i = 1:N
    prob(i,1) = wave(i,1)*conj(wave(i,1));
end

check = 0;
for i = 1:N
    check = check + prob(i,1);
end

domain = 1:N;
plot(domain, prob);
end
