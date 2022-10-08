test = Prob(2.048)

% M=2, G=0.838
% M=3, G=1.248
% M=4, G=1.523
% M=5, G=1.733
% M=6, G=1.903
% M=7, G=2.048

function scatter = Prob(g)
L = 50;
t = 14;
sz = L^2;
alpha = 1.5;
k0 = pi/2;
H = zeros(sz,sz);
tau = -1;
for i = 1:sz
    for j = 1:sz
        if j == i-1 || i == j-1 || j == i-L || i == j-L 
            H(i,j) = tau;
        end
    end
end
i = L:L:sz-L;
j = L+1:L:sz-(L-1);
for k = 1:numel(i)
    H(i(k),j(k)) = 0;
    H(j(k),i(k)) = 0;
end
i = 1:L:sz;
j = L:L:sz;
for k = 1:numel(i)
    H(i(k),j(k)) = tau;
    H(j(k),i(k)) = tau;
end
i = 1:L;
j = sz-(L-1):sz;
for k = 1:numel(i)
    H(i(k),j(k)) = tau;
    H(j(k),i(k)) = tau;
end
H(1225,1225) = g;
[V,~] = eig(H);
E = eig(H);

%t0 Gaussians
t0 = zeros(sz,1);
xs = 10;
ys = 10;
for k = 1:sz
    xcoord = mod(k,L);
    if xcoord == 0
        xcoord = 50;
    end
    ycoord = ceil(k/L);
    t0(k) = exp(-0.25*(xcoord-xs)^2/alpha^2)*exp(-0.25*(ycoord-ys)^2/alpha^2)*exp(1i*k0*(ycoord-ys))*exp(1i*k0*(xcoord-xs)) / (14.1372^0.5);
end

wave = zeros(sz,1);
for j = 1:sz
    wave = wave + V(:,j)*dot(V(:,j),t0)*exp(-E(j)*1i*t);
end

prob = zeros(sz,1);
for i = 1:sz
    prob(i) = wave(i)*conj(wave(i));
end

test3 = 0;
for i = 1:25
    for j = 1:25
        test3 = test3 + prob(50*(i-1)+j,1);
    end
end
scatter = test3*4;
end