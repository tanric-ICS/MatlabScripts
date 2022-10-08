% Purpose: Given one sheet with dimensions LxL, model the movement of a 2D
% quantum packet as it encounters a potential difference
L = 50;
alpha = 1.5;
k0 = pi/2;
sz = L^2;
tau = -1;
H = zeros(sz,sz);
t = 500;
g = 0.5;
% Creating the Hamiltonian matrix
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

% Initial quantum packet
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

% Creation of figure
figure;
axis([0,L,0,L,0,0.2]);
x = linspace(1,L,L);
y = linspace(1,L,L);
[X,Y]= meshgrid(x);
temp = zeros(L,L);
s = surf(X,Y,temp);
colorbar

% Updating the figure as t increases
for k = 0:t
    test = 0;
    axis([0,L,0,L,0,0.1]);
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    zlabel('Probability');
    z = Prob(sz,V,E,t0,k);
    for i = 1:25
        for j = 1:25
            test = test + z(50*(i-1)+j);
        end
    end
    % Calculating the scattering of the particle after it encounters the
    % potential
    scattering = test*4
    z = reshape(z,L,L);
    s.ZData = z;
    % keeping track of the time (k)
    k
    pause(1)
end

% This function calculates the evolution of the Gaussian over time
function prob = Prob(sz,V,E,t0,t)
wave = zeros(sz,1);
for j = 1:sz
    wave = wave + V(:,j)*dot(V(:,j),t0)*exp(-E(j)*1i*t);
end

prob = zeros(sz,1);
for i = 1:sz
    prob(i) = wave(i)*conj(wave(i));
end
end