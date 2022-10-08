% Purpose: investigate the specific wave functions of the resulting waves
% when a quantum wave hits a potential - comparing negative potentials with
% positive potentials

n = 60; % number of points on the circle
alpha = 3;
k0 = pi/2;

% calls two functions with equal but opposite potentials
y = CircleDynamicsPotential(n,11.7,alpha,k0,-4);
u1 = y(15)
u2 = y(45)
z = CircleDynamicsPotential(n,11.7,alpha,k0,4);
u3 = z(15)
u4 = z(45)

% this function determines probability
function wave = CircleDynamicsPotential(N,t,alpha,k0,p)
% creating Hamiltonian matrix
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

% Setting the potential
H(ceil(N/2),ceil(N/2)) = p;

% finding eigenvalues and eigenvectors
[V,~] = eig(H);
E = eig(H);

% Gaussian Function for t=0
t0 = zeros(N,1);
xs = ceil(20);
for k = 1:N
    t0(k) = (1/(2*pi*alpha^2)^0.25)*exp(-0.25*(k-xs)^2/alpha^2)*exp(1i*k0*(k-xs));
end

% using the formula to generate the wave function at time t
wave = zeros(N,1);
for j = 1:N
    wave = wave + V(:,j)*dot(V(:,j),t0)*exp(-E(j)*1i*t);
end

% calculating the probabilities at each point using the wave function
prob = zeros(N,1);
for i = 1:N
    prob(i) = wave(i)*conj(wave(i));
end

% optional phase matrix that can be returned
phase = zeros(N,1);
for i = 1:N
    phase(i) = angle(wave(i));
end
end