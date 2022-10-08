% Purpose: model the behaviour of a quantum wave in a circle with a
% potential applied at the centre

n = 60; % number of points on the circle
t = 500; % time
alpha = 3;
k0 = pi/2;
p = 50;

% creates an animated plot of the probability at each point over time
h = animatedline;
axis([1,n,0,0.3]);
x = linspace(1,n,n);
title('Dynamics of a Quantum Particle in a Circle with a Potential Over Time');
xlabel('Location');
ylabel('Probability');
for k = 0:0.1:t
    reflect = 0;
    [y,~] = CircleDynamicsGaussian(n,k,alpha,k0,p);
    addpoints(h,x,y);
    drawnow;
    clearpoints(h);
    for i = 1:ceil(n/2)-1
        reflect = reflect + y(i);
    end
    reflect
    k
end

% this function determines probability
function [prob, phase] = CircleDynamicsGaussian(N,t,alpha,k0,p)
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
H(ceil(N/2),ceil(N/2)) = p;

% finding eigenvalues and eigenvectors
[V,~] = eig(H);
E = eig(H);

% Gaussian Function for t=0
t0 = zeros(N,1);
xs = ceil(N/3);
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
phase = zeros(N,1);
for i = 1:N
    phase(i) = angle(wave(i));
end
end