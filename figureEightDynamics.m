n = 60; % number of points on one of the circles
t = 500; % time
alpha = 3;
k0 = pi/4;

% creates an animated plot of the probability at each point over time
h = animatedline;
axis([1,2*n-1,0,0.3]);
x = linspace(1,2*n-1,2*n-1);
title('Dynamics of a Quantum Particle in a Figure Eight Over Time');
xlabel('Location');
ylabel('Probability');
for k = 0:0.1:t
    y = FigureEightDynamics(n,k,alpha,k0);
    addpoints(h,x,y);
    drawnow;
    clearpoints(h);
end

% this function determines probability
function prob = FigureEightDynamics(N,t,alpha,k0)
% creating Hamiltonian matrix
H = zeros(2*N-1);
tau = -1;
for i = 1:2*N-1
    for j = 1:2*N-1
        if abs(i-j) == 2
            H(i,j) = tau;
        end
    end
end
H(1,2) = tau;
H(2,1) = tau;
H(N-1,N+1) = 0;
H(N-1,N) = tau;
H(N,[N-1,N+1]) = tau;
H(N+1,N-1) = 0;
H(N+1,N) = tau;
H(2*N-2,2*N-1) = tau;
H(2*N-1,2*N-2) = tau;

% finding eigenvalues and eigenvectors
[V,~] = eig(H);
E = eig(H);

% Gaussian Function for t=0
t0 = zeros(2*N-1,1);
xs = N/2;
for k = 1:2*N-1
    t0(k) = (1/(2*pi*alpha^2)^0.25)*exp(-0.25*(k-xs)^2/alpha^2)*exp(1i*k0*(k-xs));
end

% using the formula to generate the wave function at time t
wave = zeros(2*N-1,1);
for j = 1:2*N-1
    wave = wave + V(:,j)*dot(V(:,j),t0)*exp(-E(j)*1i*t);
end

% calculating the probabilities at each point using the wave function
prob = zeros(2*N-1,1);
for i = 1:2*N-1
    prob(i) = wave(i)*conj(wave(i));
end
end