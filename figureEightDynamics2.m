% Purpose: model a quantum wave in a figure eight geometry using two
% figures to represent each loop

n = 60; % number of points on one of the circles
t = 100; % time
alpha = 3;
k0 = pi/2;

% creating two figures with labels and axes
figure;
subplot(1,2,1);
axis([1,n,0,0.3]);
axis1 = gca;
xlabel('Location');
ylabel('Probability');
subplot(1,2,2);
axis([n,2*n-1,0,0.3]);
axis2 = gca;
xlabel('Location');
ylabel('Probability');

% creates two animated plots of the probability at each point over time
h = animatedline(axis1);
x1 = linspace(1,n,n);
l = animatedline(axis2);
x2 = linspace(n,2*n-1,n);

% displays animation for each loop
for k = 0:0.1:t
    y1 = FigureEightDynamics(n,k,alpha,k0);
    y1 = y1(1:n);
    addpoints(h,x1,y1);
    drawnow;
    clearpoints(l);
    y2 = FigureEightDynamics(n,k,alpha,k0);
    y2 = y2(n:2*n-1);
    addpoints(l,x2,y2);
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
        if abs(i-j) == 1
            H(i,j) = tau;
        end
    end
end
H(1,N) = tau;
H(N,1) = tau;
H(2*N-1,N) = tau;
H(N,2*N-1) = tau;

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