% Purpose: plot the amount of a wave packet that reflects as the potential
% increases

n = 60;
t = 12;
alpha = 3;
k0 = pi/2;
p = 50;

figure;
axis([0,p+2,0,1]);
x = linspace(0,p,p+1);
y = zeros(p+1);
for i = 1:p+1
    y(i) = Reflect(n,t,alpha,k0,i-1);
end
y = y(:,1);
xx = 0:.001:p;
yy = spline(x,y,xx);
plot(x,y,'o',xx,yy)
xlabel('Potential');
ylabel('Reflectance');
hold on

function reflected = Reflect(N,t,alpha,k0,p)
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
reflected = 0;
for i = 1:ceil(N/2)-1
    reflected = reflected + prob(i);
end
end