% Purpose: Combine the two plots PotentialPlot and OpenWiresPlot into one
% graph

n1 = 40; % Open wires graph
m1 = 50; 
t1 = 17;
alpha = 3; 
k0 = pi/2; 
n2 = 60; % potential graph
t2 = 12;
p = 50;

figure;
axis([0,m1+1,0,1]);
x1 = linspace(2,m1,m1-1);
y1 = zeros(m1-1);
for i = 2:m1
    y1(i-1) = Reflect1(n1,i,alpha,k0,t1);
end
plot(x1,y1,'r')
xlabel('Number of Wires/Potential');
ylabel('Reflectance');
hold on

x2 = linspace(0,p,p+1);
y2 = zeros(p+1);
for i = 1:p+1
    y2(i) = Reflect2(n2,t2,alpha,k0,i-1);
end
y2 = y2(:,1);
xx = 0:.001:p;
yy = spline(x,y,xx);
plot(x2,y2)


function reflected = Reflect1(n,m,alpha,k0,t)
% creating Hamiltonian matrix
sz = n*m-(m-1);
tau = -1;
H = zeros(sz);
for i = 1:sz
    for j = 1:sz
        if j == i-1 || i == j-1
            H(i,j) = tau;
        end
    end
end
for i = 1:m
    H(sz,(n-1)*i-(n-2)) = tau;
    H((n-1)*i-(n-2),sz) = tau;
end
H(sz,[1,sz-1]) = 0;
H([1,sz-1],sz) = 0;
H(sz,n-1) = tau;
H(n-1,sz) = tau;
for i = 1:m-1
    H(i*n-i,i*n-(i-1)) = 0;
    H(i*n-(i-1), i*n-i) = 0;
end

% finding eigenvalues and eigenvectors
[V,~] = eig(H);
E = eig(H);

% Gaussian Function for t=0
t0 = zeros(sz,1);
xs = n/2;
for k = 1:sz
    t0(k) = (1/(2*pi*alpha^2)^0.25)*exp(-0.25*(k-xs)^2/alpha^2)*exp(1i*k0*(k-xs));
end

% using the given formula to generate the wave function at time t
wave = zeros(sz,1);
for j = 1:sz
    wave = wave + V(:,j)*dot(V(:,j),t0)*exp(-E(j)*1i*t);
end

% calculating the probabilities at each point using the wave function
prob = zeros(sz,1);
for i = 1:sz
    prob(i) = wave(i)*conj(wave(i));
end
reflected = 0;
for i = 1:n-1
    reflected = reflected + prob(i);
end
end

function reflected = Reflect2(N,t,alpha,k0,p)
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