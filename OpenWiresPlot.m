% Purpose: plot the reflectance of a wave packet as the number of wires
% increases

n = 40; % number of positions in a wire 
m = 80; % maximum number of wires
t = 17;
alpha = 3; % width of Gaussian function
k0 = pi/2; % changes the spread and initial velocity of the packet

figure;
axis([0,m+1,0,1]);
x = linspace(2,m,m-1);
y = zeros(m-1);
for i = 1:m-1
    y(i) = Reflect(n,i+1,alpha,k0,t);
end
plot(x,y)
xlabel('Number of Wires');
ylabel('Reflectance');
hold on


function reflected = Reflect(n,m,alpha,k0,t)
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