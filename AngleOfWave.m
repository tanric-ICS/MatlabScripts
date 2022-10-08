n = 40; % number of positions in a loop 
m = 4; % number of loops
t = 15.5; % length of animation
alpha = 3; % width of Gaussian function
k0 = pi/2; % changes the spread and initial velocity of the packet

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
    H(sz,[(n-1)*i-(n-2),(n-1)*i]) = tau;
    H([(n-1)*i-(n-2),(n-1)*i],sz) = tau;
end
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

x = zeros(n-1,m);
figure;
for i = 1:m
    subplot(1,m,i);
    
    x(:,i) = linspace((n-1)*i-(n-2),(n-1)*i,n-1);
    x = x(:,i);
    [y,a] = nLoopsProb(sz,V,E,t0,t);
    y = y((n-1)*i-(n-2):(n-1)*i);
    a = a((n-1)*i-(n-2):(n-1)*i);
    for ii = 1:length(x)
        if a(ii) < 0
            c = 'b*';
        else
            c ='r*';
        end
        plot(x(ii),y(ii),c)
        hold on
    end
    xlabel('Location');
    ylabel('Probability');
    axis([(n-1)*i-(n-2),(n-1)*i,0,0.3]);
end



function [prob,phase] = nLoopsProb(sz,V,E,t0,t)
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
phase = zeros(sz,1);
for i = 1:sz
    phase(i) = angle(wave(i));
end
end
