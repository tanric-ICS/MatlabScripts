% Purpose: Given M sheets with dimensions LxL, model the movement of a
% 2D Gaussian quantum packet as it approaches a node that exists on all
% sheets

M = 2; % Number of Sheets
L = 50; % Dimensions
sz = M*(L^2)-(M-1); % Size of matrix
pos = 1225; % Position of node on the matrix
t = 20000; % Ending time
alpha = 1.5;
k0 = pi/2;
H = zeros(sz,sz);
tau = -1;
% Creating the Hamiltonion
for m = 1:M
    if m == 1
        t = reshape(1:L^2,L,L);
    else
        t = L^2*(m-1)-(m-3):L^2*m-(m-1);
        t = [t(1:pos-1),pos,t(pos:end)];
        t = reshape(t,L,L);
    end
    for l = 1:L
        i = t(:,l);
        i2 = t(l,:);
        if l ~= L
            j = t(:,l+1);
            j2 = t(l+1,:);
        else
            j = t(:,1);
            j2 = t(1,:);
        end
        for k = 1:L
            H(i(k),j(k)) = tau;
            H(j(k),i(k)) = tau;
            H(i2(k),j2(k)) = tau;
            H(j2(k),i2(k)) = tau;
        end
    end
end
% Finding eigenvectors and eigenvalues
[V,~] = eig(H);
E = eig(H);

% Creating initial Gaussian
t0 = zeros(sz,1);
xs = 10;
ys = 10;
for k = 1:sz
    xcoord = mod(k,L);
    if xcoord == 0
        xcoord = L;
    end
    ycoord = ceil(k/L);
    t0(k) = exp(-0.25*(xcoord-xs)^2/alpha^2)*exp(-0.25*(ycoord-ys)^2/alpha^2)*exp(1i*k0*(ycoord-ys))*exp(1i*k0*(xcoord-xs)) / (14.1372^0.5);
end

% Creating M subplots for each sheet
figure;
for i = 1:M
    subplot(1,M,i);
    x = linspace(1,L,L);
    y = linspace(1,L,L);
    [X,Y]= meshgrid(x);
    temp = zeros(L,L);
    s(i) = surf(X,Y,temp);
    axis([0,L,0,L,0,0.2]);
end

% Iterating through each sheet and changing data for probability as time
% increases
for k = 0:t
    temp = Prob(sz,V,E,t0,k);
    z = zeros(L^2,M);
    z(:,1) = temp(1:L^2);
    % Separating the vector into each sheet
    for i = 2:M
        z(1:pos-1,i) = temp(((L^2)-1)*(i-1)+2:((L^2)-1)*(i-1)+pos);
        z(pos,i) = temp(pos);
        z(pos+1:L^2,i) = temp(((L^2)-1)*(i-1)+pos+1:((L^2)-1)*i+1);
    end
    % Tracking probability on Sheet B (2nd sheet)
    test = 0;
    for i = 1:2500
        test = test + z(i,2);
    end
    scatteringB = test
    % Tracking the scattering of the particle on Sheet A (1st sheet)
    test3 = 0;
    for i = 1:25
        for j = 1:25
            test3 = test3 + z(50*(i-1)+j,1);
        end
    end
    scatteringA = test3*4
    % Changing probability data
    for j = 1:M
        z1 = z(:,j);
        z1 = reshape(z1,L,L);
        s(j).ZData = z1;
        axis([0,L,0,L,0,0.1]);
    end
    % Tracking time
    k
    pause(1)
end

% Function used to calculate probability as time progresses (movement of
% particle)
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