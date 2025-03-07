
clc; clear; close all;

m = 0;
w0 = 1e-3;
pts = 100;
phi = 0;
a = 1;
clim = [0 1];

x = linspace(-a*w0,a*w0,pts);
y = linspace(-a*w0,a*w0,pts);
z = 0;
[X, Y] = meshgrid(x,y);
R = sqrt(X.^2 + Y.^2);

lambda = 660e-9;
k = 2*pi/lambda;
kx = k*cos(pi/2);
ky = k*cos(pi/2);
kz = k*0;

E = (R ./ w0).^abs(m) .* exp(-R.^2/w0^2) .* exp(1i*m*phi);

E1 = E;
E2 = E .* exp(1i.*(kx.*X + ky.*Y + kz.*z));

I = abs(E1 + E2).^2;

figure;
imagesc(I)


%%

clc; clear; close all;


w0 = 1e-3;
lambda = 663e-9;
k = 2*pi/lambda;
N = 256;
x = linspace(-3*w0,3*w0,N);
y = x;
[X, Y] = meshgrid(x,y);

GB = @(X,Y,w) exp(-(X.^2 + Y.^2)/w^2);
phi = 1*pi/180;
theta = 0;
z = 0;

E1 = GB(X,Y,w0);
E2 = GB(X,Y,w0) .* exp(1i*(k.*X.*cos(theta).*sin(phi) + k.*Y.*sin(phi).*sin(theta) + k.*z.*cos(phi)));
E = E1 + E2;

% figure;
% imagesc(abs(E).^2);

mirror = exp(1i*pi);
E3 = 0.5 * E1 .* mirror;
E4 = 0.5 * E1 .* mirror;

a = 0.4;
bs_reflect = exp(1i * a);

E5 = 0.5 * E2 .* mirror .* bs_reflect .* mirror;
E6 = 0.5 * E2 .* mirror .* bs_reflect;


foco = -0.1;
lente = exp(1i*(k/(2*foco)) * (X.^2 + Y.^2));

I1 = abs(E3+E5.*lente).^2;
I2 = abs(E4+E6).^2;

figure; 
imagesc(I1);
figure; 
imagesc(I2);





