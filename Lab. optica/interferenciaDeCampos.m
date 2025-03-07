
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

figure;
imagesc(abs(E).^2);


