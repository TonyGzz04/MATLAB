
clc; close all; clear;

w0 = 1.3e-3;
lambda = 633e-9;

zmax = 100;
% Ventana numérica
N = 2^10;


x = linspace(-20*w0, 20*w0, N);
y = x;
z = linspace(0, zmax, N);
[X, Y] = meshgrid(x, y);
dz = 0.0001;


kt = sqrt(2*pi./X.^2 + 2*pi./Y.^2);
k = 2*pi/lambda;
kr = 500;

U0 = 1;
% U = U0 * exp(-(X.^2+Y.^2) ./ w0^2);
U = U0 * ( sqrt((X + 0.0003).^2 + Y.^2)<0.1 + sqrt((X - 0.0003).^2 + Y.^2)<0.1 );

U1 = U0 * sqrt((X + 0.003).^2 + Y.^2)<0.001;
U2 = U0 * sqrt((X - 0.003).^2 + Y.^2)<0.001;
U = U1+U2;




figure; colormap('gray');
for z = 0:dz:0.6
    res = ifft2(fft2(U) .* exp(-1i*0.5*kt.^2./k .* z));
    imagesc(abs(res).^2);
    pause(0.1);
end
