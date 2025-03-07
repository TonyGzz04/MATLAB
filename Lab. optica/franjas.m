
clc; clear; close all;

dominio = 2*pi * 65;
pts = 1000; 

x = linspace(0, dominio, 1920);
y = linspace(0, dominio, 1080);
[Y, X] = meshgrid(x, y);



Z = sin(X+3*pi/2); 
Z_norm = (Z - min(Z(:))) / (max(Z(:)) - min(Z(:)));
Z_img = uint8(Z_norm .* 255);
imwrite(Z_img, 'franjas4.png');

figure;

imshow('franjas4.png');

%%

Z = sin(X+pi); 
Z_norm = (Z - min(Z(:))) / (max(Z(:)) - min(Z(:)));
Z_img = uint8(Z_norm .* 255);
imwrite(Z_img, 'franjas.png');

figure;
imshow('franjas.png');

%%

clc; clear; close all;

% Read the image
imgRGB = imread('1.png');

% Check the size of the image
disp(size(imgRGB)); % Should show [a, b, 3] for an RGB image

% Convert to grayscale
imgGray = rgb2gray(imgRGB);

% Display the grayscale image
figure;
imshow(imgGray);
title('Converted Grayscale Image');

% Display the size of the grayscale image
disp(size(imgGray)); % Should show [a, b]
