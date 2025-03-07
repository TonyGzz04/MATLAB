
clc; clear; close all;

s01 = rgb2gray(im2double(imread("DSC_3937.JPG")));
s02 = rgb2gray(im2double(imread("DSC_3938.JPG")));
s21 = rgb2gray(im2double(imread("DSC_3939.JPG")));
s22 = rgb2gray(im2double(imread("DSC_3940.JPG")));
s31 = rgb2gray(im2double(imread("DSC_3941.JPG")));
s32 = rgb2gray(im2double(imread("DSC_3942.JPG")));

S0 = s01 + s02;
S1 = s01 - s02;
S2 = s21 - s22;
S3 = s31 - s32;

S = [mean(S0,'all'), mean(S1,'all'), mean(S2,'all'), mean(S3,'all')]; 
disp(S ./ mean(S0,'all'));

figure;
subplot(2,2,1); imshow(S0);
subplot(2,2,2); imshow(S1);
subplot(2,2,3); imshow(S2);
subplot(2,2,4); imshow(S3);











