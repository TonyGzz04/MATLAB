

clc; clear; close all;

cien_franjas = false;


if cien_franjas == true
    img1 = imread('100_0.JPG');   % Read the image
    img2 = imread('100_1.JPG');
    img3 = imread('100_2.JPG');
    img4 = imread('100_3.JPG');
else
    img1 = imread('65_0.JPG');   % Read the image
    img2 = imread('65_1.JPG');
    img3 = imread('65_2.JPG');
    img4 = imread('65_3.JPG');
    
    simg1 = imread('s65_0.JPG');
    simg2 = imread('s65_1.JPG');
    simg3 = imread('s65_2.JPG');
    simg4 = imread('s65_3.JPG');
end


 
img1 = im2double(rgb2gray(img1)); % Convert to grayscale
img2 = im2double(rgb2gray(img2));
img3 = im2double(rgb2gray(img3));
img4 = im2double(rgb2gray(img4));

simg1 = im2double(rgb2gray(simg1));
simg2 = im2double(rgb2gray(simg2));
simg3 = im2double(rgb2gray(simg3));
simg4 = im2double(rgb2gray(simg4));

alpha = pi/2;

% phi = atan2( -(img2 - img4) , (img1 - img3));
% phi = atan((img1-img3)./(sqrt(2)*(img2-img4)) - 1);

%phi = atan2( (img1 - img3 - sqrt(2).*(img2-img4)) , (sqrt(2).*(img2-img4)) );

phi = atan2( -(img2 - img4) , (img1 - img3));
sphi = atan2( -(simg2 - simg4) , (simg1 - simg3));

figure;
imshow(phi);  


dphi1 = img2 - img1;
dphi2 = img3 - img2;
dphi3 = img4 - img3;

sdphi1 = simg2 - simg1;
sdphi2 = simg3 - simg2;
sdphi3 = simg4 - simg3;

w1 = atan2(sin(dphi1) , cos(dphi1));
w2 = atan2(sin(dphi2) , cos(dphi2));
w3 = atan2(sin(dphi3) , cos(dphi3));

sw1 = atan2(sin(sdphi1) , cos(sdphi1));
sw2 = atan2(sin(sdphi2) , cos(sdphi2));
sw3 = atan2(sin(sdphi3) , cos(sdphi3));


% w = [w1; w2; w3];

phi = (phi - min(phi(:))) / (max(phi(:)) - min(phi(:))) * 255;
sphi = (sphi - min(sphi(:))) / (max(sphi(:)) - min(sphi(:))) * 255;
phi8 = uint8(phi);
sphi8 = uint8(sphi);

figure;
imshow(phi8 - sphi8);  

final = phi - sphi;

% q = unwrap(phi, [], 2);
% plot3(0:size(q,1)-1, q(:,1), q(:,2), '.');


q = unwrap_phase(final);
% ssq = unwrap_phase(sphi);

q8 = uint8(q);
% ssq8 = uint8(ssq);

figure;
imshow(q8)






%%

denoised_unwrap(q)


% Called functions in execution: phase_wrap, TV_min.
function denoised_unwrapped_phase = denoised_unwrap(wrapped_phase)
% Calculating denoised derivatives (components of gradient) of unwrapped phase by using TV (Total Variation) minization algorithm [1][2]
gradient_x_wrapped_phase = diff(wrapped_phase,1,2); gradient_y_wrapped_phase = diff(wrapped_phase,1,1); 
gradient_x_unwrapped_phase = phase_wrap(gradient_x_wrapped_phase); gradient_y_unwrapped_phase = phase_wrap(gradient_y_wrapped_phase); % Implementation of Eq.(3) and (4) of Ref.[1].
%figure(11); subplot(2,1,1); imagesc(gradient_x_wrapped_phase); colormap(gray); colorbar;
subplot(2,1,2); imagesc(gradient_x_unwrapped_phase); colormap(gray); colorbar;
denoised_gradient_x_unwrapped_phase = TV_min(gradient_x_unwrapped_phase); denoised_gradient_y_unwrapped_phase = TV_min(gradient_y_unwrapped_phase); % Applying TV-minimization algorithm to derivatives.
%figure(12); subplot(2,1,1); surf(gradient_x_unwrapped_phase(:,:)); subplot(2,1,2); surf(denoised_gradient_x_unwrapped_phase(:,:));
%Integration of denoised gradient to obtain denoised-unwrapped phase
gradient_x = denoised_gradient_x_unwrapped_phase; gradient_y = denoised_gradient_y_unwrapped_phase;
denoised_unwrapped_phase = zeros(size(wrapped_phase));
denoised_unwrapped_phase(1,1) = wrapped_phase(1,1);
for j = 2 : size(denoised_unwrapped_phase,2)
    denoised_unwrapped_phase(1,j) = denoised_unwrapped_phase(1,j-1) + gradient_x(1,j-1); % Integration along 1st row.
end
for i = 2 : size(denoised_unwrapped_phase,1)
    denoised_unwrapped_phase(i,1) = denoised_unwrapped_phase(i-1,1) + gradient_y(i-1,1); % Integration along 1st column.
end
for i = 2 : size(denoised_unwrapped_phase,1) % Integration for rest matrix body.
    for j = 2: size(denoised_unwrapped_phase,2)
        denoised_unwrapped_phase(i,j) = denoised_unwrapped_phase(i,j-1) + gradient_x(i,j-1);
        denoised_unwrapped_phase(i,j) = denoised_unwrapped_phase(i,j) + (denoised_unwrapped_phase(i-1,j) + gradient_y(i-1,j));
        denoised_unwrapped_phase(i,j) = denoised_unwrapped_phase(i,j)/2; % This mean-integration gives more reliable result than using derivative of only one direction.
    end
end
end