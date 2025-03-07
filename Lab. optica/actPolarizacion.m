
clc; clear; close all;

m = 1;      % carga topológica (tiene que ver con la fase del haz)
w_0 = 7.5;    % del Gaussiano

clim = [0 0.2];

% E = @(r,phi) (r ./ w_0).^abs(m) .* exp(-r.^2/w_0^2) .* exp(1i*m*phi);

a = 10;
grid_size = 2^10;
x = linspace(-a, a, grid_size);
y = linspace(-a, a, grid_size);
[X, Y] = meshgrid(x,y);
R = sqrt(X.^2 + Y.^2);

phi = 0; % si se grafica la intensidad no importa

E = (R ./ w_0).^abs(m) .* exp(-R.^2/w_0^2) .* exp(1i*m*phi);

figure;
subplot(3,2,1); imagesc(x,y,abs(E).^2);
colorbar;
caxis(clim);
title('Intensidad inicial');


% Posibles polarizaciones iniciales
H = [1; 0];
V = [0; 1]; 
D = 1/sqrt(2) .* [1; 1];
A = 1/sqrt(2) .* [1; -1];

pol_in = D;

Ex = E .* pol_in(1,1);
Ey = E .* pol_in(2,1);

E = sqrt(Ex.^2 + Ey.^2);

subplot(3,2,3); imagesc(x,y,abs(Ex).^2);
subplot(3,2,4); imagesc(x,y,abs(Ey).^2);
colorbar;
caxis(clim);
title('Intensidad con polarización inicial');


% Pasar por polarizador a cierto ángulo teta
theta = pi/6;
rot = [cos(theta)^2, cos(theta)*sin(theta); 
       cos(theta)*sin(theta), sin(theta)^2];

Ex = Ex.*rot(1,1) + Ey.*rot(2,1);
Ey = Ex.*rot(1,2) + Ey.*rot(2,2);
E = sqrt(Ex.^2 + Ey.^2);

subplot(3,2,5); imagesc(x,y,abs(Ex).^2);
subplot(3,2,6); imagesc(x,y,abs(Ey).^2);
colorbar;
caxis(clim);
title('Intensidad con polarización del polarizador');


I_0 = sum(abs(E).^2,'all');

figure; hold on;
for theta = 0:pi/180:2*pi
    rot = [cos(theta)^2, cos(theta)*sin(theta); 
       cos(theta)*sin(theta), sin(theta)^2];
    Ex_pol = Ex.*rot(1,1) + Ey.*rot(2,1);
    Ey_pol = Ex.*rot(1,2) + Ey.*rot(2,2);
    E = sqrt(Ex_pol.^2 + Ey_pol.^2);
    I = abs(E).^2;

    I = sum(I,'all');
    plot(theta, I, '.','Color','b');
end

fplot(@(theta) I_0.*cos(theta).^2, [0 2*pi]);
hold off;


%% Ahora con retardador de lambda/4

close all; clc;

E = (R ./ w_0).^abs(m) .* exp(-R.^2/w_0^2) .* exp(1i*m*phi);

figure;
subplot(3,2,1); imagesc(x,y,abs(E).^2);
colorbar;
caxis(clim);
title('Intensidad inicial');


% Posibles polarizaciones iniciales
H = [1; 0];
V = [0; 1]; 
D = 1/sqrt(2) .* [1; 1];
A = 1/sqrt(2) .* [1; -1];

pol_in = D;

Ex = E .* pol_in(1,1);
Ey = E .* pol_in(2,1);

E = sqrt(Ex.^2 + Ey.^2);

subplot(3,2,3); imagesc(x,y,abs(Ex).^2);
subplot(3,2,4); imagesc(x,y,abs(Ey).^2);
colorbar;
caxis(clim);
title('Intensidad con polarización inicial');

LR_fourth = @(theta) 1/sqrt(2) .* [1-1i*cos(2*theta), -2i*cos(theta)*sin(theta);
                       -2i*cos(theta)*sin(theta), 1+1i*cos(2*theta)];
rot_pol = @(theta) [cos(theta).^2, cos(theta)*sin(theta); 
       cos(theta)*sin(theta), sin(theta).^2];

%{
figure; hold on;
for theta = 0:pi/180:2*pi
    
    LR = LR_fourth(theta);
    Ex_ret = Ex.*LR(1,1) + Ey.*LR(2,1);
    Ey_ret = Ex.*LR(1,2) + Ey.*LR(2,2);
    E = sqrt(Ex_ret.^2 + Ey_ret.^2);
    I = abs(E).^2;
    I = sum(I,'all');

    rot = rot_pol(theta);
    Ex_pol = Ex_ret.*rot(1,1) + Ey_ret.*rot(2,1);
    Ey_pol = Ex_ret.*rot(1,2) + Ey_ret.*rot(2,2);
    E = sqrt(Ex_pol.^2 + Ey_pol.^2);
    I_pol = abs(E).^2;
    I_pol = sum(I_pol,'all');
    
    % subplot(1,2,1); hold on;
    plot3(theta, I, '.','Color','b');
    plot3(theta, I_pol, '.','Color','g');

end
hold off;
%}

figure; hold on;
for theta_pol = 0:10*pi/180:2*pi
    for theta_ret = 0:10*pi/180:2*pi
        LR = LR_fourth(theta_ret);
        Ex_ret = Ex.*LR(1,1) + Ey.*LR(2,1);
        Ey_ret = Ex.*LR(1,2) + Ey.*LR(2,2);
        % E_ret = sqrt(Ex_ret.^2 + Ey_ret.^2);
        % I_ret = abs(E).^2;
        % I_ret = sum(I_ret,'all');

        rot = rot_pol(theta_pol);
        Ex_pol = Ex_ret.*rot(1,1) + Ey_ret.*rot(2,1);
        Ey_pol = Ex_ret.*rot(1,2) + Ey_ret.*rot(2,2);
        E_pol = sqrt(Ex_pol.^2 + Ey_pol.^2);
        I_pol = abs(E_pol).^2;
        I_pol = sum(I_pol,'all');

        plot3(theta_pol, theta_ret, I_pol, '.', 'Color', 'b');

    end
    fprintf('Faltan %f.2 pi\n', theta_pol/pi);
end
xlabel('Theta polarizador');
ylabel('Theta retardador');
zlabel('Intensidad final');
hold off;
        







