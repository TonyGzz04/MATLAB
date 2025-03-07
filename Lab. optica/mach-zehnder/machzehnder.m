
clc; clear; close all;

% Reset default text properties
set(0, 'DefaultTextFontSize', 12);
set(0, 'DefaultTextFontWeight', 'normal');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');

% Reset default axes properties
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultAxesFontWeight', 'normal');

% Reset default figure properties
set(0, 'DefaultFigureColor', get(0, 'DefaultFigureColor'));
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

p1 = im2double(rgb2gray(imread("p1.JPG")));
p2 = im2double(rgb2gray(imread("p2.JPG")));

figure;
subplot(2,2,1);
imagesc(p1);
title('\textbf{Imagen tomada desde el punto D}', 'FontSize', 17);
xlabel('Pixeles', 'FontSize', 17); 
ylabel('Pixeles', 'FontSize', 17);
subplot(2,2,2);
imagesc(p2);
title('\textbf{Imagen tomada desde el punto C}', 'FontSize', 17);
xlabel('Pixeles', 'FontSize', 17); 
ylabel('Pixeles', 'FontSize', 17);

subplot(2,2,3);
plot(p1(3072/2,:));
title('\textbf{Corte transversal del patr\''on de interferencia del punto D}', 'FontSize', 17);
xlabel('Pixeles', 'FontSize', 17); 
ylabel('Intensidad normalizada', 'FontSize', 17);
subplot(2,2,4);
plot(p2(3072/2,:));
title('\textbf{Corte transversal del patr\''on de interferencia del punto C}', 'FontSize', 17);
xlabel('Pixeles', 'FontSize', 17); 
ylabel('Intensidad normalizada', 'FontSize', 17);


% fplot(@(t) cos(t/150+3*pi/4).^2);

print(gcf, 'fig1.svg', '-dsvg');

%%

close all;

% Ventana numerica
N = 5568; 
w0 = 3e-2;    % [m]
xmax = 0.015; 
lambda = 633E-9;   % [m]
xs = xmax*(2/N)*(-N/2:N/2-1);
ys = xmax*(2/N)*(-N/2:N/2-1);
[Xs, Ys] = meshgrid(xs, ys);

GB = @(X,Y,w) exp(-(X.^2 + Y.^2)/w^2); 
BS = (1/sqrt(2))*[1 1i; 1i 1];

% semilla
U = 0.79 .* GB(Xs, Ys, w0);
umax = max(abs(U(:)).^2);

% Primer divisor
E3 = BS(1,1)*U + BS(1,2)*0;
E4 = BS(2,1)*U + BS(2,2)*0;

% Espejos
E5 = exp(1i*pi)*E4;
E6 = exp(1i*pi)*E3;

% Placa de vidrio
placa = 0.0;     % (2*pi/lambda)*n*d
E5(N/2:end, :) = exp(1i*placa) * E5(N/2:end, :);

% Vector de propagacion -- up
k = 2*pi/lambda;
thzu = 0.009*(pi/180);
thTu = 0.0;
kxu = k * sin(thzu) * cos(thTu);
kyu = k * sin(thzu) * sin(thTu);
kzu = k * cos(thzu);

% Vector de propagacion -- down
k = 2*pi/lambda;
thzd = -thzu;
thTd = 0.0;
kxd = k * sin(thzd) * cos(thTd);
kyd = k * sin(thzd) * sin(thTd);
kzd = k * cos(thzd);

z=0;
E5 = exp(1i*(kxu*Xs + kyu*Ys + kzu*z))  .* E5;
E6 = exp(1i*(kxd*Xs + kyd*Ys + kzd*z))  .* E6;

% Lente delgada
foco = -1;   % [m]
lente = (k/(2*foco)) * (Xs.^2 + Ys.^2);
%E5 = exp(1i*lente) .* E5;

% Segundo divisor
E7 = BS(1,1)*E5 + BS(1,2)*E6;
E8 = BS(2,1)*E5 + BS(2,2)*E6;

I7 = abs(E7).^2;

I = zeros(size(I7));
a = 215;

for i = 1:size(I7,2)-a
    I(:,i) = I7(:,i+a);
end
I(:,size(I7,2)-a-1:end) = 0;

I(:,1:515) = 0;
I(:,5035:end) = 0;


figure;
imagesc(I); ylim([928 4640])
colormap('hot');

a = N/2;

figure; hold on;
plot(I(3072/2,:), '--');
plot(p1(3072/2,:));
xlim([0 size(I,1)]); 
title('\textbf{Corte transversal de la imagen desde el Puerto 1}', 'FontSize', 17); 
xlabel('Pixeles', 'FontSize', 17); ylabel('Intensidad normalizada', 'FontSize', 17);
legend('Te\''orica', 'Experimental', 'FontSize', 14);
hold off;

print(gcf, 'p1.svg', '-dsvg');

error = p1(3072/2,:) - I(3072/2,:);
vari = (error).^2 ./ size(p1,2);

%% FILTRO

clear; clc; close all;

nor = im2double(rgb2gray(imread("normal1.JPG")));
f1 = im2double(rgb2gray(imread("filtroEnA.JPG")));
f2 = im2double(rgb2gray(imread("filtroEnB.JPG")));

nmax = max(nor,[],'all');
fmax = max(f2,[],'all');

disp(abs(nmax-fmax)/nmax)

maxf1 = max(f1,[],'all');
minf1 = f1(2530, 3075);

disp((maxf1 - minf1) / (maxf1 + minf1));



