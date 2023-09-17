%���ӻ�����
clear all; 
% close all; 
clc;
%% 
DATA = open('data_Va_So_Va_A.mat');
%%
F = DATA.F; 
% h = DATA.h/2; 
h = DATA.dem/2; 
% h = 0.01/2;
cp_ps = DATA.cp_ps; ki_ps = DATA.ki_ps; 
%%
fh = F*(2*h); fnd = F(end);
%Ƶɢ�����������ȱ����ߣ�1:1E3
fh = fh/1E3; fnd = fnd/1E3; cp_ps = cp_ps/1E3; %cr = cr/1E3; ct = ct/1E3; 
%dispersion curve
figure(1001);
hold on;
plot(fh,cp_ps,'b');
% axis([0 floor(fnd*(2*h)) 0 20]);
title('\fontsize{12}Lamb��Ƶɢ����');
xlabel('\fontsize{12}\fontname{����}Ƶ���\fontname{Times New Roman}\rm(MHz-mm)'); 
ylabel('\fontsize{12}\fontname{����}���ٶ�\fontname{Times New Roman}\rm(km/s)');
hold off;
%attenuation curve
figure(1002);
hold on;
plot(fh,ki_ps,'b');
% axis([0 floor(fnd*(2*h)) 0 50]);
title('\fontsize{12}Lamb��˥������');
xlabel('\fontsize{12}\fontname{����}Ƶ���\fontname{Times New Roman}\rm(MHz-mm)'); 
ylabel('\fontsize{12}\fontname{����}˥��ϵ��\fontname{Times New Roman}\rm(Np/m)');
hold off;