%可视化处理
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
%频散曲线坐标缩比比例尺：1:1E3
fh = fh/1E3; fnd = fnd/1E3; cp_ps = cp_ps/1E3; %cr = cr/1E3; ct = ct/1E3; 
%dispersion curve
figure(1001);
hold on;
plot(fh,cp_ps,'b');
% axis([0 floor(fnd*(2*h)) 0 20]);
title('\fontsize{12}Lamb波频散曲线');
xlabel('\fontsize{12}\fontname{宋体}频厚积\fontname{Times New Roman}\rm(MHz-mm)'); 
ylabel('\fontsize{12}\fontname{宋体}相速度\fontname{Times New Roman}\rm(km/s)');
hold off;
%attenuation curve
figure(1002);
hold on;
plot(fh,ki_ps,'b');
% axis([0 floor(fnd*(2*h)) 0 50]);
title('\fontsize{12}Lamb波衰减曲线');
xlabel('\fontsize{12}\fontname{宋体}频厚积\fontname{Times New Roman}\rm(MHz-mm)'); 
ylabel('\fontsize{12}\fontname{宋体}衰减系数\fontname{Times New Roman}\rm(Np/m)');
hold off;