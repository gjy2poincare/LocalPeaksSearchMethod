%---------------------------------------------------------------------------------------------
%                �ֲ��������� -- LPS (Local-Peaks Search) Method
% --------------------------------------------------------------------------------------------
% ��Program�� Local-Peaks Search Method��LPS Method�� - Dispersion equation solving method.
%         for the following Models:
%         (1) Free or Fluid-loaded, Single- or Double-Layer, Elastic or Viscoelastic Plates;
%         (2) Free or Fluid-loaded, Single- or Double-Layer, Elastic or Viscoelastic Cylindrical Shells.
% ��Developer�� Jiayuan Gong, jygong33@qq.com
% ��Notes�� Anyone who optimizes the program, please share it as a PR.
%---------------------------------------------------------------------------------------------
%%
clear all; close all; diary off; clc; tic;
diary diary.log;
disp('Initializing...');
szCheck = questdlg('To check roots manually?', 'Choice', 'Yes', 'No', 'Yes'); pause(0.1);
%% %�����趨��Parameters Setting��
global szFun szMethod szBC szMode nMode_cs szChkCho
szChkCho = 'Yes';
%------------------------------------------------------------------------------
%(ճ)���Բ��ϲ������ǰ����,�����ڵ���ǰ�
%��(Visco-) Elastic Material Parameters: Shell & Plate, Also Used for Single Layer��
% rowem = 7.84E3; Eem0 = 21E10; ytaem = 0; sigmaem = 0.28;
rowem = 1.09E3; Eem0 = 3E7; ytaem = 0.249; sigmaem = 0.49;
% rowem = 1.1E3; Eem0 = 1.4E8; ytaem = 0.23; sigmaem = 0.49;
Eem = Eem0*(1-1i*ytaem);
lamdaem = Eem*sigmaem/((1+sigmaem)*(1-2*sigmaem)); miuem = Eem/(2*(1+sigmaem));
clem = sqrt((lamdaem+2*miuem)/rowem); ctem = sqrt(miuem/rowem);
crem = RayleighWave(clem, ctem);
%------------------------------------------------------------------------------
%(ճ)���Բ��ϲ��������ǲ����,ֻ����˫��ǰ�
%��(Visco-) Elastic Material Parameters: Shell & Plate, Only Used for Double Layer��
rowvm = 1.09E3; Evm0 = 3E7; ytavm = 0.249; sigmavm = 0.49;
% rowvm = 1.1E3; Evm0 = 1.4E8; ytavm = 0.23; sigmavm = 0.49;
% rowvm = 7.84E3; Evm0 = 21E10; ytavm = 0; sigmavm = 0.28;
Evm = Evm0*(1-1i*ytavm);
lamdavm = Evm*sigmavm/((1+sigmavm)*(1-2*sigmavm)); miuvm = Evm/(2*(1+sigmavm));
clvm = sqrt((lamdavm+2*miuvm)/rowvm); ctvm = sqrt(miuvm/rowvm);
crvm = RayleighWave(clvm, ctvm);
%------------------------------------------------------------------------------
%������ʲ�����Fluid Media Parameters��
row1 = 1E3; c1 = 1.5E3;
row2 = 1E3; c2 = 1.5E3;
%------------------------------------------------------------------------------
%�ߴ������Dimension Parameters��
%(������)��ߴ硾(Coated) Plate Dimension��
dem = 0.05; dvm = 0.01; %���
%(������)Բ���ǳߴ硾(Coated) Cylindrical Shell Dimension��
a_inner = 0.01; b_middle = 0.02; c_outer = 0.03; %�뾶
%--------------------------------------------------------------------------
%����ģ�͡�Computation Model��
szFun = 'PlanarPlate'; %'CylindricalShell';
%�߽�������Boundary Condition�� %Va: vacuum, So: solid, Fl: fluid
szBC = 'Va-So-Va'; %'Va-So-Va', 'Va-So-Fl', 'Fl-So-Fl';
%'Va-So-So-Va', 'Va-So-So-Fl', 'Fl-So-So-Fl'
%ģ̬ѡ��Mode Choice��
szMode = 'A/T';% ; 'S/L', 'A/T', 'F'
nMode_cs = 1; %����ģ̬������Model Level of Cylindrical Shell��
%--------------------------------------------------------------------------
%���������Computation Parameters��
err = 1E-20; %���ȡ�Precision��
kur = 1E-3; %��ȡ�Kurtosis��
szMethod = 'Mixed'; %'Muller', 'DomainRefine', 'Mixed' %��ֵ������Roots-Finding Method��
F = 1E2:1E2:10E3; %Ƶ�ʷ�Χ��Frequency Range��
cpa = 1E1; dcp = 1E0; cpb = 5.5E2; %���ٶȷ�Χ��Phase Velocity Range��
kia = -1E0; dki = 1E0; kib = 5E2; %�����鲿��Χ��Imaginary of Wavenumber Range��
%------------------------------------------------------------------------------
%�������桾Parameters Save��
save('data.mat', 'row1', 'c1', 'row2', 'c2');
save('data.mat', 'rowvm', 'Evm0', 'ytavm', 'sigmavm', 'rowem', 'Eem0', 'ytaem', 'sigmaem', '-append');
save('data.mat', 'lamdavm', 'miuvm', 'clvm', 'ctvm', 'lamdaem', 'miuem', 'clem', 'ctem', '-append');
save('data.mat', 'dem', 'dvm', 'a_inner', 'b_middle', 'c_outer', '-append');
save('data.mat', 'F', 'szFun', 'szBC', 'szMode', '-append');
%% %Ƶɢ���������Roots Find of Dispersion Equation��
%���ϲ�����Material Parameters��
MPM = [row1, c1, rowvm, lamdavm, miuvm, dvm, rowem, lamdaem, miuem, dem, row2, c2,...
    a_inner, b_middle, c_outer];
%�������ٶȡ�Singular Phase Velocities��
SVM = [clem, ctem, crem, clvm, ctvm, crvm];
%����������Search Parameters��
SPM = [cpa, dcp, cpb, kia, dki, kib];
save('data.mat', 'MPM', 'SVM', 'SPM', '-append');
%%
nL = 1E1;
cp_ps0 = zeros(length(F), nL); cp_ps0(:,:) = NaN; ki_ps0 = cp_ps0;
cp_ps1 = cp_ps0; ki_ps1 = cp_ps0; cp_ps2 = cp_ps0;
ki_ps2 = cp_ps0; cp_ps3 = cp_ps0; ki_ps3 = cp_ps0;
MSW = zeros(length(F), 6); MSW(:, :) = NaN;
n = 1;
disp(['toc: ', num2str(toc), 's']);
disp('*******************************************************');
for f = F
    diary on;
    disp(['f = ', num2str(f), 'Hz']);
    w = 2*pi*f;
    %% %ȫ��������Global Search��
    [cp1, ki1]= GlobalDomainSearch(f, MPM, SPM, kur, err);
    %% %�ֲ�������Local Search��
    [cp2, ki2, SW]= LocalDomainSearch(f, MPM, SPM, SVM, kur, err);
    MSW(n, 1:length(SW)) = w./real(SW)+1i*imag(SW); %���첨����ʵ��Ϊ���٣��鲿Ϊ˥��ϵ��
    %��Singular Wavenumbers(SW)��real->velocity, imag.->attenuation��
    %% %��ֵ�ϲ���Roots Merge��
    [cp_0, ki_0, cp_1, ki_1, cp_2, ki_2] = RootsMerge(f, SW, cp1, ki1, cp2, ki2);
    cp_ps0(n, 1:length(cp_0)) = cp_0; ki_ps0(n, 1:1:length(ki_0)) = ki_0;
    cp_ps1(n, 1:1:length(cp_1)) = cp_1; ki_ps1(n, 1:1:length(ki_1)) = ki_1;
    cp_ps2(n, 1:1:length(cp_2)) = cp_2; ki_ps2(n, 1:1:length(ki_2)) = ki_2;
    save('data.mat', 'cp_ps0', 'ki_ps0', 'cp_ps1', 'ki_ps1', 'cp_ps2', 'ki_ps2', 'MSW', '-append');
    %% %��ֵУ�顾Roots Check��
    if strcmp(szCheck, 'Yes') && ~strcmp(szChkCho, 'Never')
        [cp_3, ki_3] = RootsCheck(f, MPM, cp_0, ki_0);
        cp_ps3(n,1:1:length(cp_3)) = cp_3; ki_ps3(n, 1:1:length(ki_3)) = ki_3;
        save('data.mat', 'cp_ps3', 'ki_ps3', '-append');
    end
    %% %����Post Process��
    disp('*******************************************************');
    diary off; n = n+1;
end
disp('done!');