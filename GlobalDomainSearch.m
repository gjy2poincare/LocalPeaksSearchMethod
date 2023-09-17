function [cp, ki]= GlobalDomainSearch(f, MPM, SPM, kur, err)
%ȫ�������еľֲ���������Local-Peaks Search in Global Domain��
cpa = SPM(1); dcp = SPM(2); cpb = SPM(3); kia = SPM(4); dki = SPM(5); kib = SPM(6);
%% 
disp('��Global Domain Search��');
%% Ƶɢ������ɢ����Dispersion Funtion Discretize��
disp('=��Function Discretizing...');
[MCP, MKI, DISP] = DispFunDiscretize(f, MPM, SPM);
disp(['toc: ', num2str(toc), 's']);
%% �ֲ���������Local-Peaks Search��
disp('=��Local-Peaks Searching...');
[cpx, kix] = LocalPeaksSearch(MCP, MKI, DISP, kur); disp(['toc: ', num2str(toc), 's']);
%% ��ֵ��ȡ��Roots Find��
disp('=��Local-Roots Computing...');
[cp, ki] = RootsFind(f, MPM, cpx, kix, dcp, dki, err);
disp(['toc: ', num2str(toc), 's']);
%%
cp = cp(~isnan(cp)); ki = ki(~isnan(ki));
if isempty(cp) || isempty(ki), warning('No Roots Found!'); end
end
