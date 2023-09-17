function [cp, ki]= GlobalDomainSearch(f, MPM, SPM, kur, err)
%全局区域中的局部峰搜索【Local-Peaks Search in Global Domain】
cpa = SPM(1); dcp = SPM(2); cpb = SPM(3); kia = SPM(4); dki = SPM(5); kib = SPM(6);
%% 
disp('【Global Domain Search】');
%% 频散函数离散化【Dispersion Funtion Discretize】
disp('=》Function Discretizing...');
[MCP, MKI, DISP] = DispFunDiscretize(f, MPM, SPM);
disp(['toc: ', num2str(toc), 's']);
%% 局部峰搜索【Local-Peaks Search】
disp('=》Local-Peaks Searching...');
[cpx, kix] = LocalPeaksSearch(MCP, MKI, DISP, kur); disp(['toc: ', num2str(toc), 's']);
%% 根值求取【Roots Find】
disp('=》Local-Roots Computing...');
[cp, ki] = RootsFind(f, MPM, cpx, kix, dcp, dki, err);
disp(['toc: ', num2str(toc), 's']);
%%
cp = cp(~isnan(cp)); ki = ki(~isnan(ki));
if isempty(cp) || isempty(ki), warning('No Roots Found!'); end
end
