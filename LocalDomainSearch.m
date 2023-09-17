function [cp, ki, SW]= LocalDomainSearch(f, MPM, SPM, SVM, kur, err)
%局部区域中的局部峰搜索【Local-Peaks Search in Local Domain】
cpa = SPM(1); dcp = SPM(2); cpb = SPM(3); kia = SPM(4); dki = SPM(5); kib = SPM(6);
clem = SVM(1); ctem = SVM(2); crem = SVM(3); clvm = SVM(4); ctvm = SVM(5); crvm = SVM(6);
global szBC 
%
w = 2*pi*f;  
%频散方程的奇异波数【Singular Wavenumbers of Dispersion Equation】   
if strcmp(szBC, 'Va-So-Va') || strcmp(szBC, 'Va-So-Fl') || strcmp(szBC, 'Fl-So-Fl')
    SW = w./[clem, ctem, crem];
elseif strcmp(szBC, 'Va-So-So-Va') || strcmp(szBC, 'Va-So-So-Fl')|| strcmp(szBC, 'Fl-So-So-Fl');
    SW = w./[clem, ctem, crem, clvm, ctvm, crvm];
end
nSW = length(SW);
%%
disp('----------------------------------------');
disp('【Local Domain Search】');
cnt = 1;
cp = zeros(1,2E1); cp( :)= NaN; ki = cp;
for n = 1:nSW    
    kx = SW(n);  cp0 = w/real(kx); ki0 = imag(kx);
    disp(['>>Singular Wavenumber: ','cp0 = ', num2str(cp0), '; ki0 = ', num2str(ki0)]);
    cpa = cp0-dcp; cpb = cp0+dcp; kia = ki0-dki; kib = ki0+dki;
    dcp2 = min(dcp/20, 1); dki2 = min(dki/20, 0.1);
    SPM = [cpa, dcp2, cpb, kia, dki2, kib];
    %% 频散函数离散化【Dispersion Funtion Discretize】
    disp('=》Function Discretizing...');
    [MCP, MKI, DISP] = DispFunDiscretize(f, MPM, SPM);
    disp(['toc: ', num2str(toc), 's']);
    %% 局部峰搜索【Local-Peaks Search】
    disp('=》Local-Peaks Searching...');
    [cpx, kix] = LocalPeaksSearch(MCP, MKI, DISP, kur); 
    disp(['toc: ', num2str(toc), 's']);
    %% 根值求取【Roots Find】
    disp('=》Local-Roots Computing...');
    [cp2, ki2] = RootsFind(f, MPM, cpx, kix, dcp2, dki2, err); 
    cp2 = cp2(~isnan(cp2)); ki2 = ki2(~isnan(ki2));
    if isempty(cp2) || isempty(ki2)
        disp('>>No Roots Find!'); disp(['toc: ', num2str(toc), 's']);
        continue; 
    end
    disp(['toc: ', num2str(toc), 's']);
    for ii = 1:length(cp2)
        cp(cnt) = cp2(ii); ki(cnt) = ki2(ii);
        cnt = cnt+1;
    end
end
end