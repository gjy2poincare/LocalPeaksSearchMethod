function [cp_0, ki_0, cp_1, ki_1, cp_2, ki_2] = RootsMerge(f, SW, cp1, ki1, cp2, ki2)
%根值融合【Roots Merge】
w = 2*pi*f; nSW = length(SW); 
%%
disp('----------------------------------------');
disp('【Roots Merge】');
%% %根值排序【Roots Sort】
cp1 = cp1(~isnan(abs(cp1))); ki1 = ki1(~isnan(abs(ki1)));
cp2 = cp2(~isnan(abs(cp2))); ki2 = ki2(~isnan(abs(ki2)));
cp_0 = [cp1 cp2]; ki_0 = [ki1 ki2];
for ii = 1:length(cp_0)
    for jj=ii+1:length(cp_0)
        if cp_0(ii)>cp_0(jj)
            temp1 = cp_0(ii);  temp2 = ki_0(ii);
            cp_0(ii) = cp_0(jj); ki_0(ii) = ki_0(jj);
            cp_0(jj) = temp1;  ki_0(jj) = temp2;
        end
    end
end
%% %重根剔除【Kill Same Roots】
cp_1 = cp_0; ki_1 = ki_0;
for ii = 1:length(cp_1)
    for jj=ii+1:length(cp_1)
        if abs(cp_1(ii)-cp_1(jj))<1E-0 && abs(ki_1(ii)-ki_1(jj))<1E-1
            cp_1(ii) = NaN; ki_1(ii) = NaN;
        end
    end
end
cp_1 = cp_1(~isnan(cp_1)); ki_1 = ki_1(~isnan(ki_1));
%% %奇异根剔除【Singular Roots Cancel: Seems not work well】
cp_2 = cp_1; ki_2 = ki_1;
for n = 1:nSW
    if n==3 || n==6, continue; end % 排除Rayleigh波【Exclude Rayleigh Wave】
    kn = SW(n);  cpn = w/real(kn); kin = imag(kn);
    for ii = 1:length(cp_2) 
        if abs(cp_2(ii)-cpn)<1E-0 && abs(ki_2(ii)-kin)<1E-1
            cp_2(ii) = NaN; ki_2(ii) = NaN;
        end
    end
end
cp_2 = cp_2(~isnan(cp_2)); ki_2 = ki_2(~isnan(ki_2));
%% %结果显示【Results Output】
if isempty(cp_0) || isempty(ki_0), warning('No Roots Found!'); end
disp('Raw Results: ');
if length(cp_0)>=1
    for ii = 1:length(cp_0)
        cpx = cp_0(ii); kix = ki_0(ii);
        disp(['>>cp = ', num2str(cpx), '  ki = ', num2str(kix)]);
    end
else
    warning('No Roots Found!');
end
disp(['toc: ', num2str(toc), 's']);
end