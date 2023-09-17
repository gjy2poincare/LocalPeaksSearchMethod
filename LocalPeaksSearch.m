function [cp, ki] = LocalPeaksSearch(MCP, MKI, DISP, kur)
%局部峰值搜索【Peaks Detection Method】
%%
FD = 1./(1+abs(DISP));
wr = 3; wc = 3; %局部窗宽度：奇数【Local Window Width: Odd Number】
[M,N] = size(FD); 
if N==1; wc = 1; end
n = 1; cp = []; ki = [];
for ii = 1:M
    if ii+wr-1>M; break; end
    for jj = 1:N
        if jj+wc-1>N; break; end
        fdw = FD(ii:ii+wr-1,jj:jj+wc-1); %局部窗【Local Window Width】
        if wc==1
            [fdw_max,Ir] = max(fdw); Ic = 1;
        else
            [~,Ir] = max(max(fdw')); [fdw_max,Ic] = max(max(fdw));
        end
        if Ir==round((wr+1)/2) && Ic==round((wc+1)/2)
            %峰值提取【Peaks Extracting】
            cp2 = MCP(ii+Ir-1); ki2 = MKI(jj+Ic-1); %峰值【Peak Value】
            kurtosis = fdw_max/((sum(sum(fdw))-fdw_max)/(wr*wc-1))-1;
            if kurtosis<kur, continue; end
            cp(n) = cp2;  ki(n) = ki2;
            n = n+1;
        end
    end
end
end