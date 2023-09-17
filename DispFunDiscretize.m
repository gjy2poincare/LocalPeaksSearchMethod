function [MCP, MKI, DISP] = DispFunDiscretize(f, MPM, SPM) 
%ÆµÉ¢º¯ÊýÀëÉ¢»¯¡¾Dispersion Equation Discretization¡¿
cpa = SPM(1); dcp = SPM(2); cpb = SPM(3); kia = SPM(4); dki = SPM(5); kib = SPM(6);
%
w = 2*pi*f;
%%
MCP = cpa:dcp:cpb; MKI = (kia-dki):dki:(kib+dki);
DISP = zeros(length(MCP),length(MKI));
for ii = 1:length(MCP)
    for jj = 1:length(MKI)
        k = w/MCP(ii)+1i*MKI(jj);
        DISP(ii,jj) = DispFunction(f, MPM, k);
    end
end
end