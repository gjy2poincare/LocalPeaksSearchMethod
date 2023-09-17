function [cp, ki] = RootsFind(f, MPM, cpx, kix, dcp, dki, err)
%根值计算：以搜索到的局部峰值为初值【Roots Computation: Based on Local Peaks Searched】
%%
global szMethod
switch szMethod     
    case 'DomainRefine'
        [cp, ki] = RootsFind_DR(f, MPM, cpx, kix, dcp, dki, err);
    case 'Muller'
        [cp, ki] = RootsFind_ML(f, MPM, cpx, kix, dcp, dki, err);
    case 'Mixed'
        [cp, ki, szSuFa] = RootsFind_ML(f, MPM, cpx, kix, dcp, dki, err);        
        if strcmp(szSuFa, 'Fail')
            disp('Muller Failed!');
            [cp, ki] = RootsFind_DR(f, MPM, cpx, kix, dcp, dki, err);
        end        
end
end
%% ###########################################################################################
function [cp, ki] = RootsFind_DR(f, MPM, cpx, kix, dcp, dki, err)
%方程求根：函数零点计算【Roots-Finding of Equation: i.e. Function's Zeros Computation】
w = 2*pi*f; dcp2 = min(dcp/10,1);
%%
cp = NaN; ki = NaN; 
nc = 1;
for n = 1:length(cpx)
    cp2 = cpx(n); ki2 = kix(n);
    %%
    %区域细化，根值求取【Domain Refine, Roots Find】
    dcpr = dcp; dkir = dki; scale = 2;
    while dcpr>err || dkir>err
        cpa = cp2-scale*dcpr; cpb = cp2+scale*dcpr; 
        kia = ki2-scale*dkir; kib = ki2+scale*dkir;
        dcpr = dcpr/10; dkir = dkir/10;
        MCP = cpa:dcpr:cpb; MKI = kia:dkir:kib;
        DISP = zeros(length(MCP),length(MKI));
        for ii = 1:length(MCP)
            for jj = 1:length(MKI)
                k = w/MCP(ii)+1i*MKI(jj);
                DISP(ii,jj) = DispFunction(f, MPM, k);
            end
        end
        FD = 1./(1+abs(DISP));
        [~,pos_cp] = max(max(FD')); [~,pos_ki] = max(max(FD));
        cp2 = MCP(pos_cp); ki2 = MKI(pos_ki);
    end
    disp(['>>cp = ', num2str(cp2), '  ki = ', num2str(ki2)]);
    cp(nc) = cp2; ki(nc) = ki2;
    nc = nc+1;
end
end
%% %###########################################################################################
function [cp, ki, szSuFa] = RootsFind_ML(f, MPM, cpx, kix, dcp, dki, err)
%零点计算:Muller法
szSuFa = 'Success';
w = 2*pi*f;
N = length(cpx);
kk = 1; cp = NaN; ki = NaN;
for n = 1:N
    cp0 = cpx(n); ki0 = kix(n);
    scale = 1.2;
    cpmin = cp0-scale*dcp; cpmax = cp0+scale*dcp; 
    kimin = ki0-scale*dki; kimax = ki0+scale*dki;
    %% %区域细化【Domain Refine】
    dcpr = min(1,dcp/10); dkir = min(0.1,dki/10);
    MCP = cpmin:dcpr:cpmax; MKI = kimin:dkir:kimax;
    DISP = zeros(length(MCP),length(MKI));
    for ii = 1:length(MCP)
        for jj = 1:length(MKI)
            k = w/MCP(ii)+1i*MKI(jj);
            DISP(ii,jj) = DispFunction(f, MPM, k);
        end
    end
    FD = 1./(1+abs(DISP));
    [~,pos_cp] = max(max(FD')); [~,pos_ki] = max(max(FD));
    cp2 = MCP(pos_cp); ki2 = MKI(pos_ki);
    scale = 2;
    cpmin = cp2-scale*dcpr; cpmax = cp2+scale*dcpr; 
    kimin = ki2-scale*dkir; kimax = ki2+scale*dkir;
    %%
    PR = [cpmin,cpmax,kimin,kimax];
    [cpt, kit, szSuFa] = RootsFind_ML2(f, PR, MPM, err);
    if strcmp(szSuFa,'Fail'), return; end
    disp(['>>cp = ', num2str(cpt),';  ki = ', num2str(kit)]);
    cp(kk) = cpt; ki(kk) = kit; kk = kk+1;
end
end
%----------------------------------------------------------------------------------------------------
function [cp, ki, szSuFa] = RootsFind_ML2(f, PR, MPM, err)
cpmin = PR(1); cpmax = PR(2); kimin = PR(3); kimax = PR(4);
%
w = 2*pi*f;
krmin = w/cpmax; krmax = w/cpmin; dkr = krmax-krmin; dki = kimax-kimin;
%%
szSuFa = 'Success';
cp = NaN; ki = NaN; scale = 1.2;
%常规计算【Normal Computation】
N = 4; m = 1; SV = zeros(N-1, N-1);
for ii = 1:N-1
    for jj = 1:N-1
        SV(m) = (krmin+dkr*jj/N)+1i*(kimin+dki*ii/N);
        m = m+1;
    end
end
M = length(SV); tm0= toc;
for ii=1:M
    for jj = 1:M
        for kk = 1:M
            if ii==jj || jj==kk || ii==kk; continue; end
            z0 = SV(ii); z1 = SV(jj); z2 = SV(kk);
            %Muller's Method
            [zn, szSuFa] = muller(f, MPM, z0, z1, z2, err);
            if((toc-tm0)>10), szSuFa = 'Fail'; return; end
            if strcmp(szSuFa, 'Fail'), continue; end
            cpn = w/real(zn); kin = imag(zn);
            if ~isnan(cpn) && ~isnan(kin) && cpn>=cpmin/scale && cpn<=cpmax*scale...
                    && kin>=kimin/scale && kin<=kimax*scale
                cp = cpn; ki = kin; return;
            end
        end
    end
end
%%
%随机计算【Ranom Computing】
tm0 = toc;
if ~isnan(cp)
    return;
else
    N = 1E2;
    dkm = min(dkr,dki);
    for kk = 1:N
        while 1
            pt = rand(1,6);
            p0r = pt(1); p0i = pt(2); p1r = pt(3); p1i = pt(4); p2r = pt(5); p2i = pt(6);
            z0 = (krmin+dkr*p0r)+1i*(kimin+dki*p0i);
            z1 = (krmin+dkr*p1r)+1i*(kimin+dki*p1i);
            z2 = (krmin+dkr*p2r)+1i*(kimin+dki*p2i);
            if abs(z0-z1)>=dkm/1E1 && abs(z1-z2)>=dkm/1E1 && abs(z2-z0)>=dkm/1E1; break; end
        end
        %Muller's Method
        [zn, szSuFa] = muller(f, MPM, z0, z1, z2, err);
        if((toc-tm0)>10), szSuFa = 'Fail'; return; end
        if strcmp(szSuFa, 'Fail'), continue; end
        cpn = w/real(zn); kin = imag(zn);
        if ~isnan(cpn) && ~isnan(kin) && cpn>=cpmin/scale && cpn<=cpmax*scale...
                && kin>=kimin/scale && kin<=kimax*scale
            cp = cpn; ki = kin; return;
        end
    end
end
szSuFa = 'Fail';
end
%----------------------------------------------------------------------------------------------
function [zn, szSuFa] = muller(f, MPM, z0, z1, z2, err)
%Muller's Method
%%
if abs(z0-z1)<err || abs(z1-z2)<err || abs(z2-z0)<err || isnan(z0) || isnan(z1) || isnan(z2)
    szSuFa = 'Fail'; return;
end
n = 1; N = 1E3; delta = inf; z3 = NaN; szSuFa = 'Success';
f0 = DispFunction(f,MPM,z0); f1 = DispFunction(f,MPM,z1); f2 = DispFunction(f,MPM,z2);
q = (z2-z1)/(z1-z0+eps);
while n<=N
    p = q+1;
    a = q^2*f0-p*q*f1+q*f2; b = q^2*f0-p^2*f1+(p+q)*f2; c = p*f2;
    temp1 = (b+sqrt(b^2-4*a*c)); temp2 = (b-sqrt(b^2-4*a*c));
    if abs(temp1)>=abs(temp2)
        temp = temp1;
    else
        temp = temp2;
    end
    q3 = -2*c/(temp+eps);
    z3 = z2+q3*(z2-z1);
    f3 = DispFunction(f,MPM,z3);
    if abs(z3)<1
        delta = abs(z3-z2);
    elseif abs(z3)>=1
        delta = abs((z3-z2)/z3);
    end
    if delta<=err && n>=1E1
        zn= z3; return;
    else
        z0 = z1; z1 = z2; z2 = z3; f0 = f1; f1 = f2; f2 = f3; q = q3;
        n = n+1; continue;
    end
end
zn = z3; szSuFa = 'Fail';
end