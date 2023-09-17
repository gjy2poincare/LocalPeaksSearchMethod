function fdisp = DispFunction(f, MPM, k)
%频散函数计算【Disperion Function Computing】
global szFun
switch szFun
    case 'PlanarPlate'
        fdisp = DispFunction_PP(f, MPM, k); 
    case 'CylindricalShell'
        fdisp = DispFunction_CS(f, MPM, k);
end
end
%#############################################################################################
function fdisp = DispFunction_PP(f, MPM, k)
%频散方程
row1 = MPM(1); c1 = MPM(2); rowvm = MPM(3); lamdavm = MPM(4); miuvm = MPM(5); dvm = MPM(6);
rowem = MPM(7); lamdaem = MPM(8); miuem = MPM(9); dem = MPM(10); row2 = MPM(11); c2 = MPM(12);
a_inner = MPM(13); b_middle = MPM(14); c_outer = MPM(15);
global szBC szMode
%
w = 2*pi*f; d = dvm; h = dem/2;
clvm = sqrt((lamdavm+2*miuvm)/rowvm); ctvm = sqrt(miuvm/rowvm);
clem = sqrt((lamdaem+2*miuem)/rowem); ctem = sqrt(miuem/rowem);
%%
switch szBC
    case 'Va-So-Va'
        %自由薄板【Free Plate】
        p = sqrt((w/clem)^2-k^2); q = sqrt((w/ctem)^2-k^2);
        if strcmp(szMode, 'S/L')
            %对称模态频散函数【Symmetrical Mode】
            fdisp = (q^2-k^2)^2*sin(q*h)*cos(p*h)+4*k^2*p*q*sin(p*h)*cos(q*h);
            % fdisp = (q^2-k^2)^2*sin(q*h)*cos(p*h)/(4*k^2*p*q)+sin(p*h)*cos(q*h);
        elseif strcmp(szMode, 'A/T')
            %反对称模态频散函数【Asymmetrical Mode】
            fdisp = (q^2-k^2)^2*sin(p*h)*cos(q*h)+4*k^2*p*q*sin(q*h)*cos(p*h);
            % fdisp = (q^2-k^2)^2*sin(p*h)*cos(q*h)/(4*k^2*p*q)+sin(q*h)*cos(p*h);
        end
    case 'Va-So-Fl'
        %一面流体一面自由薄板【One-side Fluid One-side Free Plate】
        p = sqrt((w/clem)^2-k^2); q = sqrt((w/ctem)^2-k^2); yta1 = sqrt((w/c1)^2-k^2);
        %
        M(1,1) = 2*1i*k*p*sin(p*h);
        M(1,4) = (q^2-k^2)*sin(q*h);
        %
        M(2,2) = 2*1i*k*p*cos(p*h);
        M(2,3) = (q^2-k^2)*cos(q*h);
        %
        M(3,1) = (q^2-k^2)*cos(p*h);
        M(3,4) = 2*1i*k*q*cos(q*h);
        M(3,5) = row1*w^2;
        %
        M(4,2) = (q^2-k^2)*sin(p*h);
        M(4,3) = 2*1i*k*q*sin(q*h);
        M(4,5) = -row1*w^2;
        %
        M(5,1) = p*sin(p*h);
        M(5,2) = p*cos(p*h);
        M(5,3) = 1i*k*cos(q*h);
        M(5,4) = 1i*k*sin(q*h);
        M(5,5) = -2*1i*miuem*yta1;
        %
        fdisp = det(M);
    case 'Fl-So-Fl'
        %流体中薄板【Plate in Fluid】
        p = sqrt((w/clem)^2-k^2); q = sqrt((w/ctem)^2-k^2); yta = sqrt((w/c1)^2-k^2);
        if strcmp(szMode, 'S/L')
            %对称模态频散方程【Symmetrical Modes】
            fdisp = (q^2-k^2)^2*sin(q*h)*cos(p*h)+4*k^2*p*q*sin(p*h)*cos(q*h)...
                +row1*w^2*p*((q^2+k^2)*sin(p*h)*sin(q*h))/(1i*miuem*yta);
%             fdisp = (q^2-k^2)^2*sin(q*h)*cos(p*h)/(4*k^2*p*q)+sin(p*h)*cos(q*h)...
%                 +row1*w^2*p*((q^2+k^2)*sin(p*h)*sin(q*h))/((4*k^2*p*q)*(1i*miuem*yta));
        elseif strcmp(szMode, 'A/T')
            %反对称模态频散方程【Asymmetric Modes】
            fdisp = (q^2-k^2)^2*sin(p*h)*cos(q*h)+4*k^2*p*q*sin(q*h)*cos(p*h)...
                -row1*w^2*p*((q^2+k^2)*cos(p*h)*cos(q*h))/(1i*miuem*yta);
%             fdisp = (q^2-k^2)^2*sin(p*h)*cos(q*h)/(4*k^2*p*q)+sin(q*h)*cos(p*h)...
%                 -row1*w^2*p*((q^2+k^2)*cos(p*h)*cos(q*h))/((4*k^2*p*q)*(1i*miuem*yta));
        end
    case 'Va-So-So-Va'
        %自由覆盖粘弹性层薄板【Free Coated Plate】
        alphaem = sqrt((w/clem)^2-k^2); betaem = sqrt((w/ctem)^2-k^2);
        alphavm = sqrt((w/clvm)^2-k^2); betavm = sqrt((w/ctvm)^2-k^2);
        
        M(1,5) = -2*1i*k*alphavm*sin(alphavm*d);
        M(1,6) = 2*1i*k*alphavm*cos(alphavm*d);
        M(1,7) = (betavm^2-k^2)*cos(betavm*d);
        M(1,8) = (betavm^2-k^2)*sin(betavm*d);
        %
        M(2,5) = -(betavm^2-k^2)*cos(alphavm*d);
        M(2,6) = -(betavm^2-k^2)*sin(alphavm*d);
        M(2,7) = -2*1i*k*betavm*sin(betavm*d);
        M(2,8) = 2*1i*k*betavm*cos(betavm*d);
        %
        M(3,2) = 4*1i*k*alphaem*cos(alphaem*h);
        M(3,3) = 2*(betaem^2-k^2)*cos(betaem*h);
        M(3,6) = -2*1i*k*alphavm*miuvm/miuem;
        M(3,7) = -(betavm^2-k^2)*miuvm/miuem;
        %
        M(4,1) = -2*(betaem^2-k^2)*cos(alphaem*h);
        M(4,4) = 4*1i*k*betaem*cos(betaem*h);
        M(4,5) = (betavm^2-k^2)*miuvm/miuem;
        M(4,8) = -2*1i*k*betavm*miuvm/miuem;
        %
        M(5,1) = 1i*k*cos(alphaem*h);
        M(5,2) = 1i*k*sin(alphaem*h);
        M(5,3) = betaem*sin(betaem*h);
        M(5,4) = -betaem*cos(betaem*h);
        M(5,5) = -1i*k;
        M(5,8) = betavm;
        %
        M(6,1) = -alphaem*sin(alphaem*h);
        M(6,2) = alphaem*cos(alphaem*h);
        M(6,3) = 1i*k*cos(betaem*h);
        M(6,4) = 1i*k*sin(betaem*h);
        M(6,6) = -alphavm;
        M(6,7) = -1i*k;
        %
        M(7,1) = 2*1i*k*alphaem*sin(alphaem*h);
        M(7,2) = 2*1i*k*alphaem*cos(alphaem*h);
        M(7,3) = (betaem^2-k^2)*cos(betaem*h);
        M(7,4) = -(betaem^2-k^2)*sin(betaem*h);
        %
        M(8,1) = -(betaem^2-k^2)*cos(alphaem*h);
        M(8,2) = (betaem^2-k^2)*sin(alphaem*h);
        M(8,3) = 2*1i*k*betaem*sin(betaem*h);
        M(8,4) = 2*1i*k*betaem*cos(betaem*h);
        %
        fdisp = det(M);
    case 'Va-So-So-Fl'
        %一面流体一面自由覆盖粘弹性层薄板【One-side Fluid One-side Free Coated Plate】
        alphaem = sqrt((w/clem)^2-k^2); betaem = sqrt((w/ctem)^2-k^2); 
        yta1 = sqrt((w/c1)^2-k^2);
        alphavm = sqrt((w/clvm)^2-k^2); betavm = sqrt((w/ctvm)^2-k^2);
        
        M(1,5) = -2*1i*k*alphavm*sin(alphavm*d);
        M(1,6) = 2*1i*k*alphavm*cos(alphavm*d);
        M(1,7) = (betavm^2-k^2)*cos(betavm*d);
        M(1,8) = (betavm^2-k^2)*sin(betavm*d);
        %
        M(2,5) = -(betavm^2-k^2)*cos(alphavm*d);
        M(2,6) = -(betavm^2-k^2)*sin(alphavm*d);
        M(2,7) = -2*1i*k*betavm*sin(betavm*d);
        M(2,8) = 2*1i*k*betavm*cos(betavm*d);
        M(2,9) = row1*w^2/miuvm;
        %
        M(3,5) = -alphavm*sin(alphavm*d);
        M(3,6) = alphavm*cos(alphavm*d);
        M(3,7) = 1i*k*cos(betavm*d);
        M(3,8) = 1i*k*sin(betavm*d);
        M(3,9) = -1i*yta1;
        %
        M(4,2) = 4*1i*k*alphaem*cos(alphaem*h);
        M(4,3) = 2*(betaem^2-k^2)*cos(betaem*h);
        M(4,6) = -2*1i*k*alphavm*miuvm/miuem;
        M(4,7) = -(betavm^2-k^2)*miuvm/miuem;
        %
        M(5,1) = -2*(betaem^2-k^2)*cos(alphaem*h);
        M(5,4) = 4*1i*k*betaem*cos(betaem*h);
        M(5,5) = (betavm^2-k^2)*miuvm/miuem;
        M(5,8) = -2*1i*k*betavm*miuvm/miuem;
        %
        M(6,1) = 1i*k*cos(alphaem*h);
        M(6,2) = 1i*k*sin(alphaem*h);
        M(6,3) = betaem*sin(betaem*h);
        M(6,4) = -betaem*cos(betaem*h);
        M(6,5) = -1i*k;
        M(6,8) = betavm;
        %
        M(7,1) = -alphaem*sin(alphaem*h);
        M(7,2) = alphaem*cos(alphaem*h);
        M(7,3) = 1i*k*cos(betaem*h);
        M(7,4) = 1i*k*sin(betaem*h);
        M(7,6) = -alphavm;
        M(7,7) = -1i*k;
        %
        M(8,1) = 2*1i*k*alphaem*sin(alphaem*h);
        M(8,2) = 2*1i*k*alphaem*cos(alphaem*h);
        M(8,3) = (betaem^2-k^2)*cos(betaem*h);
        M(8,4) = -(betaem^2-k^2)*sin(betaem*h);
        %
        M(9,1) = -(betaem^2-k^2)*cos(alphaem*h);
        M(9,2) = (betaem^2-k^2)*sin(alphaem*h);
        M(9,3) = 2*1i*k*betaem*sin(betaem*h);
        M(9,4) = 2*1i*k*betaem*cos(betaem*h);
        %
        fdisp = det(M);
    case 'Fl-So-So-Fl'
        %流体中覆盖粘弹性层薄板【Coated Plate in Fluid】
        alphaem = sqrt((w/clem)^2-k^2); betaem = sqrt((w/ctem)^2-k^2);
        yta1 = sqrt((w/c1)^2-k^2);
        alphavm = sqrt((w/clvm)^2-k^2); betavm = sqrt((w/ctvm)^2-k^2); 
        yta2 = sqrt((w/c2)^2-k^2);
        
        M = zeros(10);
        
        M(1,5) = -2*1i*k*alphavm*sin(alphavm*d);
        M(1,6) = 2*1i*k*alphavm*cos(alphavm*d);
        M(1,7) = (betavm^2-k^2)*cos(betavm*d);
        M(1,8) = (betavm^2-k^2)*sin(betavm*d);
        %
        M(2,5) = -(betavm^2-k^2)*cos(alphavm*d);
        M(2,6) = -(betavm^2-k^2)*sin(alphavm*d);
        M(2,7) = -2*1i*k*betavm*sin(betavm*d);
        M(2,8) = 2*1i*k*betavm*cos(betavm*d);
        M(2,9) = row1*w^2/miuvm;
        %
        M(3,5) = -alphavm*sin(alphavm*d);
        M(3,6) = alphavm*cos(alphavm*d);
        M(3,7) = 1i*k*cos(betavm*d);
        M(3,8) = 1i*k*sin(betavm*d);
        M(3,9) = -1i*yta1;
        %
        M(4,2) = 4*1i*k*alphaem*cos(alphaem*h);
        M(4,3) = 2*(betaem^2-k^2)*cos(betaem*h);
        M(4,6) = -2*1i*k*alphavm*miuvm/miuem;
        M(4,7) = -(betavm^2-k^2)*miuvm/miuem;
        %
        M(5,1) = -2*(betaem^2-k^2)*cos(alphaem*h);
        M(5,4) = 4*1i*k*betaem*cos(betaem*h);
        M(5,5) = (betavm^2-k^2)*miuvm/miuem;
        M(5,8) = -2*1i*k*betavm*miuvm/miuem;
        M(5,10) = row2*w^2/miuem;
        %
        M(6,1) = 1i*k*cos(alphaem*h);
        M(6,2) = 1i*k*sin(alphaem*h);
        M(6,3) = betaem*sin(betaem*h);
        M(6,4) = -betaem*cos(betaem*h);
        M(6,5) = -1i*k;
        M(6,8) = betavm;
        %
        M(7,2) = 2*alphaem*cos(alphaem*h);
        M(7,3) = 2*1i*k*cos(betaem*h);
        M(7,6) = -alphavm;
        M(7,7) = -1i*k;
        M(7,10) = 1i*yta2;
        %
        M(8,1) = 2*1i*k*alphaem*sin(alphaem*h);
        M(8,2) = 2*1i*k*alphaem*cos(alphaem*h);
        M(8,3) = (betaem^2-k^2)*cos(betaem*h);
        M(8,4) = -(betaem^2-k^2)*sin(betaem*h);
        %
        M(9,1) = -(betaem^2-k^2)*cos(alphaem*h);
        M(9,2) = (betaem^2-k^2)*sin(alphaem*h);
        M(9,3) = 2*1i*k*betaem*sin(betaem*h);
        M(9,4) = 2*1i*k*betaem*cos(betaem*h);
        M(9,10) = row2*w^2/miuem;
        %
        M(10,1) = alphaem*sin(alphaem*h);
        M(10,2) = alphaem*cos(alphaem*h);
        M(10,3) = 1i*k*cos(betaem*h);
        M(10,4) = -1i*k*sin(betaem*h);
        M(10,10) = 1i*yta2;
        %
        fdisp = det(M);
end
%%
if isnan(abs(fdisp)) || isinf(abs(fdisp));
    fdisp = NaN;
    warning('Function Computing Failed');
end
end
%#############################################################################################
function fdisp = DispFunction_CS(f, MPM, k)
row1 = MPM(1); c1 = MPM(2); rowvm = MPM(3); lamdavm = MPM(4); miuvm = MPM(5); dvm = MPM(6);
rowem = MPM(7); lamdaem = MPM(8); miuem = MPM(9); dem = MPM(10); row2 = MPM(11); c2 = MPM(12);
a_inner = MPM(13); b_middle = MPM(14); c_outer = MPM(15);
%
w = 2*pi*f; a = a_inner; b = b_middle; c = c_outer;
clvm = sqrt((lamdavm+2*miuvm)/rowvm); ctvm = sqrt(miuvm/rowvm);
clem = sqrt((lamdaem+2*miuem)/rowem); ctem = sqrt(miuem/rowem);
%%
switch szBC
    case 'Va-So-Va'
        %自由柱壳中的类Lamb波【Lamb-like Waves in Free Cylindrical Shells】
        alpha = sqrt(w^2/clem^2-k^2); beta = sqrt(w^2/ctem^2-k^2);
        if strcmp(szMode,'S/L')
            %自由柱壳：n = 0,纵向模态【Free: n=0, Longitudinal Mode】
            n=0;
            bj_n_alfa = besselj(n,alpha*a); bj_np1_alfa = besselj(n+1,alpha*a);
            bj_n_bta = besselj(n,beta*a); bj_np1_bta = besselj(n+1,beta*a);
            by_n_alfa = bessely(n,alpha*a); by_np1_alfa = bessely(n+1,alpha*a);
            by_n_bta = bessely(n,beta*a); by_np1_bta = bessely(n+1,beta*a);
            bj_n_alfb = besselj(n,alpha*b); bj_np1_alfb = besselj(n+1,alpha*b);
            bj_n_btb = besselj(n,beta*b); bj_np1_btb = besselj(n+1,beta*b);
            by_n_alfb = bessely(n,alpha*b); by_np1_alfb = bessely(n+1,alpha*b);
            by_n_btb = bessely(n,beta*b); by_np1_btb = bessely(n+1,beta*b);
            
            d11 = -(beta^2-k^2)*a^2*bj_n_alfa+2*alpha*a*bj_np1_alfa ;
            d13 = 2*1i*k*(-beta^2*a^2*bj_n_bta+beta*a*bj_np1_bta);
            d14 = -(beta^2-k^2)*a^2*by_n_alfa+2*alpha*a*by_np1_alfa;
            d16 = 2*1i*k*(-beta^2*a^2*by_n_bta+beta*a*by_np1_bta);
            
            d31 = -2*1i*k*alpha*a*bj_np1_alfa;
            d33 = -(beta^2-k^2)*beta*a*bj_np1_bta;
            d34 = -2*1i*k*alpha*a*by_np1_alfa;
            d36 = -(beta^2-k^2)*beta*a*by_np1_bta;
            
            d41 = -(beta^2-k^2)*b^2*bj_n_alfb+2*alpha*b*bj_np1_alfb;
            d43 = 2*1i*k*(-beta^2*b^2*bj_n_btb+beta*b*bj_np1_btb);
            d44 = -(beta^2-k^2)*b^2*by_n_alfb+2*alpha*b*by_np1_alfb;
            d46 = 2*1i*k*(-beta^2*b^2*by_n_btb+beta*b*by_np1_btb);
            
            d61 = -2*1i*k*alpha*b*bj_np1_alfb;
            d63 = -(beta^2-k^2)*beta*b*bj_np1_btb;
            d64 =  -2*1i*k*alpha*b*by_np1_alfb;
            d66 = -(beta^2-k^2)*beta*b*by_np1_btb;
            D1(1,1) = d11; D1(1,2) = d13; D1(1,3) = d14; D1(1,4) = d16;
            D1(2,1) = d31; D1(2,2) = d33; D1(2,3) = d34; D1(2,4) = d36;
            D1(3,1) = d41; D1(3,2) = d43; D1(3,3) = d44; D1(3,4) = d46;
            D1(4,1) = d61; D1(4,2) = d63; D1(4,3) = d64; D1(4,4) = d66;
            fdisp = det(D1);
        elseif strcmp(szMode, 'A/T')
            %自由柱壳：n = 0,扭转模态【Free: n=0, Tortional Mode】
            d22 = -beta^2*a^2*besselj(0,beta*a)+2*beta*a*besselj(1,beta*a);
            d25 = -beta^2*a^2*bessely(0,beta*a)+2*beta*a*bessely(1,beta*a);
            d52 = -beta^2*b^2*besselj(0,beta*b)+2*beta*b*besselj(1,beta*b);
            d55 = -beta^2*b^2*bessely(0,beta*b)+2*beta*b*bessely(1,beta*b);
            D2(1,1) = d22; D2(1,2) = d25; D2(2,1) = d52; D2(2,2) = d55;
            fdisp = det(D2);
        elseif strcmp(szMode,'F')
            %自由柱壳：n = 0，1，2……【Free Shell: n=0, 1, 2, ...】
            %n = 0时，纵向和扭转模态之和；【Free: n=0, Longitudinal and Tortional Mode】
            %n = 1，2，……，弯曲模态【n=1, 2, 3, ... ,Flexural Mode】
            n =  nMode_cs;
            D = zeros(6);
            bj_n_alfa = besselj(n,alpha*a); bj_np1_alfa = besselj(n+1,alpha*a);
            bj_n_bta = besselj(n,beta*a); bj_np1_bta = besselj(n+1,beta*a);
            by_n_alfa = bessely(n,alpha*a); by_np1_alfa = bessely(n+1,alpha*a);
            by_n_bta = bessely(n,beta*a); by_np1_bta = bessely(n+1,beta*a);
            bj_n_alfb = besselj(n,alpha*b); bj_np1_alfb = besselj(n+1,alpha*b);
            bj_n_btb = besselj(n,beta*b); bj_np1_btb = besselj(n+1,beta*b);
            by_n_alfb = bessely(n,alpha*b); by_np1_alfb = bessely(n+1,alpha*b);
            by_n_btb = bessely(n,beta*b); by_np1_btb = bessely(n+1,beta*b);
            
            D(1,1) = (2*n*(n-1)-(beta^2-k^2)*a^2)*bj_n_alfa +2*alpha*a*bj_np1_alfa;
            D(1,2) = 2*n*(n-1)*bj_n_bta-2*n*beta*a*bj_np1_bta;
            D(1,3) = 2*1i*k*((n*(n-1)-beta^2*a^2)*bj_n_bta+beta*a*bj_np1_bta);
            D(1,4) = (2*n*(n-1)-(beta^2-k^2)*a^2)*by_n_alfa+2*alpha*a*by_np1_alfa;
            D(1,5) = 2*n*(n-1)*by_n_bta-2*n*beta*a*by_np1_bta;
            D(1,6) = 2*1i*k*((n*(n-1)-beta^2*a^2)*by_n_bta+beta*a*by_np1_bta);
            
            D(2,1) = -2*n*(n-1)*bj_n_alfa+2*n*alpha*a*bj_np1_alfa;
            D(2,2) = (beta^2*a^2-2*n*(n-1))*bj_n_bta-2*beta*a*bj_np1_bta;
            D(2,3) = 2*1i*k*(-n*(n-1)*bj_n_bta+n*beta*a*bj_np1_bta);
            D(2,4) = -2*n*(n-1)*by_n_alfa+2*n*alpha*a*by_np1_alfa;
            D(2,5) = (beta^2*a^2-2*n*(n-1))*by_n_bta-2*beta*a*by_np1_bta;
            D(2,6) = 2*1i*k*(-n*(n-1)*by_n_bta+n*beta*a*by_np1_bta);
            
            D(3,1) = 2*1i*k*(n*bj_n_alfa-alpha*a*bj_np1_alfa);
            D(3,2) = 1i*k*n*bj_n_bta;
            D(3,3) = (beta^2-k^2)*(n*bj_n_bta-beta*a*bj_np1_bta);
            D(3,4) = 2*1i*k*(n*by_n_alfa-alpha*a*by_np1_alfa);
            D(3,5) = 1i*k*n*by_n_bta;
            D(3,6) = (beta^2-k^2)*(n*by_n_bta-beta*a*by_np1_bta);
            
            D(4,1) = (2*n*(n-1)-(beta^2-k^2)*b^2)*bj_n_alfb+2*alpha*b*bj_np1_alfb;
            D(4,2) = 2*n*(n-1)*bj_n_btb-2*n*beta*b*bj_np1_btb;
            D(4,3) = 2*1i*k*((n*(n-1)-beta^2*b^2)*bj_n_btb+beta*b*bj_np1_btb);
            D(4,4) = (2*n*(n-1)-(beta^2-k^2)*b^2)*by_n_alfb+2*alpha*b*by_np1_alfb;
            D(4,5) = 2*n*(n-1)*by_n_btb-2*n*beta*b*by_np1_btb;
            D(4,6) = 2*1i*k*((n*(n-1)-beta^2*b^2)*by_n_btb+beta*b*by_np1_btb);
            
            D(5,1) = -2*n*(n-1)*bj_n_alfb+2*n*alpha*b*bj_np1_alfb;
            D(5,2) = (beta^2*b^2-2*n*(n-1))*bj_n_btb-2*beta*b*bj_np1_btb;
            D(5,3) = 2*1i*k*(-n*(n-1)*bj_n_btb+n*beta*b*bj_np1_btb);
            D(5,4) = -2*n*(n-1)*by_n_alfb+2*n*alpha*b*by_np1_alfb;
            D(5,5) = (beta^2*b^2-2*n*(n-1))*by_n_btb-2*beta*b*by_np1_btb;
            D(5,6) = 2*1i*k*(-n*(n-1)*by_n_btb+n*beta*b*by_np1_btb);
            
            D(6,1) = 2*1i*k*(n*bj_n_alfb-alpha*b*bj_np1_alfb);
            D(6,2) = 1i*k*n*bj_n_btb;
            D(6,3) = (beta^2-k^2)*(n*bj_n_btb-beta*b*bj_np1_btb);
            D(6,4) = 2*1i*k*(n*by_n_alfb-alpha*b*by_np1_alfb);
            D(6,5) = 1i*k*n*by_n_btb;
            D(6,6) = (beta^2-k^2)*(n*by_n_btb-beta*b*by_np1_btb);
            
            fdisp = det(D);
        end
    case 'Va-So-Fl'
        %%
        %外部充液内部自由圆柱壳【Outside Fluid Inside Free Cylindrical Shell】
        alpha = sqrt(w^2/clem^2-k^2); beta = sqrt(w^2/ctem^2-k^2);
        gamma2 = sqrt(w^2/c2^2-k^2);
        n =  nMode_cs;
        
        bj_n_alfa = besselj(n,alpha*a); bj_np1_alfa = besselj(n+1,alpha*a);
        bj_n_bta = besselj(n,beta*a); bj_np1_bta = besselj(n+1,beta*a);
        by_n_alfa = bessely(n,alpha*a); by_np1_alfa = bessely(n+1,alpha*a);
        by_n_bta = bessely(n,beta*a); by_np1_bta = bessely(n+1,beta*a);
        bj_n_alfb = besselj(n,alpha*b); bj_np1_alfb = besselj(n+1,alpha*b);
        bj_n_btb = besselj(n,beta*b); bj_np1_btb = besselj(n+1,beta*b);
        by_n_alfb = bessely(n,alpha*b); by_np1_alfb = bessely(n+1,alpha*b);
        by_n_btb = bessely(n,beta*b); by_np1_btb = bessely(n+1,beta*b);
        
        D = zeros(7);
        D(1,1) = (2*n*(n-1)-(beta^2-k^2)*a^2)*bj_n_alfa +2*alpha*a*bj_np1_alfa;
        D(1,2) = 2*n*(n-1)*bj_n_bta-2*n*beta*a*bj_np1_bta;
        D(1,3) = 2*1i*k*((n*(n-1)-beta^2*a^2)*bj_n_bta+beta*a*bj_np1_bta);
        D(1,4) = (2*n*(n-1)-(beta^2-k^2)*a^2)*by_n_alfa+2*alpha*a*by_np1_alfa;
        D(1,5) = 2*n*(n-1)*by_n_bta-2*n*beta*a*by_np1_bta;
        D(1,6) = 2*1i*k*((n*(n-1)-beta^2*a^2)*by_n_bta+beta*a*by_np1_bta);
        D(1,7) = 0;
        
        D(2,1) = -2*n*(n-1)*bj_n_alfa+2*n*alpha*a*bj_np1_alfa;
        D(2,2) = (beta^2*a^2-2*n*(n-1))*bj_n_bta-2*beta*a*bj_np1_bta;
        D(2,3) = 2*1i*k*(-n*(n-1)*bj_n_bta+n*beta*a*bj_np1_bta);
        D(2,4) = -2*n*(n-1)*by_n_alfa+2*n*alpha*a*by_np1_alfa;
        D(2,5) = (beta^2*a^2-2*n*(n-1))*by_n_bta-2*beta*a*by_np1_bta;
        D(2,6) = 2*1i*k*(-n*(n-1)*by_n_bta+n*beta*a*by_np1_bta);
        D(2,7) = 0;
        
        D(3,1) = 2*1i*k*(n*bj_n_alfa-alpha*a*bj_np1_alfa);
        D(3,2) = 1i*k*n*bj_n_bta;
        D(3,3) = (beta^2-k^2)*(n*bj_n_bta-beta*a*bj_np1_bta);
        D(3,4) = 2*1i*k*(n*by_n_alfa-alpha*a*by_np1_alfa);
        D(3,5) = 1i*k*n*by_n_bta;
        D(3,6) = (beta^2-k^2)*(n*by_n_bta-beta*a*by_np1_bta);
        D(3,7) = 0;
        
        D(4,1) = (2*n*(n-1)-(beta^2-k^2)*b^2)*bj_n_alfb+2*alpha*b*bj_np1_alfb;
        D(4,2) = 2*n*(n-1)*bj_n_btb-2*n*beta*b*bj_np1_btb;
        D(4,3) = 2*1i*k*((n*(n-1)-beta^2*b^2)*bj_n_btb+beta*b*bj_np1_btb);
        D(4,4) = (2*n*(n-1)-(beta^2-k^2)*b^2)*by_n_alfb+2*alpha*b*by_np1_alfb;
        D(4,5) = 2*n*(n-1)*by_n_btb-2*n*beta*b*by_np1_btb;
        D(4,6) = 2*1i*k*((n*(n-1)-beta^2*b^2)*by_n_btb+beta*b*by_np1_btb);
        D(4,7) = b^2*row2*w^2*besselh(n,1,gamma2*b)/miuem;
        
        D(5,1) = -2*n*(n-1)*bj_n_alfb+2*n*alpha*b*bj_np1_alfb;
        D(5,2) = (beta^2*b^2-2*n*(n-1))*bj_n_btb-2*beta*b*bj_np1_btb;
        D(5,3) = 2*1i*k*(-n*(n-1)*bj_n_btb+n*beta*b*bj_np1_btb);
        D(5,4) = -2*n*(n-1)*by_n_alfb+2*n*alpha*b*by_np1_alfb;
        D(5,5) = (beta^2*b^2-2*n*(n-1))*by_n_btb-2*beta*b*by_np1_btb;
        D(5,6) = 2*1i*k*(-n*(n-1)*by_n_btb+n*beta*b*by_np1_btb);
        D(5,7) = 0;
        
        D(6,1) = 2*1i*k*(n*bj_n_alfb-alpha*b*bj_np1_alfb);
        D(6,2) = 1i*k*n*bj_n_btb;
        D(6,3) = (beta^2-k^2)*(n*bj_n_btb-beta*b*bj_np1_btb);
        D(6,4) = 2*1i*k*(n*by_n_alfb-alpha*b*by_np1_alfb);
        D(6,5) = 1i*k*n*by_n_btb;
        D(6,6) = (beta^2-k^2)*(n*by_n_btb-beta*b*by_np1_btb);
        D(6,7) = 0;
        
        D(7,1) = n*bj_n_alfb-alpha*b*bj_np1_alfb;
        D(7,2) = n*bj_n_btb;
        D(7,3) = 1i*k*(n*bj_n_btb-beta*b*bj_np1_btb);
        D(7,4) = n*by_n_alfb-alpha*b*by_np1_alfb;
        D(7,5) = n*by_n_btb;
        D(7,6) = 1i*k*(n*by_n_btb-beta*b*by_np1_btb);
        D(7,7) = -(n*besselh(n,1,gamma2*b)-gamma2*b*besselh(n+1,1,gamma2*b));
        
        fdisp = det(D);
    case 'Fl-So-Fl'
        %外部充液内部充液圆柱壳【Outside and Inside Fluid Shell】
        alpha = sqrt(w^2/clem^2-k^2); beta = sqrt(w^2/ctem^2-k^2);
        gamma1 = sqrt(w^2/c1^2-k^2); gamma2 = sqrt(w^2/c2^2-k^2);
        n =  nMode_cs;
        
        bj_n_alfa = besselj(n,alpha*a); bj_np1_alfa = besselj(n+1,alpha*a);
        bj_n_bta = besselj(n,beta*a); bj_np1_bta = besselj(n+1,beta*a);
        by_n_alfa = bessely(n,alpha*a); by_np1_alfa = bessely(n+1,alpha*a);
        by_n_bta = bessely(n,beta*a); by_np1_bta = bessely(n+1,beta*a);
        bj_n_alfb = besselj(n,alpha*b); bj_np1_alfb = besselj(n+1,alpha*b);
        bj_n_btb = besselj(n,beta*b); bj_np1_btb = besselj(n+1,beta*b);
        by_n_alfb = bessely(n,alpha*b); by_np1_alfb = bessely(n+1,alpha*b);
        by_n_btb = bessely(n,beta*b); by_np1_btb = bessely(n+1,beta*b);
        
        D = zeros(8);
        
        D(1,1) = (2*n*(n-1)-(beta^2-k^2)*a^2)*bj_n_alfa +2*alpha*a*bj_np1_alfa;
        D(1,2) = 2*n*(n-1)*bj_n_bta-2*n*beta*a*bj_np1_bta;
        D(1,3) = 2*1i*k*((n*(n-1)-beta^2*a^2)*bj_n_bta+beta*a*bj_np1_bta);
        D(1,4) = (2*n*(n-1)-(beta^2-k^2)*a^2)*by_n_alfa+2*alpha*a*by_np1_alfa;
        D(1,5) = 2*n*(n-1)*by_n_bta-2*n*beta*a*by_np1_bta;
        D(1,6) = 2*1i*k*((n*(n-1)-beta^2*a^2)*by_n_bta+beta*a*by_np1_bta);
        D(1,7) = a^2*row1*w^2*besselj(n,gamma1*a)/miuem;
        D(1,8) = 0;
        
        D(2,1) = -2*n*(n-1)*bj_n_alfa+2*n*alpha*a*bj_np1_alfa;
        D(2,2) = (beta^2*a^2-2*n*(n-1))*bj_n_bta-2*beta*a*bj_np1_bta;
        D(2,3) = 2*1i*k*(-n*(n-1)*bj_n_bta+n*beta*a*bj_np1_bta);
        D(2,4) = -2*n*(n-1)*by_n_alfa+2*n*alpha*a*by_np1_alfa;
        D(2,5) = (beta^2*a^2-2*n*(n-1))*by_n_bta-2*beta*a*by_np1_bta;
        D(2,6) = 2*1i*k*(-n*(n-1)*by_n_bta+n*beta*a*by_np1_bta);
        D(2,7) = 0;
        D(2,8) = 0;
        
        D(3,1) = 2*1i*k*(n*bj_n_alfa-alpha*a*bj_np1_alfa);
        D(3,2) = 1i*k*n*bj_n_bta;
        D(3,3) = (beta^2-k^2)*(n*bj_n_bta-beta*a*bj_np1_bta);
        D(3,4) = 2*1i*k*(n*by_n_alfa-alpha*a*by_np1_alfa);
        D(3,5) = 1i*k*n*by_n_bta;
        D(3,6) = (beta^2-k^2)*(n*by_n_bta-beta*a*by_np1_bta);
        D(3,7) = 0;
        D(3,8) = 0;
        
        D(4,1) = (2*n*(n-1)-(beta^2-k^2)*b^2)*bj_n_alfb+2*alpha*b*bj_np1_alfb;
        D(4,2) = 2*n*(n-1)*bj_n_btb-2*n*beta*b*bj_np1_btb;
        D(4,3) = 2*1i*k*((n*(n-1)-beta^2*b^2)*bj_n_btb+beta*b*bj_np1_btb);
        D(4,4) = (2*n*(n-1)-(beta^2-k^2)*b^2)*by_n_alfb+2*alpha*b*by_np1_alfb;
        D(4,5) = 2*n*(n-1)*by_n_btb-2*n*beta*b*by_np1_btb;
        D(4,6) = 2*1i*k*((n*(n-1)-beta^2*b^2)*by_n_btb+beta*b*by_np1_btb);
        D(4,7) = 0;
        D(4,8) = b^2*row2*w^2*besselh(n,1,gamma2*b)/miuem;
        
        D(5,1) = -2*n*(n-1)*bj_n_alfb+2*n*alpha*b*bj_np1_alfb;
        D(5,2) = (beta^2*b^2-2*n*(n-1))*bj_n_btb-2*beta*b*bj_np1_btb;
        D(5,3) = 2*1i*k*(-n*(n-1)*bj_n_btb+n*beta*b*bj_np1_btb);
        D(5,4) = -2*n*(n-1)*by_n_alfb+2*n*alpha*b*by_np1_alfb;
        D(5,5) = (beta^2*b^2-2*n*(n-1))*by_n_btb-2*beta*b*by_np1_btb;
        D(5,6) = 2*1i*k*(-n*(n-1)*by_n_btb+n*beta*b*by_np1_btb);
        D(5,7) = 0;
        D(5,8) = 0;
        
        D(6,1) = 2*1i*k*(n*bj_n_alfb-alpha*b*bj_np1_alfb);
        D(6,2) = 1i*k*n*bj_n_btb;
        D(6,3) = (beta^2-k^2)*(n*bj_n_btb-beta*b*bj_np1_btb);
        D(6,4) = 2*1i*k*(n*by_n_alfb-alpha*b*by_np1_alfb);
        D(6,5) = 1i*k*n*by_n_btb;
        D(6,6) = (beta^2-k^2)*(n*by_n_btb-beta*b*by_np1_btb);
        D(6,7) = 0;
        D(6,8) = 0;
        
        D(7,1) = n*bj_n_alfa-alpha*a*bj_np1_alfa;
        D(7,2) = n*bj_n_bta;
        D(7,3) = 1i*k*(n*bj_n_bta-beta*a*bj_np1_bta);
        D(7,4) = n*by_n_alfa-alpha*a*by_np1_alfa;
        D(7,5) = n*by_n_bta;
        D(7,6) = 1i*k*(n*by_n_bta-beta*a*by_np1_bta);
        D(7,7) = -(n*besselj(n,gamma1*a)-gamma1*a*besselj(n+1,gamma1*a));
        D(7,8) = 0;
        
        D(8,1) = n*bj_n_alfb-alpha*b*bj_np1_alfb;
        D(8,2) = n*bj_n_btb;
        D(8,3) = 1i*k*(n*bj_n_btb-beta*b*bj_np1_btb);
        D(8,4) = n*by_n_alfb-alpha*b*by_np1_alfb;
        D(8,5) = n*by_n_btb;
        D(8,6) = 1i*k*(n*by_n_btb-beta*b*by_np1_btb);
        D(8,7) = 0;
        D(8,8) = -(n*besselh(n,1,gamma2*b)-gamma2*b*besselh(n+1,1,gamma2*b));
        
        fdisp = det(D);
    case 'Va-So-So-Va'
        %外部自由内部自由覆阻尼圆柱壳【Outside and Inside Free Coated Shell】
        alphaem = sqrt(w^2/clem^2-k^2); betaem = sqrt(w^2/ctem^2-k^2);
        alphavm = sqrt(w^2/clvm^2-k^2); betavm = sqrt(w^2/ctvm^2-k^2);
        n =  nMode_cs;
        
        bj_n_alfema = besselj(n,alphaem*a); bj_np1_alfema = besselj(n+1,alphaem*a);
        bj_n_btema = besselj(n,betaem*a); bj_np1_btema = besselj(n+1,betaem*a);
        by_n_alfema = bessely(n,alphaem*a); by_np1_alfema = bessely(n+1,alphaem*a);
        by_n_btema = bessely(n,betaem*a); by_np1_btema = bessely(n+1,betaem*a);
        bj_n_alfemb = besselj(n,alphaem*b); bj_np1_alfemb = besselj(n+1,alphaem*b);
        bj_n_btemb = besselj(n,betaem*b); bj_np1_btemb = besselj(n+1,betaem*b);
        by_n_alfemb = bessely(n,alphaem*b); by_np1_alfemb = bessely(n+1,alphaem*b);
        by_n_btemb = bessely(n,betaem*b); by_np1_btemb = bessely(n+1,betaem*b);
        
        bj_n_alfvmb = besselj(n,alphavm*b); bj_np1_alfvmb = besselj(n+1,alphavm*b);
        bj_n_btvmb = besselj(n,betavm*b); bj_np1_btvmb = besselj(n+1,betavm*b);
        by_n_alfvmb = bessely(n,alphavm*b); by_np1_alfvmb = bessely(n+1,alphavm*b);
        by_n_btvmb = bessely(n,betavm*b); by_np1_btvmb = bessely(n+1,betavm*b);
        bj_n_alfvmc = besselj(n,alphavm*c); bj_np1_alfvmc = besselj(n+1,alphavm*c);
        bj_n_btvmc = besselj(n,betavm*c); bj_np1_btvmc = besselj(n+1,betavm*c);
        by_n_alfvmc = bessely(n,alphavm*c); by_np1_alfvmc = bessely(n+1,alphavm*c);
        by_n_btvmc = bessely(n,betavm*c); by_np1_btvmc = bessely(n+1,betavm*c);
        
        D = zeros(12);
        
        D(1,1) = (2*n*(n-1)-(betaem^2-k^2)*a^2)*bj_n_alfema +2*alphaem*a*bj_np1_alfema;
        D(1,2) = 2*n*(n-1)*bj_n_btema-2*n*betaem*a*bj_np1_btema;
        D(1,3) = 2*1i*k*((n*(n-1)-betaem^2*a^2)*bj_n_btema+betaem*a*bj_np1_btema);
        D(1,4) = (2*n*(n-1)-(betaem^2-k^2)*a^2)*by_n_alfema+2*alphaem*a*by_np1_alfema;
        D(1,5) = 2*n*(n-1)*by_n_btema-2*n*betaem*a*by_np1_btema;
        D(1,6) = 2*1i*k*((n*(n-1)-betaem^2*a^2)*by_n_btema+betaem*a*by_np1_btema);
        D(1,7) = 0; D(1,8) = 0; D(1,9) = 0; D(1,10) = 0; D(1,11) = 0; D(1,12) = 0;
        
        D(2,1) = -2*n*(n-1)*bj_n_alfema+2*n*alphaem*a*bj_np1_alfema;
        D(2,2) = (betaem^2*a^2-2*n*(n-1))*bj_n_btema-2*betaem*a*bj_np1_btema;
        D(2,3) = 2*1i*k*(-n*(n-1)*bj_n_btema+n*betaem*a*bj_np1_btema);
        D(2,4) = -2*n*(n-1)*by_n_alfema+2*n*alphaem*a*by_np1_alfema;
        D(2,5) = (betaem^2*a^2-2*n*(n-1))*by_n_btema-2*betaem*a*by_np1_btema;
        D(2,6) = 2*1i*k*(-n*(n-1)*by_n_btema+n*betaem*a*by_np1_btema);
        D(2,7) = 0; D(2,8) = 0; D(2,9) = 0; D(2,10) = 0; D(2,11) = 0; D(2,12) = 0;
        
        D(3,1) = 2*1i*k*(n*bj_n_alfema-alphaem*a*bj_np1_alfema);
        D(3,2) = 1i*k*n*bj_n_btema;
        D(3,3) = (betaem^2-k^2)*(n*bj_n_btema-betaem*a*bj_np1_btema);
        D(3,4) = 2*1i*k*(n*by_n_alfema-alphaem*a*by_np1_alfema);
        D(3,5) = 1i*k*n*by_n_btema;
        D(3,6) = (betaem^2-k^2)*(n*by_n_btema-betaem*a*by_np1_btema);
        D(3,7) = 0; D(3,8) = 0; D(3,9) = 0; D(3,10) = 0; D(3,11) = 0; D(3,12) = 0;
        
        D(4,1) = (2*n*(n-1)-(betaem^2-k^2)*b^2)*bj_n_alfemb+2*alphaem*b*bj_np1_alfemb;
        D(4,2) = 2*n*(n-1)*bj_n_btemb-2*n*betaem*b*bj_np1_btemb;
        D(4,3) = 2*1i*k*((n*(n-1)-betaem^2*b^2)*bj_n_btemb+betaem*b*bj_np1_btemb);
        D(4,4) = (2*n*(n-1)-(betaem^2-k^2)*b^2)*by_n_alfemb+2*alphaem*b*by_np1_alfemb;
        D(1,5) = 2*n*(n-1)*by_n_btemb-2*n*betaem*b*by_np1_btemb;
        D(4,6) = 2*1i*k*((n*(n-1)-betaem^2*b^2)*by_n_btemb+betaem*b*by_np1_btemb);
        D(4,7) = -miuvm/miuem*((2*n*(n-1)-(betavm^2-k^2)*b^2)*bj_n_alfvmb+2*alphavm*b*bj_np1_alfvmb);
        D(4,8) = -miuvm/miuem*(2*n*(n-1)*bj_n_btvmb-2*n*betavm*b*bj_np1_btvmb);
        D(4,9) = -miuvm/miuem*(2*1i*k*((n*(n-1)-betavm^2*b^2)*bj_n_btvmb+betavm*b*bj_np1_btvmb));
        D(4,10) = -miuvm/miuem*((2*n*(n-1)-(betavm^2-k^2)*b^2)*by_n_alfvmb+2*alphavm*b*by_np1_alfvmb);
        D(4,11) = -miuvm/miuem*(2*n*(n-1)*by_n_btvmb-2*n*betavm*b*by_np1_btvmb);
        D(4,12) = -miuvm/miuem*(2*1i*k*((n*(n-1)-betavm^2*b^2)*by_n_btvmb+betavm*b*by_np1_btvmb));
        
        D(5,1) = -2*n*(n-1)*bj_n_alfemb+2*n*alphaem*b*bj_np1_alfemb;
        D(5,2) = (betaem^2*b^2-2*n*(n-1))*bj_n_btemb-2*betaem*b*bj_np1_btemb;
        D(5,3) = 2*1i*k*(-n*(n-1)*bj_n_btemb+n*betaem*b*bj_np1_btemb);
        D(5,4) = -2*n*(n-1)*by_n_alfemb+2*n*alphaem*b*by_np1_alfemb;
        D(5,5) = (betaem^2*b^2-2*n*(n-1))*by_n_btemb-2*betaem*b*by_np1_btemb;
        D(5,6) = 2*1i*k*(-n*(n-1)*by_n_btemb+n*betaem*b*by_np1_btemb);
        D(5,7) = -miuvm/miuem*(-2*n*(n-1)*bj_n_alfvmb+2*n*alphavm*b*bj_np1_alfvmb);
        D(5,8) = -miuvm/miuem*((betavm^2*b^2-2*n*(n-1))*bj_n_btvmb-2*betavm*b*bj_np1_btvmb);
        D(5,9) = -miuvm/miuem*(2*1i*k*(-n*(n-1)*bj_n_btvmb+n*betavm*b*bj_np1_btvmb));
        D(5,10) = -miuvm/miuem*(-2*n*(n-1)*by_n_alfvmb+2*n*alphavm*b*by_np1_alfvmb);
        D(5,11) = -miuvm/miuem*((betavm^2*b^2-2*n*(n-1))*by_n_btvmb-2*betavm*b*by_np1_btvmb);
        D(5,12) = -miuvm/miuem*(2*1i*k*(-n*(n-1)*by_n_btvmb+n*betavm*b*by_np1_btvmb));
        
        D(6,1) = 2*1i*k*(n*bj_n_alfemb-alphaem*b*bj_np1_alfemb);
        D(6,2) = 1i*k*n*bj_n_btemb;
        D(6,3) = (betaem^2-k^2)*(n*bj_n_btemb-betaem*b*bj_np1_btemb);
        D(6,4) = 2*1i*k*(n*by_n_alfemb-alphaem*b*by_np1_alfemb);
        D(6,5) = 1i*k*n*by_n_btemb;
        D(6,6) = (betaem^2-k^2)*(n*by_n_btemb-betaem*b*by_np1_btemb);
        D(6,7) = -miuvm/miuem*(2*1i*k*(n*bj_n_alfvmb-alphavm*b*bj_np1_alfvmb));
        D(6,8) = -miuvm/miuem*(1i*k*n*bj_n_btvmb);
        D(6,9) = -miuvm/miuem*((betavm^2-k^2)*(n*bj_n_btvmb-betavm*b*bj_np1_btvmb));
        D(6,10) = -miuvm/miuem*(2*1i*k*(n*by_n_alfvmb-alphavm*b*by_np1_alfvmb));
        D(6,11) = -miuvm/miuem*(1i*k*n*by_n_btvmb);
        D(6,12) = -miuvm/miuem*((betavm^2-k^2)*(n*by_n_btvmb-betavm*b*by_np1_btvmb));
        
        D(7,1) = n*bj_n_alfemb-alphaem*b*bj_np1_alfemb;
        D(7,2) = n*bj_n_btemb;
        D(7,3) = 1i*k*(n*bj_n_btemb-betaem*b*bj_np1_btemb);
        D(7,4) = n*by_n_alfemb-alphaem*b*by_np1_alfemb;
        D(7,5) = n*by_n_btemb;
        D(7,6) = 1i*k*(n*by_n_btemb-betaem*b*by_np1_btemb);
        D(7,7) = -(n*bj_n_alfvmb-alphavm*b*bj_np1_alfvmb);
        D(7,8) = -(n*bj_n_btvmb);
        D(7,9) = -1i*k*(n*bj_n_btvmb-betavm*b*bj_np1_btvmb);
        D(7,10) = -(n*by_n_alfvmb-alphavm*b*by_np1_alfvmb);
        D(7,11) = -(n*by_n_btvmb);
        D(7,12) = -1i*k*(n*by_n_btvmb-betavm*b*by_np1_btvmb);
        
        D(8,1) = n*bj_n_alfemb;
        D(8,2) = n*bj_n_btemb-betaem*b*bj_np1_btemb;
        D(8,3) = 1i*k*n*bj_n_btemb;
        D(8,4) = n*by_n_alfemb;
        D(8,5) = n*by_n_btemb-betaem*b*by_np1_btemb;
        D(8,6) = 1i*k*n*by_n_btemb;
        D(8,7) = -(n*bj_n_alfvmb);
        D(8,8) = -(n*bj_n_btvmb-betavm*b*bj_np1_btvmb);
        D(8,9) = -(1i*k*n*bj_n_btvmb);
        D(8,10) = -(n*by_n_alfvmb);
        D(8,11) = -(n*by_n_btvmb-betavm*b*by_np1_btvmb);
        D(8,12) = -(1i*k*n*by_n_btvmb);
        
        D(9,1) = 1i*k*bj_n_alfemb;
        D(9,2) = 0;
        D(9,3) = betaem^2*bj_n_btemb;
        D(9,4) = 1i*k*by_n_alfemb;
        D(9,5) = 0;
        D(9,6) = betaem^2*by_n_btemb;
        D(9,7) = -(1i*k*bj_n_alfvmb);
        D(9,8) = 0;
        D(9,9) = -(betavm^2*bj_n_btvmb);
        D(9,10) = -(1i*k*by_n_alfvmb);
        D(9,11) = 0;
        D(9,12) = -(betavm^2*by_n_btvmb);
        
        D(10,1) = 0; D(10,2) = 0; D(10,3) = 0; D(10,4) = 0; D(10,5) = 0; D(10,6) = 0;
        D(10,7) = (2*n*(n-1)-(betavm^2-k^2)*c^2)*bj_n_alfvmc+2*alphavm*c*bj_np1_alfvmc;
        D(10,8) = 2*n*(n-1)*bj_n_btvmc-2*n*betavm*c*bj_np1_btvmc;
        D(10,9) = 2*1i*k*((n*(n-1)-betavm^2*c^2)*bj_n_btvmc+betavm*c*bj_np1_btvmc);
        D(10,10) = (2*n*(n-1)-(betavm^2-k^2)*c^2)*by_n_alfvmc+2*alphavm*c*by_np1_alfvmc;
        D(10,11) = 2*n*(n-1)*by_n_btvmc-2*n*betavm*c*by_np1_btvmc;
        D(10,12) = 2*1i*k*((n*(n-1)-betavm^2*c^2)*by_n_btvmc+betavm*c*by_np1_btvmc);
        
        D(11,1) = 0; D(11,2) = 0; D(11,3) = 0; D(11,4) = 0; D(11,5) = 0; D(11,6) = 0;
        D(11,7) = -2*n*(n-1)*bj_n_alfvmc+2*n*alphavm*c*bj_np1_alfvmc;
        D(11,8) = (betavm^2*c^2-2*n*(n-1))*bj_n_btvmc-2*betavm*c*bj_np1_btvmc;
        D(11,9) = 2*1i*k*(-n*(n-1)*bj_n_btvmc+n*betavm*c*bj_np1_btvmc);
        D(11,10) = -2*n*(n-1)*by_n_alfvmc+2*n*alphavm*c*by_np1_alfvmc;
        D(11,11) = (betavm^2*c^2-2*n*(n-1))*by_n_btvmc-2*betavm*c*by_np1_btvmc;
        D(11,12) = 2*1i*k*(-n*(n-1)*by_n_btvmc+n*betavm*c*by_np1_btvmc);
        
        D(12,1) = 0; D(12,2) = 0; D(12,3) = 0; D(12,4) = 0; D(12,5) = 0; D(12,6) = 0;
        D(12,7) = 2*1i*k*(n*bj_n_alfvmc-alphavm*c*bj_np1_alfvmc);
        D(12,8) = 1i*k*n*bj_n_btvmc;
        D(12,9) = (betavm^2-k^2)*(n*bj_n_btvmc-betavm*c*bj_np1_btvmc);
        D(12,10) = 2*1i*k*(n*by_n_alfvmc-alphavm*c*by_np1_alfvmc);
        D(12,11) = 1i*k*n*by_n_btvmc;
        D(12,12) = (betavm^2-k^2)*(n*by_n_btvmc-betavm*c*by_np1_btvmc);
        
        fdisp = det(D);
    case 'Va-So-So-Fl'
        %外部自由内部充液覆阻尼圆柱壳【Outside Free Inside Fluid Coated Shell】
        alphaem = sqrt(w^2/clem^2-k^2); betaem = sqrt(w^2/ctem^2-k^2);
        alphavm = sqrt(w^2/clvm^2-k^2); betavm = sqrt(w^2/ctvm^2-k^2);
        gamma2 = sqrt(w^2/c2^2-k^2);
        n =  nMode_cs;
        
        bj_n_alfema = besselj(n,alphaem*a); bj_np1_alfema = besselj(n+1,alphaem*a);
        bj_n_btema = besselj(n,betaem*a); bj_np1_btema = besselj(n+1,betaem*a);
        by_n_alfema = bessely(n,alphaem*a); by_np1_alfema = bessely(n+1,alphaem*a);
        by_n_btema = bessely(n,betaem*a); by_np1_btema = bessely(n+1,betaem*a);
        bj_n_alfemb = besselj(n,alphaem*b); bj_np1_alfemb = besselj(n+1,alphaem*b);
        bj_n_btemb = besselj(n,betaem*b); bj_np1_btemb = besselj(n+1,betaem*b);
        by_n_alfemb = bessely(n,alphaem*b); by_np1_alfemb = bessely(n+1,alphaem*b);
        by_n_btemb = bessely(n,betaem*b); by_np1_btemb = bessely(n+1,betaem*b);
        
        bj_n_alfvmb = besselj(n,alphavm*b); bj_np1_alfvmb = besselj(n+1,alphavm*b);
        bj_n_btvmb = besselj(n,betavm*b); bj_np1_btvmb = besselj(n+1,betavm*b);
        by_n_alfvmb = bessely(n,alphavm*b); by_np1_alfvmb = bessely(n+1,alphavm*b);
        by_n_btvmb = bessely(n,betavm*b); by_np1_btvmb = bessely(n+1,betavm*b);
        bj_n_alfvmc = besselj(n,alphavm*c); bj_np1_alfvmc = besselj(n+1,alphavm*c);
        bj_n_btvmc = besselj(n,betavm*c); bj_np1_btvmc = besselj(n+1,betavm*c);
        by_n_alfvmc = bessely(n,alphavm*c); by_np1_alfvmc = bessely(n+1,alphavm*c);
        by_n_btvmc = bessely(n,betavm*c); by_np1_btvmc = bessely(n+1,betavm*c);
        
        D = zeros(13);
        
        D(1,1) = (2*n*(n-1)-(betaem^2-k^2)*a^2)*bj_n_alfema +2*alphaem*a*bj_np1_alfema;
        D(1,2) = 2*n*(n-1)*bj_n_btema-2*n*betaem*a*bj_np1_btema;
        D(1,3) = 2*1i*k*((n*(n-1)-betaem^2*a^2)*bj_n_btema+betaem*a*bj_np1_btema);
        D(1,4) = (2*n*(n-1)-(betaem^2-k^2)*a^2)*by_n_alfema+2*alphaem*a*by_np1_alfema;
        D(1,5) = 2*n*(n-1)*by_n_btema-2*n*betaem*a*by_np1_btema;
        D(1,6) = 2*1i*k*((n*(n-1)-betaem^2*a^2)*by_n_btema+betaem*a*by_np1_btema);
        D(1,7) = 0; D(1,8) = 0; D(1,9) = 0; D(1,10) = 0; D(1,11) = 0; D(1,12) = 0; D(1,13) = 0;
        
        D(2,1) = -2*n*(n-1)*bj_n_alfema+2*n*alphaem*a*bj_np1_alfema;
        D(2,2) = (betaem^2*a^2-2*n*(n-1))*bj_n_btema-2*betaem*a*bj_np1_btema;
        D(2,3) = 2*1i*k*(-n*(n-1)*bj_n_btema+n*betaem*a*bj_np1_btema);
        D(2,4) = -2*n*(n-1)*by_n_alfema+2*n*alphaem*a*by_np1_alfema;
        D(2,5) = (betaem^2*a^2-2*n*(n-1))*by_n_btema-2*betaem*a*by_np1_btema;
        D(2,6) = 2*1i*k*(-n*(n-1)*by_n_btema+n*betaem*a*by_np1_btema);
        D(2,7) = 0; D(2,8) = 0; D(2,9) = 0; D(2,10) = 0; D(2,11) = 0; D(2,12) = 0; D(2,13) = 0;
        
        D(3,1) = 2*1i*k*(n*bj_n_alfema-alphaem*a*bj_np1_alfema);
        D(3,2) = 1i*k*n*bj_n_btema;
        D(3,3) = (betaem^2-k^2)*(n*bj_n_btema-betaem*a*bj_np1_btema);
        D(3,4) = 2*1i*k*(n*by_n_alfema-alphaem*a*by_np1_alfema);
        D(3,5) = 1i*k*n*by_n_btema;
        D(3,6) = (betaem^2-k^2)*(n*by_n_btema-betaem*a*by_np1_btema);
        D(3,7) = 0; D(3,8) = 0; D(3,9) = 0; D(3,10) = 0; D(3,11) = 0; D(3,12) = 0; D(3,13) = 0;
        
        D(4,1) = (2*n*(n-1)-(betaem^2-k^2)*b^2)*bj_n_alfemb+2*alphaem*b*bj_np1_alfemb;
        D(4,2) = 2*n*(n-1)*bj_n_btemb-2*n*betaem*b*bj_np1_btemb;
        D(4,3) = 2*1i*k*((n*(n-1)-betaem^2*b^2)*bj_n_btemb+betaem*b*bj_np1_btemb);
        D(4,4) = (2*n*(n-1)-(betaem^2-k^2)*b^2)*by_n_alfemb+2*alphaem*b*by_np1_alfemb;
        D(1,5) = 2*n*(n-1)*by_n_btemb-2*n*betaem*b*by_np1_btemb;
        D(4,6) = 2*1i*k*((n*(n-1)-betaem^2*b^2)*by_n_btemb+betaem*b*by_np1_btemb);
        D(4,7) = -miuvm/miuem*((2*n*(n-1)-(betavm^2-k^2)*b^2)*bj_n_alfvmb +2*alphavm*b*bj_np1_alfvmb);
        D(4,8) = -miuvm/miuem*(2*n*(n-1)*bj_n_btvmb-2*n*betavm*b*bj_np1_btvmb);
        D(4,9) = -miuvm/miuem*(2*1i*k*((n*(n-1)-betavm^2*b^2)*bj_n_btvmb+betavm*b*bj_np1_btvmb));
        D(4,10) = -miuvm/miuem*((2*n*(n-1)-(betavm^2-k^2)*b^2)*by_n_alfvmb+2*alphavm*b*by_np1_alfvmb);
        D(4,11) = -miuvm/miuem*(2*n*(n-1)*by_n_btvmb-2*n*betavm*b*by_np1_btvmb);
        D(4,12) = -miuvm/miuem*(2*1i*k*((n*(n-1)-betavm^2*b^2)*by_n_btvmb+betavm*b*by_np1_btvmb));
        D(4,13) = 0;
        
        D(5,1) = -2*n*(n-1)*bj_n_alfemb+2*n*alphaem*b*bj_np1_alfemb;
        D(5,2) = (betaem^2*b^2-2*n*(n-1))*bj_n_btemb-2*betaem*b*bj_np1_btemb;
        D(5,3) = 2*1i*k*(-n*(n-1)*bj_n_btemb+n*betaem*b*bj_np1_btemb);
        D(5,4) = -2*n*(n-1)*by_n_alfemb+2*n*alphaem*b*by_np1_alfemb;
        D(5,5) = (betaem^2*b^2-2*n*(n-1))*by_n_btemb-2*betaem*b*by_np1_btemb;
        D(5,6) = 2*1i*k*(-n*(n-1)*by_n_btemb+n*betaem*b*by_np1_btemb);
        D(5,7) = -miuvm/miuem*(-2*n*(n-1)*bj_n_alfvmb+2*n*alphavm*b*bj_np1_alfvmb);
        D(5,8) = -miuvm/miuem*((betavm^2*b^2-2*n*(n-1))*bj_n_btvmb-2*betavm*b*bj_np1_btvmb);
        D(5,9) = -miuvm/miuem*(2*1i*k*(-n*(n-1)*bj_n_btvmb+n*betavm*b*bj_np1_btvmb));
        D(5,10) = -miuvm/miuem*(-2*n*(n-1)*by_n_alfvmb+2*n*alphavm*b*by_np1_alfvmb);
        D(5,11) = -miuvm/miuem*((betavm^2*b^2-2*n*(n-1))*by_n_btvmb-2*betavm*b*by_np1_btvmb);
        D(5,12) = -miuvm/miuem*(2*1i*k*(-n*(n-1)*by_n_btvmb+n*betavm*b*by_np1_btvmb));
        D(5,13) = 0;
        
        D(6,1) = 2*1i*k*(n*bj_n_alfemb-alphaem*b*bj_np1_alfemb);
        D(6,2) = 1i*k*n*bj_n_btemb;
        D(6,3) = (betaem^2-k^2)*(n*bj_n_btemb-betaem*b*bj_np1_btemb);
        D(6,4) = 2*1i*k*(n*by_n_alfemb-alphaem*b*by_np1_alfemb);
        D(6,5) = 1i*k*n*by_n_btemb;
        D(6,6) = (betaem^2-k^2)*(n*by_n_btemb-betaem*b*by_np1_btemb);
        D(6,7) = -miuvm/miuem*(2*1i*k*(n*bj_n_alfvmb-alphavm*b*bj_np1_alfvmb));
        D(6,8) = -miuvm/miuem*(1i*k*n*bj_n_btvmb);
        D(6,9) = -miuvm/miuem*((betavm^2-k^2)*(n*bj_n_btvmb-betavm*b*bj_np1_btvmb));
        D(6,10) = -miuvm/miuem*(2*1i*k*(n*by_n_alfvmb-alphavm*b*by_np1_alfvmb));
        D(6,11) = -miuvm/miuem*(1i*k*n*by_n_btvmb);
        D(6,12) = -miuvm/miuem*((betavm^2-k^2)*(n*by_n_btvmb-betavm*b*by_np1_btvmb));
        D(6,13) = 0;
        
        D(7,1) = n*bj_n_alfemb-alphaem*b*bj_np1_alfemb;
        D(7,2) = n*bj_n_btemb;
        D(7,3) = 1i*k*(n*bj_n_btemb-betaem*b*bj_np1_btemb);
        D(7,4) = n*by_n_alfemb-alphaem*b*by_np1_alfemb;
        D(7,5) = n*by_n_btemb;
        D(7,6) = 1i*k*(n*by_n_btemb-betaem*b*by_np1_btemb);
        D(7,7) = -(n*bj_n_alfvmb-alphavm*b*bj_np1_alfvmb);
        D(7,8) = -(n*bj_n_btvmb);
        D(7,9) = -1i*k*(n*bj_n_btvmb-betavm*b*bj_np1_btvmb);
        D(7,10) = -(n*by_n_alfvmb-alphavm*b*by_np1_alfvmb);
        D(7,11) = -(n*by_n_btvmb);
        D(7,12) = -1i*k*(n*by_n_btvmb-betavm*b*by_np1_btvmb);
        D(7,13) = 0;
        
        D(8,1) = n*bj_n_alfemb;
        D(8,2) = n*bj_n_btemb-betaem*b*bj_np1_btemb;
        D(8,3) = 1i*k*n*bj_n_btemb;
        D(8,4) = n*by_n_alfemb;
        D(8,5) = n*by_n_btemb-betaem*b*by_np1_btemb;
        D(8,6) = 1i*k*n*by_n_btemb;
        D(8,7) = -(n*bj_n_alfvmb);
        D(8,8) = -(n*bj_n_btvmb-betavm*b*bj_np1_btvmb);
        D(8,9) = -(1i*k*n*bj_n_btvmb);
        D(8,10) = -(n*by_n_alfvmb);
        D(8,11) = -(n*by_n_btvmb-betavm*b*by_np1_btvmb);
        D(8,12) = -(1i*k*n*by_n_btvmb);
        D(8,13) = 0;
        
        D(9,1) = 1i*k*bj_n_alfemb;
        D(9,2) = 0;
        D(9,3) = betaem^2*bj_n_btemb;
        D(9,4) = 1i*k*by_n_alfemb;
        D(9,5) = 0;
        D(9,6) = betaem^2*by_n_btemb;
        D(9,7) = -(1i*k*bj_n_alfvmb);
        D(9,8) = 0;
        D(9,9) = -(betavm^2*bj_n_btvmb);
        D(9,10) = -(1i*k*by_n_alfvmb);
        D(9,11) = 0;
        D(9,12) = -(betavm^2*by_n_btvmb);
        D(9,13) = 0;
        
        D(10,1) = 0; D(10,2) = 0; D(10,3) = 0; D(10,4) = 0; D(10,5) = 0; D(10,6) = 0;
        D(10,7) = (2*n*(n-1)-(betavm^2-k^2)*c^2)*bj_n_alfvmc +2*alphavm*c*bj_np1_alfvmc;
        D(10,8) = 2*n*(n-1)*bj_n_btvmc-2*n*betavm*c*bj_np1_btvmc;
        D(10,9) = 2*1i*k*((n*(n-1)-betavm^2*c^2)*bj_n_btvmc+betavm*c*bj_np1_btvmc);
        D(10,10) = (2*n*(n-1)-(betavm^2-k^2)*c^2)*by_n_alfvmc+2*alphavm*c*by_np1_alfvmc;
        D(10,11) = 2*n*(n-1)*by_n_btvmc-2*n*betavm*c*by_np1_btvmc;
        D(10,12) = 2*1i*k*((n*(n-1)-betavm^2*c^2)*by_n_btvmc+betavm*c*by_np1_btvmc);
        D(10,13) = c^2*row2*w^2*besselh(n,1,gamma2*c)/miuvm;
        
        D(11,1) = 0; D(11,2) = 0; D(11,3) = 0; D(11,4) = 0; D(11,5) = 0; D(11,6) = 0;
        D(11,7) = -2*n*(n-1)*bj_n_alfvmc+2*n*alphavm*c*bj_np1_alfvmc;
        D(11,8) = (betavm^2*c^2-2*n*(n-1))*bj_n_btvmc-2*betavm*c*bj_np1_btvmc;
        D(11,9) = 2*1i*k*(-n*(n-1)*bj_n_btvmc+n*betavm*c*bj_np1_btvmc);
        D(11,10) = -2*n*(n-1)*by_n_alfvmc+2*n*alphavm*c*by_np1_alfvmc;
        D(11,11) = (betavm^2*c^2-2*n*(n-1))*by_n_btvmc-2*betavm*c*by_np1_btvmc;
        D(11,12) = 2*1i*k*(-n*(n-1)*by_n_btvmc+n*betavm*c*by_np1_btvmc);
        D(11,13) = 0;
        
        D(12,1) = 0; D(12,2) = 0; D(12,3) = 0; D(12,4) = 0; D(12,5) = 0; D(12,6) = 0;
        D(12,7) = 2*1i*k*(n*bj_n_alfvmc-alphavm*c*bj_np1_alfvmc);
        D(12,8) = 1i*k*n*bj_n_btvmc;
        D(12,9) = (betavm^2-k^2)*(n*bj_n_btvmc-betavm*c*bj_np1_btvmc);
        D(12,10) = 2*1i*k*(n*by_n_alfvmc-alphavm*c*by_np1_alfvmc);
        D(12,11) = 1i*k*n*by_n_btvmc;
        D(12,12) = (betavm^2-k^2)*(n*by_n_btvmc-betavm*c*by_np1_btvmc);
        D(11,13) = 0;
        
        D(13,1) = 0; D(13,2) = 0; D(13,3) = 0; D(13,4) = 0; D(13,5) = 0; D(13,6) = 0;
        D(13,7) = n*bj_n_alfvmc-alphavm*c*bj_np1_alfvmc;
        D(13,8) = n*bj_n_btvmc;
        D(13,9) = 1i*k*(n*bj_n_btvmc-betavm*c*bj_np1_btvmc);
        D(13,10) = n*by_n_alfvmc-alphavm*c*by_np1_alfvmc;
        D(13,11) = n*by_n_btvmc;
        D(13,12) = 1i*k*(n*by_n_btvmc-betavm*c*by_np1_btvmc);
        D(13,13) = -(n*besselh(n,1,gamma2*c)-gamma2*c*besselh(n+1,1,gamma2*c));
        
        fdisp = det(D);
    case 'Fl-So-So-Fl'
        %外部充液内部充液覆阻尼圆柱壳【Outside and Inside Fluid Coated Shell】
        alphaem = sqrt(w^2/clem^2-k^2); betaem = sqrt(w^2/ctem^2-k^2);
        alphavm = sqrt(w^2/clvm^2-k^2); betavm = sqrt(w^2/ctvm^2-k^2);
        gamma1 = sqrt(w^2/c1^2-k^2); gamma2 = sqrt(w^2/c2^2-k^2);
        n =  nMode_cs;
        
        bj_n_alfema = besselj(n,alphaem*a); bj_np1_alfema = besselj(n+1,alphaem*a);
        bj_n_btema = besselj(n,betaem*a); bj_np1_btema = besselj(n+1,betaem*a);
        by_n_alfema = bessely(n,alphaem*a); by_np1_alfema = bessely(n+1,alphaem*a);
        by_n_btema = bessely(n,betaem*a); by_np1_btema = bessely(n+1,betaem*a);
        bj_n_alfemb = besselj(n,alphaem*b); bj_np1_alfemb = besselj(n+1,alphaem*b);
        bj_n_btemb = besselj(n,betaem*b); bj_np1_btemb = besselj(n+1,betaem*b);
        by_n_alfemb = bessely(n,alphaem*b); by_np1_alfemb = bessely(n+1,alphaem*b);
        by_n_btemb = bessely(n,betaem*b); by_np1_btemb = bessely(n+1,betaem*b);
        
        bj_n_alfvmb = besselj(n,alphavm*b); bj_np1_alfvmb = besselj(n+1,alphavm*b);
        bj_n_btvmb = besselj(n,betavm*b); bj_np1_btvmb = besselj(n+1,betavm*b);
        by_n_alfvmb = bessely(n,alphavm*b); by_np1_alfvmb = bessely(n+1,alphavm*b);
        by_n_btvmb = bessely(n,betavm*b); by_np1_btvmb = bessely(n+1,betavm*b);
        bj_n_alfvmc = besselj(n,alphavm*c); bj_np1_alfvmc = besselj(n+1,alphavm*c);
        bj_n_btvmc = besselj(n,betavm*c); bj_np1_btvmc = besselj(n+1,betavm*c);
        by_n_alfvmc = bessely(n,alphavm*c); by_np1_alfvmc = bessely(n+1,alphavm*c);
        by_n_btvmc = bessely(n,betavm*c); by_np1_btvmc = bessely(n+1,betavm*c);
        
        D = zeros(14);
        
        D(1,1) = (2*n*(n-1)-(betaem^2-k^2)*a^2)*bj_n_alfema +2*alphaem*a*bj_np1_alfema;
        D(1,2) = 2*n*(n-1)*bj_n_btema-2*n*betaem*a*bj_np1_btema;
        D(1,3) = 2*1i*k*((n*(n-1)-betaem^2*a^2)*bj_n_btema+betaem*a*bj_np1_btema);
        D(1,4) = (2*n*(n-1)-(betaem^2-k^2)*a^2)*by_n_alfema+2*alphaem*a*by_np1_alfema;
        D(1,5) = 2*n*(n-1)*by_n_btema-2*n*betaem*a*by_np1_btema;
        D(1,6) = 2*1i*k*((n*(n-1)-betaem^2*a^2)*by_n_btema+betaem*a*by_np1_btema);
        D(1,7) = 0; D(1,8) = 0; D(1,9) = 0; D(1,10) = 0; D(1,11) = 0; D(1,12) = 0;
        D(1,13) = a^2*row1*w^2*besselj(n,gamma1*a)/miuem;
        D(1,14) = 0;
        
        D(2,1) = -2*n*(n-1)*bj_n_alfema+2*n*alphaem*a*bj_np1_alfema;
        D(2,2) = (betaem^2*a^2-2*n*(n-1))*bj_n_btema-2*betaem*a*bj_np1_btema;
        D(2,3) = 2*1i*k*(-n*(n-1)*bj_n_btema+n*betaem*a*bj_np1_btema);
        D(2,4) = -2*n*(n-1)*by_n_alfema+2*n*alphaem*a*by_np1_alfema;
        D(2,5) = (betaem^2*a^2-2*n*(n-1))*by_n_btema-2*betaem*a*by_np1_btema;
        D(2,6) = 2*1i*k*(-n*(n-1)*by_n_btema+n*betaem*a*by_np1_btema);
        D(2,7) = 0; D(2,8) = 0; D(2,9) = 0; D(2,10) = 0; D(2,11) = 0;
        D(2,12) = 0; D(2,13) = 0; D(2,14) = 0;
        
        D(3,1) = 2*1i*k*(n*bj_n_alfema-alphaem*a*bj_np1_alfema);
        D(3,2) = 1i*k*n*bj_n_btema;
        D(3,3) = (betaem^2-k^2)*(n*bj_n_btema-betaem*a*bj_np1_btema);
        D(3,4) = 2*1i*k*(n*by_n_alfema-alphaem*a*by_np1_alfema);
        D(3,5) = 1i*k*n*by_n_btema;
        D(3,6) = (betaem^2-k^2)*(n*by_n_btema-betaem*a*by_np1_btema);
        D(3,7) = 0; D(3,8) = 0; D(3,9) = 0; D(3,10) = 0; D(3,11) = 0;
        D(3,12) = 0; D(3,13) = 0; D(3,14) = 0;
        
        D(4,1) = (2*n*(n-1)-(betaem^2-k^2)*b^2)*bj_n_alfemb+2*alphaem*b*bj_np1_alfemb;
        D(4,2) = 2*n*(n-1)*bj_n_btemb-2*n*betaem*b*bj_np1_btemb;
        D(4,3) = 2*1i*k*((n*(n-1)-betaem^2*b^2)*bj_n_btemb+betaem*b*bj_np1_btemb);
        D(4,4) = (2*n*(n-1)-(betaem^2-k^2)*b^2)*by_n_alfemb+2*alphaem*b*by_np1_alfemb;
        D(1,5) = 2*n*(n-1)*by_n_btemb-2*n*betaem*b*by_np1_btemb;
        D(4,6) = 2*1i*k*((n*(n-1)-betaem^2*b^2)*by_n_btemb+betaem*b*by_np1_btemb);
        D(4,7) = -miuvm/miuem*((2*n*(n-1)-(betavm^2-k^2)*b^2)*bj_n_alfvmb +2*alphavm*b*bj_np1_alfvmb);
        D(4,8) = -miuvm/miuem*(2*n*(n-1)*bj_n_btvmb-2*n*betavm*b*bj_np1_btvmb);
        D(4,9) = -miuvm/miuem*(2*1i*k*((n*(n-1)-betavm^2*b^2)*bj_n_btvmb+betavm*b*bj_np1_btvmb));
        D(4,10) = -miuvm/miuem*((2*n*(n-1)-(betavm^2-k^2)*b^2)*by_n_alfvmb+2*alphavm*b*by_np1_alfvmb);
        D(4,11) = -miuvm/miuem*(2*n*(n-1)*by_n_btvmb-2*n*betavm*b*by_np1_btvmb);
        D(4,12) = -miuvm/miuem*(2*1i*k*((n*(n-1)-betavm^2*b^2)*by_n_btvmb+betavm*b*by_np1_btvmb));
        D(4,13) = 0; D(4,14) = 0;
        
        D(5,1) = -2*n*(n-1)*bj_n_alfemb+2*n*alphaem*b*bj_np1_alfemb;
        D(5,2) = (betaem^2*b^2-2*n*(n-1))*bj_n_btemb-2*betaem*b*bj_np1_btemb;
        D(5,3) = 2*1i*k*(-n*(n-1)*bj_n_btemb+n*betaem*b*bj_np1_btemb);
        D(5,4) = -2*n*(n-1)*by_n_alfemb+2*n*alphaem*b*by_np1_alfemb;
        D(5,5) = (betaem^2*b^2-2*n*(n-1))*by_n_btemb-2*betaem*b*by_np1_btemb;
        D(5,6) = 2*1i*k*(-n*(n-1)*by_n_btemb+n*betaem*b*by_np1_btemb);
        D(5,7) = -miuvm/miuem*(-2*n*(n-1)*bj_n_alfvmb+2*n*alphavm*b*bj_np1_alfvmb);
        D(5,8) = -miuvm/miuem*((betavm^2*b^2-2*n*(n-1))*bj_n_btvmb-2*betavm*b*bj_np1_btvmb);
        D(5,9) = -miuvm/miuem*(2*1i*k*(-n*(n-1)*bj_n_btvmb+n*betavm*b*bj_np1_btvmb));
        D(5,10) = -miuvm/miuem*(-2*n*(n-1)*by_n_alfvmb+2*n*alphavm*b*by_np1_alfvmb);
        D(5,11) = -miuvm/miuem*((betavm^2*b^2-2*n*(n-1))*by_n_btvmb-2*betavm*b*by_np1_btvmb);
        D(5,12) = -miuvm/miuem*(2*1i*k*(-n*(n-1)*by_n_btvmb+n*betavm*b*by_np1_btvmb));
        D(5,13) = 0; D(5,14) = 0;
        
        D(6,1) = 2*1i*k*(n*bj_n_alfemb-alphaem*b*bj_np1_alfemb);
        D(6,2) = 1i*k*n*bj_n_btemb;
        D(6,3) = (betaem^2-k^2)*(n*bj_n_btemb-betaem*b*bj_np1_btemb);
        D(6,4) = 2*1i*k*(n*by_n_alfemb-alphaem*b*by_np1_alfemb);
        D(6,5) = 1i*k*n*by_n_btemb;
        D(6,6) = (betaem^2-k^2)*(n*by_n_btemb-betaem*b*by_np1_btemb);
        D(6,7) = -miuvm/miuem*(2*1i*k*(n*bj_n_alfvmb-alphavm*b*bj_np1_alfvmb));
        D(6,8) = -miuvm/miuem*(1i*k*n*bj_n_btvmb);
        D(6,9) = -miuvm/miuem*((betavm^2-k^2)*(n*bj_n_btvmb-betavm*b*bj_np1_btvmb));
        D(6,10) = -miuvm/miuem*(2*1i*k*(n*by_n_alfvmb-alphavm*b*by_np1_alfvmb));
        D(6,11) = -miuvm/miuem*(1i*k*n*by_n_btvmb);
        D(6,12) = -miuvm/miuem*((betavm^2-k^2)*(n*by_n_btvmb-betavm*b*by_np1_btvmb));
        D(6,13) = 0; D(6,14) = 0;
        
        D(7,1) = n*bj_n_alfemb-alphaem*b*bj_np1_alfemb;
        D(7,2) = n*bj_n_btemb;
        D(7,3) = 1i*k*(n*bj_n_btemb-betaem*b*bj_np1_btemb);
        D(7,4) = n*by_n_alfemb-alphaem*b*by_np1_alfemb;
        D(7,5) = n*by_n_btemb;
        D(7,6) = 1i*k*(n*by_n_btemb-betaem*b*by_np1_btemb);
        D(7,7) = -(n*bj_n_alfvmb-alphavm*b*bj_np1_alfvmb);
        D(7,8) = -(n*bj_n_btvmb);
        D(7,9) = -1i*k*(n*bj_n_btvmb-betavm*b*bj_np1_btvmb);
        D(7,10) = -(n*by_n_alfvmb-alphavm*b*by_np1_alfvmb);
        D(7,11) = -(n*by_n_btvmb);
        D(7,12) = -1i*k*(n*by_n_btvmb-betavm*b*by_np1_btvmb);
        D(7,13) = 0; D(7,14) = 0;
        
        D(8,1) = n*bj_n_alfemb;
        D(8,2) = n*bj_n_btemb-betaem*b*bj_np1_btemb;
        D(8,3) = 1i*k*n*bj_n_btemb;
        D(8,4) = n*by_n_alfemb;
        D(8,5) = n*by_n_btemb-betaem*b*by_np1_btemb;
        D(8,6) = 1i*k*n*by_n_btemb;
        D(8,7) = -(n*bj_n_alfvmb);
        D(8,8) = -(n*bj_n_btvmb-betavm*b*bj_np1_btvmb);
        D(8,9) = -(1i*k*n*bj_n_btvmb);
        D(8,10) = -(n*by_n_alfvmb);
        D(8,11) = -(n*by_n_btvmb-betavm*b*by_np1_btvmb);
        D(8,12) = -(1i*k*n*by_n_btvmb);
        D(8,13) = 0; D(8,14) = 0;
        
        D(9,1) = 1i*k*bj_n_alfemb;
        D(9,2) = 0;
        D(9,3) = betaem^2*bj_n_btemb;
        D(9,4) = 1i*k*by_n_alfemb;
        D(9,5) = 0;
        D(9,6) = betaem^2*by_n_btemb;
        D(9,7) = -(1i*k*bj_n_alfvmb);
        D(9,8) = 0;
        D(9,9) = -(betavm^2*bj_n_btvmb);
        D(9,10) = -(1i*k*by_n_alfvmb);
        D(9,11) = 0;
        D(9,12) = -(betavm^2*by_n_btvmb);
        D(9,13) = 0; D(9,14) = 0;
        
        D(10,1) = 0; D(10,2) = 0; D(10,3) = 0; D(10,4) = 0; D(10,5) = 0; D(10,6) = 0;
        D(10,7) = (2*n*(n-1)-(betavm^2-k^2)*c^2)*bj_n_alfvmc +2*alphavm*c*bj_np1_alfvmc;
        D(10,8) = 2*n*(n-1)*bj_n_btvmc-2*n*betavm*c*bj_np1_btvmc;
        D(10,9) = 2*1i*k*((n*(n-1)-betavm^2*c^2)*bj_n_btvmc+betavm*c*bj_np1_btvmc);
        D(10,10) = (2*n*(n-1)-(betavm^2-k^2)*c^2)*by_n_alfvmc+2*alphavm*c*by_np1_alfvmc;
        D(10,11) = 2*n*(n-1)*by_n_btvmc-2*n*betavm*c*by_np1_btvmc;
        D(10,12) = 2*1i*k*((n*(n-1)-betavm^2*c^2)*by_n_btvmc+betavm*c*by_np1_btvmc);
        D(10,13) = 0;
        D(10,14) = c^2*row2*w^2*besselh(n,1,gamma2*c)/miuvm;
        
        D(11,1) = 0; D(11,2) = 0; D(11,3) = 0; D(11,4) = 0; D(11,5) = 0; D(11,6) = 0;
        D(11,7) = -2*n*(n-1)*bj_n_alfvmc+2*n*alphavm*c*bj_np1_alfvmc;
        D(11,8) = (betavm^2*c^2-2*n*(n-1))*bj_n_btvmc-2*betavm*c*bj_np1_btvmc;
        D(11,9) = 2*1i*k*(-n*(n-1)*bj_n_btvmc+n*betavm*c*bj_np1_btvmc);
        D(11,10) = -2*n*(n-1)*by_n_alfvmc+2*n*alphavm*c*by_np1_alfvmc;
        D(11,11) = (betavm^2*c^2-2*n*(n-1))*by_n_btvmc-2*betavm*c*by_np1_btvmc;
        D(11,12) = 2*1i*k*(-n*(n-1)*by_n_btvmc+n*betavm*c*by_np1_btvmc);
        D(11,13) = 0; D(11,14) = 0;
        
        D(12,1) = 0; D(12,2) = 0; D(12,3) = 0; D(12,4) = 0; D(12,5) = 0; D(12,6) = 0;
        D(12,7) = 2*1i*k*(n*bj_n_alfvmc-alphavm*c*bj_np1_alfvmc);
        D(12,8) = 1i*k*n*bj_n_btvmc;
        D(12,9) = (betavm^2-k^2)*(n*bj_n_btvmc-betavm*c*bj_np1_btvmc);
        D(12,10) = 2*1i*k*(n*by_n_alfvmc-alphavm*c*by_np1_alfvmc);
        D(12,11) = 1i*k*n*by_n_btvmc;
        D(12,12) = (betavm^2-k^2)*(n*by_n_btvmc-betavm*c*by_np1_btvmc);
        D(11,13) = 0; D(11,14) = 0;
        
        D(13,1) = n*bj_n_alfema-alphaem*a*bj_np1_alfema;
        D(13,2) = n*bj_n_btema;
        D(13,3) = 1i*k*(n*bj_n_btema-betaem*a*bj_np1_btema);
        D(13,4) = n*by_n_alfema-alphaem*a*by_np1_alfema;
        D(13,5) = n*by_n_btema;
        D(13,6) = 1i*k*(n*by_n_btema-betaem*a*by_np1_btema);
        D(3,7) = 0; D(3,8) = 0; D(3,9) = 0; D(3,10) = 0; D(3,11) = 0; D(3,12) = 0;
        D(3,13) = -(n*besselj(n,gamma1*a)-gamma1*a*besselj(n+1,gamma1*a));
        D(3,14) = 0;
        
        D(14,1) = 0; D(14,2) = 0; D(14,4) = 0; D(14,4) = 0; D(14,5) = 0; D(14,6) = 0;
        D(14,7) = n*bj_n_alfvmc-alphavm*c*bj_np1_alfvmc;
        D(14,8) = n*bj_n_btvmc;
        D(14,9) = 1i*k*(n*bj_n_btvmc-betavm*c*bj_np1_btvmc);
        D(14,10) = n*by_n_alfvmc-alphavm*c*by_np1_alfvmc;
        D(14,11) = n*by_n_btvmc;
        D(14,12) = 1i*k*(n*by_n_btvmc-betavm*c*by_np1_btvmc);
        D(13,14) = 0;
        D(14,14) = -(n*besselh(n,1,gamma2*c)-gamma2*c*besselh(n+1,1,gamma2*c));
        
        fdisp = det(D);
end
%%
if isnan(abs(fdisp)) || isinf(abs(fdisp));
    fdisp = NaN;
    %     warning('Function Computing Failed!');
end
end