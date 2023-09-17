function cr = RayleighWave(cl, ct)
%Rayleigh Wave
% clear; clc;
% cl = 5.85E3; %×Ý²¨²¨ËÙ£ºm/s
% ct = 3.23E3; %ºá²¨²¨ËÙ£ºm/s
%%
% if abs(imag(cl))<1 && abs(imag(ct))<1
cl = abs(cl); ct = abs(ct);
m = (ct/cl)^2;
b = -8; c = 8*(3-2*m); d = 16*(m-1);
p = c-b^2/3; q = d-b*c/3+2*b^3/27;
u = -q/2+sqrt(p^3/27+q^2/4); absu = abs(u); argu = angle(u);
v = -q/2-sqrt(p^3/27+q^2/4); absv = abs(v); argv = angle(v);
Am = (absu)^(1/3)*exp(1i*[argu argu+2*pi argu+4*pi]/3); Bm = (absv)^(1/3)*exp(1i*[argv argv+2*pi argv+4*pi]/3); omg = -1/2+1i*sqrt(3)/2;
for ii = 1:3
    for jj = 1:3
        if abs(Am(ii)*Bm(jj)+p/3)<=1E-6
            A = Am(ii); B = Bm(jj);
        end
    end
end
x1 = A+B-b/3; x2 = omg*A+omg^2*B-b/3; x3 = omg^2*A+omg*B-b/3;
x = [x1 x2 x3];
delta = x(real(x)<1);
cr =  real(ct*sqrt(delta));
% else
%%
%     m = (ct/cl)^2;
%     b = -8; c = 8*(3-2*m); d = 16*(m-1);
%     fun = inline(solve('x^3 + b*x^2 + c*x + d'));
%     x = fun(b,c,d);
%     cr = ct*sqrt(x);
end 