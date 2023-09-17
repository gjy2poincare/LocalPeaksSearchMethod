function [cp_3, ki_3] = RootsCheck(f, MPM, cp_0, ki_0) 
%根值校验【Roots Check】
w = 2*pi*f; N = length(cp_0);
dcp = 1; dki = 1; scale = 10; cp_3 = cp_0; cp_3(:,:) = NaN; ki_3 = cp_3;
%%  
global szChkCho
szChkCho = questdlg('Begin to check roots?', 'Roots Check', 'Yes', 'No', 'Never', 'Yes'); pause(0.1);
if ~strcmp(szChkCho, 'Yes'), pause(0.1); return; end
disp('----------------------------------------');
disp('【Peaks Manual Check】');
ii = 1;
for n = 1:N
    cp0 = cp_0(n); ki0 = ki_0(n);    
    [CP,KI] = meshgrid(max(cp0-scale*dcp, 1):dcp:cp0+scale*dcp,...
        max(ki0-scale*dki, -1):dki:ki0+scale*dki);
    [nr,nc] = size(CP); DISP = zeros(nr,nc); DISP(:,:) = NaN;
    for ii = 1:nr
        for jj = 1:nc
            k = w/CP(ii,jj)+1i*KI(ii,jj);
            DISP(ii,jj) = DispFunction(f, MPM, k);
        end
    end
    FD = 1./(1+abs(DISP));
    figure(5784137); surf(CP,KI,20*log10(FD)); title('根值检验【Roots Check】');
    %人工校验【Manual Check】
    choice = questdlg('Is this a root?--(Type "Ctrl+C" to change current figure)'); pause(0.1);
    switch choice
        case 'Yes'
            cp_3(ii) = cp0; ki_3(ii) = ki0; ii=ii+1;
            disp(['>>cp = ', num2str(cp0), '; ki = ', num2str(ki0)]);
        case 'No'
            continue;
        otherwise
            close(5784137); pause(0.1); cp_3 = cp_0; ki_3 = ki_0; return;
    end    
end
close(5784137); pause(0.1);
end