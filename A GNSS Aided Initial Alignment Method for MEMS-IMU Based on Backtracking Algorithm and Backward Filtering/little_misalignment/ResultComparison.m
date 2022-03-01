clc; clear; close all;
load alpha_LM.mat;
load alpha_LMBT.mat;
load alpha_NL.mat;
load alpha_NLBT.mat;

N = length(alpha_LM(:,:,1));
RMS_LM = zeros(N,4);
RMS_LMBT = zeros(N,4);
RMS_NL = zeros(N,4);
RMS_NLBT = zeros(N,4);
for k=1:3
    tmp = reshape(alpha_LM(:,k,:),N,30);
    RMS_LM(:,k)   = rms(tmp,2);
    tmp = reshape(alpha_LMBT(:,k,:),N,30);
    RMS_LMBT(:,k) = rms(tmp,2);
    tmp = reshape(alpha_NL(:,k,:),N,30);
    RMS_NL(:,k)   = rms(tmp,2);
    tmp = reshape(alpha_NLBT(:,k,:),N,30);
    RMS_NLBT(:,k) = rms(tmp,2);
end

errMat = [RMS_LM(:,1:3) RMS_LMBT(:,1:3) RMS_NL(:,1:3) RMS_NLBT(:,1:3) alpha_LM(:,4,1)];

myfigure('x',4);
p=plot(errMat(:,13), errMat(:,1:3:10), 'linewidth',1.5); grid on;
p(2).LineStyle='--'; p(3).LineStyle='-.'; p(4).LineStyle=':';
xlabel('\itt \rm/ s'); ylabel('\it\phi_x \rm/ (arcmin)'); set(gca,'gridlinestyle','--');
legend('LM','LMBT','NM','NMBT'); set(gca,'gridalpha',0.65);
myfigure('y',4);
p=plot(errMat(:,13), errMat(:,2:3:11), 'linewidth',1.5); grid on;
p(2).LineStyle='--'; p(3).LineStyle='-.'; p(4).LineStyle=':';
xlabel('\itt \rm/ s'); ylabel('\it\phi_y \rm/ (arcmin)'); set(gca,'gridlinestyle','--');
legend('LM','LMBT','NM','NMBT'); set(gca,'gridalpha',0.65);
myfigure('z',4);
p=plot(errMat(:,13), errMat(:,3:3:12), 'linewidth',1.5); grid on;
p(2).LineStyle='--'; p(3).LineStyle='-.'; p(4).LineStyle=':';
xlabel('\itt \rm/ s'); ylabel('\it\phi_z \rm/ (arcmin)'); set(gca,'gridlinestyle','--');
legend('LM','LMBT','NM','NMBT'); set(gca,'gridalpha',0.65);


