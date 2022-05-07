% Unscented Kalman filter(UKF) simulation with large misalignment angles.
% See also  test_align_some_methods, test_align_ekf.
% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 05/09/2013
clear;
glvs;
load('trj_align_600s.mat');
psinstypedef('test_align_ukf_def');
[nn, ts, nts] = nnts(4, trj.ts);
avp0 = trj.avp0; qnb0 = a2qua(avp0(1:3)');
eth = earth(avp0(7:9), avp0(4:6));
imuerr = imuerrset(100, 1000, 0.2, 20);
% imu = imuadderr(trj.imu, imuerr);
imu = trj.imu;
imup = imu;
imur = flip(imup, 1);
imur(:,1:3) = -imur(:,1:3);
vpp = trj.avp(:,4:9);
vpr = flip(vpp, 1);
vpr(:,1:3) = -vpr(:,1:3);
afa = [85; 85; 135]*glv.deg;    % large misalignment angles
qnb = qaddafa(qnb0,afa);
avp0(1:3)=q2att(qnb);
kf = kfinit(nts, imuerr); kf.s = 1.01; % forgetting factor
kf.coef_fb = 0.2;
ins = insinit(avp0, ts);
len = length(imup); [res, xkpk] = prealloc(fix(len/nn), 7, 2*kf.n+6+1);

it = 2*1+1;
timebar(nn, len*it, 'UKF align simulation.');
dir = 0;

for i = 1:it
    dir = mod(dir+1,2);         % 1 正向  0 逆向
    
    if dir == 1
        imu = imup;             % 正向数据
        vpgnss = vpp;
        sign_wie = 1;
    else
        imu = imur;             % 逆向数据
        vpgnss = vpr;
        sign_wie = -1;
    end
    if i>1
        ins.vn = -ins.vn;        % 速度取反
        kf.xk([4:6,7:9]) = -kf.xk([4:6,7:9]);
    end
    
    ki=1;
    for k=1:nn:len-nn+1
        k1 = k+nn-1;
        wvm = imu(k:k1,1:6);  t = imu(k1,7);
        [phim, dvbm] = cnscl(wvm);
        phim = phim-ins.eb*nts; dvbm = dvbm-ins.db*nts;  % calibration
        dvn = qmulv(ins.qnb,dvbm); ins.vn = ins.vn + dvn + eth.gn*nts;
        ins.qnb = qupdt(ins.qnb, phim-qmulv(qconj(ins.qnb),sign_wie.*eth.wnie*nts));
        res(ki,:) = [q2att(ins.qnb); ins.vn; t]';
        if mod(k1,100)~=0
            kf.px = [phim; dvbm; [q2att(ins.qnb); ins.vn; avp0(7:9)]; nts];
            kf = ukf(kf);
        else
            kf = ukf(kf, ins.vn-vpgnss(k1,1:3)');  % UKF filter(包含一次时间更新)
            % Feedback
            xfb_tmp = kf.coef_fb.*kf.xk;
            ins.qnb = qdelphi(ins.qnb, xfb_tmp(1:3));
            ins.vn = ins.vn - xfb_tmp(4:6);
            kf.xk(1:6) = kf.xk(1:6) - xfb_tmp(1:6);
            ins.eb = ins.eb + xfb_tmp(7:9);
            ins.db = ins.db + xfb_tmp(10:12);
            kf.xk(7:12) = kf.xk(7:12) - xfb_tmp(7:12);
        end
        xkpk(ki,:) = [kf.xk; diag(kf.Pxk); ins.eb; ins.db; t]';
        ki = ki+1;
        timebar;
    end
%     ins.qnb = qdelphi(ins.qnb, kf.xk(1:3));
%     ins.vn = ins.vn - kf.xk(4:6);
%     kf.xk(1:6) = zeros(6,1);
%     kf.xk(1:6) = kf.xk(1:6) - xfb_tmp(1:6);
end

% 


