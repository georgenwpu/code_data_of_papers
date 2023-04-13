%%==============================
% 不同半锥角下的圆锥运动
%%==============================
clc; clear; close all;
glvs;
afa = 3*10.^(-2:1)*glv.deg;     % half-apex angle
f = 2;  w = 2*pi*f;             % frequency

ss = 4;
[nn, ts, nts] = nnts(ss, 0.01);  % subsamples & sampling interval
T = 2;                           % simulation time
len = T/ts; res = zeros(fix(len/nn), 3*10+1);

N = length(afa);
wq = cell(N,1);
for k=1:N
    [wm, qr] = conesimu(afa(k), f, ts, T);
    wq{k}.wm = wm;
    wq{k}.qr = qr;
end

%% calculation of Maneuver Intensity Parameter
len = length(wm);
sigma_d1 = zeros(len,2);
sigma_d2 = zeros(len,2);
alpha_t = zeros(len,2);
jj = 1;
for k=1:N
    wm = wq{k}.wm;
    for ii = ss:ss:length(wm)
        sigma_ii = norm(cross(wm(ii-ss+1, :), wm(ii-ss+2,:))) + ...
                   norm(cross(wm(ii-ss+1, :), wm(ii-ss+3,:))) + ...
                   norm(cross(wm(ii-ss+1, :), wm(ii-ss+4,:)));
        sigma_d1(jj,:) = [sigma_ii, jj*nts];
        sigma_ii = norm(cross(wm(ii-ss+1, :), wm(ii-ss+4,:)));
        sigma_d2(jj,:) = [sigma_ii, jj*nts];
        alpha_t(jj,:) = [afa(k), jj*nts];
        jj =  jj+1;
    end
end

figure;  colororder({glv.Lcolor{1}, glv.Lcolor{2}}); mki = 10:10:200;
yyaxis left;
p(1)=xlogyEx(sigma_d1(:,2), sigma_d1(:,1),'-o', 1,mki); hold on;
p(2)=xlogyEx(sigma_d2(:,2), sigma_d2(:,1),'-^', 1,mki);
xylabelsEx('time / s','MIP / (rad^2/s^2)');
yyaxis right;
p(3)=xlogyEx(alpha_t(:,2), alpha_t(:,1)./glv.deg,'-*', 2,mki); set(gca, 'yscale','log');
xylabelsEx('time / s','Half-apex angle / (deg)');
l=legend(p, '\sigma_1','\sigma_2','\alpha');
l.Orientation='horizontal'; l.NumColumns=1;
l.FontSize=12; l.Location='best';

%% Calculation of the NCE
global LPBCoe LPICoe LPDCoe;
[LPBCoe, LPICoe, LPDCoe] = CalLP_OID;   % 生成Legender多项式函数
kudt = 1;       % 更新次数
res_q = zeros(floor(len),4*7+1);       % 姿态更新结果
err_phi = zeros(floor(len),3*7+1);     % 姿态误差
for condition = 1:N
    wm = wq{condition}.wm;
    qr = wq{condition}.qr;
    q = repmat(qr(1,:)',[1,7]);       % 姿态更新四元数
    
    for k=ss:ss:len
        k0 = k-ss+1;  k1 = k;
        wmk = wm(k0:k1,1:3);            % 角增量输入
        qrk = qr(k1+1,1:4);             % 参考四元数
        dq = GetdQ_QTS(wmk, ts, 2);     % 4次展开
        q(:,1) = qmul(q(:,1), dq);
        dq = GetdQ_QTS(wmk, ts, 4);     % 6次展开
        q(:,2) = qmul(q(:,2), dq);
        dq = GetdQ_QTS(wmk, ts, 6);     % 8次展开
        q(:,3) = qmul(q(:,3), dq);
        dq = GetdQ_QTS(wmk, ts, 8);    % 10次展开
        q(:,4) = qmul(q(:,4), dq);
        dq = GetdQ_QTS(wmk, ts, 10);    % 10次展开
        q(:,5) = qmul(q(:,5), dq);
        res_q(kudt,:) = [reshape(q,1,28) kudt*nts];
        for kArg = 1:5
            err_phi(kudt,kArg*3-2:kArg*3) = abs(qq2phi(q(:,kArg),qrk));
        end
        err_phi(kudt,end) = kudt*nts;
        kudt = kudt+1;
    end
end

figure;
xlogyEx(err_phi(:,end), err_phi(:,1:3:13)./glv.sec);
xylabelsEx('time / s','attitude error / (^\prime^\prime)');
legendEx({'QHN-2','QHN-4','QHN-6','QHN-8','QHN-10'});

opt_res = [err_phi(1:50,1); err_phi(51:100,4); err_phi(101:150,7); ...
           err_phi(151:200,10)];
figure;
xlogyEx(err_phi(:,end), err_phi(:,1:3:13)./glv.sec); hold on;
l = xlogyEx(err_phi(:,end), opt_res./glv.sec, '-.'); l.LineWidth = 5;
xylabelsEx('time / s','attitude error / (^\prime^\prime)');
legendEx({'QHN-2','QHN-4','QHN-6','QHN-8','QHN-10','QHN-MIA'});

