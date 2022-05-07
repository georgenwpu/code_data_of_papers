function [ Xk ] = f_ukf( Xk_1, para )
%f_state 非线性状态更新函数
    imu = para(1:6); avp = para(7:15); ts = para(16);
    ak_1 = Xk_1(1:3);
    dvk_1 = Xk_1(4:6);
    ebk_1 = Xk_1(7:9);
    dbk_1 = Xk_1(10:12);
    Cw = a2Cw(ak_1);
    Cn_n = a2Cn_n(ak_1);
    eth = earth(avp(7:9), avp(4:6));
    wnin = eth.wnin;
    Cnb = a2mat(avp(1:3));
    fsf = imu(4:6)./ts;
    
    dak = (Cw^-1)*((eye(3)-Cn_n)*wnin - Cnb*ebk_1);
    ak = ak_1 + dak*ts;
    ddvk = (eye(3)-Cn_n')*Cnb*fsf + Cn_n'*Cnb*dbk_1;
    dvk = dvk_1 + ddvk*ts;
    ebk = ebk_1;
    dbk = dbk_1;

    Xk = [ak; dvk; ebk; dbk];
    
end

