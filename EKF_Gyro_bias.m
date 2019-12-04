function [Q2,Bais2,Pk2]=EKF_Gyro_bias(Q1,Bais1,ImuData,t,Vm,Pk1)
%EKF filter to Gyro atitude with Accelerate & Magnetic
% the Bias of the Acc、Gyro、Mag are set as State.
% but the Bias of Acc and Mag are castrated.
% author Zhang Xin


para_Pk_Q=0.00001;
para_Pk_Ba=0;             % castrate Bias of Acc
para_PK_Bg=0.00001;
para_PK_Bm=0;             % castrate Bias of Mag


if isempty(Pk1)  
    Pk1=diag([[1,1,1,1]*para_Pk_Q,[1,1,1]*para_Pk_Ba,[1,1,1]*para_PK_Bg,[1,1,1]*para_PK_Bm]);
    
end

para_Qk_Q=0.000001;
para_Qk_Ba=0;
para_Qk_Bg=0.000001;
para_Qk_Bm=0;

Qk=diag([[1,1,1,1]*para_Qk_Q,[1,1,1]*para_Qk_Ba,[1,1,1]*para_Qk_Bg,[1,1,1]*para_Qk_Bm]);

wx=ImuData(1,5)*t;
wy=ImuData(1,6)*t;
wz=ImuData(1,7)*t;

Ak=[ 1    , -wx/2 , -wy/2 , -wz/2 , 0, 0, 0,  Q1(2)*t/2,  Q1(3)*t/2,  Q1(4)*t/2, 0, 0, 0 ;...
     wx/2 ,   1   ,  wz/2 , -wy/2 , 0, 0, 0, -Q1(1)*t/2,  Q1(4)*t/2, -Q1(3)*t/2, 0, 0, 0 ;...
     wy/2 , -wz/2 ,   1   ,  wx/2 , 0, 0, 0, -Q1(4)*t/2, -Q1(1)*t/2,  Q1(2)*t/2, 0, 0, 0 ;...
     wz/2 ,  wy/2 , -wx/2 ,   1   , 0, 0, 0,  Q1(3)*t/2, -Q1(2)*t/2, -Q1(1)*t/2, 0, 0, 0 ;...
     0    ,    0  ,   0   ,   0   , 1, 0, 0,      0  ,        0   ,         0  , 0, 0, 0 ;... 
     0    ,    0  ,   0   ,   0   , 0, 1, 0,      0  ,        0   ,         0  , 0, 0, 0 ;...
     0    ,    0  ,   0   ,   0   , 0, 0, 1,      0  ,        0   ,         0  , 0, 0, 0 ;...
     0    ,    0  ,   0   ,   0   , 0, 0, 0,      1  ,        0   ,         0  , 0, 0, 0 ;... 
     0    ,    0  ,   0   ,   0   , 0, 0, 0,      0  ,        1   ,         0  , 0, 0, 0 ;...
     0    ,    0  ,   0   ,   0   , 0, 0, 0,      0  ,        0   ,         1  , 0, 0, 0 ;...
     0    ,    0  ,   0   ,   0   , 0, 0, 0,      0  ,        0   ,         0  , 1, 0, 0 ;... 
     0    ,    0  ,   0   ,   0   , 0, 0, 0,      0  ,        0   ,         0  , 0, 1, 0 ;...
     0    ,    0  ,   0   ,   0   , 0, 0, 0,      0  ,        0   ,         0  , 0, 0, 1  ];


S_temp=Ak*[Q1';Bais1'];

Qp=S_temp(1:4)/norm(S_temp(1:4));
Biasp=S_temp(5:13);

P_k=Ak*Pk1*Ak'+Qk;

norm_a=norm(ImuData(1,2:4)'-Biasp(1:3));
norm_g=norm(ImuData(1,5:7)'-Biasp(4:6));


if abs(norm_a-9.8)<2 && norm_g< 2
    
    R2=[2*Qp(2)*Qp(3)+2*Qp(1)*Qp(4);...
        2*Qp(1)^2+2*Qp(3)^2-1;...
        2*Qp(3)*Qp(4)-2*Qp(1)*Qp(2)];

    R3=[2*Qp(2)*Qp(4)-2*Qp(1)*Qp(3);...
        2*Qp(3)*Qp(4)+2*Qp(1)*Qp(2);...
        2*Qp(1)^2+2*Qp(4)^2-1 ];   

    J2=2*[ Qp(4), Qp(3), Qp(2) , Qp(1);...
           Qp(1),-Qp(2), Qp(3) ,-Qp(4);...
          -Qp(2),-Qp(1), Qp(4) , Qp(3)] ;   

    J3=2*[-Qp(3), Qp(4),-Qp(1) , Qp(2);...
           Qp(2), Qp(1), Qp(4) , Qp(3);...
           Qp(1),-Qp(2),-Qp(3) , Qp(4)] ;   
         
     d_Acc_d_Q=J2;
     
     d_Mag_d_Q=Vm(2)*J2+Vm(3)*J3;
    
     d_Q_d_Bg=-[-Q1(2)*t/2, -Q1(3)*t/2, -Q1(4)*t/2 ;...
                 Q1(1)*t/2, -Q1(4)*t/2,  Q1(3)*t/2 ;...
                 Q1(4)*t/2,  Q1(1)*t/2, -Q1(2)*t/2 ;...
                -Q1(3)*t/2,  Q1(2)*t/2,  Q1(1)*t/2 ];
           
     Acc_p=R3;
     Mag_p=Vm(2)*R2+Vm(3)*R3; 
    
     Acc_sensor=(ImuData(1,2:4)'-Biasp(1:3))/norm_a;
     norm_m=norm((ImuData(1,8:10)'-Biasp(7:9)));
     Mag_sensor=(ImuData(1,8:10)'-Biasp(7:9))/norm_m;
     
     d_Acc_d_Ba=-eye(3)/norm_a;
     
     d_Acc_d_Bm=-zeros(3,3);
     
     d_Mag_d_Bm=-eye(3)/norm_m;
     
     d_Mag_d_Ba=-zeros(3,3);
    
    Hk=[d_Acc_d_Q,d_Acc_d_Ba,d_Acc_d_Q*d_Q_d_Bg,d_Acc_d_Bm;...
        d_Mag_d_Q,d_Mag_d_Ba,d_Mag_d_Q*d_Q_d_Bg,d_Mag_d_Bm];
    
    para_Rk_acc=0.1;
    para_Rk_mag=0.2;
    Rk=diag([[1,1,1]*para_Rk_acc,[1,1,1]*para_Rk_mag]);
    
    Kk=P_k*Hk'*inv(Hk*P_k*Hk'+Rk);
    
    S_temp2=[Qp;Biasp]+ Kk*[Acc_sensor-Acc_p;Mag_sensor-Mag_p];
    
    Q2=S_temp2(1:4)/norm(S_temp2(1:4));
    
    Bais2=S_temp2(5:13);
    
   Pk2=(eye(13)-Kk*Hk)*P_k;
    
else
    
    Q2=Qp;
    Bais2=Biasp;
    Pk2=P_k;
    
end


if Q2(1)<0
    Q2=-Q2;
end
Q2=Q2';
Bais2=Bais2';

end
