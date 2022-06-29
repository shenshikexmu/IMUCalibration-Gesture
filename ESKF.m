function [Q2,Bg2,delta_X2,Pk2]=ESKF(Q1,Bg1,delta_X1,ImuData,t,Vm,Pk1)

% ESKF( Error State Kalman Filter)
% derivation <Quaternion Kinematics for the Error-state Kalman Filter>
% author Zhang xin 


if isempty(Pk1)
    Para_delta_theta0=0.00001;
    Para_delta_Bg0=0.00001;
    Pk1=diag([[1,1,1]*Para_delta_theta0,[1,1,1]*Para_delta_Bg0]);
end

norm_a=norm(ImuData(1,2:4));
%norm_g=norm(ImuData(1,5:7)-Bg1);
wx=(ImuData(1,5)-Bg1(1))*t;
wy=(ImuData(1,6)-Bg1(2))*t;
wz=(ImuData(1,7)-Bg1(3))*t;

Qp=[ 1    , -wx/2 , -wy/2 , -wz/2  ;...
     wx/2 ,   1   ,  wz/2 , -wy/2  ;...
     wy/2 , -wz/2 ,   1   ,  wx/2  ;...
     wz/2 ,  wy/2 , -wx/2 ,   1   ]*Q1';
Qp=Qp/norm(Qp);
 
delta_theta=norm([wx,wy,wz]);

if delta_theta==0
    R_u_delta_theta=eye(3);
else
    u= [wx,wy,wz]'/delta_theta;
    R_u_delta_theta=eye(3)-Skew_symmetric(u)*sin(delta_theta)+Skew_symmetric(u)*Skew_symmetric(u)*(1-cos(delta_theta));
end

R_u_delta_theta=eye(3)-Skew_symmetric(u)*sin(delta_theta)+Skew_symmetric(u)*Skew_symmetric(u)*(1-cos(delta_theta));

F_delta_X=[R_u_delta_theta, -eye(3)*t;...
           zeros(3,3),       eye(3) ];
Para_delta_theta=0.000001;
Para_delta_Bg=0.000001;
Q_delta_X=diag([[1,1,1]*Para_delta_theta,[1,1,1]*Para_delta_Bg]);
       
delta_Xp = F_delta_X*delta_X1';     

P_k=F_delta_X*Pk1*F_delta_X'+Q_delta_X;

if abs(norm_a-9.8)<2 %&& norm_g< 2
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
       
    h_acc= R3;
    h_mag = Vm(2)*R2+Vm(3)*R3;
    
    H_acc = [J3,zeros(3,3)];
    H_mag = [Vm(2)*J2+Vm(3)*J3,zeros(3,3)];
    
    Q_delta_theta= 1/2*[ -Qp(2), -Qp(3), -Qp(4) ;...
                          Qp(1), -Qp(4),  Qp(3) ;...
                          Qp(4),  Qp(1), -Qp(2) ;...
                         -Qp(3),  Qp(2),  Qp(1) ];
     
    X_delta_x=[Q_delta_theta, zeros(4,3);...
                zeros(3,3),      eye(3)];
    
    Hk= [H_acc;H_mag]* X_delta_x;     
    
    para_Rk_acc=0.1;
    para_Rk_mag=0.2;
    Rk=diag([[1,1,1]*para_Rk_acc,[1,1,1]*para_Rk_mag]);
    
    Kk=P_k*Hk'*inv(Hk*P_k*Hk'+Rk);
    
    delta_X_hat=Kk*([ImuData(1,2:4)'/norm_a;ImuData(1,8:10)'/norm(ImuData(1,8:10))]-[h_acc;h_mag]);
    
    Q2=quaternProd(Qp, [1;delta_X_hat(1:3)/2])';
    Q2=Q2/norm(Q2);
    Bg2=Bg1+delta_X_hat(4:6)';
    
    Pk2_=(eye(6)-Kk*Hk)*P_k;
    
    delta_X2=delta_Xp'+delta_X_hat';
    
    G=[eye(3)-Skew_symmetric(-delta_X2(1:3)/2),zeros(3,3);...
        zeros(3,3)   ,                 eye(3)];
    
    Pk2=G*Pk2_*G';
    
else
    Q2=Qp';
    Bg2=Bg1;
    delta_X2=delta_Xp';
    Pk2=P_k;
    
end
    delta_X2(1:3)=[0,0,0];
    
    if Q2(1)<0
        Q2=-Q2;
    end


end


function S=Skew_symmetric(u)
%Skew_Operator

S=[  0  , -u(3) ,  u(2) ;...
    u(3),   0   , -u(1) ;...
   -u(2),  u(1),    0 ];

end

function ab = quaternProd(a, b)
    ab(1,1) = a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4);
    ab(2,1) = a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3);
    ab(3,1) = a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2);
    ab(4,1) = a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1);
    if ab(1)<0
        ab=-ab;
    end
end


