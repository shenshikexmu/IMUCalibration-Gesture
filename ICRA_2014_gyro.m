function [Tg,Kg]=ICRA_2014_gyro(rotation,a0)


if nargin<2 
   a0=[0,0,0,0,0,0,0.001,0.001,0.001];
end


% options=optimset('Algorithm','Levenberg-Marquardt',...
% 'Display','iter','TolX',1e-4,'MaxIter',10);
options=optimset('Algorithm','Levenberg-Marquardt',...
'Display','iter','TolX',1e-5,'MaxIter',18);
a=lsqnonlin(@rotation_gyro,a0,[],[],options,rotation);

Tg=[   1  , -a(1) ,  a(2) ;...
     a(3) ,   1   , -a(4) ;...
    -a(5) ,  a(6) ,    1];
Kg=[ a(7) ,  0   ,  0  ;...
     0    , a(8) ,  0  ;...
     0    ,  0   , a(9)];

end


function E=rotation_gyro(a,rotation)

Tg=[   1  , -a(1) ,  a(2) ;...
     a(3) ,   1   , -a(4) ;...
    -a(5) ,  a(6) ,    1];
Kg=[ a(7) ,  0   ,  0  ;...
     0    , a(8) ,  0  ;...
     0    ,  0   , a(9)];
 Ta=rotation{end-3};
 Ka=rotation{end-2};
 Ba=rotation{end-1};
 Bg=rotation{end};
 for i=1:size(rotation,1)-4
    data=rotation{i};
    g_start=Ta*Ka*(data(1,2:4)'+Ba);
    g_end=Ta*Ka*(data(end,2:4)'+Ba);
    Q(:,1)=[1;0;0;0];
    for j=2:size(data,1)
        gyro0=Tg*Kg*(data(j-1,5:7)'+Bg);
        gyro1=Tg*Kg*(data(j,5:7)'+Bg);
        dt=(data(j,1)-data(j-1,1));
        Q(:,j)=attitude_update_RK4(Q(:,j-1),dt,gyro0,gyro1);
        
    end
    R = quatern2rotMat(Q(:,j)');
    g_end_compute=R*g_start;
    E((i-1)*3+1,1)=g_end(1)-g_end_compute(1);
    E((i-1)*3+2,1)=g_end(2)-g_end_compute(2);
    E((i-1)*3+3,1)=g_end(3)-g_end_compute(3);
 end
%E(size(rotation,1)-3:size(rotation,1),1)=[0;0;0;0];


end


function [Qk_plus1]=attitude_update_RK4(Qk,dt,gyro0,gyro1)
% RK4
% conference: A Robust and Easy to implement method for imu
% calibration without External Equipments

q_1=Qk;
k1=(1/2)*omegaMatrix(gyro0)*q_1;
q_2=Qk+dt*(1/2)*k1;
k2=(1/2)*omegaMatrix((1/2)*(gyro0+gyro1))*q_2;
q_3=Qk+dt*(1/2)*k2;
k3=(1/2)*omegaMatrix((1/2)*(gyro0+gyro1))*q_3;
q_4=Qk+dt*k3;
k4=(1/2)*omegaMatrix(gyro1)*q_4;
Qk_plus1=Qk+dt*(k1/6+k2/3+k3/3+k4/6);
Qk_plus1=Qk_plus1/norm(Qk_plus1);

end

function [omega]=omegaMatrix(data)

% wx=data(1)*pi/180;
% wy=data(2)*pi/180;
% wz=data(3)*pi/180;
wx=data(1);
wy=data(2);
wz=data(3);

omega=[0  , -wx , -wy , -wz ;...
       wx ,  0  ,  wz , -wy ;...
       wy , -wz ,  0  ,  wx ;...
       wz ,  wy , -wx ,  0   ];

end

function R = quatern2rotMat(q)
    [rows cols] = size(q);
    R = zeros(3,3, rows);
    R(1,1,:) = 2.*q(:,1).^2-1+2.*q(:,2).^2;
    R(1,2,:) = 2.*(q(:,2).*q(:,3)+q(:,1).*q(:,4));
    R(1,3,:) = 2.*(q(:,2).*q(:,4)-q(:,1).*q(:,3));
    R(2,1,:) = 2.*(q(:,2).*q(:,3)-q(:,1).*q(:,4));
    R(2,2,:) = 2.*q(:,1).^2-1+2.*q(:,3).^2;
    R(2,3,:) = 2.*(q(:,3).*q(:,4)+q(:,1).*q(:,2));
    R(3,1,:) = 2.*(q(:,2).*q(:,4)+q(:,1).*q(:,3));
    R(3,2,:) = 2.*(q(:,3).*q(:,4)-q(:,1).*q(:,2));
    R(3,3,:) = 2.*q(:,1).^2-1+2.*q(:,4).^2;
end