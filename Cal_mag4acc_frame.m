function [Tm2a,Bm,Vm,mag_strength]=Cal_mag4acc_frame(rotation,fix_point,Tg,Kg,a0)
% author  Zhang Xin

if nargin==4
   a0=[1,1,1,1,1,1,1,1,0,0,0,0];
end
n=size(rotation,1);
rotation{n+1}=Tg;
rotation{n+2}=Kg;

%  options=optimset('TolX',1e-6,'Algorithm','Levenberg-Marquardt',...
%  'Display','iter');
% [a]=lsqnonlin(@mag_in_diff_gesture,a0,[],[],options,rotation);

fprintf('Calibration Magnetometer:\n');
TolX=1e-6;
TolFun=0;
MaxIter=inf;
a=Optimize_my_LM(@mag_in_diff_gesture,a0,rotation,TolX,TolFun,MaxIter);

Tm2a=[a(1)   , a(2),  a(3);...
    a(4) ,  a(5)   , a(6);...
    a(7) ,  a(8),   -1];

Bm=[a(9);a(10);a(11)];

mag_strength=a(12);

Ba=rotation{end-3};
Ka=rotation{end-4};
Ta=rotation{end-5};
for i=1:size(fix_point,1)
    
    Acc(:,i)=Ta*Ka*(fix_point(i,1:3)'+Ba);
    Mag(:,i)=Tm2a*(fix_point(i,7:9)'+Bm);     
    SS(i)=Acc(:,i)'*Mag(:,i)/(norm(Acc(:,i))*norm(Mag(:,i)));
end

Mag_z=mean(SS);
Mag_y=sqrt(1-Mag_z^2);
Vm=[0,Mag_y,Mag_z];

end


function E=mag_in_diff_gesture(a,rotation)

Tm2a=[a(1)   , a(2),  a(3);...
    a(4) ,  a(5)   , a(6);...
    a(7) ,  a(8),   -1];

Bm=[a(9);a(10);a(11)];

mag_strength=a(12);

Kg=rotation{end};
Tg=rotation{end-1};
Bg=rotation{end-2};
% Ba=rotation{end-3};
% Ka=rotation{end-4};
% Ta=rotation{end-5};

for i=1:size(rotation,1)-6
    data=rotation{i};
    mag_start=Tm2a*(data(1,8:10)'+Bm);
    mag_end=Tm2a*(data(end,8:10)'+Bm);
    Q(:,1)=[1;0;0;0]; 
    for j=2:size(data,1)
        gyro0=Tg*Kg*(data(j-1,5:7)'+Bg);
        gyro1=Tg*Kg*(data(j,5:7)'+Bg);
        dt=(data(j,1)-data(j-1,1));
        Q(:,j)=attitude_update_RK4(Q(:,j-1),dt,gyro0,gyro1); 
    end
    R = quatern2rotMat(Q(:,j)');
    mag_end_compute=R*mag_start;
    E((i-1)*4+1,1)=mag_end(1)-mag_end_compute(1);
    E((i-1)*4+2,1)=mag_end(2)-mag_end_compute(2);
    E((i-1)*4+3,1)=mag_end(3)-mag_end_compute(3);
    E((i-1)*4+4,1)=mag_strength-norm(mag_start);

end

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
