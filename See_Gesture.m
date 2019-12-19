function See_Gesture( data,Ta,Ka,Ba,Tg,Kg,Bg,Tm2a,Bm,Vm,Set_Bias_Gyro)
% show the gesture of sensor, compare different algorithm
% data is raw data ;  
% Vm is mag vector in world frame
% author  Zhang Xin

m=size(data,1);

%====================================================
t(1)=0;
Pk_EKF=[];
eInt=[];
Pk_EKF_bias=[];
Bias_Ekf=zeros(1,9);
Pk_ESKF=[];
delta_X=zeros(1,6);
Bg_ESKF=zeros(1,3);

for i=1:m
    
     data(i,2:4)=(Ta*Ka*(data(i,2:4)'+Ba))';
     data(i,5:7)=(Tg*Kg*(data(i,5:7)'+Bg))';
     data(i,8:10)=(Tm2a*(data(i,8:10)'+Bm))';
    
     norm_a(i)=norm(data(i,2:4));
     norm_g(i)=norm(data(i,5:7));
     
     if i>1
         t(i)=data(i,1)-data(i-1,1);
     end
    if norm_g(i)<0.0873    %5*pi/180
        q(i,:)=[1,0,0,0];
    else
        q(i,:) = axisAngle2quatern(data(i,5:7)/norm_g(i), norm_g(i)*t(i));
    end
    if i==1 
        Q(i,:)  = accMeg2qRichard(data(i,:));  % only gyro
        
        Q_RK4(i,:)=Q(i,:);              % only gyro in using RK4 updata
        
        Q_random= randn(1,4);
        Q_random= Q_random/norm(Q_random);
        QfuseHL(i,:)=Q_random;            % high low pass filter
        QfuseEKF(i,:)=Q_random;           % EKF filter
        QfuseMahony(i,:)=Q_random;        % Mahony filter 
        QfuseEKF_bias(i,:)=Q_random;      % EKF filter to gyro bias
        QfuseESKF(i,:)=Q_random;          % ESKF filter
       
    else
        Q(i,:)=quaternProd(Q(i-1,:),q(i,:));    %Q(i-1,:)*q(i,:)
               
        Q_RK4(i,:)=attitude_update_RK4(Q_RK4(i-1,:)',t(i),data(i-1,5:7)',data(i,5:7)')';
        
        QfuseHL(i,:)=HighLowPassFilter(QfuseHL(i-1,:),data(i,:),t(i));

        [QfuseEKF(i,:),Pk_EKF]=EkfFilter(QfuseEKF(i-1,:),data(i,:),t(i),Vm,Pk_EKF);
        
        [QfuseMahony(i,:),eInt]=MahonyFilter(QfuseMahony(i-1,:),data(i,:),t(i),Vm,eInt);  
        
        data_bias=[0,0,0,0,Set_Bias_Gyro,0,0,0];
        
        [QfuseEKF_bias(i,:),Bias_Ekf(i,:),Pk_EKF_bias]=EKF_Gyro_bias(QfuseEKF_bias(i-1,:),Bias_Ekf(i-1,:),data(i,:)+data_bias,t(i),Vm,Pk_EKF_bias);
       
        [QfuseESKF(i,:),Bg_ESKF(i,:),delta_X(i,:),Pk_ESKF]=ESKF(QfuseESKF(i-1,:),Bg_ESKF(i-1,:),delta_X(i-1,:),data(i,:)+data_bias,t(i),Vm,Pk_ESKF);
        
    end
end

figure(1)
p1(1)=subplot(6,1,1);
plot(1:m,QfuseHL(:,1),'r',1:m,QfuseEKF(:,1),'g', 1:m,QfuseMahony(:,1),'b',1:m,Q(:,1),'m',1:m,Q_RK4(:,1),'c');
ylim([-1,1]);
legend('HLP','EKF','Mahony','only gyro','gyro in RK4');

p1(2)=subplot(6,1,2);
plot(1:m,QfuseHL(:,2),'r',1:m,QfuseEKF(:,2),'g', 1:m,QfuseMahony(:,2),'b',1:m,Q(:,2),'m',1:m,Q_RK4(:,2),'c');
ylim([-1,1]);

p1(3)=subplot(6,1,3);
plot(1:m,QfuseHL(:,3),'r',1:m,QfuseEKF(:,3),'g', 1:m,QfuseMahony(:,3),'b',1:m,Q(:,3),'m',1:m,Q_RK4(:,3),'c');
ylim([-1,1]);

p1(4)=subplot(6,1,4);
plot(1:m,QfuseHL(:,4),'r',1:m,QfuseEKF(:,4),'g', 1:m,QfuseMahony(:,4),'b',1:m,Q(:,4),'m',1:m,Q_RK4(:,4),'c');
ylim([-1,1]);

p1(5)=subplot(6,1,5);
plot(1:m,norm_g);

p1(6)=subplot(6,1,6);
plot(1:m,norm_a);
linkaxes(p1,'x');


figure('NumberTitle', 'off', 'Name', 'sensor attitude ');
set(gcf, 'doublebuffer', 'on');
% writerObj=VideoWriter('J:\out.avi');
% open(writerObj);
for i = 1:40:m       %show sensor attitude
   
    R0=quatern2rotMat(Q(i,:));                 % only gyro attitude
    R1=quatern2rotMat(Q_RK4(i,:));             % only gyro attitude in using RK4 updata
    R2=accMag2rotMat(data(i,:));               % acc & mag  attitude
    R3=quatern2rotMat(QfuseHL(i,:));           % high low pass filter attitude
    R4=quatern2rotMat(QfuseEKF(i,:));          % EKF filter attitude
    R5=quatern2rotMat(QfuseMahony(i,:));       % Mahony filter attitude
    R6=quatern2rotMat(QfuseEKF_bias(i,:));     % EKF for bias filter attitude
    R7=quatern2rotMat(QfuseESKF(i,:));         % ESKF filter attitude

    accW=accWorldframe(R5,data(i,:));   %acceleration in world coordinata
    
    % x                  y                  z    axis coordinate
    r0=R0(1,:);        g0=R0(2,:);       b0=R0(3,:);
    r1=R1(1,:);        g1=R1(2,:);       b1=R1(3,:);
    r2=R2(1,:);        g2=R2(2,:);       b2=R2(3,:);
    r3=R3(1,:);        g3=R3(2,:);       b3=R3(3,:);
    r4=R4(1,:);        g4=R4(2,:);       b4=R4(3,:);
    r5=R5(1,:);        g5=R5(2,:);       b5=R5(3,:);
    r6=R6(1,:);        g6=R6(2,:);       b6=R6(3,:);
    r7=R7(1,:);        g7=R7(2,:);       b7=R7(3,:);
    
    plot3([-6,r0(1)-6],[0,r0(2)],[0,r0(3)],'r',...
        [-6,g0(1)-6],[0,g0(2)],[0,g0(3)],'g',...
        [-6,b0(1)-6],[0,b0(2)],[0,b0(3)],'b',...
        ...
        [-3,r1(1)-3],[0,r1(2)],[0,r1(3)],'r',...
        [-3,g1(1)-3],[0,g1(2)],[0,g1(3)],'g',...
        [-3,b1(1)-3],[0,b1(2)],[0,b1(3)],'b',...
        ...
        [0,r2(1)],[0,r2(2)],[0,r2(3)],'r',...
        [0,g2(1)],[0,g2(2)],[0,g2(3)],'g',...
        [0,b2(1)],[0,b2(2)],[0,b2(3)],'b',...
        ...
        [3,r3(1)+3],[0,r3(2)],[0,r3(3)],'r',...
        [3,g3(1)+3],[0,g3(2)],[0,g3(3)],'g',...
        [3,b3(1)+3],[0,b3(2)],[0,b3(3)],'b',...
        ...
        [6,r4(1)+6],[0,r4(2)],[0,r4(3)],'r',...
        [6,g4(1)+6],[0,g4(2)],[0,g4(3)],'g',...
        [6,b4(1)+6],[0,b4(2)],[0,b4(3)],'b',...
        ...
        [9,r5(1)+9],[0,r5(2)],[0,r5(3)],'r',...
        [9,g5(1)+9],[0,g5(2)],[0,g5(3)],'g',...
        [9,b5(1)+9],[0,b5(2)],[0,b5(3)],'b',...
        ...
        [12,r6(1)+12],[0,r6(2)],[0,r6(3)],'r',...
        [12,g6(1)+12],[0,g6(2)],[0,g6(3)],'g',...
        [12,b6(1)+12],[0,b6(2)],[0,b6(3)],'b',...
        ...
        [15,r7(1)+15],[0,r7(2)],[0,r7(3)],'r',...
        [15,g7(1)+15],[0,g7(2)],[0,g7(3)],'g',...
        [15,b7(1)+15],[0,b7(2)],[0,b7(3)],'b',...
        ...
        [18,accW(1)+18],[0,0],[0,0],'r',...
        [18,18],[0,accW(2)],[0,0],'g',...
        [18,18],[0,0],[0,accW(3)],'b');
    axis equal
    set(gca,'XLim',[-8 21]);
    set(gca,'YLim',[-2.5 2.5]);
    set(gca,'ZLim',[-2.5 2.5]);
    
    xlabel('X');  
    ylabel('Y');  
    zlabel('Z');  
    
    title(['i=' num2str(i)]);
    text(-9,0,4,'only gyro ');
    text(-6.5,0,4,'gyro in RK4 ');
    text(-3,0,4,'acc & mag ');
    text(1,0,4,'highlow pass ');
    text(4,0,4,'EKF ');
    text(7,0,4,'Mahony ');
    text(10,0,4,'EKF_{bias}');
    text(13,0,4,'ESKF');
    text(16,0,4,'acc in world ');
    
    text(9,0,-1.8,'Bg_{set}=');
    text(9,0,-3,num2str(Set_Bias_Gyro));
    text(12,0,-1.8,'Bg_{ekf}=' );
    text(12,0,-3, num2str(roundn(Bias_Ekf(i,4:6)',-4)));
    text(15,0,-1.8,'Bg_{eskf}=');
    text(15,0,-3, num2str(roundn(Bg_ESKF(i,:)',-4)));
    drawnow
    
%     %%%%
%     frame=getframe;
%     writeVideo(writerObj,frame);
end
%close(writerObj);



% end
end


function q = axisAngle2quatern(axis, angle)
    q0 = cos(angle./2);
    q1 = axis(:,1)*sin(angle./2);
    q2 = axis(:,2)*sin(angle./2);
    q3 = axis(:,3)*sin(angle./2); 
    q = [q0 q1 q2 q3];
end

function ab = quaternProd(a, b)
    ab(1) = a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4);
    ab(2) = a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3);
    ab(3) = a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2);
    ab(4) = a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1);
    if ab(1)<0
        ab=-ab;
    end
end


function R = quatern2rotMat(q)

    R(1,1) = q(1)^2+q(2)^2-q(3)^2-q(4)^2;
    R(1,2) = 2*(q(2)*q(3)+q(1)*q(4));
    R(1,3) = 2*(q(2)*q(4)-q(1)*q(3));
    R(2,1) = 2*(q(2)*q(3)-q(1)*q(4));
    R(2,2) = q(1)^2+q(3)^2-q(2)^2-q(4)^2;
    R(2,3) = 2*(q(3)*q(4)+q(1)*q(2));
    R(3,1) = 2*(q(2)*q(4)+q(1)*q(3));
    R(3,2) = 2*(q(3)*q(4)-q(1)*q(2));
    R(3,3) = q(1)^2+q(4)^2-q(2)^2-q(3)^2;
end



function q = accMeg2qRichard(data)


vX=cross(data(1,8:10),data(1,2:4));
vX=vX/norm(vX);
vY=cross(data(1,2:4),vX);
vY=vY/norm(vY);

qX = qUtoV(vX,[1,0,0]);

y= qMultiVec(vY, qX);
qY = qUtoV(y,[0,1,0]);

qx=[-qX(1),qX(2:4)];
qy=[-qY(1),qY(2:4)];

q =qMultiQ(qx,qy);
q=[q(1),-q(2:4)];
if q(1)<0
    q=-q;
end
end


function [qq]=qMultiQ(p,q)   %p*q
qq=[...
        p(1) * q(1) - p(2) * q(2) - p(3) * q(3) - p(4) * q(4)...
       ,p(2) * q(1) + p(1) * q(2) - p(4) * q(3) + p(3) * q(4)...
       ,p(3) * q(1) + p(4) * q(2) + p(1) * q(3) - p(2) * q(4)...
       ,p(4) * q(1) - p(3) * q(2) + p(2) * q(3) + p(1) * q(4)  ];

end

function q = qUtoV(u, v)        %two vetor rotation to quaternions
nu = u/norm(u);
nv = v/norm(v);

if (u*v' == -1)
    q = [0, [1,0,0]];
else
    half = (nu + nv)/norm(nu + nv);
    q = [nu*half',cross(nu, half)];
end
end

function [vector]=qMultiVec(vec,q)  %sensor frame to world frame
x = q(2);
y = q(3);
z = q(4);
w = q(1);

vecx = vec(1);
vecy = vec(2);
vecz = vec(3);

x_ =  w * vecx  +  y * vecz  -  z * vecy;
y_ =  w * vecy  +  z * vecx  -  x * vecz;
z_ =  w * vecz  +  x * vecy  -  y * vecx;
w_ = -x * vecx  -  y * vecy  -  z * vecz;

vector = [x_ * w  +  w_ * -x  +  y_ * -z  -  z_ * -y ...
    , y_ * w  +  w_ * -y  +  z_ * -x  -  x_ * -z ...
    , z_ * w  +  w_ * -z  +  x_ * -y  -  y_ * -x ...
    ];

end



function [R]=accMag2rotMat(data)

VerticalX=cross(data(1,8:10),data(1,2:4));
VerticalX=VerticalX/norm(VerticalX);
VerticalY=cross(data(1,2:4),VerticalX);
VerticalY=VerticalY/norm(VerticalY);
VerticalZ=data(1,2:4)/norm(data(1,2:4));
R=[VerticalX',VerticalY',VerticalZ'];   

end

function [accW]=accWorldframe(R,data)

accW=(R'*data(1,2:4)'-[0;0;9.8])/9.8;


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

if Qk_plus1(1)<0
    Qk_plus1=-Qk_plus1;
end

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
