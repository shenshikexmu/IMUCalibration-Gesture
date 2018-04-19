function See_Gesture( data,Ta,Ka,Ba,Tg,Kg,Bg,Tm2a,Bm,Vm)
% show the gesture of sensor, compare different algorithm
% data is raw data ;  
% Vm is mag vector in world frame


m=size(data,1);

%====================================================
t(1)=0;

for i=1:m
     data(i,2:4)=(data(i,2:4)-A')./APP;
     data(i,5:7)=(data(i,5:7)-G)./[16.4,16.4,16.4];
     data(i,8:10)=(data(i,8:10)-M');
     
     norm_g(i)=norm(data(i,5:7));
     norm_a(i)=norm(data(i,2:4));
     if i>1
         t(i)=data(i,1)-data(i-1,1);
         if abs(t(i))>1000
            t(i)=10;
         end
     end
    if norm_g(i)<5
        q(i,1)=1;
        q(i,2:4)=0;
    else
        P=norm_g(i);
        q(i,1:4) = axisAngle2quatern([data(i,5)/P,data(i,6)/P,data(i,7)/P]...
            ,  (P*pi*t(i))/(1000*180));
        
    end
    if i==1
        
        Q(i,:)  = accMeg2qRichard(data(i,:),Mag_R);
        
        QfuseHL(i,:)=Q(i,:);           %high low pass filter
        QfuseEKF(i,:)=Q(i,:);          % EKF filter
        QfuseMahony(i,:)=Q(i,:);          % Mahony filter
 %       QfuseMahony1D(i,:)=Q(i,:); 
    else
        Q(i,:)=quaternProd(Q(i-1,:),q(i,:));    %Q(i-1,:)*q(i,:)
        
  %      tic;
        QfuseHL(i,:)=HighLowPassFilter(QfuseHL(i-1,:),data(i,:),t(i),Mag_R);
 %       t1=t1+toc;
        if i==2
            
            [QfuseEKF(i,:),Pk]=EkfFilter2(QfuseEKF(i-1,:),data(i,:),t(i),Mag_R,Mag_V);
            [QfuseMahony(i,:),eInt]=MahonyFilter2(QfuseMahony(i-1,:),data(i,:),t(i),Mag_R,Mag_V);
%             [QfuseEKF(i,:),Pk]=EkfFilter2(QfuseEKF(i-1,:),data(i,:),t(i),Mag_R);
%             [QfuseMahony(i,:),eInt]=MahonyFilter2(QfuseMahony(i-1,:),data(i,:),t(i),Mag_R);
        else
        %    tic;
            [QfuseEKF(i,:),Pk]=EkfFilter2(QfuseEKF(i-1,:),data(i,:),t(i),Mag_R,Mag_V,Pk);
           % [QfuseEKF(i,:),Pk]=EkfFilter2(QfuseEKF(i-1,:),data(i,:),t(i),Mag_R,Pk);
%             t2=t2+toc;
%             tic
            [QfuseMahony(i,:),eInt]=MahonyFilter2(QfuseMahony(i-1,:),data(i,:),t(i),Mag_R,Mag_V,eInt);
        %    [QfuseMahony(i,:),eInt]=MahonyFilter2(QfuseMahony(i-1,:),data(i,:),t(i),Mag_R,eInt);
%             t3=t3+toc;          
        end

    end
 
end
figure(3)
p1(1)=subplot(5,1,1);
plot(1:m,QfuseHL(:,1),'r',1:m,QfuseEKF(:,1),'g', 1:m,QfuseMahony(:,1),'b',1:m,Q(:,1),'m');
ylim([-1,1]);
legend('HLP','EKF','Mahony','original');

p1(2)=subplot(5,1,2);
plot(1:m,QfuseHL(:,2),'r',1:m,QfuseEKF(:,2),'g', 1:m,QfuseMahony(:,2),'b',1:m,Q(:,2),'m');
ylim([-1,1]);

p1(3)=subplot(5,1,3);
plot(1:m,QfuseHL(:,3),'r',1:m,QfuseEKF(:,3),'g', 1:m,QfuseMahony(:,3),'b',1:m,Q(:,3),'m');
ylim([-1,1]);

p1(4)=subplot(5,1,4);
plot(1:m,QfuseHL(:,4),'r',1:m,QfuseEKF(:,4),'g', 1:m,QfuseMahony(:,4),'b',1:m,Q(:,4),'m');
ylim([-1,1]);

p1(5)=subplot(5,1,5);
plot(1:m,norm_g);
linkaxes(p1,'x');

% figure(4)
% plot(1:m,norm_a);


figure('NumberTitle', 'off', 'Name', 'sensor attitude ');
set(gcf, 'doublebuffer', 'on');
% writerObj=VideoWriter('J:\out.avi');
% open(writerObj);
for i = 1:20:m       %show sensor attitude

    R=quatern2rotMat(Q(i,:));             % gyro attitude
    RR=accMag2rotMat(data(i,:),Mag_R);    % acc & mag  attitude

    RRR=quatern2rotMat(QfuseHL(i,:));     % high low pass filter attitude
    RRRR=quatern2rotMat(QfuseEKF(i,:));   % EKF filter attitude
    RRRRR=quatern2rotMat(QfuseMahony(i,:));   % Mahony filter attitude

   accW=accWorldframe(RRRRR,data(i,:));
   
   
    r1=R(1,:);                            % x axis coordinate
    g1=R(2,:);
    b1=R(3,:);
    r2=RR(1,:);        
    g2=RR(2,:);
    b2=RR(3,:);
    r3=RRR(1,:);        
    g3=RRR(2,:);
    b3=RRR(3,:);
    r4=RRRR(1,:);        
    g4=RRRR(2,:);
    b4=RRRR(3,:);
    r5=RRRRR(1,:);        
    g5=RRRRR(2,:);
    b5=RRRRR(3,:);
%     r6=RRRRRR(1,:);        
%     g6=RRRRRR(2,:);
%     b6=RRRRRR(3,:);

    
    plot3([-2,r1(1)-2],[0,r1(2)],[0,r1(3)],'r',...
        [-2,g1(1)-2],[0,g1(2)],[0,g1(3)],'g',...
        [-2,b1(1)-2],[0,b1(2)],[0,b1(3)],'b',...
        ...
        [0,r2(1)],[0,r2(2)],[0,r2(3)],'r',...
        [0,g2(1)],[0,g2(2)],[0,g2(3)],'g',...
        [0,b2(1)],[0,b2(2)],[0,b2(3)],'b',...
        ...
        [2,r3(1)+2],[0,r3(2)],[0,r3(3)],'r',...
        [2,g3(1)+2],[0,g3(2)],[0,g3(3)],'g',...
        [2,b3(1)+2],[0,b3(2)],[0,b3(3)],'b',...
        ...
        [4,r4(1)+4],[0,r4(2)],[0,r4(3)],'r',...
        [4,g4(1)+4],[0,g4(2)],[0,g4(3)],'g',...
        [4,b4(1)+4],[0,b4(2)],[0,b4(3)],'b',...
        ...
        [6,r5(1)+6],[0,r5(2)],[0,r5(3)],'r',...
        [6,g5(1)+6],[0,g5(2)],[0,g5(3)],'g',...
        [6,b5(1)+6],[0,b5(2)],[0,b5(3)],'b',...
        ...
        [10,accW(1)+10],[0,0],[0,0],'r',...
        [10,10],[0,accW(2)],[0,0],'g',...
        [10,10],[0,0],[0,accW(3)],'b');
    axis equal
    set(gca,'XLim',[-3.5 11.5]);
    set(gca,'YLim',[-1.5 1.5]);
    set(gca,'ZLim',[-1.5 1.5]);
    
    xlabel('X');  
    ylabel('Y');  
    zlabel('Z');  
    title(['i=' num2str(i)]) ;
 
    
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



function q = accMeg2qRichard(data,Mag_R)


vX=cross(data(1,8:10)*Mag_R',data(1,2:4));
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
end
    half = (nu + nv)/norm(nu + nv);
    q = [nu*half',cross(nu, half)];
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



function [R]=accMag2rotMat(data,Mag_R)

VerticalX=cross(data(1,8:10)*Mag_R',data(1,2:4));
VerticalX=VerticalX/norm(VerticalX);
VerticalY=cross(data(1,2:4),VerticalX);
VerticalY=VerticalY/norm(VerticalY);
VerticalZ=data(1,2:4)/norm(data(1,2:4));
R=[VerticalX',VerticalY',VerticalZ'];   %��ͨ�Ӽ��ٶ���ת����
% VerticalZ'����ֱ����
end

function [accW]=accWorldframe(R,data)
%��ʾ���������ϵ�µĴ������ܵļ��ٶ\EF\BF?
accW=R'*data(1,2:4)'-[0;0;1];


end

