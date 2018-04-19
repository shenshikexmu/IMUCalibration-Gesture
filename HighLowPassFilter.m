function [Qfuse2]=HighLowPassFilter(Qfuse1,ImuData,t)
%high and low pass filter to Gyro atitude with Accelerate & Magnetic
% author  Zhang Xin

a=0.1;
norm_g=norm(ImuData(1,5:7));
norm_a=norm(ImuData(1,2:4));
if norm_g<0.0873    %5*pi/180    
    q=[1,0,0,0];
else
    q=axisAngle2quatern(ImuData(5:7)/norm_g,  norm_g*t); 
    
end
if abs(norm_a-9.8)<2 && norm_g< 2
    
    Qtemp1= accMeg2qRichard(ImuData);            
    Qtemp2=quaternProd(Qfuse1,q);        
    if Qtemp1*Qtemp2'>0       
        Qfuse2=(1-a)*Qtemp2+a*Qtemp1;       
    else      
        Qfuse2=(1-a)*Qtemp2-a*Qtemp1;          
    end    
    Qfuse2=Qfuse2/norm(Qfuse2);      
else   
    Qfuse2=quaternProd(Qfuse1,q); 
end
if Qfuse2(1)<0
    Qfuse2=-Qfuse2;
end
end


function q = axisAngle2quatern(axis, angle)
    q0 = cos(angle./2);
    q1 = axis(:,1)*sin(angle./2);
    q2 = axis(:,2)*sin(angle./2);
    q3 = axis(:,3)*sin(angle./2); 
    q = [q0 q1 q2 q3];
end

function ab = quaternProd(a, b)
    ab(:,1) = a(:,1).*b(:,1)-a(:,2).*b(:,2)-a(:,3).*b(:,3)-a(:,4).*b(:,4);
    ab(:,2) = a(:,1).*b(:,2)+a(:,2).*b(:,1)+a(:,3).*b(:,4)-a(:,4).*b(:,3);
    ab(:,3) = a(:,1).*b(:,3)-a(:,2).*b(:,4)+a(:,3).*b(:,1)+a(:,4).*b(:,2);
    ab(:,4) = a(:,1).*b(:,4)+a(:,2).*b(:,3)-a(:,3).*b(:,2)+a(:,4).*b(:,1);
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