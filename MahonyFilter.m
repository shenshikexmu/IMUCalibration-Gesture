function [Q2,eInt]=MahonyFilter(Q1,ImuData,t,Vm,eInt)
% Complementary filter to Gyro atitude with Accelerate & Magnetic
% conference <Nonlinear Complementery Filters on the Special Orthogonal Group>
% inspired by    http://blog.csdn.net/luoshi006/article/details/51513580
% author  Zhang Xin

if isempty(eInt)
   eInt=[0,0,0];
end

Kp=0.0015;  %0.175;
Ki=0.000001; %     0.001;

wx=ImuData(1,5)*t;
wy=ImuData(1,6)*t;
wz=ImuData(1,7)*t;

Ak=[ 1    , -wx/2 , -wy/2 , -wz/2 ;...
     wx/2 ,   1   ,  wz/2 , -wy/2 ;...
     wy/2 , -wz/2 ,   1   ,  wx/2 ;...
     wz/2 ,  wy/2 , -wx/2 ,   1   ];
Qp=Ak*Q1';     %Q_predict

%Qp=Qp/norm(Qp);
acc=ImuData(1,2:4);
norm_a=norm(acc);
norm_g=norm(ImuData(1,5:7));

if abs(norm_a-9.8)<2 && norm_g< 2          % when acc and mag desirable
    R2=[2*Qp(2)*Qp(3)+2*Qp(1)*Qp(4);...
        2*Qp(1)^2+2*Qp(3)^2-1;...
        2*Qp(3)*Qp(4)-2*Qp(1)*Qp(2)];
    R3=[2*Qp(2)*Qp(4)-2*Qp(1)*Qp(3);...
        2*Qp(3)*Qp(4)+2*Qp(1)*Qp(2);...
        2*Qp(1)^2+2*Qp(4)^2-1 ];   

    h1=R3;                           %acc prediction from attitude in sensor fixed frame
    h2=Vm(2)*R2+Vm(3)*R3;      %mag prediction from attitude in sensor fixed frame
    e1=cross(acc/norm_a,h1');
    
    mag=ImuData(1,8:10);
    mag=mag/norm(mag);               

    e2=cross(mag,h2');                % rotation difference bewteen real mag & prediction
    e=e1+e2;
else 
    e=[0,0,0];
end

eInt = eInt + Ki * e;                  %  I in the PI control
g = [wx,wy,wz] + Kp * e + eInt;           % PI control


Akk=[ 1    , -g(1)/2 , -g(2)/2 , -g(3)/2 ;...
     g(1)/2 ,   1   ,  g(3)/2 , -g(2)/2 ;...
     g(2)/2 , -g(3)/2 ,   1   ,  g(1)/2 ;...
     g(3)/2 ,  g(2)/2 , -g(1)/2 ,   1   ];
Q2=Akk*Q1';

Q2=Q2'/norm(Q2);

if Q2(1)<0
    Q2=-Q2;
end
    

end




