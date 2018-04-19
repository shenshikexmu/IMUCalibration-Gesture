function [Ta,Ka,Ba,axis]=ICRA2014_acc(fix_point,a0)
% conference: A Robust and Easy to implement method for imu
% calibration without External Equipments
% author  Zhang Xin

if nargin<2 
   a0=[0,0,0,0.0048,0.0048,0.0048,0,0,0];
end

B=fix_point(:,1:3)-ones(size(fix_point,1),1)*(mean(fix_point(:,1:3),1));
 options=optimset('TolX',1e-6,'Algorithm','Levenberg-Marquardt',...
 'Display','iter');
% options=optimset('Algorithm','Levenberg-Marquardt',...
% 'Display','iter');
%a=lsqnonlin(@elliposoid_acc,a0,[],[],options,B(:,1:3),zeros(size(B,1),1));
[a,resnorm]=lsqnonlin(@elliposoid_acc,a0,[],[],options,B(:,1:3));


Ta=[1 , -a(1),  a(2);...
    0 ,  1   , -a(3);...
    0 ,  0   ,   1];
Ka=[a(4) ,  0   ,  0;...
    0    , a(5) ,  0;...
    0    ,  0   , a(6)];
Ba=-mean(fix_point(:,1:3),1)'+[a(7);a(8);a(9)];
axis=[1/a(4);1/a(5);1/a(6)];
E=elliposoid_acc(a, B(:,1:3) );

end

function E=elliposoid_acc(a, x )

Ta=[1 , -a(1),  a(2);...
    0 ,  1   , -a(3);...
    0 ,  0   ,   1];
Ka=[a(4) ,  0   ,  0;...
    0    , a(5) ,  0;...
    0    ,  0   , a(6)];
Ba=[a(7);a(8);a(9)];

for i=1:size(x,1)
    E(i,1)=(9.8-(norm(Ta*Ka*(x(i,1:3)'+Ba))));
end


% E=0;
% for i=1:size(x,1)
%     E=E+(9.8^2-norm(Ta*Ka*(x(i,1:3)'+Ba))^2)^2;
% end
end