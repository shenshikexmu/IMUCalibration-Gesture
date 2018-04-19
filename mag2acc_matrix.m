function [Tm2a,Bm,Vm]=mag2acc_matrix(fix_point,Ta,Ka,Ba,a0)

% author  Zhang Xin

if nargin==4 
   a0=[1,1,1,1,1,1,1,1,0,0,0];
end
B{1}=fix_point;
B{2}=Ta;
B{3}=Ka;
B{4}=Ba;
 options=optimset('TolX',1e-6,'Algorithm','Levenberg-Marquardt',...
 'Display','iter');
% options=optimset('Algorithm','Levenberg-Marquardt',...
% 'Display','iter');
%a=lsqnonlin(@elliposoid_acc,a0,[],[],options,B(:,1:3),zeros(size(B,1),1));
[a]=lsqnonlin(@elliposoid_H_acc,a0,[],[],options,B);

Tm2a=[a(1)   , a(2),  a(3);...
    a(4) ,  a(5)   , a(6);...
    a(7) ,  a(8),   -1];

Bm=[a(9);a(10);a(11)];
for i=1:size(fix_point,1)
    
    Acc(:,i)=Ta*Ka*(fix_point(i,1:3)'+Ba);
    Mag(:,i)=Tm2a*(fix_point(i,7:9)'+Bm);     
    SS(i)=Acc(:,i)'*Mag(:,i)/(norm(Acc(:,i))*norm(Mag(:,i)));
end

Mag_z=mean(SS);
Mag_y=sqrt(1-Mag_z^2);
Vm=[0,Mag_y,Mag_z];

end

function E=elliposoid_H_acc(a, x )


B=x{1};
Ta=x{2};
Ka=x{3};
Ba=x{4};
m=size(B,1);
Tm2a=[a(1)   , a(2),  a(3);...
    a(4) ,  a(5)   , a(6);...
    a(7) ,  a(8),   -1];

Bm=[a(9);a(10);a(11)];

for i=1:m
    Acc(:,i)=Ta*Ka*(B(i,1:3)'+Ba);
    Mag(:,i)=Tm2a*(B(i,7:9)'+Bm);    
       
end
for i=1:m
    if i==1
        E(i)=Acc(:,i)'*Mag(:,i)-Acc(:,m)'*Mag(:,m); 
        E(i+m)=norm(Mag(:,i))-norm(Mag(:,m));
    else
        E(i)=Acc(:,i)'*Mag(:,i)-Acc(:,i-1)'*Mag(:,i-1);
        E(i+m)=norm(Mag(:,i))-norm(Mag(:,i-1));
    end
end

end