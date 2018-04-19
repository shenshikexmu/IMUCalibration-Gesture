function [Ta,Ka,Ba,Tg,Kg,Bg,Tm2a,Bm,Vm]=ImuCalibration_Gesture(data)



[~,fix_point,rotation]=FindFixData(data,30);

[Ta,Ka,Ba]=ICRA2014_acc(fix_point);

Bg=-mean(fix_point(:,4:6),1)';

n=size(rotation,1);

rotation{n+1}=Ta;
rotation{n+2}=Ka;
rotation{n+3}=Ba;
rotation{n+4}=Bg;

[Tg,Kg]=ICRA_2014_gyro(rotation);


[Tm2a,Bm,Vm]=mag2acc_matrix(fix_point,Ta,Ka,Ba);



end