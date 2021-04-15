function  [a,resnorm]=Optimize_my_LM2(Loss_fun,a0,data,TolX,TolFun,MaxIter)

% author  Zhang Xin


tao=1e-10;
xk=a0;

v=2;


Jacobi=Get_Jacobi(Loss_fun,xk,data);
Ek=Loss_fun(xk,data);

g=Jacobi'*Ek;

found=logical(norm(g)<=TolFun);
mou=tao*max(diag(Jacobi'*Jacobi));

k=0;
while (~found && k<MaxIter+1)
    
    %delta_x=-(Jacobi'*Jacobi+mou*eye(size(a0,2)))\Jacobi'*Ek;    
    
    delta_x=-[Jacobi;sqrt(mou)*eye(size(a0,2))]\[Ek;zeros(size(a0,2),1)];      

    
    if (norm(delta_x)<=TolX*(norm(xk)+TolX))
        found=true;

    else
        xk_new=xk+delta_x';
        Ek=Loss_fun(xk,data);
        Ek_new=Loss_fun(xk_new,data);
        L0=delta_x'*Jacobi'*Ek;
        L_delta=delta_x'*Jacobi'*Jacobi*delta_x;
        rho=(Ek'*Ek-Ek_new'*Ek_new)/(-L0-L_delta);
        
        
        if rho>0
            
            fprintf('Iterations: %d, Residual: %d, Step: %d \n',k,Ek'*Ek,norm(delta_x));
            k=k+1;
            
            xk=xk_new;
            Jacobi=Get_Jacobi(Loss_fun,xk,data);
            Ek=Loss_fun(xk,data);
            g=Jacobi'*Ek;
            found=(norm(g)<=TolFun);
            mou=mou*max([1/3,1-(2*rho-1)^3]);
            v=2;
        else
            mou=mou*v;
            v=2*v;
           
        end
  
    end

end


xk=xk+delta_x';
Ek=Loss_fun(xk,data);
fprintf('Iterations: %d, Residual: %d, Step: %d \n',k,Ek'*Ek,norm(delta_x));



a=xk;

resnorm=Ek'*Ek;




end


function Jacobi=Get_Jacobi(Loss_fun,xk,data)

scale=1e-4;

%Ek=Loss_fun(xk,data);

for i=1:length(xk)
    x_temp1=xk;
    x_temp2=xk;
    if abs(x_temp1(i))>scale
  
        delta=x_temp1(i)*scale;
    else
        delta=scale;
    end
    x_temp1(i)=x_temp1(i)+delta;
    x_temp2(i)=x_temp2(i)-delta;

    E_temp1=Loss_fun(x_temp1,data);

    E_temp2=Loss_fun(x_temp2,data);

    Jacobi(:,i)=(E_temp1-E_temp2)/delta/2;

end

end



