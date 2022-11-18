% Initialization of variables

n = 48; % Number of nodes
T_i=(25)*ones(n,n); % Initial condition
T=zeros(n,n);
T_old=T_i;


l=0.1; % length in meters
h=0.25; % height in meters
x=linspace(0,l,n);
y=linspace(0,h,n);
dx=x(2)-x(1);
dy=y(2)-y(1);
dt=power(10,-3);


%%% Thermal properties
h=10;
Ta=30;
%k=15;
rho =7500;
Fo=alpha*(dt/dx^2);
k=14.6+(1.27*10^(-2)*T_old)
cp=450+(0.28*T_old)-(2.91*10^(-4)*T_old^(-2))+(1.34*10^(-7)*T_old^(-3))
alpha=k/(rho*cp)
%alpha=3.3*10^(-6); 
tol=1*power(10,-5); % Tolerance Value
t=0;
n_iter=0; % Number of iterations
error=10; % Error value


if (Fo>=10^(-3) && Fo<=10^(-6))

while(error>tol)
    for i=1:n
    
       
          %  Middle nodes
        if(2<=i && i<=n-1)
            for j=2:n-1
                k=14.6+(1.27*10^(-2)*T_old(i,j));
                cp=450+(0.28*T_old(i,j))-(2.91*10^(-4)*(T_old(i,j))^(-2))+(1.34*10^(-7)*(T_old(i,j))^(-3));
                alpha=k/(rho*cp);
                t1=alpha*dt*(T_old(i,j+1)+T(i,j-1)-2*T_old(i,j))/(dx^2);
                t2=alpha*dt*(T_old(i+1,j)+T(i-1,j)-2*T_old(i,j))/(dy^2);
                t3=alpha*dt*(8.32*exp(((dx*i)^2))-((dy*j)^3))/k;
                T(i,j)=T_old(i,j)+t1+t2+t3;
            end
        end
          % Left wall
        
            for i=1:n
                T(i,1)=60;
            end
    
       
          % Right wall
        
            for i=1:n
                T(i,n)=T(i,n-1);
            end
        
       
          % Bottom wall
       
        
        for j=2:n-1
           k = 14.6+(1.27*10^(-2)*T_old(i,j));
           T(n,j)=T(n-1,j)+((power(10,3)*50*(dx*i)*dy/k));
        end
       
          % Top wall
       
       
        for j=2:n-1
            k = 14.6+(1.27*10^(-2)*T_old(i,j));
            w=1-(h*dy)/k;
            T(1,j)=T_old(2,j)-(h*dy*Ta/k);
            T(1,j)=T(1,j)/w;
        end
    end
   
    
    n_iter=n_iter+1;
    t=n_iter*dt;
    error=max(max(abs(T-T_old)));
    T_old=T;
    
     if (t== 0.407675e+03)
        figure(1)
        [C,u]=contourf(x,y,T);
        clabel(C,u)
        colorbar;
        colormap(jet);
        set(gca,'ydir','reverse');
        xlabel('X axis')
        ylabel('Y axis')
     end
     if(t==0.81535e+03)
        figure(2)
        [C,u]=contourf(x,y,T);
        clabel(C,u)
        colorbar;
        colormap(jet);
        set(gca,'ydir','reverse');
        xlabel('X axis')
        ylabel('Y axis')
     end
     if(t==1.223025e+03)
        figure(3)
        [C,u]=contourf(x,y,T);
        clabel(C,u)
        colorbar;
        colormap(jet);
        set(gca,'ydir','reverse');
        xlabel('X axis')
        ylabel('Y axis')
     end
     end
    
     if(error<tol)
            disp('STEADYSTATE REACHED')
            disp('Time to reach= ')
            disp(t)
            disp('MAXIMUM TEMPERATURE IS= ')
            disp(max(max(T)))
     end
end
     
     figure(4)
     plot(T(:,n/4),y)
     ylabel('HEIGHT')
     xlabel('TEMPERATURE')
     
     figure(5)
     plot(T(:,(3*n)/4),y)
     ylabel('HEIGHT')
     xlabel('TEMPERATURE')
     
     figure(6)
     plot(x,T(n/4,:))
     xlabel('LENGTH')
     ylabel('TEMPERATURE')
     
     figure(7)
     plot(x,T(n/12,:))
     xlabel('LENGTH')
     ylabel('TEMPERATURE')

     figure(8)
     [C,u]=contourf(x,y,T);
     clabel(C,u)
     colorbar;
     colormap(jet);
     set(gca,'ydir','reverse');
     xlabel('X axis')
     ylabel('Y axis')           