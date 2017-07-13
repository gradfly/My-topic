clc;
close all;
clear;
%%
%*************************************
%          Flexible Wings  
%          without  Control
%*************************************
nx=12;                                    %��չС��
nt=5*10^4;                                %ʱ���

L=2;                                      %��չ����L
tmax=100;                                 
Ttr=50;

dx=L/(nx-1);   
dt=tmax/(nt-1);

%parameters
m=10;                                    %�����ܶ�
EIb=3;                                   %�˶��ն�
xec=0.25;                                %���ĵ����ļ������ĵľ���
xac=0.05;                                %��������ѧ���ĵ����ļ������ĵľ���
eta=0.022;                               %Kelvin-Voigt����ϵ��
Ip=1.5;                                  %�������ļ����Ծ�
GJ=0.2;                                  %Ťת�ն�



w_free=zeros(nx,nt);              %λ�ƾ��� ������12*��5*10^4����0����
th_free=zeros(nx,nt);             %ת���ǶȾ��� ����12*��5*10^4����0����
Fb=w_free;                        %Fb�����ţ�ͬ��
F=zeros(nt,1);                    %��������ľ����˶���������5*10^4��*1�е��о���           
M=F;                              %��������ľ���Ťת��������5*10^4��*1�е��о���             

w_3D_free=zeros(Ttr,nx);          %        ����50*12�ľ���
th_3D_free=w_3D_free;             %        ����50*12�ľ���


% disturbance
for i=1:nx
    for j=1:nt
        Fb(i,j)=( 1+sin(1*(j-1)*pi)+3*cos((j-1)*3*pi) )*(i-1)*dx/1;
    end
end


% initial condition
for i=1:nx
    w_free(i,1)=(i-1)*dx/L;        %    w(x,0)=x/L
    th_free(i,1)=pi/2*(i-1)*dx/L;  %    th(x,0)=pi*x/(2L)
end

w_free(:,2)=w_free(:,1);     %\dot{w(x,0)}=0
th_free(:,2)=th_free(:,1);    %\dot{th(x,0)}=0


% main cycle
for j=2:nt-1
    w_free(1,j+1)=0;                  %w(0,t)=0
    w_free(2,j+1)=w_free(1,j+1);      %w'(0,t)=0
    th_free(1,j+1)=0;                %th(0,t)=0
    
    
    % when i=2
    ddw=(w_free(2,j+1)-2*w_free(2,j)+w_free(2,j-1))/dt^2;               %    \ddot{w(2,t)}  
    dthxx=((th_free(3,j)-2*th_free(2,j)+th_free(1,j))-(th_free(3,j-1)-2*th_free(2,j-1)+th_free(1,j-1)))/dt/dx^2;      %   \dot{th(2,t)}''
    thxx=(th_free(3,j)-2*th_free(2,j)+th_free(1,j))/dx^2;               %  \th''(2,t) 
    S3=-xac*Fb(2,j)+GJ*thxx+eta*GJ*dthxx+m*xec*ddw;                          %I_p \ddot{th(2,t)}
    th_free(2,j+1)=2*th_free(2,j)-th_free(2,j-1)+S3*dt^2/Ip;                 %��ʵ����\ddot{th(2,t)}�ı��Σ�2��ʽ���õ�th_free(2, j+1)
    
    
    for i=3:nx-2
        
        dwxxxx=( ( w_free(i+2,j)-4*w_free(i+1,j)+6*w_free(i,j)-4*w_free(i-1,j)+w_free(i-2,j) )-( w_free(i+2,j-1)-4*w_free(i+1,j-1)+6*w_free(i,j-1)-4*w_free(i-1,j-1)+w_free(i-2,j-1) ))/dt/dx^4;  %\dot{w}''''(x,t)
        wxxxx=( w_free(i+2,j)-4*w_free(i+1,j)+6*w_free(i,j)-4*w_free(i-1,j)+w_free(i-2,j) )/dx^4;      %w''''(x,t)
        dthxx=( ( (th_free(i+1,j)-2*th_free(i,j)+th_free(i-1,j))/dx^2 )-( (th_free(i+1,j-1)-2*th_free(i,j-1)+th_free(i-1,j-1))/dx^2 ) )/dt;  %\dot{th(x,t)}''
        thxx=(th_free(i+1,j)-2*th_free(i,j)+th_free(i-1,j))/dx^2 ;        %   th''(x,t)
        
        e=Fb(i,j)-EIb*wxxxx-eta*EIb*dwxxxx;                           %e
        f=-xac*Fb(i,j)+GJ*thxx+eta*GJ*dthxx;                          %f
        
        a=m;
        b=-m*xec;
        c=Ip;
        d=-m*xec;

       
      S1=(c*e-b*f)/(a*c-b*d)  ;      %ddw=(ce-bf)/(ac-bd)                            
      S2=(a*f-d*e)/(a*c-b*d)  ;      %ddth=(ed-af)/(bd-ac),ԭ����������(a*f-b*e)/(a*c-b*d)
       
       
       w_free(i,j+1)=2*w_free(i,j)-w_free(i,j-1)+S1*dt^2;     %central differences method
       th_free(i,j+1)=2*th_free(i,j)-th_free(i,j-1)+S2*dt^2;
    end
    
    
    w_free(nx,j+1)=( F(j-1)*dx^3/EIb+w_free(nx-3,j+1)-1.5*w_free(nx-2,j+1) )/-0.5;  %�߽����������޲��
    w_free(nx-1,j+1)=( w_free(nx,j+1)+w_free(nx-2,j+1) )/2;                         %��ֵ
    
    
    
    % when i=nx-1
    
    ddw=( w_free(nx-1,j+1)-2*w_free(nx-1,j)+w_free(nx-1,j-1) )/dt^2;    %\ddot{w(nx-1,t)}
    dthxx=( ( th_free(nx,j)-2*th_free(nx-1,j)+th_free(nx-2,j) )-( th_free(nx,j-1)-2*th_free(nx-1,j-1)+th_free(nx-2,j-1) ) )/dt/dx^2;   %\dot{th(nx-1,t)}''
    thxx=( th_free(nx,j)-2*th_free(nx-1,j)+th_free(nx-2,j) )/dx^2;       % th(nx-1,t)''���Ը���-xac*Fb����th_free(nx-1,j+1)
    
    S4=-xac*Fb(nx-1,j)+GJ*thxx+eta*GJ*dthxx+m*xec*ddw;         
    
    th_free(nx-1,j+1)=2*th_free(nx-1,j)-th_free(nx-1,j-1)+S4*dt^2/Ip;   %  th(nx-1,t)
    
    % when i=nx    
    th_free(nx,j+1)=th_free(nx-1,j+1)+M(j-1)*dx/GJ;    %�߽�����            
    

       if mod(j,nt/Ttr)==0
           w_3D_free(1+j*Ttr/nt,:)=w_free(:,j)';
           th_3D_free(1+j*Ttr/nt,:)=th_free(:,j)';
       end
end
    

w_3D_free(1,:)=w_free(:,1)';
th_3D_free(1,:)=th_free(:,1)';



%%
%*************************************
%          Flexible Wings  
%           with Control
%*************************************



%patameters of controllers

alpha=50;
beta=1;
k1=1;
k2=1.5*10^-2;


%create matrixes to save data
w_control=zeros(nx,nt);
th_control=zeros(nx,nt);
F=zeros(nt,1);
M=F;
U=F;
V=M;

w_3D_control=zeros(Ttr,nx);
th_3D_control=w_3D_control; 


% initial condition
for i=1:nx
    w_control(i,1)=(i-1)*dx/L;
    th_control(i,1)=pi/2*(i-1)*dx/L;
end

w_control(:,2)=w_control(:,1);
th_control(:,2)=th_control(:,1);


% main cycle
for j=2:nt-1
    w_control(1,j+1)=0;
    w_control(2,j+1)=w_control(1,j+1);
    th_control(1,j+1)=0;
    
    
    % when i=2
    ddw=( w_control(2,j+1)-2*w_control(2,j)+w_control(2,j-1) )/dt^2;
    dthxx=( ( th_control(3,j)-2*th_control(2,j)+th_control(1,j) )-( th_control(3,j-1)-2*th_control(2,j-1)+th_control(1,j-1) ) )/dt/dx^2;
    thxx=( th_control(3,j)-2*th_control(2,j)+th_control(1,j) )/dx^2;
    
    S3=-xac*Fb(2,j)+GJ*thxx+eta*GJ*dthxx+m*xec*ddw;
    th_control(2,j+1)=2*th_control(2,j)-th_control(2,j-1)+S3*dt^2/Ip;    
    
    
    for i=3:nx-2
        
        dwxxxx=( ( w_control(i+2,j)-4*w_control(i+1,j)+6*w_control(i,j)-4*w_control(i-1,j)+w_control(i-2,j) )-( w_control(i+2,j-1)-4*w_control(i+1,j-1)+6*w_control(i,j-1)-4*w_control(i-1,j-1)+w_control(i-2,j-1) ))/dt/dx^4;
        wxxxx=( w_control(i+2,j)-4*w_control(i+1,j)+6*w_control(i,j)-4*w_control(i-1,j)+w_control(i-2,j) )/dx^4;
        dthxx=( ( (th_control(i+1,j)-2*th_control(i,j)+th_control(i-1,j))/dx^2 )-( (th_control(i+1,j-1)-2*th_control(i,j-1)+th_control(i-1,j-1))/dx^2 ) )/dt;
        thxx=(th_control(i+1,j)-2*th_control(i,j)+th_control(i-1,j))/dx^2 ;
        
        e=Fb(i,j)-EIb*wxxxx-eta*EIb*dwxxxx;%e
        f=-xac*Fb(i,j)+GJ*thxx+eta*GJ*dthxx;%f
        
        a=m;
        b=-m*xec;
        c=Ip;
        d=-m*xec;

       
        S1=(c*e-b*f)/(a*c-b*d)  ;                                        %S1=(a*f-c*e)/(a*d-b*c);  
        S2=(a*f-d*e)/(a*c-b*d)  ;                                        %S2=-( b*f-d*e )/(a*d-b*c);
       
       
       w_control(i,j+1)=2*w_control(i,j)-w_control(i,j-1)+S1*dt^2;
       th_control(i,j+1)=2*th_control(i,j)-th_control(i,j-1)+S2*dt^2;
    end  
   
    
    
    w_control(nx,j+1)=( F(j-1)*dx^3/EIb+w_control(nx-3,j+1)-1.5*w_control(nx-2,j+1) )/-0.5;
    w_control(nx-1,j+1)=( w_control(nx,j+1)+w_control(nx-2,j+1) )/2;
    
    
    
    % when i=nx-1
    
    ddw=( w_control(nx-1,j+1)-2*w_control(nx-1,j)+w_control(nx-1,j-1) )/dt^2;
    dthxx=( ( th_control(nx,j)-2*th_control(nx-1,j)+th_control(nx-2,j) )-( th_control(nx,j-1)-2*th_control(nx-1,j-1)+th_control(nx-2,j-1) ) )/dt/dx^2;
    thxx=( th_control(nx,j)-2*th_control(nx-1,j)+th_control(nx-2,j) )/dx^2;
    
    S4=-xac*Fb(nx-1,j)+GJ*thxx+eta*GJ*dthxx+m*xec*ddw;
    
    th_control(nx-1,j+1)=2*th_control(nx-1,j)-th_control(nx-1,j-1)+S4*dt^2/Ip;
%     th(nx-1,j+1)=th(nx-2,j+1);
    
    % when i=nx
    
    th_control(nx,j+1)=th_control(nx-1,j+1)+M(j-1)*dx/GJ;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ϺͲ��ӿ�����һ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    % compute controllers
    
    U(j)=k1*( alpha*w_control(nx,j+1)+beta*( w_control(nx,j+1)-w_control(nx,j) )/dt );
    V(j)=-k2*( alpha*th_control(nx,j+1)+beta*( th_control(nx,j+1)-th_control(nx,j) )/dt );
    
    
    F(j)=( U(j)*dt+eta*F(j-1) )/( dt+eta );
    M(j)=( V(j)*dt+eta*M(j-1) )/( dt+eta );

          
       if mod(j,nt/Ttr)==0
           w_3D_control(1+j*Ttr/nt,:)=w_control(:,j)';
           th_3D_control(1+j*Ttr/nt,:)=th_control(:,j)';
       end
end
    

w_3D_control(1,:)=w_control(:,1)';
th_3D_control(1,:)=th_control(:,1)';


x1=linspace(0,L,nx);
t_tr=linspace(0,tmax,Ttr);  

figure(1);
surf(x1,t_tr,w_3D_free);view([60 30]);
xlabel('x [m]');ylabel('t [s]');zlabel('w(x,t) [m]');
title('Bending displacement w(x,t) without control');
    
figure(2);
surf(x1,t_tr,th_3D_free);view([60 30]);
xlabel('x [m]');ylabel('t [s]');zlabel('\theta(x,t) [rad]');     
title('Twist displacement \theta(x,t) without control');
  
figure(3);
surf(x1,t_tr,w_3D_control);view([60 30]);
xlabel('x [m]');ylabel('t [s]');zlabel('w(x,t) [m]');
title('Bending displacement w(x,t) with control');
    
figure(4);
surf(x1,t_tr,th_3D_control);view([60 30]);
xlabel('x [m]');ylabel('t [s]');zlabel('\theta(x,t) [rad]]');     
title('Twist displacement \theta(x,t) with control');       

figure(5);
subplot(211)
plot(linspace(0,tmax,nt),F);
xlabel('t [s]');ylabel('F(t) [N]');title('Tip force control input F(t)');
subplot(212)
plot(linspace(0,tmax,nt),M);
xlabel('t [s]');ylabel('M(t)  [Nm]');title('Twisting moment control input M(t)');
 
figure(6);
plot(linspace(0,tmax,nt),w_free(nx/3,:),'r:');
hold on
plot(linspace(0,tmax,nt),w_control(nx/3,:),'b');
xlabel('t [s]');ylabel('w(L/3,t) [m]');
legend('w(L/3,t) without control','w(L/3,t) with control');
 
figure(7);
plot(linspace(0,tmax,nt),w_free(nx,:),'r:');
hold on
plot(linspace(0,tmax,nt),w_control(nx,:),'b');
xlabel('t [s]');ylabel('w(2L/3,t) [m]');
legend('w(2L/3,t) without control','w(2L/3,t) with control');
  
figure(8);
plot(linspace(0,tmax,nt),w_free(nx,:),'r:');
hold on
plot(linspace(0,tmax,nt),w_control(nx,:),'b');
xlabel('t [s]');ylabel('w(L,t) [m]');
legend('w(L,t) without control','w(L,t) with control');
  
figure(9);
plot(linspace(0,tmax,nt),th_free(nx,:),'r:');
hold on
plot(linspace(0,tmax,nt),th_control(nx,:),'b');
xlabel('t [s]');ylabel('\theta(L/3,t) [rad]');
legend('\theta(L/3,t) without control','\theta(L/3,t) with control');
 
 figure(10);
plot(linspace(0,tmax,nt),th_free(nx,:),'r:');
hold on
plot(linspace(0,tmax,nt),th_control(nx,:),'b');
xlabel('t [s]');ylabel('\theta(2L/3,t) [rad]');
legend('\theta(2L/3,t) without control','\theta(2L/3,t) with control');
 
figure(11);
plot(linspace(0,tmax,nt),th_free(nx,:),'r:');
hold on
plot(linspace(0,tmax,nt),th_control(nx,:),'b');
xlabel('t [s]');ylabel('\theta(L,t) [rad]');
legend('\theta(L,t) without control','\theta(L,t) with control');
 

