% 欧拉伯努利梁振动仿真
clc; 
close all;
clear;
tmax=7;
L=1;
dx=0.05;
dt=dx^2/5;
nx=round(L/dx+1);
nt=round(tmax/dt+1);
Ttr=round(0.1/dt);
EI=5;
T=6;
rho=1;
M=5;

u0=zeros(nx,1);u1=u0;u2=u1;
wl_free=zeros(nt,1);
w_free=zeros(Ttr,nx);
%% 调试
% for j=1:nt   
%     d(j)=( 2+sin(5*(j-1)*dt*pi) )/10^-1;
% end
d=zeros(1,nt);
%%
% initial condition
for i=3:nx
    u0(i)=(i-2)*dx/L/10;
    u1(i)=u0(i);
end
%revise for drawing
w_free(1,:)=u0;
wl_free([1,2])=repmat(u0(nx),2,1);
for j=3:nt
      u2(1)=0;
      u2(2)=0;
    for i=3:nx-2
        wxx=( u1(i+1)-2*u1(i)+u1(i-1) )/dx^2;
        wxxxx=( u1(i+2)-4*u1(i+1)+6*u1(i)-4*u1(i-1)+u1(i-2) )/dx^4;
%         u2(i)=2*u1(i)-u0(i)+ ( -EI*wxxxx+T*wxx )*dt^2/rho;
        u2(i)=2*u1(i)-u0(i)+ ( -EI*wxxxx )*dt^2/rho;
    end
    wxxxl=( -u1(nx)+2*u1(nx-1)-u1(nx-2) )/dx^3;
    wxl=( u1(nx)-u1(nx-2) )/dx/2;
    
%     u2(nx)=2*u1(nx)-u0(nx)+( d(j)+EI*wxxxl-T*wxl )*dt^2/M;
    u2(nx)=2*u1(nx)-u0(nx)+( d(j)+EI*wxxxl )*dt^2/M;
    u2(nx-1)=( u2(nx)+u2(nx-2) )/2;
    
    if mod(j-1,round((nt-1)/Ttr))==0
        w_free(1+round((j-1)/(nt-1)*Ttr),:)=u2;
    end
    
    wl_free(j)=u2(nx);
    
    u0=u1;
    u1=u2;
end

x=linspace(0,L,nx);t_tr1=linspace(0,tmax,Ttr+1);t_tr2=linspace(0,tmax,nt);
figure;
plot(t_tr2,d);
xlabel('t [s]');ylabel('d(t) [N]');
title('Boundary disturbance');
figure;
surf(x,t_tr1,w_free);view([65 30]);
title('Euler-Bernoulli beam without control');
xlabel('x [m]');ylabel('t [s]');zlabel('w(x,t) [m]');
%% 动画显示
figure
filename='wave1D.gif';
h1=plot(x,w_free(1,:),'linewidth',1);
axis([0 L min(min(w_free)) max(max(w_free))])
set(gcf,'renderer','zbuffer','position',[423 184 520 300])
set(gca,'xcolor','k','ycolor','k')
set(gcf,'color','w')
grid off
im=frame2im(getframe(gcf));
[A,map]=rgb2ind(im,256);
imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',tmax/Ttr);
for i=1:Ttr
    set(h1,'XData',x,'YData',w_free(i,:));
    title(['t=',num2str(t_tr1(i)),' s'])
    drawnow;pause(tmax/Ttr)%grid off
    im=frame2im(getframe(gcf));
    [A,map]=rgb2ind(im,256);
    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',tmax/Ttr);
end
