clc;
clear;
close all
k=10;L=1;
T=10;              %仿真时间终点
dx=0.05;           %空间离散步长
dt=dx^2/5;         %时间离散步长
x=0:dx:L;          %生成空间网格点
t=0:dt:T;          %生成时间网格点
t_presice=0.01;    %时间显示精度
t_cnt=round(t_presice/dt);
n=length(x)-1;     %空间网格的数目
m=length(t)-1;     %时间仿真计算次数

u0 = sin(x)/10;    %原系统的初始位置状态
u1 = 1+cos(x);     %原系统的初始速度状态

%%
%中心差分求解初值-边值问题
%w_tt+w_xxxx=0,
%w(0,t)=w_x(0,t)=0,w_xx(L,t)=0,
%w_xxx(L,t)=kw_xt(L,t),
%w(x,0)=u0,w'(x,0)=u1.
%% 构造矩阵
c = ones(n-1,1);
A = spdiags([c -4*c 6*c -4*c c],-1:3,n-1,n+1);
A(1,2)=7;A(n-1,n:n+1)=[5,-2];
c = ones(m+1,1);
B = spdiags([c -2*c],[-1,0],m+1,m+1);
w=zeros(n+1,m+1);w(w==0)=nan;w(1,:)=0;
w(:,1)=u0;w(:,2)=u0+dt*u1;w(1,2)=0;
r=dt/dx^2;
%% 动态显示
figure
surf(x,t(1:10*t_cnt:end),w(:,1:10*t_cnt:end)');
set(gcf,'renderer','zbuffer','position',[0 184 500 300])
xlabel x,ylabel t,zlabel('w'),view([65 30]);
h=gca;
figure
filename='wave1D1.gif';
h1=plot(x,w(:,1),'linewidth',1);
h2=gca;
axis([0 L min(min(w)) max(max(w))])
set(gcf,'renderer','zbuffer','position',[524 184 500 300])
set(gcf,'color','w')
grid off
h3=gcf;
im=frame2im(getframe(h3));
[png,map]=rgb2ind(im,256);
imwrite(png,map,filename,'gif','LoopCount',Inf,'DelayTime',t_presice);

%%
for j=2:m
    for i=2:n
        w(i,j+1)=-r^2*A(i-1,:)*w(:,j)-B(j,:)*w(i,:)';
    end
    w(n+1,j+1)=r/k*(-w(n-1,j)+2*w(n,j)-w(n+1,j))+w(n,j+1)+w(n+1,j)-w(n,j);
    if mod(j,t_cnt)==1
        axis(h2,[0 L min(min(w)) max(max(w))])
        set(h1,'XData',x,'YData',w(:,j));
        title(h2,['t=',num2str((j-1)*dt,'%.2f'),' s'])
        drawnow;%pause(t_presice)
        im=frame2im(getframe(h3));
        [png,map]=rgb2ind(im,256);
        imwrite(png,map,filename,'gif','WriteMode','append','DelayTime',t_presice);
    end
    if mod(j,10*t_cnt)==1
        surf(h,x,t(1:10*t_cnt:end),w(:,1:10*t_cnt:end)');
        xlabel x,ylabel t,zlabel('w'),view(h,[65 30]);
    end
end
