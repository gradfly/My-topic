clc;
clear;
close all
L=1;               %������
M=5;EI=5;rho=1;
T=7;               %����ʱ���յ�
dx=0.05;           %�ռ���ɢ����
dt=dx^2/EI*rho;    %ʱ����ɢ����
t_presice=0.01;    %ʱ����ʾ����
t_cnt=round(t_presice/dt);
x=0:dx:L;          %���ɿռ������
t=0:dt:T;          %����ʱ�������
n=length(x)-1;     %�ռ��������Ŀ
m=length(t)-1;     %ʱ�����������

u0 = x/L/10;       %ԭϵͳ�ĳ�ʼλ��״̬
u1 = 0;            %ԭϵͳ�ĳ�ʼ�ٶ�״̬

%%
%���Ĳ������ֵ-��ֵ����
%rho*w_tt+EI*w_xxxx=0,
%w(0,t)=w_x(0,t)=0,w_xx(L,t)=0,
%M*w_tt(L,t)-EI*w_xxx(L,t)=0,
%w(x,0)=u0,w'(x,0)=u1.
%% �������
c = ones(n-1,1);
A = spdiags([c -4*c 6*c -4*c c],-1:3,n-1,n+1);
A(1,2)=7;A(n-1,n:n+1)=[5,-2];
c = ones(m+1,1);
B = spdiags([c -2*c],[-1,0],m+1,m+1);
w=zeros(n+1,m+1);
w(:,1)=u0;w(:,2)=u0+dt*u1;
r1=EI/dx^4/rho*dt^2;
r2=EI*(dt)^2/M/(dx)^3;
for j=2:m
    for i=2:n
        w(i,j+1)=-r1*A(i-1,:)*w(:,j)-B(j,:)*w(i,:)';
    end
    w(n+1,j+1)=r2*(-w(n-1,j)+2*w(n,j)-w(n+1,j))+2*w(n+1,j)-w(n+1,j-1);j/m
end
surf(x,t(1:t_cnt:end),w(:,1:t_cnt:end)');
xlabel x,ylabel t,zlabel('w'),view([65 30]);

%% ������ʾ
figure
filename='wave1D_11.gif';
h1=plot(x,w(:,1),'linewidth',1);
axis([0 L min(min(w)) max(max(w))])
set(gcf,'renderer','zbuffer','position',[423 184 520 300])
set(gca,'xcolor','k','ycolor','k')
set(gcf,'color','w')
grid off
im=frame2im(getframe(gcf));
[A,map]=rgb2ind(im,256);
imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',t_presice);
for j=1:t_cnt:m
    set(h1,'XData',x,'YData',w(:,j));
    title(['t=',num2str((j-1)*dt,'%.2f'),' s'])
    drawnow;%pause(t_presice)
    im=frame2im(getframe(gcf));
    [A,map]=rgb2ind(im,256);
    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',t_presice);
end