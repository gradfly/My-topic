clc;
clear;
close all
k=6;
L=1;                %������
M=5;EI=5;rho=1;
T=10;               %����ʱ���յ�
dx=0.05;            %�ռ���ɢ����
dt=dx^2/5;          %ʱ����ɢ����
t_presice=0.01;     %ʱ����ʾ����
t_cnt=round(t_presice/dt);
x=0:dx:L;           %���ɿռ������
t=0:dt:T;          %����ʱ�������
n=length(x)-1;     %�ռ��������Ŀ
m=length(t)-1;     %ʱ�����������

u0 = x/L;          %ԭϵͳ�ĳ�ʼλ��״̬
u0 = [0,u0(1:end-1)]/50;
u1 = 0;            %ԭϵͳ�ĳ�ʼ�ٶ�״̬

%%
%���Ĳ������ֵ-��ֵ����
%w_tt+w_xxxx=0,
%w(0,t)=w_x(0,t)=0,
%w_xx(L,t)=0,w_xxx(L,t)=kw_tt(L,t),
%w(x,0)=u0,w'(x,0)=u1.
%�������
c = ones(n,1);
A = spdiags([c -4*c 6*c -4*c c],-2:2,n,n);
A(1)=7;A(n-1,n-1:n)=[5,-2];A(n,n-2:n)=[2,-4,2-4*M*dx^3/EI/dt^2];
c = ones(m+1,1);
B = spdiags([c -2*c],[-1,0],m+1,m+1);
w=zeros(n+1,m+1);
w(:,1)=u0;w(:,2)=u0+u1*dx;
r1=EI/dx^4/rho*dt^2;
for j=2:m
    for i=2:n+1
        w(i,j+1)=-r1*A(i-1,:)*w(2:end,j)-B(j,:)*w(i,:)';
        if i==n+1
            w(i,j+1)=w(i,j+1)-r1*(-4*w(n,j)+w(n-1,j));
        end
    end
    j/m
end
surf(x,t(1:t_cnt:end),w(:,1:t_cnt:end)');
xlabel x,ylabel t,zlabel('w'),view([65 30]);

%%
%������ʾ
figure
filename='wave1D.gif';
h1=plot(x,w(:,1),'linewidth',1);
axis([0 L min(min(w)) max(max(w))])
set(gcf,'renderer','zbuffer','position',[423 184 520 300])
set(gca,'xcolor','k','ycolor','k')
set(gcf,'color','w')
grid off
im=frame2im(getframe(gcf));
[A,map]=rgb2ind(im,256);
% pause(1)
imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',t_presice);
for j=1:t_cnt:m
    set(h1,'XData',x,'YData',w(:,j));
    title(['t=',num2str((j-1)*dt,'%.2f'),' s'])
    drawnow;%pause(0.1)%grid off
    im=frame2im(getframe(gcf));
    [A,map]=rgb2ind(im,256);
    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',t_presice);
end