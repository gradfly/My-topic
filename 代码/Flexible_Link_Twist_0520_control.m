clc;
close all;
clear ;
%%
%*************************************
%        Rigid Flexible Wings
%          without  Control
%*************************************
nx = 12;                           %翼展小段
nt = 5 * 10^4;                     %时间

L2 = 2;                             %翼展长度
tmax = 100;
Ttr = 50;

dx = L2 / ( nx - 1 );
dt = tmax / ( nt - 1 ); 

%flexible-link parameters
rho = 10;                           %质量密度
EIb = 3;                           %扑动刚度
xec = 0.25;                        %质心到翅膀的剪切中心的距离
xac = 0.05;                        %空气动力学中心到翅膀的剪切中心的距离
eta = 0.022;                       %Kelvin-Voigt阻尼系数
Ih = 20;
Ip = 1.5;                          %机翼截面的极惯性距
GJ = 0.2;                          %扭转刚度

w_free = zeros(nx, nt);            %形变矩阵，产生12*（5*10^4）的0矩阵
y_free = zeros(nx, nt);            %位移矩阵，产生12* (5*10^4) 的0矩阵
phi_free = zeros(nx, nt);          %扭转角度矩阵，产生12*（5*10^4）的0矩阵
Fb = w_free;                       %Fb（干扰）同上
th2_free = zeros(nt, 1);           %flexible-link转动角度矩阵，产生 (5*10^4) *1行的列矩阵
F2 = zeros(nt,1);                  %控制输入的矩阵（rigid-link 扑动）产生 (5*10^4) *1行的列矩阵
M2 = F2;                           %控制输入的矩阵 (flexible-link 扭转) 产生 (5*10^4) *1行的列矩阵

w_3D_free = zeros(Ttr, nx);        %产生50*12的矩阵
y_3D_free = w_3D_free;             %产生50*12的矩阵
phi_3D_free = w_3D_free;           %产生50*12的矩阵

% disturbance
for i = 1 : nx
    for j = 1 : nt
        Fb(i, j) = ( 1 + sin(1 * (j - 1) * pi) + 3 * cos((j - 1) * 3 * pi) ) * (i - 1) * dx;
    end
end

% flexible link
th2_free(1) = pi /10;                                     % th2(0) = pi /10
th2_free(2) = th2_free(1);
y_free( 1, :) = 0;
for i = 1 : nx
     w_free(i, 1) = ( (i - 1) * dx ) /L2 /20;            % w(x, t) = x / L2    %??
     w_free(i, 2) = w_free(i, 1);                  % \dot{w(x, 0)} = 0
     y_free(i, 1) = w_free(i, 1) + ( (i - 1) * dx ) * th2_free(1);
     y_free(i, 2) = w_free(i, 2) + ( (i - 1) * dx ) * th2_free(2);
     phi_free(i, 1) = (pi / 2) * ( i - 1 ) * dx /L2;%pi /2 * ( i - 9 ) * dx /L2 /2;  % phi(x, t) = (pi /2) * x /L2  %??
end
%w_free(:, 2) = w_free(:, 1);       % \dot{w(x, 0)} = 0
phi_free(:, 2) = phi_free(:, 1);   % \dot{phi(x, 0)} = 0

% main cycle
for j = 2 : ( nt - 1 )
    
    w_free(1, j + 1) = 0;
    w_free(2, j + 1) = w_free(1, j + 1);
    %phi_free(1,j + 1) = 0;
    
    % when L = L1 + dx
    
    wxx = ( w_free(3, j) - w_free(2, j) ) /dx^2;
    dwxx = ( ( w_free(3, j) - w_free(2, j) ) - ( w_free(3, j - 1) - w_free(2, j - 1) ) ) /dt /(dx^2);
    ddth2 = - (F2(j) + EIb * wxx + eta * dwxx) / Ih;
    %ddth2 = F2(j) / Ih;
    th2_free(j + 1) = 2 * th2_free(j) - th2_free(j - 1) + ddth2 * dt^2;
    %y_free( 9 + 1, j + 1) = y_free( 9, j + 1 ) + th2_free(j + 1) * dx;
    y_free(1, j + 1) = w_free(1, j + 1);
    y_free(2, j + 1) = w_free(2, j + 1) + dx * th2_free(j + 1);
    
    %phi_free(2, j + 1) = ( -M2(j - 1) * dx * dt - GJ * dt * phi_free(2, j) + GJ * dt * phi_free(1, j) + eta * GJ * phi_free(1, j + 1) + eta * GJ * phi_free(2, j) - eta * GJ * phi_free(1, j) ) /(eta * GJ);
    
    ddw = ( w_free(2, j + 1) - 2 * w_free(2, j) + w_free(2, j - 1) ) /dt^2;    % \ddot{w(L1 + dx, t)}
    dphixx = ( ( phi_free(2 + 1, j) - 2 * phi_free(2, j) + phi_free(2 - 1, j) ) - ( phi_free(2 + 1, j - 1) - 2 * phi_free(2, j - 1) + phi_free(2 - 1, j - 1) ) ) /dt /(dx^2);  % \dot{phi(L1 + dx, t)}''
    phixx = ( phi_free(2 + 1, j) - 2 * phi_free(2, j) + phi_free(2 - 1, j) ) /dx^2;        % \phi''(L1 + dx, t)
    
    S3 = -xac * Fb(2, j) + GJ * phixx + eta * GJ * dphixx + rho * xec * ddw;     % I_p \ddot{phi(L1 + dx, t)}
    phi_free(2, j + 1) = 2 * phi_free(2, j) - phi_free(2, j - 1) + S3 * dt^2 /Ip;  %其实就是\ddot{th(L1 + dx, t)}的变形
    
    %phi_free(2 - 1, j + 1) = ( M2(j) * dx * dt + GJ * dt * phi_free(2, j) - GJ * dt * phi_free(2 - 1, j) + eta * GJ * phi_free(2, j + 1) - eta * GJ * phi_free(2, j) + eta * GJ * phi_free(2 - 1, j) ) /(eta * GJ);
    
%     %phi_free(9 + 1, j + 1) = ( -M2(j) * dx * dt - GJ * dt * phi_free(9 + 1, j) + GJ * dt * phi_free(9, j) + eta * GJ * phi_free(9, j + 1) + eta * GJ * phi_free(9 + 1, j) - eta * GJ * phi_free(9, j) ) /(eta * GJ);
    
    for i = 3 : nx - 2
        
        wxx = ( w_free(3, j) - w_free(2, j) ) /dx^2;
        dwxx = ( ( w_free(3, j) - w_free(2, j) ) - ( w_free(3, j - 1) - w_free(2, j - 1) ) ) /dt /(dx^2);
        %ddth2 = - (F2(j - 1) + EIb * wxx + eta * dwxx) / Ih;
        ddth2 = F2(j) / Ih;
        th2_free(j + 1) = 2 * th2_free(j) - th2_free(j - 1) + ddth2 * dt^2;
        
        dwxxxx = ( ( w_free(i + 2, j) - 4 * w_free(i + 1, j) + 6 * w_free(i, j) - 4 * w_free(i - 1, j) + w_free(i - 2, j) ) - ( w_free(i + 2, j - 1) - 4 * w_free(i + 1, j - 1) + 6 * w_free(i, j - 1) - 4 * w_free(i - 1, j - 1) + w_free(i - 2, j - 1) ) ) /dt /dx^4;% \dot{w}''''(x,t)
        wxxxx = ( w_free(i + 2, j) - 4 * w_free(i + 1, j) + 6 * w_free(i, j) - 4 * w_free(i - 1, j) + w_free(i - 2, j) ) /dx^4;% w''''(x,t)
        dphixx = ( ( (phi_free(i + 1, j) - 2 * phi_free(i, j) + phi_free(i - 1, j) ) /dx^2 ) - ( (phi_free(i + 1, j - 1) - 2 * phi_free(i, j - 1) + phi_free(i - 1, j - 1) ) /dx^2 ) ) /dt;% \dot{phi(x,t)}''
        phixx = ( phi_free(i + 1, j) - 2 * phi_free(i, j) + phi_free(i - 1, j) ) /dx^2 ;% phi''(x,t)
        
        e = Fb(i, j) - EIb * wxxxx - eta * EIb * dwxxxx - rho * ( (i - 1) * dx ) * ddth2;
        f = -xac * Fb(i, j) + GJ * phixx + eta * GJ * dphixx;
        
        a = rho;
        b = -rho * xec;
        c = Ip;
        d = -rho * xec;
        
        S1 = (c * e - b * f) /(a * c - b * d); % \ddot{w}(x, t)
        S2 = (a * f - b * e) /(a * c - b * d); % \ddot{phi}(x, t)
        
        %ddy = S1 + ( (i - 1) * dx - 8 * dx ) * ddth2 + 8 * dx * ddth1;
        
        w_free(i, j + 1) = 2 * w_free(i, j) - w_free(i, j - 1) + S1 * dt^2;
        %y_free(i, j + 1) = 2 * y_free(i, j) - y_free(i, j - 1) + ddy * dt^2;
        y_free(i, j + 1) = w_free(i, j + 1) + ( (i - 1) * dx ) * th2_free(j + 1);
        phi_free(i, j + 1) = 2 * phi_free(i, j) - phi_free(i, j - 1) + S2 * dt^2;
    end
    
    %w_free(nx, j + 1) = ( F2(j - 1) * dx^3 * dt /EIb - dt * w_free(nx - 3, j + 1) + 1.5 * dt * w_free(nx - 2, j + 1) - eta * w_free(nx - 3, j + 1) + 1.5 * eta * w_free(nx - 2, j + 1) + 0.5 * eta * w_free(nx, j) + eta * w_free(nx - 3, j) - 1.5 * eta * w_free(nx - 2, j) ) /(0.5 * dt + 0.5 * eta);
%    %w_free(nx, j + 1) = ( F2(j - 1) * dx^3 /EIb + w_free(nx - 3, j + 1) - 1.5 * w_free(nx - 2, j + 1) ) /-0.5;   %？
    
    w_free(nx, j + 1) = ( 1.5 * dt * w_free(nx - 2, j + 1) - dt * w_free(nx - 3, j + 1) + 1.5 * eta * w_free(nx - 1, j + 1) - eta * w_free(nx - 3, j + 1) + 0.5 * eta * w_free(nx, j) - 1.5 * eta * w_free(nx - 2, j) + eta * w_free(nx - 3, j) ) /(0.5 * (dt + eta));
    w_free(nx - 1, j + 1) = ( w_free(nx, j + 1) + w_free(nx - 2, j + 1) ) /2;
    
    y_free(nx, j + 1) = w_free(nx, j + 1) + ( (nx - 1) * dx ) * th2_free(j + 1);
    y_free(nx - 1, j + 1) = w_free(nx - 1, j + 1) + ( (nx - 2) * dx ) * th2_free(j + 1);
    
    % when i=nx-1
    
    ddw = ( w_free(nx - 1, j + 1) - 2 * w_free(nx - 1, j) + w_free(nx - 1, j - 1) ) /dt^2;     % \ddot{w(nx-1,t)}
    dphixx = ( ( phi_free(nx, j) - 2 * phi_free(nx - 1, j) + phi_free(nx - 2, j) ) - ( phi_free(nx, j - 1) - 2 * phi_free(nx - 1, j - 1) + phi_free(nx - 2, j - 1) ) ) /dt /dx^2;    % \dot{phi(nx-1,t)}''
    phixx = ( phi_free(nx, j) - 2 * phi_free(nx - 1, j) + phi_free(nx - 2, j) ) /dx^2; % phi(nx-1,t)''
    
    S4 = -xac * Fb(nx - 1, j) + GJ * phixx + eta * GJ * dphixx + rho * xec * ddw;
    
    phi_free(nx - 1, j + 1) = 2 * phi_free(nx - 1, j) - phi_free(nx - 1, j - 1) + S4 * dt^2 /Ip; % phi(nx-1,t)
    %phi(nx - 1, j + 1) = phi(nx - 2, j + 1);
    
    % when i=nx
    
    %phi_free(nx, j + 1) = M2(j) * dx * dt /GJ * (dt + eta) + (eta /GJ * (dt + eta)) * (phi_free(nx, j) - phi_free(nx - 1, j)) + phi_free(nx - 1, j + 1);
    %phi_free(nx, j + 1) = phi_free(nx - 1, j + 1) + M2(j - 1) * dx /GJ;    %？
    
    phi_free(nx, j + 1) = (-dt * phi_free(nx, j) + dt * phi_free(nx - 1, j) + eta * phi_free(nx - 1, j + 1) + eta * phi_free(nx, j) - eta * phi_free(nx - 1, j) ) /eta;
   
    if mod(j, nt /Ttr) == 0
        w_3D_free(1 + j * Ttr /nt, :) = w_free(:, j)';
        y_3D_free(1 + j * Ttr /nt, :) = y_free(:, j)';
        phi_3D_free(1 + j * Ttr /nt, :) = phi_free(:, j)';
    end
end

w_3D_free(1,:)=w_free(:,1)';
y_3D_free(1,:)=y_free(:,1)';
phi_3D_free(1,:)=phi_free(:,1)';



%%
%*************************************
%       Rigid Flexible Wings
%           with Control
%*************************************

%patameters of controllers  

alpha = 9; %10^(-60)
beta = 0.01;   %10^(-80)

p = 0.5;  %66
p1 = 0.4; %28
p2 = 5; %36

q = 10^(-60); %60

th2_d = pi /12;

%create matrixes to save data
w_control = zeros(nx, nt);
y_control = zeros(nx, nt);
phi_control = zeros(nx, nt);
th2_control = zeros(nt, 1);
Fb = w_control;
F2 = zeros(nt, 1);
M2 = F2;

w_3D_control = zeros(Ttr, nx);
y_3D_control = w_3D_control; 
phi_3D_control = w_3D_control; 

% disturbance
for i = 1 : nx
    for j = 1 : nt
        Fb(i, j) = ( 1 + sin(1 * (j - 1) * pi) + 3 * cos((j - 1) * 3 * pi) ) * (i - 1) * dx;
    end
end

% flexible link
th2_control(1) = pi /10;                                     % th2(0) = pi /10
th2_control(2) = th2_control(1);
y_free( 1, 1) = 0;
for i = 1 : nx
     w_control(i, 1) = ( (i - 1) * dx ) /L2;            % w(x, t) = x / L2    %??
     w_control(i, 2) = w_control(i, 1);                  % \dot{w(x, 0)} = 0
     y_control(i, 1) = w_control(i, 1) + ( (i - 1) * dx ) * th2_control(1);
     y_control(i, 2) = w_control(i, 2) + ( (i - 1) * dx ) * th2_control(2);
     phi_control(i, 1) = (pi /2) * ( i - 1 ) * dx /L2;%pi /2 * ( i - 9 ) * dx /L2 /2;  % phi(x, t) = (pi /2) * x /L2  %??
end
%w_free(:, 2) = w_free(:, 1);       % \dot{w(x, 0)} = 0
phi_control(:, 2) = phi_control(:, 1);   % \dot{phi(x, 0)} = 0

% main cycle
for j = 2 : ( nt - 1 )
    
    w_control(1, j + 1) = 0;
    w_control(2, j + 1) = 0;
    %phi_control(1, j + 1) = 0;
    
    % when L = L1 + dx
    
    wxx = ( w_control(3, j) - w_control(2, j) ) /dx^2;
    dwxx = ( ( w_control(3, j) - w_control(2, j) ) - ( w_control(3, j - 1) - w_control(2, j - 1) ) ) /dt /(dx^2);
    ddth2 = - (F2(j) + EIb * wxx + eta * dwxx) / Ih;
    %ddth2 = F2(j) / Ih;
    th2_control(j + 1) = 2 * th2_control(j) - th2_control(j - 1) + ddth2 * dt^2;
    %y_control( 9 + 1, j + 1) = y_control( 9, j + 1 ) + th2_control(j + 1) * dx;
    y_control(1, j + 1) = w_control(1, j + 1);
    y_control(2, j + 1) = w_control(2, j + 1) + dx * th2_control(j + 1);
    
    %phi_control(2, j + 1) = ( -M2(j - 1) * dx * dt - GJ * dt * phi_control(2, j) + GJ * dt * phi_control(1, j) + eta * GJ * phi_control(1, j + 1) + eta * GJ * phi_control(2, j) - eta * GJ * phi_control(1, j) ) /(eta * GJ);
    
    ddw = ( w_control(2, j + 1) - 2 * w_control(2, j) + w_control(2, j - 1) ) /dt^2;    % \ddot{w(L1 + dx, t)}
    dphixx = ( ( phi_control(2 + 1, j) - 2 * phi_control(2, j) + phi_control(2 - 1, j) ) - ( phi_control(2 + 1, j - 1) - 2 * phi_control(2, j - 1) + phi_control(2 - 1, j - 1) ) ) /dt /(dx^2);  % \dot{phi(L1 + dx, t)}''
    phixx = ( phi_control(2 + 1, j) - 2 * phi_control(2, j) + phi_control(2 - 1, j) ) /dx^2;        % \phi''(L1 + dx, t)
    
    S3 = -xac * Fb(2, j) + GJ * phixx + eta * GJ * dphixx + rho * xec * ddw;     % I_p \ddot{phi(L1 + dx, t)}
    phi_control(2, j + 1) = 2 * phi_control(2, j) - phi_control(2, j - 1) + S3 * dt^2 /Ip;  %其实就是\ddot{th(L1 + dx, t)}的变形
    
    phi_control(2 - 1, j + 1) = ( M2(j) * dx * dt + GJ * dt * phi_control(2, j) - GJ * dt * phi_control(2 - 1, j) + eta * GJ * phi_control(2, j + 1) - eta * GJ * phi_control(2, j) + eta * GJ * phi_control(2 - 1, j) ) /(eta * GJ);
    
%     %phi_control(9 + 1, j + 1) = ( -M2(j) * dx * dt - GJ * dt * phi_control(9 + 1, j) + GJ * dt * phi_control(9, j) + eta * GJ * phi_control(9, j + 1) + eta * GJ * phi_control(9 + 1, j) - eta * GJ * phi_control(9, j) ) /(eta * GJ);
    
    for i = 3 : nx - 2
        
        wxx = ( w_control(3, j) - w_control(2, j) ) /dx^2;
        dwxx = ( ( w_control(3, j) - w_control(2, j) ) - ( w_control(3, j - 1) - w_control(2, j - 1) ) ) /dt /(dx^2);
        ddth2 = - (F2(j - 1) + EIb * wxx + eta * dwxx) / Ih;
        %ddth2 = F2(j) / Ih;
        th2_control(j + 1) = 2 * th2_control(j) - th2_control(j - 1) + ddth2 * dt^2;
        
        dwxxxx = ( ( w_control(i + 2, j) - 4 * w_control(i + 1, j) + 6 * w_control(i, j) - 4 * w_control(i - 1, j) + w_control(i - 2, j) ) - ( w_control(i + 2, j - 1) - 4 * w_control(i + 1, j - 1) + 6 * w_control(i, j - 1) - 4 * w_control(i - 1, j - 1) + w_control(i - 2, j - 1) ) ) /dt /dx^4;% \dot{w}''''(x,t)
        wxxxx = ( w_control(i + 2, j) - 4 * w_control(i + 1, j) + 6 * w_control(i, j) - 4 * w_control(i - 1, j) + w_control(i - 2, j) ) /dx^4;% w''''(x,t)
        dphixx = ( ( (phi_control(i + 1, j) - 2 * phi_control(i, j) + phi_control(i - 1, j) ) /dx^2 ) - ( (phi_control(i + 1, j - 1) - 2 * phi_control(i, j - 1) + phi_control(i - 1, j - 1) ) /dx^2 ) ) /dt;% \dot{phi(x,t)}''
        phixx = ( phi_control(i + 1, j) - 2 * phi_control(i, j) + phi_control(i - 1, j) ) /dx^2 ;% phi''(x,t)
        
        e = Fb(i, j) - EIb * wxxxx - eta * EIb * dwxxxx - rho * ( (i - 1) * dx ) * ddth2;
        f = -xac * Fb(i, j) + GJ * phixx + eta * GJ * dphixx;
        
        a = rho;
        b = -rho * xec;
        c = Ip;
        d = -rho * xec;
        
        S1 = (c * e - b * f) /(a * c - b * d); % \ddot{w}(x, t)
        S2 = (a * f - b * e) /(a * c - b * d); % \ddot{phi}(x, t)
        
        %ddy = S1 + ( (i - 1) * dx - 8 * dx ) * ddth2 + 8 * dx * ddth1;
        
        w_control(i, j + 1) = 2 * w_control(i, j) - w_control(i, j - 1) + S1 * dt^2;
        %y_control(i, j + 1) = 2 * y_control(i, j) - y_control(i, j - 1) + ddy * dt^2;
        y_control(i, j + 1) = w_control(i, j + 1) + ( (i - 1) * dx ) * th2_control(j + 1);
        phi_control(i, j + 1) = 2 * phi_control(i, j) - phi_control(i, j - 1) + S2 * dt^2;
    end
    
    %w_control(nx, j + 1) = ( F2(j - 1) * dx^3 * dt /EIb - dt * w_control(nx - 3, j + 1) + 1.5 * dt * w_control(nx - 2, j + 1) - eta * w_control(nx - 3, j + 1) + 1.5 * eta * w_control(nx - 2, j + 1) + 0.5 * eta * w_control(nx, j) + eta * w_control(nx - 3, j) - 1.5 * eta * w_control(nx - 2, j) ) /(0.5 * dt + 0.5 * eta);
%    %w_control(nx, j + 1) = ( F2(j - 1) * dx^3 /EIb + w_control(nx - 3, j + 1) - 1.5 * w_control(nx - 2, j + 1) ) /-0.5;   %？
    
    w_control(nx, j + 1) = ( 1.5 * dt * w_control(nx - 2, j + 1) - dt * w_control(nx - 3, j + 1) + 1.5 * eta * w_control(nx - 1, j + 1) - eta * w_control(nx - 3, j + 1) + 0.5 * eta * w_control(nx, j) - 1.5 * eta * w_control(nx - 2, j) + eta * w_control(nx - 3, j) ) /(0.5 * (dt + eta));
    w_control(nx - 1, j + 1) = ( w_control(nx, j + 1) + w_control(nx - 2, j + 1) ) /2;
    
    y_control(nx, j + 1) = w_control(nx, j + 1) + ( (nx - 1) * dx ) * th2_control(j + 1);
    y_control(nx - 1, j + 1) = w_control(nx - 1, j + 1) + ( (nx - 2) * dx ) * th2_control(j + 1);
    
    % when i=nx-1
    
    ddw = ( w_control(nx - 1, j + 1) - 2 * w_control(nx - 1, j) + w_control(nx - 1, j - 1) ) /dt^2;     % \ddot{w(nx-1,t)}
    dphixx = ( ( phi_control(nx, j) - 2 * phi_control(nx - 1, j) + phi_control(nx - 2, j) ) - ( phi_control(nx, j - 1) - 2 * phi_control(nx - 1, j - 1) + phi_control(nx - 2, j - 1) ) ) /dt /dx^2;    % \dot{phi(nx-1,t)}''
    phixx = ( phi_control(nx, j) - 2 * phi_control(nx - 1, j) + phi_control(nx - 2, j) ) /dx^2; % phi(nx-1,t)''
    
    S4 = -xac * Fb(nx - 1, j) + GJ * phixx + eta * GJ * dphixx + rho * xec * ddw;
    
    phi_control(nx - 1, j + 1) = 2 * phi_control(nx - 1, j) - phi_control(nx - 1, j - 1) + S4 * dt^2 /Ip; % phi(nx-1,t)
    %phi(nx - 1, j + 1) = phi(nx - 2, j + 1);
    
    % when i=nx
    
    %phi_control(nx, j + 1) = M2(j) * dx * dt /GJ * (dt + eta) + (eta /GJ * (dt + eta)) * (phi_control(nx, j) - phi_control(nx - 1, j)) + phi_control(nx - 1, j + 1);
    %phi_control(nx, j + 1) = phi_control(nx - 1, j + 1) + M2(j - 1) * dx /GJ;    %？
    
    phi_control(nx, j + 1) = (-dt * phi_control(nx, j) + dt * phi_control(nx - 1, j) + eta * phi_control(nx - 1, j + 1) + eta * phi_control(nx, j) - eta * phi_control(nx - 1, j) ) /eta;
   
    dth2 = ( th2_control(j + 1) - th2_control(j) ) /dt;
    u2a = -dth2 - ( th2_control(j + 1) - th2_d );
    wxx = ( w_control(3, j + 1) - w_control(2, j + 1) ) /dx^2;
    dwxx = ( ( w_control(3, j + 1) - w_control(2, j + 1) ) - ( w_control(3, j) - w_control(2, j) ) ) /dt /(dx^2);
    F2(j + 1) = -Ih * dth2 + p * u2a - p1 * dth2 - p2 * ( th2_control(j + 1) - th2_d ) - EIb * wxx - eta * EIb * dwxx;
    
    %dphi = ( phi_control(nx, j + 1) - phi_control(nx, j) ) /dt;
    
    dphi = ( phi_control(1, j + 1) - phi_control(1, j) ) /dt;
    M2(j + 1) = -q * ( alpha * dphi + beta * phi_control(1, j + 1) );
    
    if mod(j,nt/Ttr)==0
        w_3D_control(1+j*Ttr/nt,:)=w_control(:,j)';
        phi_3D_control(1+j*Ttr/nt,:)=phi_control(:,j)';
    end
end

w_3D_control(1,:)=w_control(:,1)';
phi_3D_control(1,:)=phi_control(:,1)';



%%
%*************************************
%         Rigid Flexible Wings  
%                Plot
%*************************************
x1 = linspace(0, L2, nx);
t_tr = linspace(0, tmax, Ttr);

figure(1);
surf(x1, t_tr, y_3D_free); view([60 30]);
xlabel('x [m]'); ylabel('t [s]'); zlabel('y(x,t) [m]');
title('Bending displacement y(x,t) without control');

figure(2);
surf(x1, t_tr, w_3D_free); view([60 30]);
xlabel('x [m]'); ylabel('t [s]'); zlabel('w(x,t) [m]');
title('Bending displacement w(x,t) without control');

figure(3);
surf(x1, t_tr, phi_3D_free); view([60 30]);
xlabel('x [m]'); ylabel('t [s]'); zlabel('\phi(x,t) [rad]');
title('Twist displacement \phi(x,t) without control');

figure(4);
plot(linspace(0, tmax, nt), th2_free);
xlabel('t [s]');ylabel('\theta_2 [rad]');zlabel('\theta_2(x,t) [rad]');
title('Twist displacement \theta_2(t) without control');

figure(5);
surf(x1, t_tr, y_3D_control); view([60 30]);
xlabel('x [m]'); ylabel('t [s]'); zlabel('y(x,t) [m]');
title('Bending displacement y(x,t) with control');

figure(6);
surf(x1, t_tr, w_3D_control); view([60 30]);
xlabel('x [m]'); ylabel('t [s]'); zlabel('w(x,t) [m]');
title('Bending displacement w(x,t) with control');

figure(7);
surf(x1, t_tr, phi_3D_control); view([60 30]);
xlabel('x [m]'); ylabel('t [s]'); zlabel('\phi(x,t) [rad]');
title('Twist displacement \phi(x,t) with control');

figure(8);
plot(linspace(0, tmax, nt), th2_control);
xlabel('t [s]');ylabel('\theta_2 [rad]');zlabel('\theta_2(x,t) [rad]');
title('Twist displacement \theta_2(t) with control');

figure(9);
subplot(211)
plot(linspace(0, tmax, nt), F2);
xlabel('t [s]'); ylabel('F2(t) [N]'); title('Control input F2(t)');
subplot(212)
plot(linspace(0, tmax, nt), M2);
xlabel('t [s]'); ylabel('M(t)  [Nm]'); title('Control input M2(t)');

D2 = pi /12 * ones(1, nt); %constraints

figure(11);
plot(linspace(0, tmax, nt), th2_free, ':b', 'Linewidth', 1.5);
hold on
plot(linspace(0, tmax, nt), th2_control, 'k');
hold on
plot(linspace(0, tmax, nt), D2(1,:), '--r');
xlabel('t [s]'); ylabel('\theta_2(t) [m]');
legend('\theta_2(t) without control', '\theta_2(t) with control');
