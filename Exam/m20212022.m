A = [89 12 -0.713 0;21 -45 42 5;32 96 0 29;5 -9.54 54 2.14]
B = [A(2,2) A(2,4);A(3,1) A(4,1)]
C = [fliplr(B(1,:)); B(2,:)]%fliplr(A) flipud(A) 左右对调 上下对调
D = diag(diag(C))
E = inv(D)

%ode45解法
tspan = [0 1]; % 从t=0到t=5  
y0 = 1; % 初始条件y(0)=1  
[t, y] = ode45(@ODE2, tspan, y0);  
plot(t, y, '-o');  
xlabel('时间 t');  
ylabel('函数值 y');  
title('ODE45函数的数值解');  
grid on;
%符号变量
syms y(t)
eqn = diff(y,t,1) == (t^2+1)*y;
cond = [y(0)==1];
ySol(t) = dsolve(eqn,cond)

x_s = 0:0.03:1
y_s = exp((x_s.*(x_s.^2 + 3))./3)
plot(x_s,y_s,'-o')
title('符号变量求解');
xlabel('Time t');
ylabel('Solution y');

%simulink plot 
plot(tout,simout)
xlabel('Time t');
ylabel('Simulink y');

[X, Y] = meshgrid(-5:0.1:5, -5:0.1:5);  
% 定义 Z 数据（例如，一个高斯函数）  
Z = X./(X.^4+Y.^4+1); 
% 绘制三维曲面  
surf(X, Y, Z);  

x = [0.3 0.5 0.7 0.9 1.1 1.3 1.5]
y = [0.3 0.6 0.9 1.1 1.3 1.6 1.8]
I=trapz(x,y)
P1=polyfit(x,y,3)
format short
fx = @(x) 0.3472*x.^3-1.0863*x.^2+2.2465*x-0.2856
[S,n]=quad(fx,0.3,1.5)

P1 = [5 4 3]
P2 = [1 1]
P = conv(P1,P2)
[p,q]=polyder(P1,P2)

% 定义符号变量
syms x;
p1 = 5*x^2 + 4*x + 3;
p2 = x + 1;
p = expand(p1 * p2);
disp('The expression of p(x) = p1(x) * p2(x) is:');
disp(p);
q = p1 / p2;
q_prime = expand(diff(q, x));
disp('The derivative of p1(x) / p2(x) is:');
disp(q_prime);

format short
syms x pi
f = (sqrt(pi)-sqrt(acos(x)))/(sqrt(x+1))
Lim=limit(f,x,-1,'right')
Lim = double(Lim)


A = sqrt(2)
f = 50
T = 1/f
N = 1024
M =16
xmax = 16*T
fs = N/xmax
x = (0:N-1).*(xmax/N);
y = A.*sin(2.*pi.*f.*x);
plot(x,y)
P = mean(y.^2)%计算信号的功率

noise = normrnd(0, sqrt(P/10000), 1, 1024);%生成均值为0，标准差为信号功率的10^-4(40dB  10log10P1/P2)
figure(2) % 创建一个新的图形窗口，编号为2
plot(x,noise); % 绘制噪声
var_noise = var(noise); % 计算噪声的方差(功率)
% 计算信噪比，使用10*log10来将比值转换为分贝（dB）
SNR = 10*log10(P/var_noise);

y_noise = y+noise
plot(x,y_noise)

Y = fft(y); %为加噪声的信号频谱
P2 = abs(Y/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1); %把量纲 大小还原

f = fs/N*(0:(N/2)); %只画一半 幅值恢复正确
plot(f,10*log10(P1),'LineWidth',1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

phi = pi/3
pphi = acos(1.01*cos(phi))
e_phi = pphi - phi

phi = pi/3;
pphi = acos(1.01*cos(phi));
f = 1638400;
N0 = 32768;
dt =1/f/N0;
n = (phi-pphi)/2/pi/f/dt;
n = round(n)

