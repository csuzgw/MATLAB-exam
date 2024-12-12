A = [23 10;41 -45;32 5]
B = A'
C = B(:,[1 3])
D = diag(diag(C))
E = inv(D)

f1 = @(x,y) -sqrt(x.^2+y.^2);
fsurf(f1)
hold on
f2 = @(x,y) sqrt(x.^2+y.^2);
fsurf(f2)
hold off

[X, Y] = meshgrid(-1:0.05:1, -1:0.05:1);  
Z = X.^2+Y.^2; 
figure;  
surf(X, Y, Z); 

f = @(x) -x./(1+x.^3);  
% 使用 fminbnd 函数找到 [0, 5] 区间内的最小值点  
x_min = fminbnd(f, 0, 5);  
f_min = -f(x_min);  
disp(['最大值 f(x) = ', num2str(f_min)]);
plot(0:0.1:5,-f(0:0.1:5))
%0.52913

x = [0.3 0.5 0.7 0.9 1.1 1.3 1.5]
y = [0.3 0.6 0.9 1.1 1.3 1.6 1.8]
I=trapz(x,y)
%1.3100
p=polyfit(x,y,3); %前面定义f是个函数 如然后这个地方要输入x序列和y序列 所以要写出f（x）
f = @(x)0.3472*x.^3-1.0863*x.^2+2.2465*x-0.2856
[S,n]=quad(f,0.3,1.5)
%1.3099
x = linspace(0,pi,10)
y = cos(x)
delta1_y = diff(y);
plot(x(1:9),delta1_y)

syms x
y1 = exp(-x)
taylor(y1,x,0,'Order',5)
%x^4/24 - x^3/6 + x^2/2 - x + 1

syms x a
y1 = x^2+a*x+1 ==0
d1 = solve(y1,x)
% - a/2 - ((a - 2)*(a + 2))^(1/2)/2
 %  ((a - 2)*(a + 2))^(1/2)/2 - a/2

%绘制时域波形
dataX = simout.data;
timex = simout.time;
plot(timex,dataX) %simulink type


T=0.5;
N=250;
fs=N/T;
Ts=1/fs;
t=(0:N-1)*Ts;
h=sin(2*pi.*50*t)+0.1*sin(2*pi*50*3.*t)+0.01*sin(2*pi*50*5.*t);
Y=2*fft(h)/N;
f=(0:N/2)*fs/N;
amplitude=20*log10(((abs(Y(1:N/2+1)))));
plot(f,amplitude)
xlabel('Frequency')
ylabel('Amplitude(dB)')

Fs = 600;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 120;             % Length of signal %采样点数要是信号周期的整数倍 2*pi*f/fs=数字角频率 T=2*pi/数字角频率
t = (0:L-1)*T;        % Time vector

S = sin(2*pi.*50*t)+0.1*sin(2*pi*50*3.*t)+0.01*sin(2*pi*50*5.*t);

Y = fft(S); %为加噪声的信号频谱
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1); %把量纲 大小还原

f = Fs/L*(0:(L/2)); %只画一半 幅值恢复正确
stem(f,abs(P1),'LineWidth',1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|(dB)')
