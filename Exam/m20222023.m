%1
[X, Y] = meshgrid(-5:0.1:5, -5:0.1:5);  
Z = (X+Y)./(X.^2+Y.^4+2); 
surf(X, Y, Z);  
view(3)
%2
A = [89 12 -0.713 0;
21 -45 42 5;
32 96 0 29;
5 -9.54 54 2.14
];
B = [A(2,1) A(2,2);A(3,1) A(2,4);A(4,1) A(3,4)]
C = flipud(B)
D = C'
%3
syms x
f1 = 5*x.^2+32*x+21
f2 = 29*x.^2+5*x-45
P = f1*f2
P = expand(P)
x1=fzero('fun20222023',0) %要先建立函数文件

f = @(x) -(145*x.^4 + 953*x.^3 + 544*x.^2 - 1335*x - 945); 
x_min = fminbnd(f, -2, 1);  
f_min = -f(x_min);  
disp(['最大值点x=',num2str(x_min),'最大值 f(x) = ', num2str(f_min)]);
%4
dataX = simout.data;
timex = simout.time;
plot(timex,dataX) %simulink type

tspan = [0 5];%ode45求解
y0 = [0; 1]; % 初始条件 [y(0); y'(0); y''(0)]  
[t, y] = ode45(@ODE4, tspan, y0);  
figure;  
plot(t, y(:,1), '-o', 'DisplayName', 'y(t)'); % 绘制y(t)  
hold on;  
plot(t, y(:,2), '-x', 'DisplayName', ' dy(t) '); % 绘制y''(t)  
xlabel('时间 t');  
ylabel('函数值');  
title('高阶ODE的数值解');  
legend show;  
grid on;

syms y(t)%符号变量求解
%y''+ty'+y=0 y(0)=0 y'(0)=1
D2y=diff(y,t,2);
Dy=diff(y,t,1);
res=dsolve(D2y+t*Dy+y==0,y(0)==0,Dy(0)==1) %解＝dsolve（函数，初始值）
fplot(res,[0,5])

format short
xk = [0.3 0.5 0.7 0.9 1.1 1.3 1.5]
fxk = [0.3 0.6 0.9 1.1 1.3 1.6 1.8]
X1 = [0.4 0.6 0.8 1.0 1.2 1.4]
X2 = [0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5]
LinearY1=interp1(xk,fxk,X1,'linear')
SplineY1=interp1(xk,fxk,X1,'spline')
[P1,S1]=polyfit(xk,fxk,5)
values_fit = polyval(P1, X1);  
disp(values_fit)

format short
fx = @(x) x.*log(x)  %使用函数句柄
[S,n]=quad(fx,1,exp(1))

% 已知参数

% 符号变量定义
syms Z_T_real Z_T_imag Z_A_real Z_A_imag f L C real;  % 定义实数变量
syms omega real;  % 定义角频率 omega 为实数
% 已知值
Z_T_real = 150;        % 发射机输出阻抗的实部
Z_T_imag = 75;         % 发射机输出阻抗的虚部
Z_A_real = 75;         % 天线输入阻抗的实部
Z_A_imag = 15;         % 天线输入阻抗的虚部
f = 2e9;               % 频率 (Hz)
omega = 2 * pi * f;    % 角频率 (rad/s)
Z_T = Z_T_real + 1i * Z_T_imag;  % 发射机输出阻抗
Z_A = Z_A_real + 1i * Z_A_imag;  % 天线输入阻抗
Z_M = conj(Z_A);                 % 匹配阻抗，与天线输入阻抗共轭
X1_eq = real( (( Z_T*-1j/(omega*C) )/ ( Z_T-1j/(omega*C)))+1j*omega*L ) == real(Z_M);
X2_eq = imag( (( Z_T*-1j/(omega*C) )/ ( Z_T-1j/(omega*C)))+1j*omega*L) == imag(Z_M);
solutions = solve([X1_eq, X2_eq], [L, C]);
L_result = double(solutions.L) * 1e9;   % 将电感结果转换为 nH
C_result = double(solutions.C) * 1e12; % 将电容结果转换为 pF
fprintf('符号解法结果：\n');
fprintf('匹配电感 L: %.2f nH\n', L_result);
fprintf('匹配电容 C: %.2f pF\n', C_result);

% 定义系统函数 H(z)
b = [1 -4 4]; % 分子系数
a = [1 0 0];      % 分母系数（无极点）

% 计算零点和极点
[z, p, k] = tf2zp(b, a); % 零极点计算
fprintf('零点：\n');
disp(z);
fprintf('极点：\n');
disp(p);

% 频率范围
omega = linspace(0, 2*pi, 1000); % 频率范围 [0, 2π]
% 手动计算频率响应 |H(e^jω)| 和 ∠H(e^jω)
magnitude = abs(exp(1j.*omega)-2).*abs(exp(1j.*omega)-2)./abs(exp(1j.*omega))./abs(exp(1j.*omega)) % 幅度响应
phase = angle(exp(1j.*omega)-2)+angle(exp(1j.*omega)-2)-angle(exp(1j.*omega))-angle(exp(1j.*omega));   % 相位响应
% 绘制幅度和相位响应
figure;
subplot(2, 1, 1);
plot(omega./2./pi, magnitude, 'LineWidth', 1.5);
xlabel('\omega (2pi)');
ylabel('|H(e^{j\omega})|');
title('幅度响应');
grid on;

subplot(2, 1, 2);
plot(omega./2./pi, phase./2./pi, 'LineWidth', 1.5);
xlabel('\omega (2pi)');
ylabel('\angleH(e^{j\omega}) (2pi)');
title('相位响应');
grid on;

% 使用 freqz 函数计算频率响应
[H_f, w] = freqz(b, a, 1000, 'whole'); % 'whole' 表示频率范围为 [0, 2π]
% 幅度和相位响应
magnitude_f = abs(H_f); % 幅度响应
phase_f = angle(H_f);   % 相位响应
% 绘制频率响应
figure;
subplot(2, 1, 1);
plot(w./2./pi, magnitude_f, 'LineWidth', 1.5);
xlabel('\omega (2pi)');
ylabel('|H(e^{j\omega})|');
title('使用 freqz 的幅度响应');
grid on;

subplot(2, 1, 2);
plot(w./2./pi, phase_f./2./pi, 'LineWidth', 1.5);
xlabel('\omega (2pi)');
ylabel('\angleH(e^{j\omega}) (2pi)');
title('使用 freqz 的相位响应');
grid on;

