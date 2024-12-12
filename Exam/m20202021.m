[x, y] = meshgrid(-5:0.1:5, -5:0.1:5);
z = x.^2./(x.^2+y.^4+5)
surf(x, y, z);
xlabel('x');
ylabel('y');
zlabel('z');
title('3D Surface Plot');

f = @(x) -sin(x).*x.^2;
% 使用 fminbnd 函数找到 [0, 10] 区间内的最小值点  
x_min = fminbnd(f, 0, 2*pi);  
f_min = -f(x_min);  
disp(['最大值点x=',num2str(x_min),'最大值 f(x) = ', num2str(f_min)]);
plot(0:0.1:2*pi,-f(0:0.1:2*pi));

A = [15 12 87; 32 26 55];
B = A';
C = [A(1,3) A(2,1); B(2,:)]; % 注意：这里假设A(1,3)是A的第一行第三个元素，B(2,:)是B的第二行
D = diag(diag(C));
E = inv(D);
disp('Matrix A:');
disp(A);
disp('Matrix B:');
disp(B);
disp('Matrix C:');
disp(C);
disp('Matrix D:');
disp(D);
disp('Matrix E:');
disp(E);

p1 = @(x) x^3 + 2*x^2 + x;
p2 = @(x) x + 1;
p3 = @(x) 5*x^2 +4*x+3;
p = @(x) p1(x) * p2(x) + p3(x);
syms x;
p_sym = expand(p(x));
disp('Expression of p(x):');
disp(p_sym);

x = linspace(0, pi, 10);
y = cos(x);
dy = diff(y)./diff(x);
plot(x(1:end-1), dy);
xlabel('x');
ylabel('dy/dx');
title('Numerical Differentiation');

syms x;
f = sin(x);
T = taylor(f, 'ExpansionPoint', 0, 'Order', 6);
disp('Taylor Expansion:');
disp(T);


f = 50; % 频率
A = 1; % 幅值
Fs = 750; % 采样率
t = 0:1/Fs:0.02; % 采样时间
x = A*sin(2*pi*f*t);
% 找到第一个过零点两侧的两个采样点
for i = 1:16
 if x(i)>0&&x(i+1)<0
    break
 end
    
end
t1 = t(i);
t2 = t(i+1);
y1 = A*sin(2*pi*f*t1)
y2 = A*sin(2*pi*f*t2)
% 线性插值
t_zero=(((y2-y1)/(t2-t1)*t2)-y2)*((t2-t1)/(y2-y1))
theoretical_zero = 1/(2*f);
relative_error = abs(t_zero - theoretical_zero) / theoretical_zero;
disp(['Relative Error: ', num2str(relative_error)]);


% 定义等电阻值（实部）和等电抗值（虚部）的集合
r_values = [0, 0.5, 1, 1.5, 2];
x_values = [-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2];

% 创建一个 0 到 2π 的角度数组，用于绘制圆
theta = linspace(0, 2 * pi, 100);

% 开始绘图
figure;
hold on;
axis square;
xlim([-1 1]);
ylim([-1 1]);
xlabel('实部 (\Gamma_r)');
ylabel('虚部 (\Gamma_i)');
title('等电阻圆和等电抗圆');

% 绘制等电阻圆
for r = r_values
    gamma_i = 1 / (1 + r);
    plot(gamma_i * cos(theta)+(r/(1+r)), gamma_i * sin(theta), 'b');
end

for x = x_values
    gamma_x = 1/x;
    plot(1+gamma_x* cos(theta), (1/x)+gamma_x*sin(theta), 'r');
end
hold on
% 在史密斯图上标出特定点，例如 z = 1 + j 的位置
z = 1 + 1i;
gamma_z = (z - 1) / (z + 1);
plot(real(gamma_z), imag(gamma_z), 'go', 'MarkerSize', 10, 'DisplayName', 'z = 1 + j');

% 输出该点的模和角度
mag_gamma_z = abs(gamma_z);
angle_gamma_z = angle(gamma_z);
fprintf('反射系数的模: %f\n', mag_gamma_z);
fprintf('反射系数的角度: %f 弧度\n', angle_gamma_z);
hold on
% 绘制从原点到反射系数点的线段
plot([0, real(gamma_z)], [0, imag(gamma_z)], 'k-', 'DisplayName', '反射系数向量');

% 图例
legend({'等电阻圆', '等电抗圆', 'z = 1 + j', '反射系数向量'}, 'Location', 'best');
grid on
hold off;
