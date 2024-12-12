%AABBA
%0.473
%0.021736
%0.8179
%7.7495
%0.7248
%1.27679930
%
syms x;
% 定义函数，例如 f(x) = x^2 + 3x + 2
y = x*cos(x);
% 计算函数 f 关于 x 的导数
df = diff(y, x);
df=matlabFunction(df)
ddf =diff(df, x);
ddf=matlabFunction(ddf)
dddf =diff(ddf, x);
dddf=matlabFunction(dddf)
% 显示结果
disp(ddf(5))
disp(dddf(5))

[t,T]=ode45(@vdp5,[0 2],[1 -0.5]);

fx = @(x) (abs(x)+sin(x)).*x.^2  %使用函数句柄
[S,n]=quad(fx,0,1)

f = @(x) -sin(x+1)./(x.^2+1);  
% 使用 fminbnd 函数找到 [0, 10] 区间内的最小值点  
x_min = fminbnd(f, 0, 10);  
f_min = -f(x_min);  
disp(['最大值 f(x) = ', num2str(f_min)]);

syms y(x)
%y''+2y'+2y=0 y(0)=y'(0)=1
% 符号解/解析解
D2y=diff(y,x,2);
Dy=diff(y,x,1);
res=dsolve(D2y+2*Dy+2*y==0,y(0)==1,Dy(0)==1) %解＝dsolve（函数，初始值）
a=matlabFunction(res); %符号变成一个@（x）的函数

syms n k x
eval(symsum((k+2)/k^2,k,1,50))

P = [5 0 0 0 0 -1]
roots(P)

n = 100;
A = zeros(n, n);
% 填充矩阵 A
for m = 1:n
    for n_col = 1:n
        A(m, n_col) = 1 / (m + 2 * n_col);
    end
end
eigenvalues = eig(A);
max_eigenvalue_magnitude = max(abs(eigenvalues));
% 输出结果
fprintf('矩阵 A 的特征值模的最大值为: %.16f\n', max_eigenvalue_magnitude);



syms u2(t)%符号变量求解
%y''+ty'+y=0 y(0)=0 y'(0)=1
Dy=diff(u2,t,1);
res=dsolve(Dy+2*u2-cos(t)+2*sin(t)==0,u2(0)==0) %解＝dsolve（函数，初始值）
figure
fplot(res, [0,10], 'Color', [1, 0, 0], 'LineStyle', '-.', 'LineWidth', 2, 'Marker', '*');
xlabel('Time (t)');
ylabel('Output u_2(t)');
title('系统的输出 u_2(t)');
grid on;

dataX = simout.data;
timex = simout.time;
figure('Color', 'white'); % 设置背景颜色为白色
plot(timex, dataX, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', 4, 'Marker', '*', 'MarkerSize', 6);
grid on; % 显示网格线

