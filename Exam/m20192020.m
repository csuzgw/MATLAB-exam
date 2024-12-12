x = randi([1000, 1100]);
is_prime = isprime(x);
fprintf('随机整数 x = %d, 是否为质数: %d\n', x, is_prime);

p1 = [2 4 0 5]; % p1(x) = 2x^3 + 4x^2 + 5
p2 = [1 2];     % p2(x) = x + 2
p3 = [0 0 1 2 3];     % p3(x) = x^2 +2x+ 3

% 求 p(x) = (p1 * p2) + p3
p = conv(p1, p2);
p_result = p+p3

% 显示结果
disp('p(x)的表达式: ');
disp(poly2sym(p_result));

% 求 p(x) = 0 的根
roots_p = roots(p_result);
disp('p(x) = 0 的根: ');
disp(roots_p);

[x, y] = meshgrid(-5:0.1:5, -5:0.1:5);
f = (x.^4 + y.^4) ./ (x.^2 + y.^2 + 2);

figure;
surf(x, y, f);
xlabel('x');
ylabel('y');
zlabel('f(x, y)');
title('三维曲面图形');

A = [10 5 23;89 2 -35];
B = A';

% 构成方阵 C
C = B([1, 3], :);

% 形成对角阵 D，并求其逆矩阵 E
D = diag(diag(C));
E = inv(D);

disp('矩阵 A:'); disp(A);
disp('矩阵 B:'); disp(B);
disp('矩阵 C:'); disp(C);
disp('矩阵 D:'); disp(D);
disp('矩阵 E:'); disp(E);

f = @(x) x.^5 + x.^3 + 1;%可加可不加点
x0 = -1;
x_sol = fzero(f, x0);
fprintf('方程的数值解: x = %f\n', x_sol);
est=f(-0.837620)

f = @(x) -(x + 1) ./ (x.^3 + 1);
x = linspace(0, 2, 100);
y = -f(x);

% 使用 fminbnd 函数找到 [0, 10] 区间内的最小值点  
x_min = fminbnd(f, 0, 2);  
f_min = -f(x_min);  
disp(['最大值 f(x) = ', num2str(f_min)]);

figure;
plot(x, y);
xlabel('x');
ylabel('f(x)');
title('函数曲线和最大值');

x = linspace(0, pi, 10);
y = cos(x);
dy_dx = diff(y) ./ diff(x);

% 绘图
figure;
plot(x(1:end-1), dy_dx, '-o');
xlabel('x');
ylabel('一阶数值微分');
title('cos(x)在10个点上的一阶数值微分');

syms x;
y = exp(x);
taylor_expansion = taylor(y, x, 'Order', 5);
disp('泰勒展开式:');
disp(taylor_expansion);


%(1) ZL = 0 
syms f d ZL
Z0 = 50;
ZL = 0
d = 0.1
vp = 2e8
beta = 2*pi*f/2e8
SigmaL = (ZL - Z0)./(ZL + Z0)
Sigmad = SigmaL.*exp(-2.*j.*beta.*d)
Zd = Z0.*(1+Sigmad)./(1-Sigmad);
figure(1)
fplot(abs(Zd),[0 4e9])
xlabel('f频率（Hz）')
ylabel('abs(Zd)')
title('负载短路时的Zd的幅值随频率变化')
%(2)ZL = 0 d 0--lamada
syms f d ZL dlamada
Z0 = 50;
ZL = 0
vp = 2e8
SigmaL = (ZL - Z0)./(ZL + Z0)
Sigmad = SigmaL.*exp(-2.*j.*2*pi.*(dlamada))
Zd = Z0.*(1+Sigmad)./(1-Sigmad);
figure(2)
fplot(abs(Zd),[0 1])
xlabel('d/lamada')
ylabel('abs(Zd)')
title('负载短路时的Zd的幅值随d/lamada变化')
%(3)
syms f d ZL
Z0 = 50;
ZL = 25
vp = 2e8
SigmaL = (ZL - Z0)./(ZL + Z0)
Sigmad = SigmaL.*exp(-2.*j.*(2.*pi.*f./2e8).*d)
Zd = Z0.*(1+Sigmad)./(1-Sigmad);
figure(3)
fsurf(abs(Zd),[0 4e9 0.1 0.2])
xlabel('f频率（Hz）')
ylabel('距离d（m）')
zlabel('abs(Zd)')
title('Zd的幅值随频率和距离变化')


%致远的代码
syms Zl f d
Z0=50;
vp=2e8;
% 1)
d=0.1;
sigma0=-1;
beta=2*pi*f/vp;
Td_1=sigma0*exp(-1i.*2.*beta.*d);
Zd_1=Z0.*(1+Td_1)./(1-Td_1);
figure(1)
fplot(abs(Zd_1),[0 4e9])
% 2)
Zd_2=Z0*(1+sigma0*exp(-2*1i*2*pi*d))./(1-sigma0*exp(-2*1i*2*pi*d));
figure(2)
fplot(abs(Zd_2),[0 1])
% 3)
Tl=(25-Z0)/(25+Z0);
Td_3=Tl.*exp(-1i*2*2*pi*(f./vp).*d);
Zd_3=Z0*(1+Td_3)./(1-Td_3);
figure(3)
fsurf(abs(Zd_3),[0.1 0.2 0 4e9])

