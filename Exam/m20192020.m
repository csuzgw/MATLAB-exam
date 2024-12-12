x = randi([1000, 1100]);
is_prime = isprime(x);
fprintf('������� x = %d, �Ƿ�Ϊ����: %d\n', x, is_prime);

p1 = [2 4 0 5]; % p1(x) = 2x^3 + 4x^2 + 5
p2 = [1 2];     % p2(x) = x + 2
p3 = [0 0 1 2 3];     % p3(x) = x^2 +2x+ 3

% �� p(x) = (p1 * p2) + p3
p = conv(p1, p2);
p_result = p+p3

% ��ʾ���
disp('p(x)�ı��ʽ: ');
disp(poly2sym(p_result));

% �� p(x) = 0 �ĸ�
roots_p = roots(p_result);
disp('p(x) = 0 �ĸ�: ');
disp(roots_p);

[x, y] = meshgrid(-5:0.1:5, -5:0.1:5);
f = (x.^4 + y.^4) ./ (x.^2 + y.^2 + 2);

figure;
surf(x, y, f);
xlabel('x');
ylabel('y');
zlabel('f(x, y)');
title('��ά����ͼ��');

A = [10 5 23;89 2 -35];
B = A';

% ���ɷ��� C
C = B([1, 3], :);

% �γɶԽ��� D������������� E
D = diag(diag(C));
E = inv(D);

disp('���� A:'); disp(A);
disp('���� B:'); disp(B);
disp('���� C:'); disp(C);
disp('���� D:'); disp(D);
disp('���� E:'); disp(E);

f = @(x) x.^5 + x.^3 + 1;%�ɼӿɲ��ӵ�
x0 = -1;
x_sol = fzero(f, x0);
fprintf('���̵���ֵ��: x = %f\n', x_sol);
est=f(-0.837620)

f = @(x) -(x + 1) ./ (x.^3 + 1);
x = linspace(0, 2, 100);
y = -f(x);

% ʹ�� fminbnd �����ҵ� [0, 10] �����ڵ���Сֵ��  
x_min = fminbnd(f, 0, 2);  
f_min = -f(x_min);  
disp(['���ֵ f(x) = ', num2str(f_min)]);

figure;
plot(x, y);
xlabel('x');
ylabel('f(x)');
title('�������ߺ����ֵ');

x = linspace(0, pi, 10);
y = cos(x);
dy_dx = diff(y) ./ diff(x);

% ��ͼ
figure;
plot(x(1:end-1), dy_dx, '-o');
xlabel('x');
ylabel('һ����ֵ΢��');
title('cos(x)��10�����ϵ�һ����ֵ΢��');

syms x;
y = exp(x);
taylor_expansion = taylor(y, x, 'Order', 5);
disp('̩��չ��ʽ:');
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
xlabel('fƵ�ʣ�Hz��')
ylabel('abs(Zd)')
title('���ض�·ʱ��Zd�ķ�ֵ��Ƶ�ʱ仯')
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
title('���ض�·ʱ��Zd�ķ�ֵ��d/lamada�仯')
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
xlabel('fƵ�ʣ�Hz��')
ylabel('����d��m��')
zlabel('abs(Zd)')
title('Zd�ķ�ֵ��Ƶ�ʺ;���仯')


%��Զ�Ĵ���
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

