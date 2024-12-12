[x, y] = meshgrid(-5:0.1:5, -5:0.1:5);
z = x.^2./(x.^2+y.^4+5)
surf(x, y, z);
xlabel('x');
ylabel('y');
zlabel('z');
title('3D Surface Plot');

f = @(x) -sin(x).*x.^2;
% ʹ�� fminbnd �����ҵ� [0, 10] �����ڵ���Сֵ��  
x_min = fminbnd(f, 0, 2*pi);  
f_min = -f(x_min);  
disp(['���ֵ��x=',num2str(x_min),'���ֵ f(x) = ', num2str(f_min)]);
plot(0:0.1:2*pi,-f(0:0.1:2*pi));

A = [15 12 87; 32 26 55];
B = A';
C = [A(1,3) A(2,1); B(2,:)]; % ע�⣺�������A(1,3)��A�ĵ�һ�е�����Ԫ�أ�B(2,:)��B�ĵڶ���
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


f = 50; % Ƶ��
A = 1; % ��ֵ
Fs = 750; % ������
t = 0:1/Fs:0.02; % ����ʱ��
x = A*sin(2*pi*f*t);
% �ҵ���һ����������������������
for i = 1:16
 if x(i)>0&&x(i+1)<0
    break
 end
    
end
t1 = t(i);
t2 = t(i+1);
y1 = A*sin(2*pi*f*t1)
y2 = A*sin(2*pi*f*t2)
% ���Բ�ֵ
t_zero=(((y2-y1)/(t2-t1)*t2)-y2)*((t2-t1)/(y2-y1))
theoretical_zero = 1/(2*f);
relative_error = abs(t_zero - theoretical_zero) / theoretical_zero;
disp(['Relative Error: ', num2str(relative_error)]);


% ����ȵ���ֵ��ʵ�����͵ȵ翹ֵ���鲿���ļ���
r_values = [0, 0.5, 1, 1.5, 2];
x_values = [-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2];

% ����һ�� 0 �� 2�� �ĽǶ����飬���ڻ���Բ
theta = linspace(0, 2 * pi, 100);

% ��ʼ��ͼ
figure;
hold on;
axis square;
xlim([-1 1]);
ylim([-1 1]);
xlabel('ʵ�� (\Gamma_r)');
ylabel('�鲿 (\Gamma_i)');
title('�ȵ���Բ�͵ȵ翹Բ');

% ���Ƶȵ���Բ
for r = r_values
    gamma_i = 1 / (1 + r);
    plot(gamma_i * cos(theta)+(r/(1+r)), gamma_i * sin(theta), 'b');
end

for x = x_values
    gamma_x = 1/x;
    plot(1+gamma_x* cos(theta), (1/x)+gamma_x*sin(theta), 'r');
end
hold on
% ��ʷ��˹ͼ�ϱ���ض��㣬���� z = 1 + j ��λ��
z = 1 + 1i;
gamma_z = (z - 1) / (z + 1);
plot(real(gamma_z), imag(gamma_z), 'go', 'MarkerSize', 10, 'DisplayName', 'z = 1 + j');

% ����õ��ģ�ͽǶ�
mag_gamma_z = abs(gamma_z);
angle_gamma_z = angle(gamma_z);
fprintf('����ϵ����ģ: %f\n', mag_gamma_z);
fprintf('����ϵ���ĽǶ�: %f ����\n', angle_gamma_z);
hold on
% ���ƴ�ԭ�㵽����ϵ������߶�
plot([0, real(gamma_z)], [0, imag(gamma_z)], 'k-', 'DisplayName', '����ϵ������');

% ͼ��
legend({'�ȵ���Բ', '�ȵ翹Բ', 'z = 1 + j', '����ϵ������'}, 'Location', 'best');
grid on
hold off;
