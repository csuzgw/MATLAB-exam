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
x1=fzero('fun20222023',0) %Ҫ�Ƚ��������ļ�

f = @(x) -(145*x.^4 + 953*x.^3 + 544*x.^2 - 1335*x - 945); 
x_min = fminbnd(f, -2, 1);  
f_min = -f(x_min);  
disp(['���ֵ��x=',num2str(x_min),'���ֵ f(x) = ', num2str(f_min)]);
%4
dataX = simout.data;
timex = simout.time;
plot(timex,dataX) %simulink type

tspan = [0 5];%ode45���
y0 = [0; 1]; % ��ʼ���� [y(0); y'(0); y''(0)]  
[t, y] = ode45(@ODE4, tspan, y0);  
figure;  
plot(t, y(:,1), '-o', 'DisplayName', 'y(t)'); % ����y(t)  
hold on;  
plot(t, y(:,2), '-x', 'DisplayName', ' dy(t) '); % ����y''(t)  
xlabel('ʱ�� t');  
ylabel('����ֵ');  
title('�߽�ODE����ֵ��');  
legend show;  
grid on;

syms y(t)%���ű������
%y''+ty'+y=0 y(0)=0 y'(0)=1
D2y=diff(y,t,2);
Dy=diff(y,t,1);
res=dsolve(D2y+t*Dy+y==0,y(0)==0,Dy(0)==1) %�⣽dsolve����������ʼֵ��
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
fx = @(x) x.*log(x)  %ʹ�ú������
[S,n]=quad(fx,1,exp(1))

% ��֪����

% ���ű�������
syms Z_T_real Z_T_imag Z_A_real Z_A_imag f L C real;  % ����ʵ������
syms omega real;  % �����Ƶ�� omega Ϊʵ��
% ��ֵ֪
Z_T_real = 150;        % ���������迹��ʵ��
Z_T_imag = 75;         % ���������迹���鲿
Z_A_real = 75;         % ���������迹��ʵ��
Z_A_imag = 15;         % ���������迹���鲿
f = 2e9;               % Ƶ�� (Hz)
omega = 2 * pi * f;    % ��Ƶ�� (rad/s)
Z_T = Z_T_real + 1i * Z_T_imag;  % ���������迹
Z_A = Z_A_real + 1i * Z_A_imag;  % ���������迹
Z_M = conj(Z_A);                 % ƥ���迹�������������迹����
X1_eq = real( (( Z_T*-1j/(omega*C) )/ ( Z_T-1j/(omega*C)))+1j*omega*L ) == real(Z_M);
X2_eq = imag( (( Z_T*-1j/(omega*C) )/ ( Z_T-1j/(omega*C)))+1j*omega*L) == imag(Z_M);
solutions = solve([X1_eq, X2_eq], [L, C]);
L_result = double(solutions.L) * 1e9;   % ����н��ת��Ϊ nH
C_result = double(solutions.C) * 1e12; % �����ݽ��ת��Ϊ pF
fprintf('���Žⷨ�����\n');
fprintf('ƥ���� L: %.2f nH\n', L_result);
fprintf('ƥ����� C: %.2f pF\n', C_result);

% ����ϵͳ���� H(z)
b = [1 -4 4]; % ����ϵ��
a = [1 0 0];      % ��ĸϵ�����޼��㣩

% �������ͼ���
[z, p, k] = tf2zp(b, a); % �㼫�����
fprintf('��㣺\n');
disp(z);
fprintf('���㣺\n');
disp(p);

% Ƶ�ʷ�Χ
omega = linspace(0, 2*pi, 1000); % Ƶ�ʷ�Χ [0, 2��]
% �ֶ�����Ƶ����Ӧ |H(e^j��)| �� ��H(e^j��)
magnitude = abs(exp(1j.*omega)-2).*abs(exp(1j.*omega)-2)./abs(exp(1j.*omega))./abs(exp(1j.*omega)) % ������Ӧ
phase = angle(exp(1j.*omega)-2)+angle(exp(1j.*omega)-2)-angle(exp(1j.*omega))-angle(exp(1j.*omega));   % ��λ��Ӧ
% ���Ʒ��Ⱥ���λ��Ӧ
figure;
subplot(2, 1, 1);
plot(omega./2./pi, magnitude, 'LineWidth', 1.5);
xlabel('\omega (2pi)');
ylabel('|H(e^{j\omega})|');
title('������Ӧ');
grid on;

subplot(2, 1, 2);
plot(omega./2./pi, phase./2./pi, 'LineWidth', 1.5);
xlabel('\omega (2pi)');
ylabel('\angleH(e^{j\omega}) (2pi)');
title('��λ��Ӧ');
grid on;

% ʹ�� freqz ��������Ƶ����Ӧ
[H_f, w] = freqz(b, a, 1000, 'whole'); % 'whole' ��ʾƵ�ʷ�ΧΪ [0, 2��]
% ���Ⱥ���λ��Ӧ
magnitude_f = abs(H_f); % ������Ӧ
phase_f = angle(H_f);   % ��λ��Ӧ
% ����Ƶ����Ӧ
figure;
subplot(2, 1, 1);
plot(w./2./pi, magnitude_f, 'LineWidth', 1.5);
xlabel('\omega (2pi)');
ylabel('|H(e^{j\omega})|');
title('ʹ�� freqz �ķ�����Ӧ');
grid on;

subplot(2, 1, 2);
plot(w./2./pi, phase_f./2./pi, 'LineWidth', 1.5);
xlabel('\omega (2pi)');
ylabel('\angleH(e^{j\omega}) (2pi)');
title('ʹ�� freqz ����λ��Ӧ');
grid on;

