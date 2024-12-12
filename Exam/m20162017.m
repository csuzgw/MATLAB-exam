%AABBA
%0.473
%0.021736
%0.8179
%7.7495
%0.7248
%1.27679930
%
syms x;
% ���庯�������� f(x) = x^2 + 3x + 2
y = x*cos(x);
% ���㺯�� f ���� x �ĵ���
df = diff(y, x);
df=matlabFunction(df)
ddf =diff(df, x);
ddf=matlabFunction(ddf)
dddf =diff(ddf, x);
dddf=matlabFunction(dddf)
% ��ʾ���
disp(ddf(5))
disp(dddf(5))

[t,T]=ode45(@vdp5,[0 2],[1 -0.5]);

fx = @(x) (abs(x)+sin(x)).*x.^2  %ʹ�ú������
[S,n]=quad(fx,0,1)

f = @(x) -sin(x+1)./(x.^2+1);  
% ʹ�� fminbnd �����ҵ� [0, 10] �����ڵ���Сֵ��  
x_min = fminbnd(f, 0, 10);  
f_min = -f(x_min);  
disp(['���ֵ f(x) = ', num2str(f_min)]);

syms y(x)
%y''+2y'+2y=0 y(0)=y'(0)=1
% ���Ž�/������
D2y=diff(y,x,2);
Dy=diff(y,x,1);
res=dsolve(D2y+2*Dy+2*y==0,y(0)==1,Dy(0)==1) %�⣽dsolve����������ʼֵ��
a=matlabFunction(res); %���ű��һ��@��x���ĺ���

syms n k x
eval(symsum((k+2)/k^2,k,1,50))

P = [5 0 0 0 0 -1]
roots(P)

n = 100;
A = zeros(n, n);
% ������ A
for m = 1:n
    for n_col = 1:n
        A(m, n_col) = 1 / (m + 2 * n_col);
    end
end
eigenvalues = eig(A);
max_eigenvalue_magnitude = max(abs(eigenvalues));
% ������
fprintf('���� A ������ֵģ�����ֵΪ: %.16f\n', max_eigenvalue_magnitude);



syms u2(t)%���ű������
%y''+ty'+y=0 y(0)=0 y'(0)=1
Dy=diff(u2,t,1);
res=dsolve(Dy+2*u2-cos(t)+2*sin(t)==0,u2(0)==0) %�⣽dsolve����������ʼֵ��
figure
fplot(res, [0,10], 'Color', [1, 0, 0], 'LineStyle', '-.', 'LineWidth', 2, 'Marker', '*');
xlabel('Time (t)');
ylabel('Output u_2(t)');
title('ϵͳ����� u_2(t)');
grid on;

dataX = simout.data;
timex = simout.time;
figure('Color', 'white'); % ���ñ�����ɫΪ��ɫ
plot(timex, dataX, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', 4, 'Marker', '*', 'MarkerSize', 6);
grid on; % ��ʾ������

