function dydt = vdp5(t,y)
%dxdt = 2y dydt = -0.5x
% y = y(1)
% y(1)'=y(2)
%   无论怎么说先把变量替换成y与t的一个函数 再把y换成y1 y1'换成y2
%   y2'等于y3或者等于右边那个式子
dydt = [y(2); -y(1)]; %y(1)的导数，y(2)的导数
end
