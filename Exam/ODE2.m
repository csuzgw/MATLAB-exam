% 定义一个函数，用于描述ODE  
function dydt = ODE2(t, y)  
    % 此处填写ODE的右侧表达式，例如：  
    dydt = (t^2+1)*y; % 示例方程，dy/dt = (t^2+1)*y y(0)=1 
end  
