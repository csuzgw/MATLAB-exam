% ����һ����������������ODE  
function dydt = ODE2(t, y)  
    % �˴���дODE���Ҳ���ʽ�����磺  
    dydt = (t^2+1)*y; % ʾ�����̣�dy/dt = (t^2+1)*y y(0)=1 
end  
