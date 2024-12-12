function dydt = ODE4(t, y)  
    % y1 = y, y2 = y'  
    % dydt = [y2; y3; y2*y3 - 1]  
    % y''' - y''y + 1 = 0
    dydt = [y(2);-y(1)-t*y(2);]
end
