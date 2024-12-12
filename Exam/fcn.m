function d2x = fcn(dx,x,zeta,Wn)
%#codegen

d2x = -2*zeta*Wn*dx - Wn^2*x;
