function [k] = Peri_const_mk_1d(xi,xj,deltai,voli,volj,c,delta_x)

absxi = sqrt((xj-xi)^2);
half_deltax = delta_x/2;

if abs(absxi-deltai) <= half_deltax
    scalefac = (deltai-absxi)/(2*half_deltax) + 1/2;
%     scalefac = 1;
else
    scalefac = 1;
end
% c = c*(1-absxi/deltai);

A = (xj-xi)^2;
B = [1 -1];
C = c/absxi^3*A*B;

k = C*voli*volj*scalefac;



