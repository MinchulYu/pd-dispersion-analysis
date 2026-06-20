function [k] = Peri_const_mk_2d(xi,yi,xj,yj,deltai,voli,volj,micro_modulus,delta_x)
%--------------------------------------------------------------------------
% Local stiffness contribution for 2D bond-based peridynamics.
%
% Input:
%   (xi,yi)        : coordinates of the center node
%   (xj,yj)        : coordinates of the family node
%   deltai         : horizon radius
%   voli, volj     : nodal volumes of the center and family nodes
%   micro_modulus  : bond-based PD micro-modulus
%   delta_x        : nodal spacing
%
% Output:
%   k              : 2-by-4 local stiffness contribution associated with
%                    [u_i, u_j, v_i, v_j]
%
% Note:
%   This code uses the distance-dependent influence function
%   w = 1 - |xi|/delta in Peri_const_mk_2d.m.
%   Therefore, m = 1 is generally not recommended because the nearest
%   neighbor bonds lie on the horizon boundary, where w becomes zero.
%   Use m >= 2 for the dispersion results reported in the manuscript.
%--------------------------------------------------------------------------

dx = xj - xi;
dy = yj - yi;

absxi = sqrt(dx^2 + dy^2);

if absxi < eps
    error('The bond length is zero. The center node and family node must be different.');
end

half_deltax = delta_x/2;

% Partial-volume correction near the horizon boundary
if absxi > deltai + half_deltax
    scalefac = 0;
elseif absxi >= deltai - half_deltax
    scalefac = (deltai + half_deltax - absxi)/delta_x;
else
    scalefac = 1;
end

% Distance-dependent influence function
w = 1 - absxi/deltai;
c_bond = micro_modulus*w;

A = [dx^2 dx*dy;
     dy*dx dy^2];

B = [1 -1 0  0;
     0  0 1 -1];

C = c_bond/absxi^3*A*B;

k = C*voli*volj*scalefac;
