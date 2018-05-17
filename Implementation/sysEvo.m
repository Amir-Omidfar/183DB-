function sol = sysEvo(p, dp, q, dq, theta, dtheta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Constants
[g, I_B, I_C, m_B, m_C, r_BC, r_CB, F_GC, F_GB, F_T] = getConstants();
% Unpack
q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);
dq0 = dq(1);
dq1 = dq(2);
dq2 = dq(3);
dq3 = dq(4);
theta = 0;
dtheta = 0.5;
O_z_B = quatrotate(q.',[0 0 1]);
qt = [cos(theta/2) sin(theta)*O_z_B].';
% Sym variable
ddq0 = sym('ddq0','real');
ddq1 = sym('ddq1','real');
ddq2 = sym('ddq2','real');
ddq3 = sym('ddq3','real');
syms ddpx ddpy ddpz
ddtheta = sym('ddtheta', 'real');
dR = sym(zeros(3,1));
ddR = sym(zeros(3, 1));

% Rotation Matrix
dR(1,1) = -4*(q2*dq2 + q3*dq3);
dR(1,2) = 2*(dq1*q2 + q1*dq2 + dq0*q3 + q0*dq3);
dR(1,3) = 2*(dq1*q3 + q1*dq3 - dq0*q2 - q0*dq2);
dR(2,1) = 2*(dq1*q2 + q1*dq2 - dq0*q3 - q0*dq3);
dR(2,2) = -4*(q1*dq1 + q3*dq3);
dR(2,3) = 2*(dq2*q3 + q2*dq3 + dq0*q1 + q0*dq1);
dR(3,1) = 2*(dq1*q3 + q1*dq3 + dq0*q2 + q0*dq2);
dR(3,2) = 2*(dq2*q3 + q2*dq3 - dq0*q1 - q0*dq1);
dR(3,3) = -4*(q1*dq1 + q2*dq2);

ddR(1,1) = -4*(q2*ddq2 + q3*ddq3 + dq2^2 + dq3^2);
ddR(1,2) = 2*(ddq1*q2 + 2*dq1*dq2 + q1*ddq2 + ddq0*q3 + + 2*dq0*dq3 + q0*ddq3);
ddR(1,3) = 2*(ddq1*q3 + 2*dq1*dq3 + q1*ddq3 - ddq0*q2 - q0*ddq2 - 2*dq0*dq2);
ddR(2,1) = 2*(ddq1*q2 + q1*ddq2 + 2*dq1*dq2 - ddq0*q3 - q0 * ddq3 - 2*dq0*dq3);
ddR(2,2) = -4*(q1*ddq1 + dq1*dq1 + q3*ddq3 + dq3*dq3);
ddR(2,3) = 2*(ddq2*q3+q2*ddq3 + 2*dq2*dq3 + ddq0*q1 + q0*ddq1 + 2*dq0*dq1);
ddR(3,1) = 2*(ddq1*q3 + q1*ddq3 + 2*dq1*dq3 + ddq0*q2 + q0*ddq2 + 2*dq0*dq2);
ddR(3,2) = 2*(ddq2*q3 + q2*ddq3 + 2*dq2*dq3 - ddq0*q1 - q0*ddq1 - 2*dq0*dq1);
ddR(3,3) = -4*(q1*ddq1 + dq1^2 + q2*ddq2 + dq2^2);
% Combine Sym variable
ddp = [ddpx ddpy ddpz].';
ddq = [ddq0 ddq1 ddq2 ddq3].';
dqt = [-1/2*sin(theta/2)*dtheta, 1/2*cos(theta/2)*dtheta*O_z_B+sin(theta/2)*(dR*[0 0 1].')'].';
a = -1/2 * (sin(theta/2)*ddtheta + 1/2 * cos(theta/2)*(dtheta)^2);
b = 1/2 * ((cos(theta/2)*ddtheta - 1/2 * sin(theta/2)*(dtheta)^2)*O_z_B + 1/2 * cos(theta/2) * dtheta * (dR*[0 0 1].').');
c = sin(theta/2)*(ddR*[0 0 1].').' + 1/2 * cos(theta/2)*dtheta*(dR*[0 0 1].').';
ddqt = [a b+c].';
% Force equation
eqn1 = (m_B + m_C)*ddp + m_C * ddR * r_BC == F_GC + F_GB + F_T
% Torque equation
Lterm1 = 2*qMatrix(I_B)*(quatmultiply(ddq.',quatconj(q.'))).';
Lterm2 = 2*qMatrix(I_C)*quatmultiply(quatmultiply(qt.',ddq.'), qSquare(quatmultiply(qt.',q.'))).';
Lterm3 = 2*qMatrix(I_C)*(quatmultiply((quatmultiply(ddqt.',q.')), quatconj(quatmultiply(qt.',q.')))).';
Lterm4 = [0; (cross(quatrotate(q.', r_CB.'), m_B*ddp.'-F_GB.'-quatrotate(q.',F_T.'))).'];
Rterm1 = 2*qMatrix(I_B)*(qSquare(quatmultiply(dq.',quatconj(q.')))).';
Rterm2 = 2*qMatrix(I_C)*(qSquare(quatmultiply(quatmultiply(qt.', q.') + quatmultiply(dqt.', q.'), quatconj(quatmultiply(qt.', q.'))))).';
Rterm3 = 4*qMatrix(I_C)*(quatmultiply(quatmultiply(dqt.', dq.'), quatconj(quatmultiply(qt.',q.')))).';
eqn2 = Lterm1 + Lterm2 + Lterm3 + Lterm4 == Rterm1 + Rterm2 - Rterm3
% Unit quaternion second derivative constraints
eqn3 = q0*ddq0 + q1*ddq1 + q2*ddq2 + q3*ddq3 + dq0^2 + dq1^2 + dq2^2 + dq3^2 == 0
sol = solve([eqn1, eqn2, eqn3], [ddpx, ddpy, ddpz, ddq0, ddq1, ddq2, ddq3, ddtheta])
end

