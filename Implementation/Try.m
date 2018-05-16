%% Try
p = [0 0 0].';
dp = [1 1 1].';
q = [1 0 0 0].';
dq = [0 0 0 0].';
q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);
dq0 = dq(1);
dq1 = dq(2);
dq2 = dq(3);
dq3 = dq(4);
theta = 0;
O_z_B = quatrotate(q.',[0 0 1]);
qt = [cos(theta/2) sin(theta)*O_z_B].';
% System
ddq0 = sym('ddq0','real');
ddq1 = sym('ddq1','real');
ddq2 = sym('ddq2','real');
ddq3 = sym('ddq3','real');
syms ddpx ddpy ddpz
syms ddtheta
ddp = [ddpx ddpy ddpz].';
ddq = [ddq0 ddq1 ddq2 ddq3].';
ddqt = [cos(ddtheta/2)]
% Double derivative of Rotation matrix
R = sym(zeros(3, 1));
R(1,1) = -4*(q2*ddq2 + q3*ddq3 + dq2^2 + dq3^2);
R(1,2) = 2*(ddq1*q2 + 2*dq1*dq2 + q1*ddq2 + ddq0*q3 + + 2*dq0*dq3 + q0*ddq3);
R(1,3) = 2*(ddq1*q3 + 2*dq1*dq3 + q1*ddq3 - ddq0*q2 - q0*ddq2 - 2*dq0*dq2);
R(2,1) = 2*(ddq1*q2 + q1*ddq2 + 2*dq1*dq2 - ddq0*q3 - q0 * ddq3 - 2*dq0*dq3);
R(2,2) = -4*(q1*ddq1 + dq1*dq1 + q3*ddq3 + dq3*dq3);
R(2,3) = 2*(ddq2*q3+q2*ddq3 + 2*dq2*dq3 + ddq0*q1 + q0*ddq1 + 2*dq0*dq1);
R(3,1) = 2*(ddq1*q3 + q1*ddq3 + 2*dq1*dq3 + ddq0*q2 + q0*ddq2 + 2*dq0*dq2);
R(3,2) = 2*(ddq2*q3 + q2*ddq3 + 2*dq2*dq3 - ddq0*q1 - q0*ddq1 - 2*dq0*dq1);
R(3,3) = -4*(q1*ddq1 + dq1^2 + q2*ddq2 + dq2^2);
% Force equation
eqn1 = (m_B + m_C)*ddp + m_C * R * r_BC == F_GC + F_GB + F_T
% Torque equation
Lterm1 = 2*qMatrix(I_B)*(quatmultiply(ddq.',quatconj(q.'))).' 
Lterm2 = 2*qMatrix(I_C)*quatmultiply(quatmultiply(qt.',ddq.'), qSquare(quatmultiply(qt.',q.'))).'
Lterm3 =
