function [g, I_B, I_C, m_B, m_C, r_BC, r_CB, F_GC, F_GB, F_T] = getConstants()
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
g = -9.8;
I_B = diag([1,1,1]);
I_C = diag([1,1,1]);
m_B = 10; % in grams
m_C = 8; 
r_BC = [15 0 10].';
r_CB = -r_BC;
F_GC = [0, 0, m_C * g].';
F_GB = [0, 0, m_B * g].';
F_T = [0, 0, 10].'; % IDK
end

