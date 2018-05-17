%% Try 2
p = [0 0 0].';
dp = [1 1 1].';
q = [1 0 0.021 -0.007].';
dq = [0 0 0 0].';
theta = 0;
dtheta = 0;
sol = sysEvo(p, dp, q, dq, theta, dtheta)
double(sol.ddpx)
double(sol.ddpy)
double(sol.ddpz)
double(sol.ddq0)
double(sol.ddq1)
double(sol.ddq2)
double(sol.ddq3)
double(sol.ddtheta)
