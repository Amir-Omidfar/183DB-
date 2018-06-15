%% Try 2
p = [0 0 0].';
dp = [1 1 1].';
q = [1 0 0 0].';
dq = [0 0 0 0].';
theta = 0;
dtheta = 0.1;
sol = sysEvo(p, dp, q, dq, theta, dtheta)
double(sol.ddpx)
double(sol.ddpy)
double(sol.ddpz)
double(sol.ddq0)
double(sol.ddq1)
double(sol.ddq2)
double(sol.ddq3)
double(sol.ddtheta)

%% Case 2
p = [0 0 0].';
dp = [0 0 0].';
q = [0.1 0.2 0.3 0.9273].';
dq = [0.1 0.2 0.3 -0.15098].';
theta = 0;
dtheta = 0.1;
sol = sysEvo(p, dp, q, dq, theta, dtheta)
double(sol.ddpx)
double(sol.ddpy)
double(sol.ddpz)
double(sol.ddq0)
double(sol.ddq1)
double(sol.ddq2)
double(sol.ddq3)
double(sol.ddtheta)

