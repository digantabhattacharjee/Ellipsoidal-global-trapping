%%%%%%%%%%%%%%%%
% Computes global trapping ellipsoids for given shift coordinates
%%%%%%%%%%%%%
function SDP_soln = global_TR_genloss_ell(sys_info,G,opt_specs,chi1,m)
%%%%%%%%
n = sys_info.n; 
Ix = eye(n);
O1 = zeros(n,1);
A = sys_info.A;
d = sys_info.d;
Qm = sys_info.Qm;
% m = sys_info.m;
%%%%%%%%%
pmin = opt_specs.pmin; 
pmax = opt_specs.pmax;
LMIx_tol = opt_specs.LMIx_tol;
solver_indx = opt_specs.solver_indx;
eps_E = opt_specs.eps_E;
%%%%%%%%%%%%
if size(A,2) == 2 
    Q1 = Qm{1};
    Q2 = Qm{2};
    Lm = A + 2*[m'*Q1; m'*Q2];
    cm = d + A*m + [m'*Q1*m; m'*Q2*m];
elseif size(A,2) == 3 
    Q1 = Qm{1};
    Q2 = Qm{2};
    Q3 = Qm{3};
    Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3];
    cm = d + A*m + [m'*Q1*m; m'*Q2*m; m'*Q3*m];
elseif size(A,2) == 9 
    Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
    Q4 = Qm{4}; Q5 = Qm{5}; Q6 = Qm{6};
    Q7 = Qm{7}; Q8 = Qm{8}; Q9 = Qm{9};
    Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3;...
                m'*Q4; m'*Q5; m'*Q6;...
                m'*Q7; m'*Q8; m'*Q9];
    cm = d + A*m + [m'*Q1*m; m'*Q2*m; m'*Q3*m;...
                    m'*Q4*m; m'*Q5*m; m'*Q6*m;...
                    m'*Q7*m; m'*Q8*m; m'*Q9*m];
elseif size(A,2) == 6
    Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
    Q4 = Qm{4}; Q5 = Qm{5}; Q6 = Qm{6};
    Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3;...
                m'*Q4; m'*Q5; m'*Q6];
    cm = d + A*m + [m'*Q1*m; m'*Q2*m; m'*Q3*m;...
                    m'*Q4*m; m'*Q5*m; m'*Q6*m];
end
%%%%=====
cvx_clear;
cvx_begin sdp % quiet
if solver_indx == 1
    cvx_solver sdpt3 
elseif solver_indx == 2
    cvx_solver sedumi  
end
%  
cvx_precision best 
% cvx_precision high
% Define the SDP variables
% variable eps_E %%% tolerance of the energy decrease rate: -eps_E*||y(t)||^2
% eps_E >= eps;
% eps_E = 100*eps; 
variable P(n,n) symmetric;
variable r_s;
%%%%%**************Constraints*******
P >= pmin*Ix; 
P <= pmax*Ix;
% P - pmin*Ix == semidefinite(n);
G'*P(:) == 0;
LMIx = [Lm'*P + P*Lm + chi1*P              P*cm;
           (P*cm)'                      -chi1*r_s];
%%% Desired energy decrease rate: Vdot <= -eps_E*||y(t)||^2
Q = -eps_E*Ix;
LMIx = LMIx -  [Q          O1;
                O1'        0];
%%%%%%%%%%%%%%%
LMIx <= 0;
%%%***** minimize the trapping region radius 
minimize r_s;
cvx_end
%%%%%====== Post processing
S = -P;
if r_s > 0 && eps_E > 0 && chi1 > 0 ...
   && min(eig(P)) >= 0.9999*pmin ...
   && max(eig(LMIx)) <= LMIx_tol % && max(eig(LMIx(1:n,1:n))) <= 0
    SDP_soln.chi1 = chi1;
    SDP_soln.P = P;
    SDP_soln.pmin = pmin;
    SDP_soln.Lm = Lm;
    SDP_soln.m = m;
    SDP_soln.LMIx = LMIx;
   %
    SDP_soln.radius = sqrt(r_s);
    SDP_soln.ROBB = sqrt(r_s)/sqrt(min(eig(P)));
    SDP_soln.S = S;
    SDP_soln.eps_E = eps_E;
    SDP_soln.cvxprec = cvx_precision;
    SDP_soln.cvxstatus = cvx_status;
    % SDP_soln.cvxslvtol = cvx_slvtol;
else
    SDP_soln = [];
end