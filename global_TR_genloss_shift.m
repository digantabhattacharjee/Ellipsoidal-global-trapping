%%%%%%%%%%%%%%%%
% Computes the shift coordinate perturbations for a given 
% (global) trapping ellipsoid and shift coordinates(m0)
% Requires linearizing the nonlinear matrix inequality (trapping condition) 
%%%%%%%%%%%%%
function SDP_soln = global_TR_genloss_shift(sys_info,opt_specs,chi1,m0,P)
%%%%%%%%%%%%
n = sys_info.n; 
Ix = eye(n);
O1 = zeros(n,1);
A = sys_info.A;
d = sys_info.d;
Qm = sys_info.Qm;
% m = sys_info.m;
sys_indx = sys_info.sys_indx;
%%%%%%%%%
% pmin = opt_specs.pmin; 
LMIx_tol = opt_specs.LMIx_tol;
solver_indx = opt_specs.solver_indx;
eps_E = opt_specs.eps_E;
w_deltam = opt_specs.w_deltam;
%%%%=====
cvx_clear;
cvx_begin sdp % quiet
if solver_indx == 1
    cvx_solver sdpt3 
elseif solver_indx == 2
    cvx_solver sedumi  
end
%  
cvx_precision best %high
% Define the SDP variables
% variable eps_E %%% tolerance of the energy decrease rate: -eps_E*||y(t)||^2
% eps_E >= eps;
% eps_E = 100*eps; 
% %%%%**************Constraints*******
variable Delta_m(n,1)
m = m0 + Delta_m;
variable r_s

if size(A,2) == 2 %sys_indx == 1 || sys_indx == 5 || sys_indx == 6
    Q1 = Qm{1};
    Q2 = Qm{2};
    Lm = A + 2*[m'*Q1; m'*Q2];
    % cm = d + A*m + [m'*Q1*m; m'*Q2*m];
    cm = d + A*m + [(m0'*Q1*m0 + 2*m0'*Q1*Delta_m); 
                    (m0'*Q2*m0 + 2*m0'*Q2*Delta_m)];
elseif size(A,2) == 3 %sys_indx == 2 || sys_indx == 3 || sys_indx == 7
    Q1 = Qm{1};
    Q2 = Qm{2};
    Q3 = Qm{3};
    Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3];
    % cm = d + A*m + [m'*Q1*m; m'*Q2*m; m'*Q3*m];
    cm = d + A*m + [(m0'*Q1*m0 + 2*m0'*Q1*Delta_m); 
                    (m0'*Q2*m0 + 2*m0'*Q2*Delta_m);
                    (m0'*Q3*m0 + 2*m0'*Q3*Delta_m)];
elseif size(A,2) == 9 %sys_indx == 4
    Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
    Q4 = Qm{4}; Q5 = Qm{5}; Q6 = Qm{6};
    Q7 = Qm{7}; Q8 = Qm{8}; Q9 = Qm{9};
    Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3;...
                m'*Q4; m'*Q5; m'*Q6;...
                m'*Q7; m'*Q8; m'*Q9];
    cm = d + A*m +  [(m0'*Q1*m0 + 2*m0'*Q1*Delta_m); 
                     (m0'*Q2*m0 + 2*m0'*Q2*Delta_m);
                     (m0'*Q3*m0 + 2*m0'*Q3*Delta_m);
                     (m0'*Q4*m0 + 2*m0'*Q4*Delta_m); 
                     (m0'*Q5*m0 + 2*m0'*Q5*Delta_m);
                     (m0'*Q6*m0 + 2*m0'*Q6*Delta_m);
                     (m0'*Q7*m0 + 2*m0'*Q7*Delta_m); 
                     (m0'*Q8*m0 + 2*m0'*Q8*Delta_m);
                     (m0'*Q9*m0 + 2*m0'*Q9*Delta_m)];
elseif size(A,2) == 6
    Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
    Q4 = Qm{4}; Q5 = Qm{5}; Q6 = Qm{6};
    Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3;...
                m'*Q4; m'*Q5; m'*Q6];
%     cm = d + A*m + [m'*Q1*m; m'*Q2*m; m'*Q3*m;...
%                     m'*Q4*m; m'*Q5*m; m'*Q6*m];
     cm = d + A*m + [(m0'*Q1*m0 + 2*m0'*Q1*Delta_m); 
                     (m0'*Q2*m0 + 2*m0'*Q2*Delta_m);
                     (m0'*Q3*m0 + 2*m0'*Q3*Delta_m);
                     (m0'*Q4*m0 + 2*m0'*Q4*Delta_m); 
                     (m0'*Q5*m0 + 2*m0'*Q5*Delta_m);
                     (m0'*Q6*m0 + 2*m0'*Q6*Delta_m)];
end
%%%%%%%%%%%
LMIx =[Lm'*P + P*Lm + chi1*P              P*cm;
           (P*cm)'                      -chi1*r_s];
%%% Desired energy decrease rate: Vdot <= -eps_E*||y(t)||^2
Q = -eps_E*Ix;
LMIx = LMIx -  [Q          O1;
                O1'        0];
%%%%%%%%%%%%%%%
LMIx <= 0;
%%%***** minimize the trapping region radius along with Delta m
if sys_indx == 1
    minimize r_s; 
else
    minimize r_s + w_deltam*norm(Delta_m);
end
cvx_end
%%%%%===============
% % Reformulate the LMI in terms of exact, nonlinear cm
% if size(A,2) == 2 
%     cm = d + A*m + [m'*Q1*m; m'*Q2*m];
% elseif size(A,2) == 3 
%     cm = d + A*m + [m'*Q1*m; m'*Q2*m; m'*Q3*m];
% elseif size(A,2) == 9 
%     cm = d + A*m + [m'*Q1*m; m'*Q2*m; m'*Q3*m;...
%                     m'*Q4*m; m'*Q5*m; m'*Q6*m;...
%                     m'*Q7*m; m'*Q8*m; m'*Q9*m];
% elseif size(A,2) == 6
%     cm = d + A*m + [m'*Q1*m; m'*Q2*m; m'*Q3*m;...
%                     m'*Q4*m; m'*Q5*m; m'*Q6*m];
% end
% %%%%%%%%%
% LMIx =[Lm'*P + P*Lm + chi1*P + eps_E*Ix              P*cm;
%                     (P*cm)'                      -chi1*r_s];
%%%%%%%%%%%% Post processing
if r_s > 0 && eps_E > 0 && chi1 > 0 ...
   && max(eig(LMIx)) <= LMIx_tol 
    SDP_soln.chi1 = chi1;
    % SDP_soln.chi2 = chi2;
    SDP_soln.P = P;
    % SDP_soln.pmin = pmin;
    SDP_soln.Lm = Lm;
    SDP_soln.m = m;
    SDP_soln.Deltam = Delta_m;
    SDP_soln.LMIx = LMIx;
   %
    SDP_soln.radius = sqrt(r_s);
    SDP_soln.ROBB = sqrt(r_s)/sqrt(min(eig(P)));
    SDP_soln.eps_E = eps_E;
    SDP_soln.cvxprec = cvx_precision;
    SDP_soln.cvxstatus = cvx_status;
else
    SDP_soln = [];
end