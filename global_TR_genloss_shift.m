%%%%%%%%%%%%%%%%
% Computes the shift coordinate perturbations for a given 
% (global) trapping ellipsoid and shift coordinates (m0)
% Requires linearizing the nonlinear matrix inequality (trapping condition) 
%%%%%%%%%%%%%
function SDP_soln = global_TR_genloss_shift(sys_info,opt_specs,chi1,m0,P)
%%%%%%%%%%%%
n = sys_info.n; 
A = sys_info.A;
d = sys_info.d;
Qm = sys_info.Qm;
%%%%%%%%%
eps_E = opt_specs.eps_E;
w_deltam = opt_specs.w_deltam;
%%%%=====
cvx_clear;
cvx_begin sdp % quiet
% Define the SDP variables
variable Delta_m(n,1);
m = m0 + Delta_m;
variable r_s;
%%%%%%%%%%%%%%%%%%%%%%%
if size(A,2) == 2 
    Q1 = Qm{1};
    Q2 = Qm{2};
    Lm = A + 2*[m'*Q1; m'*Q2];
    cm = d + A*m + [(m0'*Q1*m0 + 2*m0'*Q1*Delta_m); 
                    (m0'*Q2*m0 + 2*m0'*Q2*Delta_m)];
elseif size(A,2) == 3 
    Q1 = Qm{1};
    Q2 = Qm{2};
    Q3 = Qm{3};
    Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3];
    cm = d + A*m + [(m0'*Q1*m0 + 2*m0'*Q1*Delta_m); 
                    (m0'*Q2*m0 + 2*m0'*Q2*Delta_m);
                    (m0'*Q3*m0 + 2*m0'*Q3*Delta_m)];

elseif size(A,2) == 4 
    Q1 = Qm{1};
    Q2 = Qm{2};
    Q3 = Qm{3};
    Q4 = Qm{4};
    Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3; m'*Q4];
    % 
    cm = d + A*m + [(m0'*Q1*m0 + 2*m0'*Q1*Delta_m); 
                    (m0'*Q2*m0 + 2*m0'*Q2*Delta_m);
                    (m0'*Q3*m0 + 2*m0'*Q3*Delta_m);
                    (m0'*Q4*m0 + 2*m0'*Q4*Delta_m)];
elseif size(A,2) == 9 
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
    cm = d + A*m + [(m0'*Q1*m0 + 2*m0'*Q1*Delta_m); 
                     (m0'*Q2*m0 + 2*m0'*Q2*Delta_m);
                     (m0'*Q3*m0 + 2*m0'*Q3*Delta_m);
                     (m0'*Q4*m0 + 2*m0'*Q4*Delta_m); 
                     (m0'*Q5*m0 + 2*m0'*Q5*Delta_m);
                     (m0'*Q6*m0 + 2*m0'*Q6*Delta_m)];
end
%%%%%%%%%%%
LMIx =[P*Lm + Lm'*P + chi1*P                 P*cm;
           (P*cm)'                      -chi1*r_s];
LMIx = 0.5*(LMIx + LMIx');     %symmetrize 
LMIx <= -eps_E*eye(size(LMIx));
%%%***** minimize the trapping ellipsoid radius along with Delta m
minimize r_s + w_deltam*norm(Delta_m);
cvx_end
%%%%%%%%%%%% Post processing
if r_s > 0 && eps_E > 0 && chi1 > 0 && max(eig(LMIx)) < 0
    SDP_soln.chi1 = chi1;
    SDP_soln.P = P;
    SDP_soln.Lm = Lm;
    SDP_soln.m = m;
    SDP_soln.Deltam = Delta_m;
    SDP_soln.LMIx = LMIx;
    SDP_soln.radius = sqrt(r_s);
    SDP_soln.ROBB = sqrt(r_s)/sqrt(min(eig(P)));
    SDP_soln.eps_E = eps_E;
    SDP_soln.cvxprec = cvx_precision;
    SDP_soln.cvxstatus = cvx_status;
else
    SDP_soln = [];
end