%%%%%%%%%%%%%%%%
% Computes global spherical trapping regions for given shift coordinates
% using the framework in SCL et al., IJRNC 2025
% Link: https://onlinelibrary.wiley.com/doi/full/10.1002/rnc.7807
%%%%%%%%%%%%%
function SDP_soln = global_TR_SCL(m,A,Qm,d,solver_indx,LMIx_tol)
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
% 
n = size(A,1);
Ix = eye(n);
%%%%=====
cvx_clear;
cvx_begin sdp % quiet
if solver_indx == 1
    cvx_solver sdpt3 
elseif solver_indx == 2
    cvx_solver sedumi  
end
%  
cvx_precision best % cvx_precision high
Lms = 0.5*(Lm+Lm');
% Define the SDP variables
variable chi
chi >= 0;
variable r_s
% % %%%%**************Constraints*******
LMIx = [Ix + chi*Lms              (chi/2)*cm;
         (chi/2)*cm'                -r_s];
%%%%%%%%%%%%%%%
LMIx <= 0;
%%%***** minimize the trapping region radius 
minimize r_s;
cvx_end
%%%%%===============
if r_s > 0 && chi >= 0 && max(eig(LMIx)) <= LMIx_tol 
    SDP_soln.Lm = Lm;
    SDP_soln.m = m;
    SDP_soln.LMIx = LMIx;
    SDP_soln.radius = sqrt(r_s);
    SDP_soln.chi = chi;
else
    SDP_soln = [];
end