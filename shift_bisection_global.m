%%%%%%%%%%%%%%%%
% Computes shift coordinates for ellipsoidal trapping region analysis
% Employs a bisection in the ellipsoid shape and shift coordinates
%%%%%%%%%%%%%
function shift_soln = shift_bisection_global(G,sys_info,opt_specs,sys_indx,m,P)
%%%%%%%%%%%%
n = sys_info.n; 
A = sys_info.A;
d = sys_info.d;
Qm = sys_info.Qm;
Ix = eye(n);
%%%%%%%%%%%%
pmin = opt_specs.pmin; 
pmax = opt_specs.pmax;
lossless_tol = opt_specs.lossless_tol;
solver_indx = opt_specs.solver_indx;
shift_norm = opt_specs.shift_norm;
%%%%%%%%
cvx_clear;
cvx_begin sdp % quiet
if solver_indx == 1
    cvx_solver sdpt3 
elseif solver_indx == 2
    cvx_solver sedumi  
end
%  
if sys_indx == 17 % || sys_indx == 18
    cvx_precision low; 
else
    cvx_precision best;
end

% Define the SDP variables
if isempty(m) == 1
    variable m(n,1);
    if ~isempty(shift_norm)
        norm(m) <= shift_norm;
    end
    variable a;
    % SDP_indx = 1;
%%%%%%%%%%%%
elseif isempty(P) == 1
    variable P(n,n) symmetric
    P >= pmin*Ix;
    P <= pmax*Ix;
    G'*P(:) == 0;
    variable a;
    % SDP_indx = 2;
end
%%%%=====Constraints 
if size(A,2) == 2 
    Q1 = Qm{1};
    Q2 = Qm{2};
    Lm = A + 2*[m'*Q1; m'*Q2];
elseif size(A,2) == 3 
    Q1 = Qm{1};
    Q2 = Qm{2};
    Q3 = Qm{3};
    Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3];
elseif size(A,2) == 9 
    Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
    Q4 = Qm{4}; Q5 = Qm{5}; Q6 = Qm{6};
    Q7 = Qm{7}; Q8 = Qm{8}; Q9 = Qm{9};
    Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3;...
                m'*Q4; m'*Q5; m'*Q6;...
                m'*Q7; m'*Q8; m'*Q9];
elseif size(A,2) == 6
    Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
    Q4 = Qm{4}; Q5 = Qm{5}; Q6 = Qm{6};
    Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3;...
                m'*Q4; m'*Q5; m'*Q6];
end
%%%%%%%%
P*Lm + Lm'*P <= a*Ix;
%-----------------
% if SDP_indx == 1
%     minimize a + 1e-6*norm(m)
% elseif SDP_indx == 2
%     minimize a
% end
%-----------------
minimize a
cvx_end
% %%%%%============
if min(eig(P)) > 0 && norm(G'*P(:)) <= lossless_tol   
    shift_soln.a = a;
    shift_soln.m = m;
    shift_soln.P = P;
    shift_soln.Lm = Lm;
    if size(A,2) == 2 
        cm = d + A*m + [m'*Q1*m; m'*Q2*m];
    elseif size(A,2) == 3 
        cm = d + A*m + [m'*Q1*m; m'*Q2*m; m'*Q3*m];
    elseif size(A,2) == 9 
        cm = d + A*m + [m'*Q1*m; m'*Q2*m; m'*Q3*m;...
                        m'*Q4*m; m'*Q5*m; m'*Q6*m;...
                        m'*Q7*m; m'*Q8*m; m'*Q9*m];
    elseif size(A,2) == 6
        cm = d + A*m + [m'*Q1*m; m'*Q2*m; m'*Q3*m;...
                        m'*Q4*m; m'*Q5*m; m'*Q6*m];
    end
    shift_soln.cm = cm;
else
    shift_soln = [];
end
