%%%%%%%%%%%%%%%%
% Computes shift coordinates for ellipsoidal trapping region analysis
% Employs an iteration on the ellipsoid shape matrix and shift coordinates
%%%%%%%%%%%%%
function shift_soln = shift_iteration_global(G,sys_info,opt_specs,sys_indx,m,P)
%%%%%%%%%%%%
n = sys_info.n; 
A = sys_info.A;
Qm = sys_info.Qm;
Ix = eye(n);
%%%%%%%%%%%%
pmin = opt_specs.pmin; 
pmax = opt_specs.pmax;
lossless_tol = opt_specs.lossless_tol;
shift_norm = opt_specs.shift_norm;
%%%%%%%%
cvx_clear;
cvx_begin sdp % quiet
% Define the SDP variables
variable a;
if isempty(m) == 1
    variable m(n,1);
    if ~isempty(shift_norm)
        norm(m) <= shift_norm;
    end
    SDP_indx = 1;
%%%%%%%%%%%%
elseif isempty(P) == 1
    if sys_indx == 1 % enforcing the G'*vec(P) = 0 constraint on P through symbolic calculations
        P = Ix; % Symbolic calculations yield P = p*Ix ---> the lossless 
                % inner product x'*P*f(x) has to be zero for P = Ix
    elseif sys_indx == 2 % enforcing the G'*vec(P) = 0 constraint on P through symbolic calculations
        variables p11 p22; 
        P = [p11   0      0;
              0   p22     0;
              0    0     p22];
        P >= pmin*Ix;
        P <= pmax*Ix;
    %%%%%%%%%%%
    elseif sys_indx == 23 % enforcing the G'*vec(P) = 0 constraint on P through symbolic calculations
        variables p11 p14 p33 p44;
        P = [p11     0         0     p14;
              0    2*p33      0        0;
              0     0         p33     0;
              p14   0         0      p44];
        P >= pmin*Ix;
        P <= pmax*Ix;
    else  
        variable P(n,n) symmetric;
        P >= pmin*Ix;
        P <= pmax*Ix;
        G'*P(:) == 0; % enforce the G'*vec(P) = 0 contraint on P numerically
    end
    SDP_indx = 2;
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
elseif size(A,2) == 4 
    Q1 = Qm{1};
    Q2 = Qm{2};
    Q3 = Qm{3};
    Q4 = Qm{4};
    Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3; m'*Q4];
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
minimize a
%-----------------
cvx_end
%%%%%%============
if SDP_indx == 1
    shift_soln.a = a;
    shift_soln.m = m;
    shift_soln.P = P;
    shift_soln.Lm = Lm;
elseif SDP_indx == 2
    if min(eig(P)) > 0.9*pmin && max(eig(P)) < 1.1*pmax ...
       && max(abs(G'*P(:))) <= min(1e-2*min(eig(P)),lossless_tol)  
        shift_soln.a = a;
        shift_soln.m = m;
        shift_soln.P = P;
        shift_soln.Lm = Lm;
        %%%%%%%
        [~,cm] = get_shifted_ops(sys_info,m);
        shift_soln.cm = cm;
    else
        shift_soln = [];
    end
end