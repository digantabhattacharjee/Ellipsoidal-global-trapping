%%%%%%%%%%%%%%%%
% Computes global trapping ellipsoids for given shift coordinates
%%%%%%%%%%%%%
function SDP_soln = global_TR_genloss_ell(sys_info,G,opt_specs,chi1,m,P,r_s)
%%%%%%%%%%
sys_indx = sys_info.sys_indx;
pmin = opt_specs.pmin; 
pmax = opt_specs.pmax;
eps_E = opt_specs.eps_E;
lossless_tol = opt_specs.lossless_tol;
n = sys_info.n; 
Ix = eye(n);
%
[Lm,cm] = get_shifted_ops(sys_info,m);
%%%%=====
cvx_clear;
cvx_begin sdp % quiet
% Define the SDP variables
if isempty(P) == 1 && isempty(r_s) == 1
    SDP_indx = 1;
    if sys_indx == 23 % enforcing the G'*vec(P) = 0 on P through symbolic calculations
        variables p11 p14 p33 p44;
        P = [p11     0         0     p14;
              0    2*p33      0        0;
              0     0         p33     0;
              p14   0         0      p44];
        P >= pmin*Ix;
        P <= pmax*Ix;
        variable r_s;
    elseif sys_indx == 1 % enforcing the G'*vec(P) = 0 on P through symbolic calculations
        P = Ix; % Symbolic calculations yield P = p*Ix ---> the lossless
                % inner product x'*P*f(x) has to be zero for P = Ix
        variable r_s;
    elseif sys_indx == 2 % enforcing the G'*vec(P) = 0 on P through symbolic calculations
        variables p11 p22;
        P = [p11   0      0;
              0   p22     0;
              0    0     p22];
        P >= pmin*Ix;
        P <= pmax*Ix;
        %%%%%%%
        variable r_s;
    else
        variable P(n,n) symmetric;
        variable r_s;
        %%%%%**************Constraints*******
        P >= pmin*Ix; 
        P <= pmax*Ix;
        G'*P(:) == 0; % enforce the G'*vec(P) = 0 contraint on P numerically
    end
elseif isempty(chi1) == 1
    variable chi1 nonnegative; 
    SDP_indx = 0;
end
%%%%%%%%%%%%%
LMIx = [P*Lm + Lm'*P + chi1*P              P*cm;
         (P*cm)'                      -chi1*r_s];
LMIx = 0.5*(LMIx + LMIx'); %symmetrize 
LMIx <= -eps_E*eye(size(LMIx));
%%%***** minimize the trapping region radius 
if SDP_indx
    minimize r_s;
end
cvx_end
%%%%%====== Post processing
S = -P; % Generalized lossless matrix (S = P is also valid)
lam_LmP = max(eig(Lm'*P + P*Lm));
%*********************
if r_s > 0 && eps_E > 0 && chi1 > 0 ...
   && min(eig(P)) >= 0.9*pmin && max(eig(P)) < 1.1*pmax ...
   && lam_LmP < 0 && max(eig(LMIx)) < 0 ...
   && max(abs(G'*P(:))) <= lossless_tol
    SDP_soln.chi1 = chi1;
    SDP_soln.P = P;
    SDP_soln.opt_specs = opt_specs;
    SDP_soln.Lm = Lm;
    SDP_soln.m = m;
    SDP_soln.LMIx = LMIx;
    SDP_soln.radius = sqrt(r_s);
    SDP_soln.ROBB = sqrt(r_s)/sqrt(min(eig(P)));
    SDP_soln.S = S;
    SDP_soln.cvxprec = cvx_precision;
    SDP_soln.cvxstatus = cvx_status;
else
    SDP_soln = [];
end