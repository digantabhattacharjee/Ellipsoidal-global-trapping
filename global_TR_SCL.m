%%%%%%%%%%%%%%%%
% Computes global spherical trapping regions for given shift coordinates
% using the framework in SCL et al., IJRNC 2025
% Link: https://onlinelibrary.wiley.com/doi/full/10.1002/rnc.7807
%%%%%%%%%%%%%
function SDP_soln = global_TR_SCL(sys_info,LMIx_tol,m)
%%%%%%%%%%%%
n = sys_info.n;
Ix = eye(n);
[Lm,cm] = get_shifted_ops(sys_info,m);
%%%%=====
cvx_clear;
cvx_begin sdp % quiet
% 
Lms = 0.5*(Lm + Lm'); % symmetric part of the linear operator
% Define the SDP variables
variable chi nonnegative;
variable r_s;
%%%%%**************Constraints*******
LMIx = [Ix + chi*Lms      (chi/2)*cm;
         (chi/2)*cm'         -r_s];
%
LMIx = 0.5*(LMIx + LMIx'); %symmetrize 
LMIx <= 0;
%%%***** minimize the trapping region radius 
minimize r_s;
cvx_end
%%%%%===============
if r_s > 0 && chi >= 0 && max(eig(LMIx)) <= LMIx_tol 
    SDP_soln.Lm = Lm;
    SDP_soln.cm = cm;
    SDP_soln.m = m;
    SDP_soln.LMIx = LMIx;
    SDP_soln.radius = sqrt(r_s);
    SDP_soln.chi = chi;
else
    SDP_soln = [];
end