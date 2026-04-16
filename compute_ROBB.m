%%% Semidefinite program to compute ROBB (ROBB = Radius of Outer Bounding Ball)
function ROBB = compute_ROBB(opt_soln)
P_opt = opt_soln.P;
radius_opt = opt_soln.radius;
%%%%%%%%%%
S2 = P_opt; r2 = radius_opt;
S1 = eye(size(P_opt));
%%% Set containment: Ellipsoid(S2,r2) is subset to ellipsoid(S1,r1)
%%% It holds if and only if there is a chi >= 0 s.t. the LMI here holds
cvx_clear;
cvx_begin sdp % quiet
% Define the SDP variables
variable chi nonnegative
variable r1s nonnegative
%%%%%%%%%%%%%
LMI_sc = blkdiag(chi*S2-S1, r1s -chi*r2^2);
LMI_sc >= 0;
%%%%%%%%%
minimize r1s;
cvx_end
%%%%%%%%%%%%%%
if chi >= 0 && r1s >= 0
    rE_opt = sqrt(r1s); 
    ROBB = rE_opt;
else
    ROBB = [];
end