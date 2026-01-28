%%% Semidefinite program to compute ROBB (ROBB = Radius of Outer Bounding Ball)
function ROBB = compute_ROBB(opt_soln)
P_opt = opt_soln.P;
% traceP_opt = trace(P_opt);
radius_opt = opt_soln.radius;
% rE_opt_s = radius_opt^2/((1/n)*traceP_opt);
% rE_opt = sqrt(rE_opt_s);

S2 = P_opt; r2 = radius_opt;
S1 = eye(size(P_opt));
%%% Set containment: Ellipsoid(S2,r2) is a subset to the ellipsoid(S1,r1)
%%% It holds if and only if there is a chi >= 0 s.t the LMI here holds
cvx_clear;
cvx_begin sdp % quiet
cvx_solver sedumi %sdpt3 % %  
%  
cvx_precision high
% Define the SDP variables
variables chi r1s
%%%%%%%%%%%%%
LMI_sc = blkdiag(chi*S2-S1, r1s -chi*r2^2);
LMI_sc >= 0;
chi >= 0;
%%%%%%%%%
minimize r1s;
cvx_end
%%%%%%%%%%%%%%
if chi>= 0 && r1s >= 0
    rE_opt = sqrt(r1s); 
    ROBB = rE_opt;
else
    ROBB = [];
end