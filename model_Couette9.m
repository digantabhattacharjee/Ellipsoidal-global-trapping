function model = model_Couette9(Re,Lx,Lz)
% model struct object
%
% [1] A. Kalur, T. Mushtaq, P. Seiler, M. S. Hemati, Estimating Regions
% of Attraction for Transitional Flows using Quadratic Constraints

    if nargin == 0
        Lx=1.75*pi; %X-domain size
        Lz=1.5*pi; %Y-domain size
        Re=20; %Renolds number
    elseif nargin == 1
        % default parameter following implementation of [1] 
        Lx = 1.75*pi;
        Lz = 1.5*pi;
    end

    nx = 9;
    model.nx = nx;
    model.m = zeros(model.nx, 1);

    [E,L,Q_out] = Couette_AmplitudeEqns(Lx,Lz,Re);

    model.c = E;
    model.L = L;
    model.Ls = 1/2*(L+L');
    model.Q = Q_out;
    
    model.ode = @(t,x) func_xdotNL(model.c, model.L, model.Q, t, x);
end