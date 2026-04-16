%***************
%%% Ellipsoidal global trapping analysis under generalized lossless property 
% Paper reference: "Trapping Regions for Quadratic Systems with Generalized
% Lossless Nonlinearities" by Diganta Bhattacharjee, Shih-Chi Liao, Peter J. Seiler, and Maziar S. Hemati
%***************
% Subroutines needed: 
% get_system; get_6state_model; model_Couette9; Couette_AmplitudeEqns; 
% shift_iteration_global; get_gen_lossless_matrix; 
% global_TR_genloss_ell; global_TR_genloss_shift; 
% shift_SDP; global_TR_SCL; 
% get_lossless_innerprod (optional); compute_ROBB(optional)
%***************
% External packages: (1) cvx, (2) SOSTOOLS, 
% (3) Phase portrait (https://www.mathworks.com/matlabcentral/fileexchange/110785-phase-portrait-plotter-on-2d-phase-plane)
%%%%%%%%%%%%%%%%
% Demonstrated on: 
% (sys_indx = 1)  2D toy problem (SCL et al., IJRNC 2025)
% (sys_indx = 2)  Lorenz system 
% (sys_indx = 4)  9-state Couette flow model
% (sys_indx = 5)  2D toy problem
% (sys_indx = 6)  2D toy problem
% (sys_indx = 12) Lorenz83 system
% (sys_indx = 13) 6-state airfoil flow model
% (sys_indx = 14) 2D toy example
% (sys_indx = 19) Three-dimensional cylinder flow model (aka mean field model)
% (sys_indx = 21) Finance chaotic system
% (sys_indx = 22) Atmospheric oscillator model
% (sys_indx = 23) Modified Lorenz–Stenflo system
% (sys_indx = 26) Hadley cell or Hadley circulation model
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Written by DB, 2026
%%% Contribution from PJS (generalized lossless constraint)
%%% Contribution from SCL and LH (9-state model)
%%%%%%%%%%%%%%%%%
clc; 
clear all;
close all;
%%%%%%%%%%%%%%%%%% IMPORTANT %%%%%%%%%%%%%%%%%%%%%%
% Add SOSTOOLS path here 
addpath(genpath(pwd))
%%%%%%%%%%%%
sys_indx = 13; 
%%%%%%%%
sys_info = get_system(sys_indx);
sys_info.sys_indx = sys_indx;
n = sys_info.n; 
A = sys_info.A;
d = sys_info.d;
Qm = sys_info.Qm;
Ix = eye(n);
%% Specify method of analysis (spherical or ellipsoidal), shift coordinates
P_indx = 1; % 1 = Solve P through the SDP; 0 = P is fixed (Ix or scaled Ix)
% NOTE: (1) No need for multiplier grid when P is fixed (Ix or scaled Ix),
% reverts back to the setup in SCL et al., IJRNC 2025
% (2) S (generalized lossless matrix) is solved through the SDP when
% P isn't fixed; it is equal to -P (or P, the sign doesn't matter)
%%%%%%%%%%%%
lossless_tol = sqrt(eps);   % lossless constraint tolerance
pmin = 1;  
pmax = 10; % min and max eigenvalues of P for shift coordinate computation
if sys_indx == 13 
    pmin = 1e-4;
    pmax = 1;
end
%%%%%%%
if 0
    innPD = get_lossless_innerprod(n,Qm,Ix); % check lossless inner product of the nonlinearity
end
%%%%%%%% Set solver and required precision level
solver_indx = 1; %% 1 = sdpt3 (cvx default), 2 = sedumi
if solver_indx == 1
    cvx_solver sdpt3 
    cvx_precision best; % precision specified for solver
elseif solver_indx == 2
    cvx_solver sedumi  
    cvx_precision best; % precision specified for solver
end
opt_specs.pmin = pmin; 
opt_specs.pmax = pmax;
opt_specs.lossless_tol = lossless_tol;
%%%%%%
G = get_gen_lossless_matrix(n,Qm); % Generalized lossless constraint matrix
% Check if the lossless property holds *exactly* with specified 
% structure of the energy or ellipsoid shape matrix P
loss_constraint = [];
if sys_indx == 2 % enforcing the G'*vec(P) = 0 on P through symbolic calculations
    p11 = randn(1); p22 = randn(1); 
    P = [p11   0      0;
          0   p22     0;
          0    0     p22];
    loss_constraint = max(abs(G'*P(:))); % should be equal to zero
elseif sys_indx == 23    
    p11 = randn(1); p14 = randn(1); p33 = randn(1); p44 = randn(1);
    P = [p11     0         0     p14;
          0    2*p33      0        0;
          0     0         p33     0;
          p14   0         0      p44];
    loss_constraint = max(abs(G'*P(:))); % should be equal to zero
end
if ~isempty(loss_constraint) && loss_constraint ~= 0 
    fprintf("Lossless constraint not satisfied, check the form of P via symbolic calculations");
    return;
end
%%%%%%%%%%%%%%%%%%%%
if ~P_indx
    %%%%%% Compute/specify the shift coordinate
    if sys_indx ~= 1 
        P = Ix;
        if max(abs(G'*P(:))) > 1e-3 %should be ~0 if feasible
            fprintf("Nonlinearity non-lossless with respect to identity: " + ...
                    "a spherical global trapping region does not exist for this system. " + ...
                    "Update the inner product matrix choice and try again. " + ...
                    "Or, do the iteration to possibly resolve the issue");
            return;
        else
            shift_soln = shift_SDP(sys_info,P);
            m = shift_soln.m;
            Lm = shift_soln.Lm;
        end
    else % zero shift @ sys_indx == 1 (2D toy problem (SCL et al., IJRNC 2025))
        m = zeros(n,1);
        Lm = A;
    end
else  % Generalized lossless case 
    if 0 % check if S = Ix is a feasible solution for the lossless property
        S = Ix;
        max(abs(G'*S(:))) %should be ~0 if S = Ix is feasible
    end
    % Initial iteration parameters
    a_tol = 1e-6;
    max_counter = 30;
    %%%%% Specify shift coordinate norm bound
    if sys_indx == 2
        shift_norm = 35; %38;
    elseif sys_indx == 13
        shift_norm = 50; 
    elseif sys_indx == 14
        shift_norm = 1; %10; 
    elseif sys_indx == 23
        shift_norm = 30;
    else
        shift_norm = [];
    end
    opt_specs.shift_norm = shift_norm;
    %%%%%%%%%%%%%%%%%%%%%%%
    if sys_indx == 1 % zero shift @ sys_indx == 1 (2D toy problem (SCL et al., IJRNC 2025))
        m = zeros(n,1);
        Lm = A;
        Pshift = Ix;
        P = Ix;
    else
    %***** Obtain shift coordinates by iterating over ellipsoid shape P and shift coordinates m 
        a = 0; a_new = 1; counter = 0; 
        zz = randn(n,1);
        m = (1/norm(zz))*zz; 
        while norm(a - a_new) > a_tol && counter <= max_counter && a_new > 0
            shift_soln_1 = ...
            shift_iteration_global(G,sys_info,opt_specs,sys_indx,m,[]);
            if ~isempty(shift_soln_1)
                Pshift = shift_soln_1.P;
                a = shift_soln_1.a;
            else
                fprintf("Numerical issues or infeasibility, " + ...
                        "check tolerances on ellipsoid shape matrix " + ...
                        "and generalized lossless condition");
                return;
            end
            %%%%%%%
            shift_soln_2 = ...
            shift_iteration_global(G,sys_info,opt_specs,sys_indx,[],Pshift);
            if ~isempty(shift_soln_2)
                m = shift_soln_2.m;
                Lm = shift_soln_2.Lm;
                a_new = shift_soln_2.a;
            else
                fprintf("Numerical issues");
                return;
            end
          %%%%%%%%%%%%%
            counter = counter + 1; % update iteration counter
        end
        %%%%% Check the final solution from iterations 
        if a_new < 0
            shift_soln = shift_soln_2;
            shift_soln.shift_norm = shift_norm;
            m = shift_soln.m;
            Lm = shift_soln.Lm;
            P = shift_soln.P;
        else
            fprintf("Iteration didn't yield suitable shift coordinates:" + ...
                    "a global (ellipsoidal) trapping might not exist");
            return;
        end
    end
end
%% Compute trapping region (ellipsoidal or spherical)
if ~P_indx %%% This is the spherical case
    LMIx_tol = 1e-12;
    opt_soln = global_TR_SCL(sys_info,LMIx_tol,m);
    opt_soln.P = Ix;
    if ~isempty(opt_soln)
        R_opt = opt_soln.radius  % This is the optimal spherical trapping radius 
        opt_soln.ROBB = R_opt;   % needed for Vdot computations later
        Ult_bound = norm(m) + R_opt
    else
        fprintf("Infeasibility, check the choice of shift coordinates or the choice of P (if fixed)");
    end
else % Generalized lossless case (ellipsoidal trapping region)
    % Set the LMI tolerance and multiplier (chi) grid
    if sys_indx == 1 || sys_indx == 5 || sys_indx == 6 || sys_indx == 14 % 2D examples
        eps_E = 1e-6;  
        %%%%% Define the multipler grid
        Nchi = 100;
        chi_grid = logspace(-1,1,Nchi)';
    elseif sys_indx == 2  % Lorenz attractor
        eps_E = 1e-6;
        %%%%% Define the multipler grid
        Nchi = 100;
        chi_grid = logspace(-1,1,Nchi)';
     elseif sys_indx == 4  % Couette flow model
        eps_E = 1e-6;
        %%%%% Define the multipler grid
        Nchi = 100;
        chi_grid = logspace(-2,0,Nchi)';
    % elseif sys_indx == 5 || sys_indx == 6  % 2D system
    %     eps_E = 1e-6;
    %     %%%%% Define the multipler grid
    %     Nchi = 100;
    %     chi_grid = logspace(-1,1,Nchi)';
    elseif sys_indx == 12 % Lorenz83 system
        Nchi = 100;
        %%%%% Define the multipler grid
        chi_grid = logspace(-1,1,Nchi)';
        eps_E = 1e-6; 
     elseif sys_indx == 13 % Airfoil model, Heide et al., 2025 
        eps_E = 1e-6;
        %%%%% Define the multipler grid 
        Nchi = 100;
        chi_grid = logspace(-1,1,Nchi)'; 
    elseif sys_indx == 19 % Cylinder wake model
        eps_E = 1e-6;
        %%%%% Define the multipler grid
        Nchi = 100; 
        chi_grid = logspace(-1,1,Nchi)';
     elseif sys_indx == 22 % Atmospheric oscillator model 
         eps_E = 1e-6;
         %%%%% Define the multipler grid
         Nchi = 100;
         chi_grid = logspace(-2,0,Nchi)';        
     elseif sys_indx == 23 % Stochastic Lorenz–Stenflo system (modified)
        eps_E = 1e-6;
        %%%%% Define the multipler grid
        Nchi = 100; 
        chi_grid = logspace(-1,1,Nchi)';
     else
        eps_E = 1e-6; 
        %%%%% Define the multipler grid
        Nchi = 100; 
        chi_grid = logspace(-3,2,Nchi)';
    end
    opt_specs.eps_E = eps_E;
    %%*************************************************
    % NOTE: There are many ways of computing the trapping ellipsoids as 
    % the optimization is nonlinear. We choose to optimize the trapping
    % ellipsoid on a specified grid of mulpliers, take the best solution 
    % and then refine the multipler using a GEVP approach. Finally, a local
    % search is conducted to refine the choice of shift coordinates
    %%*************************************************
    trace_store = zeros(length(chi_grid),1);
    radius_store = zeros(length(chi_grid),1);
    ROBB_store = radius_store;
    RIBB_store = ROBB_store;
    max_lossless_innPD = RIBB_store;
    soln_store = cell(length(chi_grid),1);
    %%%%%%%%%%%%%%%%%%%%%%
    for ct1=1:length(chi_grid)
        chi1 = chi_grid(ct1);
        soln = global_TR_genloss_ell(sys_info,G,opt_specs,chi1,m,[],[]);   
        soln_store{ct1} = soln;
        %%%*************
        if isempty(soln)
            trace_store(ct1) = NaN;
            radius_store(ct1) = NaN;
            ROBB_store(ct1) = NaN;
            RIBB_store(ct1) = NaN;
            max_lossless_innPD(ct1,1) = NaN;
        else
            P = soln.P;
            trace_store(ct1) = trace(P);
            radius_store(ct1) = soln.radius;
            ROBB_store(ct1) = soln.ROBB;
            RIBB_store(ct1) = soln.radius/sqrt(max(eig(P)));
            max_lossless_innPD(ct1,1) = max(abs(G'*P(:)));
        end
    end
    %%%%%%%%%% Pick the optimal solution on the multiplier grid
    [~, opt_indx] = min(ROBB_store); 
    opt_grid_soln = soln_store{opt_indx};
    opt_soln = opt_grid_soln; % this is the optimal solution on the grid (without adjusting the shift coordinates)
    m = opt_soln.m;
    P = opt_soln.P;
    %*******************
    % Update the multiplier and ellipsoid radius by solving a GEVP 
    % Use a bisection on ellipsoid radius
    %*********************
    % Inputs: 
    rLB = 0; 
    rUB = 10*opt_soln.radius; 
    relTol = eps; %1e-12;
    absTol = eps; %1e-12;
    while rUB-rLB > relTol*rUB + absTol
        r = (rUB + rLB)/2;
        r_s = r^2;
        %%%%%%%%%%%
        bisection_soln = global_TR_genloss_ell(sys_info,G,opt_specs,[],m,P,r_s);
        if ~isempty(bisection_soln) % The LMI with r fixed is feasible
            rUB = r;
            opt_GEVP_soln = bisection_soln;
        else
            rLB = r;
        end
    end
    %%%%%%
    opt_soln = opt_GEVP_soln;
    %% Local search/tuning in the shift coordinates
    shift_local_search = 1;
    local_search_success = 0; % index related to local shift coordinate search
    w_deltam = 1e-6; 
    opt_specs.w_deltam = w_deltam;
    %%%%%%%%
    if ~isempty(opt_soln) 
        if shift_local_search
        % Optimize the shift coordinates with a linearized c(m)
            fprintf("Local search in shift coordinates attempted");
            opt_soln_shift1 = ...
            global_TR_genloss_shift(sys_info,opt_specs,opt_soln.chi1,opt_soln.m,opt_soln.P);
            %%%%%%%%%%%%%%%
            if ~isempty(opt_soln_shift1)
            % Optimize P and radius for the above shift coordinates
            % without any approximation in c(m)
                opt_soln_shift2 = ...
                global_TR_genloss_ell(sys_info,G,opt_specs,opt_soln.chi1,opt_soln_shift1.m,[],[]);
                if ~isempty(opt_soln_shift2) && opt_soln_shift2.ROBB < opt_soln.ROBB  
                    opt_soln = opt_soln_shift2;
                    fprintf("Local search in shift coordinates improves results:" + ...
                            "optimal solution updated");
                    local_search_success = 1;
                else
                    fprintf("Local search in shift coordinates *did not* improve results");
                end
            end
        end
        %%%%%%%%%%%%%%%%%%
        % The largest semi-major axis length equals the ROBB (= Radius of Outer Bounding Ball)
        % Similarly, the smallest principal semi-major axis length equals
        % the radius of the largest norm ball that can be fitted inside the
        % trapping ellipsoid (RIBB = Radius of Inscribed Ball)
        %%%%%%%%%%%%%%%%
        ROBB_opt = opt_soln.ROBB %opt_soln.radius/sqrt(min(eig(opt_soln.P)));
        RIBB_opt = opt_soln.radius/sqrt(max(eig(opt_soln.P)));
        %%% ROBB = compute_ROBB(opt_soln);
        Ult_bound = norm(opt_soln.m) + ROBB_opt %this is the ultimate bound in the original coordinates
    else
        ROBB_opt = NaN; Ult_bound = NaN;
        fprintf("Infeasibility, check the choice of shift coordinates or the choice of P (if fixed)");
    end
%%  
    if ~isempty(opt_soln) && exist('ROBB_store','var') ~= 0
        figure()
        loglog(chi_grid, max_lossless_innPD(:,1),'b.','LineWidth',2);
        %
        xlabel('$\chi$','Interpreter','latex');
        ylabel('$\mathrm{max}(G^\mathrm{T} \mathrm{vec}(P))$','Interpreter','latex');
        grid on;
        set(gca,'Fontsize',15);
        % xlim([min(chi_grid) max(chi_grid)]);
        title('Generalized lossless constraint');
    
        figure()
        loglog(chi_grid, ROBB_store(:,1),'bo','LineWidth',1.5);
        %
        xlabel('$\chi$','Interpreter','latex');
        ylabel('Trapping Radius (Spherical)','Interpreter','latex');
        % ylabel('$\alpha$','Interpreter','latex');
        grid on;
        set(gca,'Fontsize',15);
        % title('Radius of outer bounding norm ball');
       % xlim([min(chi_grid) max(chi_grid)]);
        
        figure()
        loglog(chi_grid, radius_store(:,1),'ro','LineWidth', 1.5);
        xlabel('$\chi$','Interpreter','latex');
        ylabel('Trapping Radius (Ellipsoidal)','Interpreter','latex');
        grid on;
        set(gca,'Fontsize',15);
        % xlim([min(chi_grid) max(chi_grid)]);
        %%%%%%%%%
    end
end
%% Check the final solution obtained for accuracy
if ~isempty(opt_soln)
    [max(abs(G'*opt_soln.P(:)))  max(eig(opt_soln.LMIx))  max(eig(opt_soln.LMIx(1:n,1:n)))  min(eig(opt_soln.P))]
    %%% Should be [~0  -ve(~0)  -ve(~0)  +ve]
    m_opt = opt_soln.m;
end
%% Check Vdot on the trapping region/ellipsoid boundary
Vdot_check_indx = 1;
if Vdot_check_indx
    if ~isempty(opt_soln)
        [Lm,cm] = get_shifted_ops(sys_info,m_opt);
        if n == 2 
            Q1 = Qm{1}; Q2 = Qm{2};
        elseif n == 3 
            Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3}; 
        elseif n == 4 
            Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3}; Q4 = Qm{4};
        elseif n == 9 
            Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
            Q4 = Qm{4}; Q5 = Qm{5}; Q6 = Qm{6};
            Q7 = Qm{7}; Q8 = Qm{8}; Q9 = Qm{9};
        elseif n == 6
            Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
            Q4 = Qm{4}; Q5 = Qm{5}; Q6 = Qm{6};
        end
        N_run = 1e6; % number of random runs for computing Vdot
        Vdot1 = zeros(N_run,1); Vdot2 = Vdot1;
        Rvec = randn(n,N_run); 
        opt_radius = opt_soln.radius;
        P_ell =  opt_soln.P/(opt_radius^2);
        P_ell = full(P_ell);
        E_ell = sqrtm(inv(P_ell));
        for ct = 1:N_run
            rvec = Rvec(:,ct); 
            rvec = (1/norm(rvec))*rvec;
            Y = E_ell*rvec; % [Y'*P_ell*Y - 1 ~ 0], i.e., Y is on the trapping ellipsoid boundary
            if n == 2
                fy = [Y'*Q1*Y; Y'*Q2*Y];
            elseif n == 3
                fy = [Y'*Q1*Y; Y'*Q2*Y; Y'*Q3*Y];
            elseif n == 4
                fy = [Y'*Q1*Y; Y'*Q2*Y; Y'*Q3*Y; Y'*Q4*Y];
            elseif n == 6
                fy = [Y'*Q1*Y; Y'*Q2*Y; Y'*Q3*Y; Y'*Q4*Y; Y'*Q5*Y; Y'*Q6*Y];
            elseif n == 9
                fy = [Y'*Q1*Y; Y'*Q2*Y; Y'*Q3*Y;...
                      Y'*Q4*Y; Y'*Q5*Y; Y'*Q6*Y;...
                      Y'*Q7*Y; Y'*Q8*Y; Y'*Q9*Y];
            end
            ydot = Lm*Y + fy + cm;
            Vdot1(ct,1) = ydot'*P_ell*Y + Y'*P_ell*ydot;                           % true derivative
            Vdot2(ct,1) = Y'*(Lm'*P_ell + P_ell*Lm)*Y + Y'*P_ell*cm + cm'*P_ell*Y; % approximate, artificially applying the lossless inner product
        end
        [max(Vdot2)   max(Vdot1)] % both should be approximately the same and negative
    end
end
%% 2D phase portrait and the associated trapping region 
phase_portrait_indx = 1; % set this to 1 for plotting phase portraits
if phase_portrait_indx
    addpath(genpath(pwd)) % Add the phase portrait plotter path here
end
if n == 2 && phase_portrait_indx
    if sys_indx == 1
        ODEfunc = @(t,x)[-x(1) - x(1)*x(2);...
                         -4*x(2) + x(1)^2  + 1];
    elseif sys_indx == 5
        ODEfunc = @(t,x)[-x(1) - x(1)*x(2);...
                         -x(2) + x(1)^2];
    elseif sys_indx == 6
        ODEfunc = @(t,x)[x(1) + x(1)*x(2);...
                         -x(2) - x(1)^2];
    elseif sys_indx == 14
        ODEfunc = @(t,x)[x(1) + 2*x(1)*x(2) - x(2)^2;...
                         -x(2) + x(1)*x(2)   - 2*x(1)^2];
    end
        % figure()
        plotpp(ODEfunc, 'xlim', [-3,3], 'ylim', [-3,3]);
        set(gca,'fontsize',15);
        xlabel('$x_1$', 'interpreter', 'latex','FontSize',15);
        ylabel('$x_2$', 'interpreter', 'latex','FontSize',15);
        hold on;
        th = (0:0.01:2*pi)';
        %
        cir_x = m_opt(1) + opt_soln.ROBB*cos(th);
        cir_y = m_opt(2) + opt_soln.ROBB*sin(th);
        plot(cir_x, cir_y, 'r','LineWidth', 5);
        axis equal;
end
%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%