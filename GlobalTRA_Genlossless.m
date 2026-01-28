%***************
%%% Ellipsoidal global trapping analysis under generalized lossless property 
%***************
% Subroutines needed: 
% get_system; get_6state_model; model_Couette9; Couette_AmplitudeEqns; 
% shift_bisection_global; get_gen_lossless_matrix; 
% global_TR_genloss_ell; global_TR_genloss_shift; 
% shift_SDP; global_TR_SCL; 
% get_lossless_innerprod (optional); compute_ROBB(optional)
%***************
% External packages: (1) cvx, (2) SOSTOOLS, 
% (3) (optional) Phase portrait (https://www.mathworks.com/matlabcentral/fileexchange/110785-phase-portrait-plotter-on-2d-phase-plane)
%%%%%%%%%%%%%%%%
% Demonstrated on: 
% (sys_indx = 1)  2D toy problem (SCL et al., IJRNC 2025); 
% (sys_indx = 2)  Lorenz system; 
% (sys_indx = 4)  9-state Couette flow model; 
% (sys_indx = 9)  Chen system;
% (sys_indx = 12) Lorenz83 system
% (sys_indx = 13) 6-state airfoil flow model 
% (sys_indx = 14) 2D toy example
% (sys_indx = 17) Three-scroll system
% (sys_indx = 19) Cylinder wake model
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Written by DB, October-November 2025
%%% Contribution from PJS (generalized lossless constraint)
%%% Contribution from SCL and LH (9-state model)
%%%%%%%%%%%%%%%%%
clc; clear;
close all;
%%%%%%%%%%%%%%%%% IMPORTANT %%%%%%%%%%%%%%%%%%%%%%
% Add SOSTOOLS path here 
% Add the phase portrait plotter path here (uncomment the
% part of the code that does the plotting, see below)
addpath(genpath(pwd))
%%%%%%%%%%
sys_indx = 1;
%
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
pmin = 1e-6; pmax = 1; %min and max eigenvalues of P for shift coordinate computation
lossless_tol = 1e-3; % lossless constraint tolerance
%%%%%%%%
solver_indx = 2; %% 1 = sdpt3, 2 = sedumi
opt_specs.solver_indx = solver_indx;
opt_specs.pmin = pmin; opt_specs.pmax = pmax;
opt_specs.lossless_tol = lossless_tol;
%%%%%%
G = get_gen_lossless_matrix(n,Qm); % Generalized lossless constraint matrix
%%%%%%%%%%%%%%%%%%%%
if ~P_indx
    %%%%%% Compute/specify the shift coordinate
    if sys_indx ~=1
        shift_indx = 1;
    else % zero shift @ sys_indx == 1 (2D toy problem (SCL et al., IJRNC 2025))
        shift_indx = 0;
    end
    if shift_indx
        P = Ix;
        if norm(G'*P(:)) > lossless_tol %should be ~0 if feasible
            fprintf("Nonlinearity non-lossless with respect to identity: " + ...
                    "a spherical global trapping region does not exist for this system. " + ...
                    "Update the inner product matrix choice and try again. " + ...
                    "Or, do a bisection to possibly resolve the issue");
            return;
        else
            shift_soln = shift_SDP(sys_info,opt_specs,P);
            m = shift_soln.m;
            Lm = shift_soln.Lm;
        end
    else
        m = zeros(n,1);
        Lm = A;
    end
else  % Generalized lossless case 
    if 0 % check if S = Ix is a feasible solution for the lossless property
        S = Ix;
        norm(G'*S(:)) %should be ~0 if feasible
    end
    % Bisection parameters
    a_tol = 1e-4;
    max_counter = 30;
    %%%%%
    % if sys_indx == 1
    %     shift_norm = 1;
    if sys_indx == 2
        shift_norm = 38;
    elseif sys_indx == 13
        shift_norm = 50; %15; %1e2; 
    elseif sys_indx == 14
        shift_norm = 10; 
    elseif sys_indx == 17
        shift_norm = 90;
    else
        shift_norm = []; %1;
    end
    opt_specs.shift_norm = shift_norm;
    %%%%%%%%%%%%%%%%%%%%%%%
    if sys_indx == 1 
        m = zeros(n,1);
        Lm = A;
    elseif sys_indx == 9
        m = [0;0;-10];
        Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
        Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3];
    % elseif sys_indx == 17
    %     m = [25.2; 13; 84.8];
    %     Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
    %     Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3];
    else
    %%%%%%%
        % solver_indx = 1; %% 1= sdpt3, 2 = sedumi
        % opt_specs.solver_indx = solver_indx;
        %%%%%%%%%%
        a = 0; a_new = 1; counter = 0; 
        rvec = randn(n,1);
        m = (1/norm(rvec))*rvec; %ones(n,1);
        %
        if sys_indx == 17 
            rvec = randn(n,1);
            m = (1/norm(rvec))*rvec + [20;10;85]; 
            % m = (1/norm(rvec))*rvec 
            % + [20.453400036972100;11.136932502126605;85.161544267953760];
        elseif sys_indx == 18
            rvec = randn(n,1);
            m = (1/norm(rvec))*rvec + [50;50;100]; 
            % m = (1/norm(rvec))*rvec 
            % + [20.453400036972100;11.136932502126605;85.161544267953760];
        end
        while norm(a - a_new) > a_tol && counter <= max_counter
            shift_soln_1 = ...
            shift_bisection_global(G,sys_info,opt_specs,sys_indx,m,[]);
            if ~isempty(shift_soln_1)
                Pshift = shift_soln_1.P;
                a = shift_soln_1.a;
            else
                fprintf("Numerical infeasibility, check tolerances");
                return;
            end
            %%%%%%%
            shift_soln_2 = ...
            shift_bisection_global(G,sys_info,opt_specs,sys_indx,[],Pshift);
            if ~isempty(shift_soln_2)
                m = shift_soln_2.m;
                Lm = shift_soln_2.Lm;
                a_new = shift_soln_2.a;
            else
                fprintf("Numerical infeasibility, check tolerances");
                return;
            end
          %%%%%%%%%%%%%
            counter = counter + 1;
           % [min(eig(Pshift))  shift_soln_2.a] % should be [+ve  -ve] for feasibility
        end
        %%%%%%%%%%% Check the solution of bisection
        if min(eig(Pshift)) > 0 &&  max(real(eig(Lm))) < 0 && a_new < 0 ...
           && max(abs(G'*shift_soln_2.P(:))) < lossless_tol
            m = shift_soln_2.m;
            Lm = shift_soln_2.Lm;
        else
            fprintf("Bisection didn't yield suitable shift coordinates:" + ...
                    "a global (ellipsoidal) trapping might not exist");
            return;
        end
    end
end
%% Compute trapping region
if P_indx
% Specify grids for chi1 and pmin (P >= pmin*Ix) 
% Set the SDP solver: solver_indx = 1 (sdpt3), 2 (sedumi)
%*****************************
eps_E = 100*eps;
pmax = 1;
 if sys_indx == 1 % 2D toy example, Liao et al. IJRNC, 2025
    Nchi = 1000;
    chi_grid = logspace(-1,1,Nchi)';
    % chi_grid = chi_grid(650); %chi_grid(575);
    pmin_grid = 1e-1; 
    % pmax = 10;
    solver_indx = 2;
    LMIx_tol = 1e-12;
    % eps_E = 1e-8; %100*eps;       
elseif sys_indx == 2  % Lorenz attractor
    Nchi = 1000;
    chi_grid = logspace(-0.5,0.5,Nchi)'; %logspace(-1,1,Nchi)'; %
    % chi_grid = chi_grid(790:851);
    pmin_grid = 1e-7; 
    solver_indx = 2;
    LMIx_tol = 1e-12;
    % eps_E = 1e-9; %100*eps;
elseif sys_indx == 9 % Chen system
    Nchi = 1000;
    chi_grid = logspace(-2,2,Nchi)';
    pmin_grid = 1e-9;  
    solver_indx = 2;
    LMIx_tol = 1e-12;
    % eps_E = 1e-12; %100*eps;
elseif sys_indx == 12 % Lorenz83 system
    Nchi = 1000;
    chi_grid = logspace(-1,1,Nchi)';
    pmin_grid = 1e-7; 
    solver_indx = 2;
    LMIx_tol = 1e-12;
    % eps_E = 1e-9; %100*eps;
 elseif sys_indx == 13  % Airfoil model, Heide et al., 2025
    Nchi = 1000;
    chi_grid = logspace(0,2,Nchi)';
    % chi_grid = chi_grid(400:600);
    pmin_grid = 1e-7;
    solver_indx = 2;
    LMIx_tol = 1e-12;
    % eps_E = 100*eps; %1e-12; %
elseif sys_indx == 14 % 2D example
    Nchi = 1000;
    chi_grid = logspace(-1,0,Nchi)';
    pmin_grid = 1e-2; 
    solver_indx = 2;
    LMIx_tol = 1e-12;
    % eps_E = 1e-9; %100*eps;
elseif sys_indx == 17 % Three-scroll chaotic attractor
    Nchi = 1000;
    chi_grid = logspace(-1,1,Nchi)'; 
    pmin_grid = 1e-7; 
    solver_indx = 2;
    LMIx_tol = 1e-12;
    % eps_E = 1e-12; %100*eps;
 else
    Nchi = 1000;
    chi_grid = logspace(-2,2,Nchi)';
    pmin_grid = 1e-6; 
    solver_indx = 2;
    LMIx_tol = 1e-12;
    % eps_E = 1e-9; %100*eps;
 end
opt_specs.pmax = pmax;
opt_specs.LMIx_tol = LMIx_tol;
opt_specs.solver_indx = solver_indx;
opt_specs.eps_E = eps_E;
%%%==============
%% Find trapping region over a grid of the multiplier chi_1
%%%============== 
trace_store = zeros(length(chi_grid),length(pmin_grid));
radius_store = zeros(length(chi_grid),length(pmin_grid));
ROBB_store = radius_store;
soln_store = cell(length(chi_grid),length(pmin_grid));
%%%%%%%%%%%%%%%%%%%%%%
for ct2 = 1: length(pmin_grid)
    pmin = pmin_grid(ct2);
    opt_specs.pmin = pmin; 
    %%%%**************
    for ct1=1:length(chi_grid)
        chi1 = chi_grid(ct1);
        soln = global_TR_genloss_ell(sys_info,G,opt_specs,chi1,m);
        soln_store{ct1,ct2} = soln;
        %%%*************
        if isempty(soln)
            trace_store(ct1,ct2) = NaN;
            radius_store(ct1,ct2) = NaN;
            ROBB_store(ct1,ct2) = NaN;
        else
            P = soln.P;
            trace_store(ct1,ct2) = trace(P);
            radius_store(ct1,ct2) = soln.radius;
            ROBB_store(ct1,ct2) = soln.ROBB;
        end
    end
%%%%%%%%%%
%% Local search in the shift coordinates
    shift_bisection = 1;
    local_search_success = 0; % index related to local shift coordinate search
    w_deltam = 1e-6;
    opt_specs.w_deltam = w_deltam;
    %
    [~, opt_indx] = min(ROBB_store(:,ct2)); % min(radius_store);
    opt_soln_grid = soln_store{opt_indx,ct2};
    opt_soln = opt_soln_grid; % this is the optimal solution on the grid (without adjusting the shift coordinates)
    if ~isempty(opt_soln) 
        if shift_bisection 
        % Optimize the shift coordinates with a linearized c(m)
            opt_soln_shift1 = ...
            global_TR_genloss_shift(sys_info,opt_specs,opt_soln.chi1,opt_soln.m,opt_soln.P);
            %%%%%%%%%%%%%%%
            if ~isempty(opt_soln_shift1)
            % Optimize P and radius for the above shift coordinates
            % without any approximation in c(m)
                opt_soln_shift2 = ...
                global_TR_genloss_ell(sys_info,G,opt_specs,opt_soln.chi1,opt_soln_shift1.m);
                if ~isempty(opt_soln_shift2) && opt_soln_shift2.ROBB < opt_soln.ROBB
                    opt_soln = opt_soln_shift2;
                    fprintf("Local search in shift coordinates improves results:" + ...
                            "optimal solution updated");
                    local_search_success = 1;
                end
            end
        end
        opt_soln_store{ct2} = opt_soln;
        ROBB_opt(ct2) = opt_soln.radius/sqrt(min(eig(opt_soln.P))); 
        % maximum 2-norm on the trapping ellipsoid or largest semi-major axis length
        % which equals the ROBB (= Radius of Outer Bounding Ball)
        %%% ROBB = compute_ROBB(opt_soln);
        Ult_bound(ct2) = norm(opt_soln.m) + ROBB_opt(ct2); %this is the ultimate bound in the original coordinates
    else
        ROBB_opt(ct2) = NaN; Ult_bound(ct2) = NaN;
        fprintf("Infeasibility, check the choice of shift coordinates or the choice of P (if fixed)");
    end
end
%% Check the final solution for accuracy
if ~isempty(opt_soln)
    [max(abs(G'*opt_soln.P(:)))  max(eig(opt_soln.LMIx))  max(eig(opt_soln.LMIx(1:n,1:n)))  min(eig(opt_soln.P))]
    %%% Should be [~0  -ve(~0)  -ve(~0)  +ve]
end
%%
if length(pmin_grid) > 1
    figure()
    loglog(pmin_grid, ROBB_opt,'ko','LineWidth',2);
    %
    xlabel('$\mathrm{p_{min}}$','Interpreter','latex');
    ylabel('$\mathrm{Optimal \ ROBB}$','Interpreter','latex');
    grid on;
    set(gca,'Fontsize',15);
    %%%%%%%%
    figure()
    loglog(pmin_grid, Ult_bound,'ko','LineWidth',2);
    %
    xlabel('$\mathrm{p_{min}}$','Interpreter','latex');
    ylabel('$\mathrm{Ultimate \ bound}$','Interpreter','latex');
    grid on;
    set(gca,'Fontsize',15);
end
%%  
ct_plot_pmax = 1;
%
if ~isempty(opt_soln)
    figure()
    loglog(chi_grid, trace_store(:,ct_plot_pmax),'b.','LineWidth',2);
    %
    xlabel('$\chi$','Interpreter','latex');
    ylabel('$\mathrm{trace}(P)$','Interpreter','latex');
    grid on;
    set(gca,'Fontsize',15);
    % xlim([min(chi_grid) max(chi_grid)]);
    % ylim ([])
    % xticks([])

    figure()
    loglog(chi_grid, ROBB_store(:,ct_plot_pmax),'bo','LineWidth',1.5);
    %
    xlabel('$\chi$','Interpreter','latex');
    ylabel('Trapping Radius (Spherical)','Interpreter','latex');
    grid on;
    set(gca,'Fontsize',15);
   % title('Radius of Outer Bounding Ball');
   % xlim([min(chi_grid) max(chi_grid)]);
    
    figure()
    % loglog(chi_grid, gam_store,'bd','LineWidth', 1.5);
    %
    % hold on;
    loglog(chi_grid, radius_store(:,ct_plot_pmax),'ro','LineWidth', 1.5);
    xlabel('$\chi$','Interpreter','latex');
    ylabel('Trapping Radius (Ellipsoidal)','Interpreter','latex');
    % ylabel('$\gamma, r$','Interpreter','latex');
    grid on;
    set(gca,'Fontsize',15);
    % xlim([min(chi_grid) max(chi_grid)]);
    % legend('$\gamma$','$r$','Interpreter','latex')
    %%%%%%%%%
end
else
    solver_indx = 2;
    LMIx_tol = 1e-12;
    opt_soln = global_TR_SCL(m,A,Qm,d,solver_indx,LMIx_tol);
    opt_soln.P = Ix;
    if ~isempty(opt_soln)
        R_opt = opt_soln.radius  % This is the optimal spherical trapping radius 
        opt_soln.ROBB = R_opt;   % needed for Vdot computations next
        Ult_bound = norm(m) + R_opt
    else
        fprintf("Infeasibility, check the choice of shift coordinates or the choice of P (if fixed)");
    end
end
%% Check Vdot on the boundary of the trapping region
if ~isempty(opt_soln)
    m_opt = opt_soln.m;
    if n == 2 
        Q1 = Qm{1}; Q2 = Qm{2};
        Lm = A + 2*[m_opt'*Q1; m_opt'*Q2];
        cm = d + A*m_opt + [m_opt'*Q1*m_opt; m_opt'*Q2*m_opt];
    elseif n == 3 
        Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
        Lm = A + 2*[m_opt'*Q1; m_opt'*Q2; m_opt'*Q3];
        cm = d + A*m_opt + [m_opt'*Q1*m_opt; m_opt'*Q2*m_opt; m_opt'*Q3*m_opt]; 
    elseif n == 9 
        Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
        Q4 = Qm{4}; Q5 = Qm{5}; Q6 = Qm{6};
        Q7 = Qm{7}; Q8 = Qm{8}; Q9 = Qm{9};
        Lm = A + 2*[m_opt'*Q1; m_opt'*Q2; m_opt'*Q3;...
                    m_opt'*Q4; m_opt'*Q5; m_opt'*Q6;...
                    m_opt'*Q7; m_opt'*Q8; m_opt'*Q9];
        cm = d + A*m_opt + [m_opt'*Q1*m_opt; m_opt'*Q2*m_opt; m_opt'*Q3*m_opt;...
                            m_opt'*Q4*m_opt; m_opt'*Q5*m_opt; m_opt'*Q6*m_opt;...
                            m_opt'*Q7*m_opt; m_opt'*Q8*m_opt; m_opt'*Q9*m_opt];
    elseif n == 6
        Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
        Q4 = Qm{4}; Q5 = Qm{5}; Q6 = Qm{6};
        Lm = A + 2*[m_opt'*Q1; m_opt'*Q2; m_opt'*Q3;...
                    m_opt'*Q4; m_opt'*Q5; m_opt'*Q6];
        cm = d + A*m_opt + [m_opt'*Q1*m_opt; m_opt'*Q2*m_opt; m_opt'*Q3*m_opt;...
                            m_opt'*Q4*m_opt; m_opt'*Q5*m_opt; m_opt'*Q6*m_opt];
    end
    N_eps = 1e6; 
    Vdot1 = zeros(N_eps,1); Vdot2 = Vdot1;
    Rvec = randn(n,N_eps); %-ones(n,N_eps) + 2*rand(n,N_eps); % 
    opt_radius = opt_soln.radius;
    P_ell =  opt_soln.P/(opt_radius^2);
    P_ell = full(P_ell);
    E_ell = sqrtm(inv(P_ell));
    for ct = 1:N_eps
        rvec = Rvec(:,ct); 
        rvec = (1/norm(rvec))*rvec;
        Y = E_ell*rvec; % [Y'*P_ell*Y - 1 ~ 0]
        if n == 2
            fy = [Y'*Q1*Y; Y'*Q2*Y];
        elseif n == 3
            fy = [Y'*Q1*Y; Y'*Q2*Y; Y'*Q3*Y];
        elseif n == 6
            fy = [Y'*Q1*Y; Y'*Q2*Y; Y'*Q3*Y; Y'*Q4*Y; Y'*Q5*Y; Y'*Q6*Y];
        elseif n == 9
            fy = [Y'*Q1*Y; Y'*Q2*Y; Y'*Q3*Y;...
                  Y'*Q4*Y; Y'*Q5*Y; Y'*Q6*Y;...
                  Y'*Q7*Y; Y'*Q8*Y; Y'*Q9*Y];
        end
        innPd = Y'*P_ell*fy;
        % if abs(innPd) < 1e-2 
        ydot = Lm*Y + fy + cm;
        % Vdot(ct,1) = ydot'*opt_soln.P*Y + Y'*opt_soln.P*ydot;
        Vdot1(ct,1) = ydot'*P_ell*Y + Y'*P_ell*ydot; % true derivative
        Vdot2(ct,1) = Y'*(Lm'*P_ell + P_ell*Lm)*Y + Y'*P_ell*cm + cm'*P_ell*Y; % approximate, artificially applying the lossless inner product
        % end
    end
    max(Vdot2) % should be <0 for set invariance
    % max(Vdot1) % This might become positive due to the inner product not being exactly zero numerically
end
%%%%%%%%%%%%% Uncomment when using the phase portrait plotter %%%%%%%%%%%
% %% 2D/3D phase portrait and the associated trapping region
% if sys_indx == 1
%     ODEfunc = @(t,x)[-x(1) - x(1)*x(2);...
%                      -4*x(2) + x(1)^2  + 1];
%     % figure()
%     plotpp(ODEfunc, 'xlim', [-3,3], 'ylim', [-3,3]);
%     set(gca,'fontsize',15);
%     xlabel('$x_1$', 'interpreter', 'latex','FontSize',15);
%     ylabel('$x_2$', 'interpreter', 'latex','FontSize',15);
%     hold on;
%     th = (0:0.01:2*pi)';
%     %
%     cir_x = m_opt(1) + opt_soln.ROBB*cos(th);
%     cir_y = m_opt(2) + opt_soln.ROBB*sin(th);
%     plot(cir_x, cir_y, 'r','LineWidth', 1.5);
%     axis equal;
% %%%%%%%%
% elseif sys_indx == 14
%     ODEfunc = @(t,x)[x(1) + 2*x(1)*x(2) - x(2)^2;...
%                      -x(2) + x(1)*x(2)   - 2*x(1)^2];
%     % figure()
%     plotpp(ODEfunc, 'xlim', [-3,3], 'ylim', [-3,3]);
%     set(gca,'fontsize',15);
%     xlabel('$x_1$', 'interpreter', 'latex','FontSize',15);
%     ylabel('$x_2$', 'interpreter', 'latex','FontSize',15);
%     hold on;
%     th = (0:0.01:2*pi)';
%     %
%     cir_x = m_opt(1) + opt_soln.ROBB*cos(th);
%     cir_y = m_opt(2) + opt_soln.ROBB*sin(th);
%     plot(cir_x, cir_y, 'r','LineWidth', 1.5);
%     axis equal;
% end