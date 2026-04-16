%%%%%%%%%%%
% Generate system-related information according to specified system index
%%%%%%%
function sys_info = get_system(sys_indx)
%%%%%%%%%%
if sys_indx == 1 %% 2D toy problem (see Liao et al., 2025)
    n = 2;
    % 
    A = [-1  0;
          0 -4];
    d = [0; 1];
    % 
    Q1 = [0 -0.5;
         -0.5  0];
    Q2 = [1 0;
          0 0];
    Qm = {Q1,Q2};
elseif sys_indx == 2 %% Lorenz attractor (see Liao et al., 2025)
    n = 3; sig = 10; alph_Lorenz = 8/3; rh = 28;
    % 
    A = [-sig    sig     0;
          rh     -1      0;
           0      0   -alph_Lorenz];
    d = zeros(n,1);
    % 
    Q1 = zeros(n,n);
    Q2 = [  0     0  -0.5;
            0     0    0;
          -0.5    0    0];
    Q3 = [ 0   0.5  0;
          0.5  0    0;
           0    0   0];
    Qm = {Q1,Q2,Q3};  
    %%%%%%%%%%%%%%%%%%%%
elseif sys_indx == 4  %% 9-state model
    model = model_Couette9;
    A = model.L;
    Qm = model.Q;
    Q1 = Qm(:,:,1); Q2 = Qm(:,:,2); Q3 = Qm(:,:,3);
    Q4 = Qm(:,:,4); Q5 = Qm(:,:,5); Q6 = Qm(:,:,6);
    Q7 = Qm(:,:,7); Q8 = Qm(:,:,8); Q9 = Qm(:,:,9);
    Qm = {Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9};
    %%%%%%%%%%%%%%%%
    d = model.c;
    n = model.nx; 
elseif sys_indx == 5 %% 2D toy problems 
    n = 2;
    % 
    A = -eye(n);
    d = zeros(n,1);
    Q1 = [0   -0.5;
          -0.5  0];
    Q2 = [1   0;
           0   0];
    Qm = {Q1,Q2};
elseif sys_indx == 6 %% 2D toy problems 
    n = 2;
    % 
    A = [1  0;
         0 -1];
    d = zeros(n,1);
    %%%%%%%%%%%%
    Q1 = [0   0.5;
          0.5  0];
    Q2 = [-1   0;
           0   0];
    %%%%%%%%%%%%%%%%
    Qm = {Q1,Q2};
elseif sys_indx == 12  %Lorenz83 (from https://www.dynamicmath.xyz/strange-attractors/)
    n = 3; a = 0.95; b = 7.91; f = 4.83; g = 4.66;
    % 
    A = diag([-a,-1,-1]);
    d = [a*f;g;0];
    % 
    Q1 = [0     0     0;
          0    -1     0;
          0     0    -1];
    
    Q2 = [0      1/2    -b/2;
          1/2      0     0;
          -b/2   0    0];
    
    Q3 = [0      b/2    1/2;
          b/2     0      0;
          1/2     0      0];
    %%%%%%%%%%
    Qm = {Q1,Q2,Q3};
elseif sys_indx == 13  %% 6-state unsteady aerodynamics model (Heide et al. 2025)
    [A, Qm, d, n] = get_6state_model;
    %%%%%%%%
elseif sys_indx == 14 %% 2D toy problem 
    n = 2;
    % 
    A = [1  0;
          0 -1];
    d = zeros(n,1);
    % 
    Q1 = [0 1;
          1 -1];
    Q2 = [-2 1/2;
          1/2 0];
    %
    Qm = {Q1,Q2};
    %%%%%%%%%%%%%%%%%%%%%%%%%
elseif sys_indx == 19 %% Three-dimensional model for cylinder flow (also called mean field model)
    % (Source: Kaptanoglu et al. 2021, Kalur et al 2023, VKC et al 2025)
    n = 3; gam = 0.1;
    % 
    A = [gam   -1  0;
         1   gam  0;
         0   0  -1];
    d = zeros(n,1);
    % 
    Q1 = [0     0     -1/2;
          0     0     0;
          -1/2     0     0];
    
    Q2 = [0       0    0;
          0       0     -1/2;
          0    -1/2     0];
    
    Q3 = [1     0        0;
          0     1        0;
          0     0        0];
    %%%%%%%%%%
    Qm = {Q1,Q2,Q3};
    %%%%%%%%%%%%%%%%%%
elseif sys_indx == 21 %% Finance chaotic system, Peng et al, POF 2025
    n = 3; alph = 0.0001; bet = 0.2; sig = 1.1;
    % 
    A = [(1/bet - alph)          0      1;
         0                    -bet      0;
         -1                        0      -sig];
    d = [0;1;0];
    % 
    Q1 = zeros(n); Q1(1,2) = 1/2; Q1(2,1) = 1/2;
    
    Q2 = zeros(n,n); Q2(1,1) = -1;
    
    Q3 = zeros(n,n);
    %%%%%%%%%%
    Qm = {Q1,Q2,Q3};
elseif sys_indx == 22 %% Atmospheric oscillator model (see Kaptanoglu et al., 2021)
    n = 3; mu1 = 0.05; mu2 = -0.01; omega = 3; sig = 1.1; kap = -2; 
    bt = -6;
    % 
    A = [mu1      0          0;
          0      mu2      omega;
           0   -omega      mu2];
    d = zeros(n,1);
    % 
    Q1 = zeros(n,n); Q1(1,2) = sig/2; Q1(2,1) = Q1(1,2);
    Q2 = [-sig     0       0;
           0       0      kap/2;
           0     kap/2     bt];
    Q3 = [0       0     0;
          0     -kap   -bt/2;
          0    -bt/2      0];
    %%%%%%%%%%%%%
    Qm = {Q1,Q2,Q3};
    %%%%%%%%%%%%%%%%
elseif sys_indx == 23
    %%% Stochastic Lorenz–Stenflo system (see Peng et al., 2025)
    %%% Nonlinearity modified
    n = 4; sig = 2; rr = 26; s = 1.5;
    bt = 0.7;
    % 
    A = [-sig    sig        0   s;
          rr      -1        0   0;
           0      0        -bt   0;
           -1     0        0    -sig];
    d = zeros(n,1);
    % 
    Q1 = zeros(n,n); 
    Q2 = Q1; Q2(1,3) = -1/2; Q2(3,1) = -1/2;
    % Q3 = Q1; Q3(1,2) = 1/2; Q3(2,1) = 1/2;
    %%%%%%% Modified nonlinearity in the third state, Q3
    Q3 = Q1; Q3(1,2) = 1; Q3(2,1) = 1;
    %%%%%
    Q4 = zeros(n,n);
    %%%%%%%%%%%%%
    Qm = {Q1,Q2,Q3,Q4};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif sys_indx == 26  %% Hadley cell or Hadley circulation model, Peng et al. 2025
    n = 3;
    a = 0.2; b = 4; f = 9; g = 1;
    A = [-a    0        0;
         0     -1       0;
         0      0      -1];
    d = [a*f;g;0];
    % 
    Q1 = [0  0   0;
          0  -1  0;
          0  0  -1];
    Q2 = [0    1/2  -b/2;
          1/2   0     0;
         -b/2   0     0];
    Q3 = [0    b/2    1/2;
          b/2   0      0;
          1/2   0      0];
    Qm = {Q1,Q2,Q3};
end
%%%%% Package info into a structure and output that
sys_info.n = n;
sys_info.d = d;
sys_info.A = A;
sys_info.Qm = Qm;

