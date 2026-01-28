%%%%%%%%%%%
% Generate system-related information according to specified system index
%%%%%%%
function sys_info = get_system(sys_indx)
%% Define system matrices, dimension and shift coordinates
if sys_indx == 1
    %%%%%%%%%% 2D toy problem (see Liao et al., 2025)
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
elseif sys_indx == 2
    %%%%%%%%%% Lorenz attractor (see Liao et al., 2025)
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
    elseif sys_indx == 4  %%%%%%%%%% 9-state model
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
    %%%%%%%%%%%%%%%%%%%% 
elseif sys_indx == 9 % Chen 
    % (from https://www.dynamicmath.xyz/strange-attractors/)
    n = 3;
    a = 5;
    b = -10;
    c = -0.38;
    % 
    A = diag([a,b,c]);
    d = [0;0;0];
    % 
    Q1 = -[0 0 0;
          0 0 0.5;
          0 0.5 0];
    Q2 = [0 0 0.5;
          0 0  0;
          0.5 0 0];
    Q3 = [0    1/6  0;
          1/6    0   0;
         0   0   0];
    %%%%%%%%%%%%%
    Qm = {Q1,Q2,Q3};
       %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif  sys_indx == 12  %Lorenz83 
    % (from https://www.dynamicmath.xyz/strange-attractors/)
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
    %
    Qm = {Q1,Q2,Q3};
        %%%%%%%%%%%%%%%%%%%%%
elseif  sys_indx == 13  %6-state model (Heide et al. 2025)
    % %%%%%%%%%%
    [A, Qm, d, n] = get_6state_model;
    %%%%%%%%
elseif sys_indx == 14
    %%%%%%%%%% 2D toy problem 
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
    % %%%%%%%%
elseif sys_indx == 17 % three-scroll chaotic attractor: 
    % from https://www.dynamicmath.xyz/strange-attractors/
    n = 3; 
    a = 32.48; b = 45.84; c = 1.18; d_p = 0.13; e = 0.57; f = 14.7; 
    % 
    A = [-a  a  0;
         b   f  0;
         0   0  c];
    d = zeros(n,1);
    % 
    Q1 = [0     0     d_p/2;
          0     0     0;
          d_p/2     0     0];

    Q2 = [0       0    -1/2;
          0       0     0;
          -1/2    0     0];

    Q3 = [-e      1/2    0;
          1/2     0      0;
          0     0      0];
    %%%%%%%%%%
    Qm = {Q1,Q2,Q3};
        %%%%%%%%       
 elseif sys_indx == 19 % Cylinder wake model
     % From: Kalur et al 2023, VKC et al 2025
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
end
%%%%% Package info into a structure and output that
sys_info.n = n;
sys_info.d = d;
sys_info.A = A;
sys_info.Qm = Qm;
