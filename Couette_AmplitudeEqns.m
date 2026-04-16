function [E,L,Q_out] = Couette_AmplitudeEqns(Lx,Lz,Re)
%% Find Derivatives of Temporal Modes for Moehlis 9-state Couette flow
% Leonid Heide,2023
%
% This model is obtained from the paper:
%   [1]: "Periodic Ornits and Chaotic Sets in a Low-Dimensional Model for Shear
%   Flows" by J.Moehlis, H. Faisst, B.Eckhardt (2005)/
%   https://sites.me.ucsb.edu/~moehlis/moehlis_papers/siads.pdf
%
% To obtain the Matrix form, the following paper was referenced:
%    [2]: "Global Stability Analysis of Fluid Flows using Sum-of-Squares" 
%    by P.Goulart, S. Chernyshenko (2011)
%    https://arxiv.org/pdf/1101.1043.pdf
%
% This function defines the amplitude equations (derivatives of temporal
% modes) for a nine-state couette flow model obtained from the reference.
%
% The amplitudes are a result of using the modes defined for this model and
% applying the ansatz:
%
%   U(x,t)=SUM_m(a_m(t)*u_m(x));
%
% where U(x,t) us the velocity field, a_m are the temporal modes (aka
% amplitudes) and u_m(x) are the spatial modes. Using this ansatz and
% performing galerkin projection yields:
%
%   da_i/dt= (beta^2/Re)*del_(i,1) - (d_i/Re)*a_i +SUM_{j,k}(N_{i,jk}a_ja_k
%
% Note that the nonlinear terms here satisfy energy equation st.
%   SUM_{ijk}(N_{i,jk}*a_i*a_j*a_k = 0
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% OUTPUTS:
%   E: a vector containing the affine terms
%   L: a matrix containing the linear terms on the diagonal
%   Q_out: a tensor containing the nonlinear terms.
%
%   The outputs define the system of ODEs as follows:
%   adot(jj)= E(jj)+L(:,jj)'*a + a'*Q(:,:,jj)'*a;
%       Where jj indicated the mode number. 
%
% INPUTS:
%   Lx: x-domain size (streamwise) s.t. 0 <=x<=Lx
%   Lz: z-domain size (spanwise) s.t. 0<=z<=Lz
%   Re: Reynolds number
%   
% NOTES:
%   -Refer to reference [1] to determine what values to pick for inputs
%   -Wall normal direction (y) is constrained to be -1<=y<=1

%% Load Variables

% Unpack the Wave Vector variables and Reynolds Number
alpha=2*pi/Lx; % spanwise
beta=pi/2; % wall normal
gamma=2*pi/Lz; % streamwise

% For simplification define the following constants:
k_ag=sqrt(alpha^2+gamma^2);
k_bg=sqrt(beta^2+gamma^2);
k_abg=sqrt(alpha^2+beta^2+gamma^2);
k_ab=sqrt(alpha^2+beta^2);

%% Initialize Matricies:

E=zeros(9,1);
L=zeros(9,9);
Q=zeros(9,9,9);
%% Find Affine Terms (E):

% Amplitude corresponding to basic profile mode
E(1)= (beta^2/Re);
E(2:9)=0;

%% Find Linear Terms (L):

% MODE 1: Amplitude corresponding to basic profile mode
L(1,1)=-(beta^2);

% MODE 2: Amplitude corresponding to the streak mode
L(2,2)=-((4*beta^2)/3 + gamma^2);

% MODE 3: Amplitude corresponding to the downstream vortex mode
L(3,3)=-k_bg^2;

% MODE 4: Amplitude corresponding to the spanwise flow mode #1
L(4,4)=-(3*alpha^2+4*beta^2)/(3);

% MODE 5: Amplitude corresponding to the spanswise flow mode #2
L(5,5)=-k_ab^2;

% MODE 6: Amplitude corresponding to normal vortex mode #1
L(6,6)=-(3*alpha^2+4*beta^2+3*gamma^2)/(3);

% MODE 7: Amplitude corresponding to normal vortex mode #2
L(7,7)=-k_abg^2;

% MODE 8: Amplitude corresponding to the 3D advection mode
L(8,8)= -k_abg^2;

% MODE 9: Amplitude corresdponing to the mode modifying the basic profile
L(9,9)=-(9*beta^2);


% Divide L by Reynolds number (as shown in [1,2]):
L=L/Re;

%% Find Quadratic Terms (Q):

% MODE 1: Amplitude corresponding to basic profile mode
Q(2,3,1)= sqrt(3/2)*beta*gamma/k_bg;
Q(6,8,1)= -sqrt(3/2)*beta*gamma/k_abg;

% MODE 2: Amplitude corresponding to the streak mode
Q(4,6,2)= 10*gamma^2/(3*sqrt(6)*k_ag);
Q(5,7,2)= -(gamma^2)/sqrt(6)/k_ag;
Q(5,8,2)= -alpha*beta*gamma/(sqrt(6)*k_ag*k_abg);
Q(1,3,2)= -sqrt(3/2)*beta*gamma/k_bg;
Q(3,9,2)= -sqrt(3/2)*beta*gamma/k_bg;

% MODE 3: Amplitude corresponding to the downstream vortex mode
Q(4,7,3)= sqrt(2/3)*(alpha*beta*gamma)/(k_ag*k_bg);
Q(5,6,3)= sqrt(2/3)*(alpha*beta*gamma)/(k_ag*k_bg);
Q(4,8,3)= (beta^2*(3*alpha^2+gamma^2)-3*gamma^2*k_ag^2)/(sqrt(6)*k_ag*k_bg*k_abg);


% MODE 4: Amplitude corresponding to the spanwise flow mode #1
Q(1,5,4)= -alpha/sqrt(6);
Q(5,9,4)= -alpha/sqrt(6);
Q(2,6,4)= -(10*alpha^2)/(3*sqrt(6)*k_ag);
Q(3,7,4)= -sqrt(3/2)*(alpha*beta*gamma)/(k_ag*k_bg);
Q(3,8,4)= -sqrt(3/2)*(alpha^2*beta^2)/(k_ag*k_bg*k_abg);


% MODE 5: Amplitude corresponding to the spanswise flow mode #2
Q(1,4,5)= alpha/sqrt(6);
Q(4,9,5)= alpha/sqrt(6);
Q(2,7,5)= (alpha^2)/sqrt(6)/k_ag;
Q(2,8,5)= -alpha*beta*gamma/(sqrt(6)*k_ag*k_abg);
Q(3,6,5)= sqrt(2/3)*(alpha*beta*gamma)/(k_ag*k_bg);

% MODE 6: Amplitude corresponding to normal vortex mode #1
Q(1,7,6)= alpha/sqrt(6);
Q(7,9,6)= alpha/sqrt(6);
Q(1,8,6)= sqrt(3/2)*(beta*gamma/k_abg);
Q(8,9,6)= sqrt(3/2)*(beta*gamma/k_abg);
Q(2,4,6)= (10/(3*sqrt(6)))*(alpha^2-gamma^2)/(k_ag);
Q(3,5,6)= -sqrt(2/3)*2*(alpha*beta*gamma)/(k_ag*k_bg);

% MODE 7: Amplitude corresponding to normal vortex mode #2
Q(1,6,7)= -alpha/sqrt(6);
Q(6,9,7)= -alpha/sqrt(6);
Q(2,5,7)= (gamma^2-alpha^2)/sqrt(6)/k_ag;
Q(3,4,7)= alpha*beta*gamma/(sqrt(6)*k_ag*k_bg);


% MODE 8: Amplitude corresponding to the 3D advection mode
Q(2,5,8)= sqrt(2/3)*(alpha*beta*gamma)/(k_ag*k_abg);
Q(3,4,8)= gamma^2*(3*alpha^2-beta^2+3*gamma^2)/(sqrt(6)*k_ag*k_bg*k_abg);

% MODE 9: Amplitude corresdponing to the mode modifying the basic profile
Q(2,3,9)= sqrt(3/2)*(beta*gamma)/k_bg;
Q(6,8,9)= -sqrt(3/2)*(beta*gamma)/k_abg;

%% Reformat Q:

% Make Q a symmetric matrix (need to divide by 2 to do this):

Q_out=zeros(9,9,9);
for i=1:9
    Q_out(:,:,i)=(1/2)*Q(:,:,i)+(1/2)*Q(:,:,i)';
end


end