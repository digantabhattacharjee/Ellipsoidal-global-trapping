%%%%%%%%%%%%%%%%
% Computes the matrix asscoiated with recasting the generalized lossless 
% condition into a set of linear equality constraints 
% for a given quadratic nonlinearity Q
% Uses custom functions from SOSTOOLS
% Written by PJS, Fall 2025
%%%%%%%%%%%%%
function V = get_gen_lossless_matrix(nx,Qm)

S = mpvar('S',[nx nx]);
x = mpvar('x',[nx,1]);
% 
if nx == 2
    f = [x'*Qm{1}*x;x'*Qm{2}*x];
elseif nx == 3
    f = [x'*Qm{1}*x;x'*Qm{2}*x;x'*Qm{3}*x];
elseif nx == 4
    f = [x'*Qm{1}*x;x'*Qm{2}*x;x'*Qm{3}*x; x'*Qm{4}*x];
elseif nx == 6
    f = [x'*Qm{1}*x;x'*Qm{2}*x;x'*Qm{3}*x;x'*Qm{4}*x;x'*Qm{5}*x;x'*Qm{6}*x];
elseif nx == 9
    f = [x'*Qm{1}*x;x'*Qm{2}*x;x'*Qm{3}*x;x'*Qm{4}*x;x'*Qm{5}*x;x'*Qm{6}*x;...
         x'*Qm{7}*x; x'*Qm{8}*x; x'*Qm{9}*x];
end

%% Collect entries of S
%  Collect p(x,S) into the form g0(S)+g(S)*h(x) where h(x) is a vector
%  of unique monomials in x. In this case g0=0.  The entries of
%  g(S) must sum to zero.
p = x'*S*f;
[g0,g,h] = collect(p,S);

%% Projection
%  Project a vector of polynomials g onto the span of the monomials
%  contained in the vector R.  In this case R=S(:) and we have
%  g = R'* V.  We can flip this around to g' = V' *R.
R = S(:);
V = poly2basis(g,R);
V = full(V);