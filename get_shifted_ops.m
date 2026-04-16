%%%*****************
% Computes the linear operator and constant vector for a given system and
% shift coordinates
%%%*********************************
function [Lm,cm] = get_shifted_ops(sys_info,m)
A = sys_info.A;
d = sys_info.d;
Qm = sys_info.Qm;
%%%%%%%%%%%%
if size(A,2) == 2 
    Q1 = Qm{1};
    Q2 = Qm{2};
    Lm = A + 2*[m'*Q1; m'*Q2];
    cm = d + A*m + [m'*Q1*m; m'*Q2*m];
elseif size(A,2) == 3 
    Q1 = Qm{1};
    Q2 = Qm{2};
    Q3 = Qm{3};
    Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3];
    cm = d + A*m + [m'*Q1*m; m'*Q2*m; m'*Q3*m];
elseif size(A,2) == 4 
    Q1 = Qm{1};
    Q2 = Qm{2};
    Q3 = Qm{3};
    Q4 = Qm{4};
    Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3; m'*Q4];
    cm = d + A*m + [m'*Q1*m; m'*Q2*m; m'*Q3*m; m'*Q4*m];
elseif size(A,2) == 9 
    Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
    Q4 = Qm{4}; Q5 = Qm{5}; Q6 = Qm{6};
    Q7 = Qm{7}; Q8 = Qm{8}; Q9 = Qm{9};
    Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3;...
                m'*Q4; m'*Q5; m'*Q6;...
                m'*Q7; m'*Q8; m'*Q9];
    cm = d + A*m + [m'*Q1*m; m'*Q2*m; m'*Q3*m;...
                    m'*Q4*m; m'*Q5*m; m'*Q6*m;...
                    m'*Q7*m; m'*Q8*m; m'*Q9*m];
elseif size(A,2) == 6
    Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
    Q4 = Qm{4}; Q5 = Qm{5}; Q6 = Qm{6};
    Lm = A + 2*[m'*Q1; m'*Q2; m'*Q3;...
                m'*Q4; m'*Q5; m'*Q6];
    cm = d + A*m + [m'*Q1*m; m'*Q2*m; m'*Q3*m;...
                    m'*Q4*m; m'*Q5*m; m'*Q6*m];
end