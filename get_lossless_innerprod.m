% Computes the lossless inner product x'*P*Q(x) 
function inner_prod = get_lossless_innerprod(n,Qm,P)
%%%%%
Nrand = 1e6;
IPs = NaN(Nrand, 1); 
for ct = 1:Nrand
    X = randn(n,1);
    if n==2
        Q1 = Qm{1};
        Q2 = Qm{2};
        fx = [X'*Q1*X; X'*Q2*X];
    elseif n==3
        Q1 = Qm{1};
        Q2 = Qm{2};
        Q3 = Qm{3};
        % 
        fx = [X'*Q1*X; X'*Q2*X; X'*Q3*X];
    elseif n==4
        Q1 = Qm{1};
        Q2 = Qm{2};
        Q3 = Qm{3};
        Q4 = Qm{4};
        % 
        fx = [X'*Q1*X; X'*Q2*X; X'*Q3*X; X'*Q4*X];
    elseif n==6
        Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
        Q4 = Qm{4}; Q5 = Qm{5}; Q6 = Qm{6};
        % 
        fx = [X'*Q1*X; X'*Q2*X; X'*Q3*X;...
              X'*Q4*X; X'*Q5*X; X'*Q6*X];
    elseif n==9
        Q1 = Qm{1}; Q2 = Qm{2}; Q3 = Qm{3};
        Q4 = Qm{4}; Q5 = Qm{5}; Q6 = Qm{6};
        Q7 = Qm{7}; Q8 = Qm{8}; Q9 = Qm{9};
        % 
        fx = [X'*Q1*X; X'*Q2*X; X'*Q3*X;...
              X'*Q4*X; X'*Q5*X; X'*Q6*X;...
              X'*Q7*X; X'*Q8*X; X'*Q9*X];
    end
    IPs(ct,1) = X'*P*fx;
end
inner_prod = mean(IPs); % return the average of the randomly computed values