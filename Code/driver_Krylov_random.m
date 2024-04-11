close all
clear all

n=20;

A=rand(n,n);
b=rand(n,1);
if norm(b,2)==0
    return
end

beta=norm(b,2);
q1=b/beta;
f='exp';
tol=1e-13; % Chosen tolerance for the Arnoldi process, due to potential numerical cancellations

[Q,H,HK1K]=Arnoldi(A,q1,tol);

[m, ~]=size(H); % m is the step of the computation reached
nR(1)=0; % initialize the generalized residual to 0

ek(1, 1)=1; % define the e_k vector of the canonical basis
e1(1, 1)=1; % define the e_1 vector of the canonical basis

err(1)=0; % error with respect to the matrix obtained using MATLAB's "expm" function

for k=1:m
    Qk=Q(:,1:k);
    Hk=H(1:k,1:k);
    if k~=1
        ek(k,1)=1; %vettore e_k della base canonica
        ek(k-1,1)=0;
        e1(k,1)=0;
    end
    nR(k)=beta*HK1K(k)*abs(ek'*expm(Hk)*e1); % Calculation of the exponential 
    % through the funm2 function, which uses the Schur-Parlette algorithm
    err(k)=norm(expm(A)*b-Qk*expm(Hk)*e1*beta,2);
end 

semilogy(1:m, nR,'b'), hold on
xlabel('k')
semilogy(1:m, err,'r')
legend('Generalized residual', 'Estimated error')