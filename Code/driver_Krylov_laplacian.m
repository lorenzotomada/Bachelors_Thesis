close all
clear all

% Delta=[0,1]
N=20;
h=1/N;
c=1e-1;
uDeltaN0=ones(N-1,1);

ADeltaN = -1*gallery('tridiag', N-1, -1, 2, -1);
T=0.1; % Time instant for evaluation
A=T*c*ADeltaN/(h^2);

b=uDeltaN0;
beta=norm(b,2);
q1=b/beta;
f='exp';
tol=1e-8; % Tolerance chosen for the Arnoldi process, to avoid numerical cancellations


[Q,H,HK1K]=Arnoldi(A,q1,tol);
[m, ~] = size(H); % m is the step of the computation reached

nR(1) = 0; % Initialize the generalized residual to 0
ek(1, 1) = 1; % Define the e_k vector of the canonical basis
e1(1, 1) = 1; % Define the e_1 vector of the canonical basis
err(1) = 0; % Error with respect to the matrix obtained using the MATLAB "expm" function

for k=1:m
    Qk=Q(:,1:k);
    Hk=H(1:k,1:k);
    if k~=1
        ek(k,1)=1;  % e_k vector of the canonical basis
        ek(k-1,1)=0;
        e1(k,1)=0;
    end
    nR(k)=beta*HK1K(k)*abs(ek'*funm2(Hk, f)*e1); % Compute the exponential
    % using the funm2 function, which utilizes the Schur-Parlette algorithm
    err(k)=norm(expm(A)*b-Qk*expm(Hk)*e1*beta,2);
end 

semilogy(1:m, nR,'b'), hold on
xlabel('k')
semilogy(1:m, err,'r')
legend('Generalized residual', 'Estimated error')