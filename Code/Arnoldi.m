function [Q, H, HK1K] = Arnoldi(A, b, tol)

% Function that implements the Arnoldi algorithm.
% Takes as input a matrix A, a vector b, with norm(b,2)=1;
% returns as output the matrices Q_k, H_k, as previously defined,
% and the matrix HK1K, whose elements are the coefficients h_{k+1,k} obtained at the k-th iteration.

n=length(b); % dimension of the problem
Q(:,1)=b; % q1=b/norm(b,2); here norm(b,2)=1

for k=1:n
    qk=Q(:,k);
    z=A*qk;

    for i=1:k
        qi=Q(:,i);
        H(i,k)=dot(qi,z);
        z=z-H(i,k)*qi;
    end

    HK1K(k)=norm(z,2);

    if HK1K(k)<=tol % if h_{k+1,k}=0
        return
    end

    H(k+1,k)=HK1K(k);
    Q(:,k+1)=z/HK1K(k);
end
end