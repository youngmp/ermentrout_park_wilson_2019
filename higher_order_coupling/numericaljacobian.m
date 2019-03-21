function[J] = numericaljacobian(myfn,point)

n = length(point);
J = zeros(n,n);
eps = 1e-6;

PM = zeros(size(J));
PP = PM;

for k = 1:n
    epsvec = zeros(size(point));
    epsvec(k) = eps;
    PP(:,k) = myfn(0,point + epsvec);
    PM(:,k) = myfn(0,point - epsvec);
end

for k = 1:n
    for j = 1:n
        J(k,j) = (PP(k,j)-PM(k,j))/(2*eps);
    end
end

