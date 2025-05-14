function [X] = RBWNNM(Y, M)
q2double = @(X) double2q(X,'inverse');

lamb1 = 200;
lamb2 = 150;
mu = 0.001;
rho = 1.5;
[n1, n2] = size(Y);

omega = find(M);  
Yomega = Y(omega);
numIter = 500;
tol = 0.0001;

%% Initialization
X = Y;
Z = X;
eta = quaternion(zeros(n1,n2),zeros(n1,n2),zeros(n1,n2),zeros(n1,n2));


for iter = 1 : numIter
  %% Update X
  tempc = 1 + mu;
  X = Z - eta./mu;
  X(omega) = (mu .* Z(omega) - eta(omega) + Yomega) ./ tempc;

  %% Update Z
 tempA = X + eta./mu;
 tempC1  =  lamb1*10*sqrt(2) * sqrt(min(n1, n2)) / mu;
 tempC2  =  lamb2*10*sqrt(2) * sqrt(min(n1, n2)) / mu;
 Z  = ComplexRBWNNM(tempA,tempC1,tempC2);

 
  %% Update Lagrange multipliers
  eta = eta + mu .* (X - Z);
  mu = mu * rho;
  %% Stop if the termination condition is met
  loss_X = (norm(X - Z, 'fro') / norm(Z, 'fro'))^2;

  if loss_X < tol
    break;
  end    
end
end