function  [X]  = ComplexRBWNNM(D,mu1,mu2)
[U, S, V]   = RBSVD(D);
[n1,n2] = size(S);
S1 = zeros(n1/2,n2/2);
S2 = zeros(n1/2,n2/2);
S1 = S(1:n1/2,1:n2/2);
S2 = S(n1/2+1:n1,n2/2+1:n2);

[SigmaZ1, ~] = ClosedWNNMRB(diag(S1), mu1, eps);
[SigmaZ2, ~] = ClosedWNNMRB(diag(S2), mu2, eps);
Z = zeros(n1,n2);
[mm1,nn1] = size(diag(SigmaZ1));
Z(1:mm1,1:nn1) = diag(SigmaZ1);
[mm2,nn2] = size(diag(SigmaZ2));
Z(n1/2+1:n1/2+mm2,n2/2+1:n2/2+nn2) = diag(SigmaZ2);
X = RBbackSVD(U,Z,V);
