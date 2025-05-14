function  [X]  = ComplexRBNMF(D,mu1,mu2)
alpha = 5;
[U, S, V]   = RBSVD(D);
[n1,n2] = size(S);
S1 = zeros(n1/2,n2/2);
S2 = zeros(n1/2,n2/2);
S1 = S(1:n1/2,1:n2/2);
S2 = S(n1/2+1:n1,n2/2+1:n2);

[SigmaZ1, ~] = ClosedRBNMF(S1,alpha, mu1);
[SigmaZ2, ~] = ClosedRBNMF(S2,alpha, mu2);
Z = zeros(n1,n2);
[mm1,nn1] = size(diag(SigmaZ1));
Z(1:mm1,1:nn1) = diag(SigmaZ1);
[mm2,nn2] = size(diag(SigmaZ2));
Z(n1/2+1:n1/2+mm2,n2/2+1:n2/2+nn2) = diag(SigmaZ2);
X = RBbackSVD(U,Z,V);
