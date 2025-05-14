function X = RBtimes(A,B)

A0 = A.w;
Ai = A.x;
Aj = A.y;
Ak = A.z;

B0 = B.w;
Bi = B.x;
Bj = B.y;
Bk = B.z;

X0 = A0*B0 - Ai*Bi + Aj*Bj - Ak*Bk;
Xi = A0*Bi + Ai*B0 + Aj*Bk + Ak*Bj;
Xj = A0*Bj - Ai*Bk + Aj*B0 - Ak*Bi;
Xk = A0*Bk + Ai*Bj + Aj*Bi + Ak*B0;
X = quaternion(X0,Xi,Xj,Xk);