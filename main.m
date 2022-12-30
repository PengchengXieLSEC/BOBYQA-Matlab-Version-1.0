%
%     Test problem for BOBYQA, the objective function being the sum of
%     the reciprocals of all pairwise distances between the points P_I,
%     I=1,2,...,M in two dimensions, where M=N/2 and where the components
%     of P_I are X(2*I-1) and X(2*I). Thus each vector X of N variables
%     defines the M points P_I. The initial X gives equally spaced points
%     on a circle. Four different choices of the pairs (N,NPT) are tried,
%     namely (10,16), (10,21), (20,26) and (20,41). Convergence to a local
%     minimum that is not global occurs in both the N=10 cases. The details
%     of the results are highly sensitive to computer rounding errors. The
%     choice IPRINT=2 provides the current X and optimal F so far whenever
%     RHO is reduced. The bound constraints of the problem require every
%     component of X to be in the interval [-1,1].
%
% IMPLICIT REAL * 8 (A - H, O - Z)
% DIMENSION X(100), XL(100), XU(100), W(500000)
X = zeros(1, 100);
XL = zeros(1, 100);
XU = zeros(1, 100);
TWOPI = 8.0e0 * atan(1.0e0);
BDL = -1.0e0;
BDU = 1.0e0;
IPRINT = 2;
MAXFUN = 500000;
RHOBEG = 1.0e-1;
RHOEND = 1.0e-6;
M = 10;
while (M <= 10)
  N = 2 * M;
  for I = 1:N
    XL(I) = BDL;
    XU(I) = BDU;
  end
  for JCASE = 2:2
    NPT = N + 6;
    if (JCASE == 2)
      NPT = 2 * N + 1;
    end
    disp(['2D output with M =', num2str(M), ...
            ',  N =', num2str(N), '  and  NPT =', num2str(NPT)])
    for J = 1:M
      TEMP = J * TWOPI / M;
      X(2 * J - 1) = cos(TEMP);
      X(2 * J) = sin(TEMP);
    end
    bobyqa(N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT, MAXFUN)
  end
  M = M + M;
end
