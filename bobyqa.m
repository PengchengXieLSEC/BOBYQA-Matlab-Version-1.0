function bobyqa (N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT, MAXFUN)
  %          IMPLICIT REAL*8 (A-H,O-Z)
  %          DIMENSION X(*),XL(*),XU(*),W(*)
  %
  %     This subroutine seeks the least value of a function of many variables,
  %     by applying a trust region method that forms quadratic models by
  %     interpolation. There is usually some freedom in the interpolation
  %     conditions, which is taken up by minimizing the Frobenius norm of
  %     the change to the second derivative of the model, beginning with the
  %     zero matrix. The values of the variables are constrained by upper and
  %     lower bounds. The arguments of the subroutine are as follows.
  %
  %     N must be set to the number of variables and must be at least two.
  %     NPT is the number of interpolation conditions. Its value must be in
  %       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
  %       recommended.
  %     Initial values of the variables must be set in X(1),X(2),...,X(N). They
  %       will be changed to the values that give the least calculated F.
  %     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper
  %       bounds, respectively, on X(I). The construction of quadratic models
  %       requires XL(I) to be strictly less than XU(I) for each I. Further,
  %       the contribution to a model from changes to the I-th variable is
  %       damaged severely by rounding errors if XU(I)-XL(I) is too small.
  %     RHOBEG and RHOEND must be set to the initial and final values of a trust
  %       region radius, so both must be positive with RHOEND no greater than
  %       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
  %       expected change to a variable, while RHOEND should indicate the
  %       accuracy that is required in the final values of the variables. An
  %       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N,
  %       is less than 2*RHOBEG.
  %     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
  %       amount of printing. Specifically, there is no output if IPRINT=0 and
  %       there is output only at the return if IPRINT=1. Otherwise, each new
  %       value of RHO is printed, with the best vector of variables so far and
  %       the corresponding value of the objective function. Further, each new
  %       value of F with its variables are output if IPRINT=3.
  %     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
  %     The array W will be used for working space. Its length must be at least
  %       (NPT+5)*(NPT+N)+3*N*(N+5)/2.
  %
  %     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
  %     F to the value of the objective function for the current values of the
  %     variables X(1),X(2),...,X(N), which are generated automatically in a
  %     way that satisfies the bounds given in XL and XU.
  %
  %     Return if the value of NPT is unacceptable.
  %
  NP = N + 1;
  if (NPT < N + 2 || NPT > ((N + 2) * NP) / 2)
    disp('Return from BOBYQA because NPT is not in  the required interval')
    return
  end
  %
  %     Partition the working space array, so that different parts of it can
  %     be treated separately during the calculation of BOBYQB. The partition
  %     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
  %     space that is taken by the last array in the argument list of BOBYQB.
  %
  SL = zeros(1, N);
  SU = zeros(1, N);
  %
  %     Return if there is insufficient space between the bounds. Modify the
  %     initial X if necessary in order to avoid conflicts between the bounds
  %     and the construction of the first quadratic model. The lower and upper
  %     bounds on moves from the updated X are set now, in the ISL and ISU
  %     partitions of W, in order to provide useful and exact information about
  %     components of X that become within distance RHOBEG from their bounds.
  %
  ZERO = 0.0e0;
  for J = 1:N
    TEMP = XU(J) - XL(J);
    if (TEMP < RHOBEG + RHOBEG)
      disp(['Return from BOBYQA because one of the', ...
            ' differences XU(I)-XL(I) is less than 2*RHOBEG.'])
      return
    end
    SL(J) = XL(J) - X(J);
    SU(J) = XU(J) - X(J);
    if (SL(J) >= -RHOBEG)
      if (SL(J) >= ZERO)
        X(J) = XL(J);
        SL(J) = ZERO;
        SU(J) = TEMP;
      else
        X(J) = XL(J) + RHOBEG;
        SL(J) = -RHOBEG;
        SU(J) = max(XU(J) - X(J), RHOBEG);
      end
    elseif (SU(J) <= RHOBEG)
      if (SU(J) <= ZERO)
        X(J) = XU(J);
        SL(J) = -TEMP;
        SU(J) = ZERO;
      else
        X(J) = XU(J) - RHOBEG;
        SL(J) = min(XL(J) - X(J), -RHOBEG);
        SU(J) = RHOBEG;
      end
    end
  end
  %
  %     Make the call of BOBYQB.
  %
  bobyqb (N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT, MAXFUN, SL, SU)
  return
end
