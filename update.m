function [BMAT, ZMAT, VLAG, W] ...
    = update(N, NPT, BMAT, ZMAT, VLAG, BETA, DENOM, KNEW, W)
  %          IMPLICIT REAL*8 (A-H,O-Z)
  %          DIMENSION BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*)
  %
  %     The arrays BMAT and ZMAT are updated, as required by the new position
  %     of the interpolation point that has the index KNEW. The vector VLAG has
  %     N+NPT components, set on entry to the first NPT and last N components
  %     of the product Hw in equation (4.11) of the Powell (2006) paper on
  %     NEWUOA. Further, BETA is set on entry to the value of the parameter
  %     with that name, and DENOM is set to the denominator of the updating
  %     formula. Elements of ZMAT may be treated as zero if their moduli are
  %     at most ZTEST. The first NDIM elements of W are used for working space.
  %
  %     Set some constants.
  %
  ONE = 1.0e0;
  ZERO = 0.0e0;
  NPTM = NPT - N - 1;
  ZTEST = ZERO;
  for K = 1:NPT
    for J = 1:NPTM
      ZTEST = max(ZTEST, abs(ZMAT(K, J)));
    end
  end
  ZTEST = 1.0e-20 * ZTEST;
  %
  %     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
  %
  for J = 2:NPTM
    if (abs(ZMAT(KNEW, J)) > ZTEST)
      TEMP = sqrt(ZMAT(KNEW, 1) ^ 2 + ZMAT(KNEW, J) ^ 2);
      TEMPA = ZMAT(KNEW, 1) / TEMP;
      TEMPB = ZMAT(KNEW, J) / TEMP;
      for I = 1:NPT
        TEMP = TEMPA * ZMAT(I, 1) + TEMPB * ZMAT(I, J);
        ZMAT(I, J) = TEMPA * ZMAT(I, J) - TEMPB * ZMAT(I, 1);
        ZMAT(I, 1) = TEMP;
      end
    end
    ZMAT(KNEW, J) = ZERO;
  end
  %
  %     Put the first NPT components of the KNEW-th column of HLAG into W,
  %     and calculate the parameters of the updating formula.
  %
  for I = 1:NPT
    W(I) = ZMAT(KNEW, 1) * ZMAT(I, 1);
  end
  ALPHA = W(KNEW);
  TAU = VLAG(KNEW);
  VLAG(KNEW) = VLAG(KNEW) - ONE;
  %
  %     Complete the updating of ZMAT.
  %
  TEMP = sqrt(DENOM);
  TEMPB = ZMAT(KNEW, 1) / TEMP;
  TEMPA = TAU / TEMP;
  for I = 1:NPT
    ZMAT(I, 1) = TEMPA * ZMAT(I, 1) - TEMPB * VLAG(I);
  end
  %
  %     Finally, update the matrix BMAT.
  %
  for J = 1:N
    JP = NPT + J;
    W(JP) = BMAT(KNEW, J);
    TEMPA = (ALPHA * VLAG(JP) - TAU * W(JP)) / DENOM;
    TEMPB = (-BETA * W(JP) - TAU * VLAG(JP)) / DENOM;
    for I = 1:JP
      BMAT(I, J) = BMAT(I, J) + TEMPA * VLAG(I) + TEMPB * W(I);
      if (I > NPT)
        BMAT(JP, I - NPT) = BMAT(I, J);
      end
    end
  end
end
