function [XNEW, XALT, ALPHA, CAUCHY, GLAG, HCOL, W] = ...
    altmov (N, NPT, XPT, XOPT, BMAT, ZMAT, SL, SU, KOPT, ...
    KNEW, ADELT, XNEW, XALT, CAUCHY, GLAG, HCOL, W)
  %          IMPLICIT REAL*8 (A-H,O-Z)
  %          DIMENSION XPT(NPT,*),XOPT(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),
  %      1     SU(*),XNEW(*),XALT(*),GLAG(*),HCOL(*),W(*)
  %
  %     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have
  %       the same meanings as the corresponding arguments of BOBYQB.
  %     KOPT is the index of the optimal interpolation point.
  %     KNEW is the index of the interpolation point that is going to be moved.
  %     ADELT is the current trust region bound.
  %     XNEW will be set to a suitable new position for the interpolation point
  %       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region
  %       bounds and it should provide a large denominator in the next call of
  %       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the
  %       straight lines through XOPT and another interpolation point.
  %     XALT also provides a large value of the modulus of the KNEW-th Lagrange
  %       function subject to the constraints that have been mentioned, its main
  %       difference from XNEW being that XALT-XOPT is a constrained version of
  %       the Cauchy step within the trust region. An exception is that XALT is
  %       not calculated if all components of GLAG (see below) are zero.
  %     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
  %     CAUCHY will be set to the square of the KNEW-th Lagrange function at
  %       the step XALT-XOPT from XOPT for the vector XALT that is returned,
  %       except that CAUCHY is set to zero if XALT is not calculated.
  %     GLAG is a working space vector of length N for the gradient of the
  %       KNEW-th Lagrange function at XOPT.
  %     HCOL is a working space vector of length NPT for the second derivative
  %       coefficients of the KNEW-th Lagrange function.
  %     W is a working space vector of length 2N that is going to hold the
  %       constrained Cauchy step from XOPT of the Lagrange function, followed
  %       by the downhill version of XALT when the uphill step is calculated.
  %
  %     Set the first NPT components of W to the leading elements of the
  %     KNEW-th column of the H matrix.
  %
  HALF = 0.5e0;
  ONE = 1.0e0;
  ZERO = 0.0e0;
  CONST = ONE + sqrt(2.0e0);
  for K = 1:NPT
    HCOL(K) = ZERO;
  end
  for J = 1:NPT - N - 1
    TEMP = ZMAT(KNEW, J);
    for K = 1:NPT
      HCOL(K) = HCOL(K) + TEMP * ZMAT(K, J);
    end
  end
  ALPHA = HCOL(KNEW);
  HA = HALF * ALPHA;
  %
  %     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
  %
  for I = 1:N
    GLAG(I) = BMAT(KNEW, I);
  end
  for K = 1:NPT
    TEMP = ZERO;
    for J = 1:N
      TEMP = TEMP + XPT(K, J) * XOPT(J);
    end
    TEMP = HCOL(K) * TEMP;
    for I = 1:N
      GLAG(I) = GLAG(I) + TEMP * XPT(K, I);
    end
  end
  %
  %     Search for a large denominator along the straight lines through XOPT
  %     and another interpolation point. SLBD and SUBD will be lower and upper
  %     bounds on the step along each of these lines in turn. PREDSQ will be
  %     set to the square of the predicted denominator for each line. PRESAV
  %     will be set to the largest admissible value of PREDSQ that occurs.
  %
  PRESAV = ZERO;
  for K = 1:NPT
    if (K == KOPT)
      continue
    end
    DDERIV = ZERO;
    DISTSQ = ZERO;
    for I = 1:N
      TEMP = XPT(K, I) - XOPT(I);
      DDERIV = DDERIV + GLAG(I) * TEMP;
      DISTSQ = DISTSQ + TEMP * TEMP;
    end
    SUBD = ADELT / sqrt(DISTSQ);
    SLBD = -SUBD;
    ILBD = 0;
    IUBD = 0;
    SUMIN = min(ONE, SUBD);
    %
    %     Revise SLBD and SUBD if necessary because of the bounds in SL and SU.
    %
    for I = 1:N
      TEMP = XPT(K, I) - XOPT(I);
      if (TEMP > ZERO)
        if (SLBD * TEMP < SL(I) - XOPT(I))
          SLBD = (SL(I) - XOPT(I)) / TEMP;
          ILBD = -I;
        end
        if (SUBD * TEMP > SU(I) - XOPT(I))
          SUBD = max(SUMIN, (SU(I) - XOPT(I)) / TEMP);
          IUBD = I;
        end
      elseif (TEMP < ZERO)
        if (SLBD * TEMP > SU(I) - XOPT(I))
          SLBD = (SU(I) - XOPT(I)) / TEMP;
          ILBD = I;
        end
        if (SUBD * TEMP < SL(I) - XOPT(I))
          SUBD = max(SUMIN, (SL(I) - XOPT(I)) / TEMP);
          IUBD = -I;
        end
      end
    end
    %
    %     Seek a large modulus of the KNEW-th Lagrange function when the index
    %     of the other interpolation point on the line through XOPT is KNEW.
    %
    if (K == KNEW)
      DIFF = DDERIV - ONE;
      STEP = SLBD;
      VLAG = SLBD * (DDERIV - SLBD * DIFF);
      ISBD = ILBD;
      TEMP = SUBD * (DDERIV - SUBD * DIFF);
      if (abs(TEMP) > abs(VLAG))
        STEP = SUBD;
        VLAG = TEMP;
        ISBD = IUBD;
      end
      TEMPD = HALF * DDERIV;
      TEMPA = TEMPD - DIFF * SLBD;
      TEMPB = TEMPD - DIFF * SUBD;
      if (TEMPA * TEMPB < ZERO)
        TEMP = TEMPD * TEMPD / DIFF;
        if (abs(TEMP) > abs(VLAG))
          STEP = TEMPD / DIFF;
          VLAG = TEMP;
          ISBD = 0;
        end
      end
      %
      %     Search along each of the other lines through XOPT and another point.
      %
    else
      STEP = SLBD;
      VLAG = SLBD * (ONE - SLBD);
      ISBD = ILBD;
      TEMP = SUBD * (ONE - SUBD);
      if (abs(TEMP) > abs(VLAG))
        STEP = SUBD;
        VLAG = TEMP;
        ISBD = IUBD;
      end
      if (SUBD > HALF)
        if (abs(VLAG) < 0.25e0)
          STEP = HALF;
          VLAG = 0.25e0;
          ISBD = 0;
        end
      end
      VLAG = VLAG * DDERIV;
    end
    %
    %     Calculate PREDSQ for the current line search and maintain PRESAV.
    %
    TEMP = STEP * (ONE - STEP) * DISTSQ;
    PREDSQ = VLAG * VLAG * (VLAG * VLAG + HA * TEMP * TEMP);
    if (PREDSQ > PRESAV)
      PRESAV = PREDSQ;
      KSAV = K;
      STPSAV = STEP;
      IBDSAV = ISBD;
    end
  end
  %
  %     Construct XNEW in a way that satisfies the bound constraints exactly.
  %
  for I = 1:N
    TEMP = XOPT(I) + STPSAV * (XPT(KSAV, I) - XOPT(I));
    XNEW(I) = max(SL(I), min(SU(I), TEMP));
  end
  if (IBDSAV < 0)
    XNEW(-IBDSAV) = SL(-IBDSAV);
  end
  if (IBDSAV > 0)
    XNEW(IBDSAV) = SU(IBDSAV);
  end
  %
  %     Prepare for the iterative method that assembles the constrained Cauchy
  %     step in W. The sum of squares of the fixed components of W is formed in
  %     WFIXSQ, and the free components of W are set to BIGSTP.
  %
  BIGSTP = ADELT + ADELT;
  IFLAG = 0;
  while (true)
    WFIXSQ = ZERO;
    GGFREE = ZERO;
    for I = 1:N
      W(I) = ZERO;
      TEMPA = min(XOPT(I) - SL(I), GLAG(I));
      TEMPB = max(XOPT(I) - SU(I), GLAG(I));
      if (TEMPA > ZERO || TEMPB < ZERO)
        W(I) = BIGSTP;
        GGFREE = GGFREE + GLAG(I) ^ 2;
      end
    end
    if (GGFREE == ZERO)
      CAUCHY = ZERO;
      return
    end
    %
    %     Investigate whether more components of W can be fixed.
    %
    while (true)
      TEMP = ADELT * ADELT - WFIXSQ;
      if (TEMP > ZERO)
        WSQSAV = WFIXSQ;
        STEP = sqrt(TEMP / GGFREE);
        GGFREE = ZERO;
        for I = 1:N
          if (W(I) == BIGSTP)
            TEMP = XOPT(I) - STEP * GLAG(I);
            if (TEMP <= SL(I))
              W(I) = SL(I) - XOPT(I);
              WFIXSQ = WFIXSQ + W(I) ^ 2;
            elseif (TEMP >= SU(I))
              W(I) = SU(I) - XOPT(I);
              WFIXSQ = WFIXSQ + W(I) ^ 2;
            else
              GGFREE = GGFREE + GLAG(I) ^ 2;
            end
          end
        end
        if (WFIXSQ <= WSQSAV || GGFREE <= ZERO)
          break
        end
      else
        break
      end
    end
    %
    %     Set the remaining free components of W and all components of XALT,
    %     except that W may be scaled later.
    %
    GW = ZERO;
    for I = 1:N
      if (W(I) == BIGSTP)
        W(I) = -STEP * GLAG(I);
        XALT(I) = max(SL(I), min(SU(I), XOPT(I) + W(I)));
      elseif (W(I) == ZERO)
        XALT(I) = XOPT(I);
      elseif (GLAG(I) > ZERO)
        XALT(I) = SL(I);
      else
        XALT(I) = SU(I);
      end
      GW = GW + GLAG(I) * W(I);
    end
    %
    %     Set CURV to the curvature of the KNEW-th Lagrange function along W.
    %     Scale W by a factor less than one if that can reduce the modulus of
    %     the Lagrange function at XOPT+W. Set CAUCHY to the final value of
    %     the square of this function.
    %
    CURV = ZERO;
    for K = 1:NPT
      TEMP = ZERO;
      for J = 1:N
        TEMP = TEMP + XPT(K, J) * W(J);
      end
      CURV = CURV + HCOL(K) * TEMP * TEMP;
    end
    if (IFLAG == 1)
      CURV = -CURV;
    end
    if (CURV > -GW && CURV < -CONST * GW)
      SCALE = -GW / CURV;
      for I = 1:N
        TEMP = XOPT(I) + SCALE * W(I);
        XALT(I) = max(SL(I), min(SU(I), TEMP));
      end
      CAUCHY = (HALF * GW * SCALE) ^ 2;
    else
      CAUCHY = (GW + HALF * CURV) ^ 2;
    end
    %
    %     If IFLAG is zero, then XALT is calculated as before after reversing
    %     the sign of GLAG. Thus two XALT vectors become available. The one that
    %     is chosen is the one that gives the larger value of CAUCHY.
    %
    if (IFLAG ~= 0)
      break
    end
    for I = 1:N
      GLAG(I) = -GLAG(I);
      W(N + I) = XALT(I);
    end
    CSAVE = CAUCHY;
    IFLAG = 1;
  end
  if (CSAVE > CAUCHY)
    for I = 1:N
      XALT(I) = W(N + I);
    end
    CAUCHY = CSAVE;
  end
end
