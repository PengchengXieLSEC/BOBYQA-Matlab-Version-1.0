function [XBASE, XPT, FVAL, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, SL, SU, NF, KOPT, ...
            VLAG, PTSAUX, PTSID, W] ...
  = rescue (N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, ...
  FVAL, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL, SU, NF, DELTA, ...
  KOPT, VLAG, PTSAUX, PTSID, W)
%          IMPLICIT REAL*8 (A-H,O-Z)
%          DIMENSION XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),XOPT(*),
%      1     GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*),
%      2     VLAG(*),PTSAUX(2,*),PTSID(*),W(*)
%
%     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
%       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as
%       the corresponding arguments of BOBYQB on the entry to RESCUE.
%     NF is maintained as the number of calls of CALFUN so far, except that
%       NF is set to -1 if the value of MAXFUN prevents further progress.
%     KOPT is maintained so that FVAL(KOPT) is the least calculated function
%       value. Its correct value must be given on entry. It is updated if a
%       new least function value is found, but the corresponding changes to
%       XOPT and GOPT have to be made later by the calling program.
%     DELTA is the current trust region radius.
%     VLAG is a working space vector that will be used for the values of the
%       provisional Lagrange functions at each of the interpolation points.
%       They are part of a product that requires VLAG to be of length NDIM.
%     PTSAUX is also a working space array. For J=1:2,...,N, PTSAUX(1,J) and
%       PTSAUX(2,J) specify the two positions of provisional interpolation
%       points when a nonzero step is taken along e_J (the J-th coordinate
%       direction) through XBASE+XOPT, as specified below. Usually these
%       steps have length DELTA, but other lengths are chosen if necessary
%       in order to satisfy the given bounds on the variables.
%     PTSID is also a working space array. It has NPT components that denote
%       provisional new positions of the original interpolation points, in
%       case changes are needed to restore the linear independence of the
%       interpolation conditions. The K-th point is a candidate for change
%       if and only if PTSID(K) is nonzero. In this case let p and q be the
%       integer parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p
%       and q are both positive, the step from XBASE+XOPT to the new K-th
%       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise
%       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or
%       p=0, respectively.
%     The first NDIM+NPT elements of the array W are used for working space.
%     The final elements of BMAT and ZMAT are set in a well-conditioned way
%       to the values that are appropriate for the new interpolation points.
%     The elements of GOPT, HQ and PQ are also revised to the values that are
%       appropriate to the final quadratic model.
%
%     Set some constants.
%
HALF = 0.5e0;
ONE = 1.0e0;
ZERO = 0.0e0;
NP = N + 1;
SFRAC = HALF / (NP);
NPTM = NPT - NP;
%
%     Shift the interpolation points so that XOPT becomes the origin, and set
%     the elements of ZMAT to zero. The value of SUMPQ is required in the
%     updating of HQ below. The squares of the distances from XOPT to the
%     other interpolation points are set at the end of W. Increments of WINC
%     may be added later to these squares to balance the consideration of
%     the choice of point that is going to become current.
%
SUMPQ = ZERO;
WINC = ZERO;
for K = 1:NPT
  DISTSQ = ZERO;
  for J = 1:N
    XPT(K, J) = XPT(K, J) - XOPT(J);
    DISTSQ = DISTSQ + XPT(K, J) ^ 2;
  end
  SUMPQ = SUMPQ + PQ(K);
  W(NDIM + K) = DISTSQ;
  WINC = max(WINC, DISTSQ);
  for J = 1:NPTM
    ZMAT(K, J) = ZERO;
  end
end
%
%     Update HQ so that HQ and PQ define the second derivatives of the model
%     after XBASE has been shifted to the trust region centre.
%
IH = 0;
for J = 1:N
  W(J) = HALF * SUMPQ * XOPT(J);
  for K = 1:NPT
    W(J) = W(J) + PQ(K) * XPT(K, J);
  end
  for I = 1:J
    IH = IH + 1;
    HQ(IH) = HQ(IH) + W(I) * XOPT(J) + W(J) * XOPT(I);
  end
end
%
%     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
%     also set the elements of PTSAUX.
%
for J = 1:N
  XBASE(J) = XBASE(J) + XOPT(J);
  SL(J) = SL(J) - XOPT(J);
  SU(J) = SU(J) - XOPT(J);
  XOPT(J) = ZERO;
  PTSAUX(1, J) = min(DELTA, SU(J));
  PTSAUX(2, J) = max(-DELTA, SL(J));
  if (PTSAUX(1, J) + PTSAUX(2, J) < ZERO)
    TEMP = PTSAUX(1, J);
    PTSAUX(1, J) = PTSAUX(2, J);
    PTSAUX(2, J) = TEMP;
  end
  if (abs(PTSAUX(2, J)) < HALF * abs(PTSAUX(1, J)))
    PTSAUX(2, J) = HALF * PTSAUX(1, J);
  end
  for I = 1:NDIM
    BMAT(I, J) = ZERO;
  end
end
FBASE = FVAL(KOPT);
%
%     Set the identifiers of the artificial interpolation points that are
%     along a coordinate direction from XOPT, and set the corresponding
%     nonzero elements of BMAT and ZMAT.
%
PTSID(1) = SFRAC;
for J = 1:N
  JP = J + 1;
  JPN = JP + N;
  PTSID(JP) = (J) + SFRAC;
  if (JPN <= NPT)
    PTSID(JPN) = J / NP + SFRAC;
    TEMP = ONE / (PTSAUX(1, J) - PTSAUX(2, J));
    BMAT(JP, J) = -TEMP + ONE / PTSAUX(1, J);
    BMAT(JPN, J) = TEMP + ONE / PTSAUX(2, J);
    BMAT(1, J) = -BMAT(JP, J) - BMAT(JPN, J);
    ZMAT(1, J) = sqrt(2.0e0) / abs(PTSAUX(1, J) * PTSAUX(2, J));
    ZMAT(JP, J) = ZMAT(1, J) * PTSAUX(2, J) * TEMP;
    ZMAT(JPN, J) = -ZMAT(1, J) * PTSAUX(1, J) * TEMP;
  else
    BMAT(1, J) = -ONE / PTSAUX(1, J);
    BMAT(JP, J) = ONE / PTSAUX(1, J);
    BMAT(J + NPT, J) = -HALF * PTSAUX(1, J) ^ 2;
  end
end
%
%     Set any remaining identifiers with their nonzero elements of ZMAT.
%
if (NPT >= N + NP)
  for K = (2 * NP):NPT
    IW = fix(((K - NP) - HALF) / (N));
    IP = K - NP - IW * N;
    IQ = IP + IW;
    if (IQ > N)
      IQ = IQ - N;
    end
    PTSID(K) = (IP) + (IQ) / (NP) + SFRAC;
    TEMP = ONE / (PTSAUX(1, IP) * PTSAUX(1, IQ));
    ZMAT(1, K - NP) = TEMP;
    ZMAT(IP + 1, K - NP) = -TEMP;
    ZMAT(IQ + 1, K - NP) = -TEMP;
    ZMAT(K, K - NP) = TEMP;
  end
end
NREM = NPT;
KOLD = 1;
KNEW = KOPT;
%
%     Reorder the provisional points in the way that exchanges PTSID(KOLD)
%     with PTSID(KNEW).
%
while (true)
  for J = 1:N
    TEMP = BMAT(KOLD, J);
    BMAT(KOLD, J) = BMAT(KNEW, J);
    BMAT(KNEW, J) = TEMP;
  end
  for J = 1:NPTM
    TEMP = ZMAT(KOLD, J);
    ZMAT(KOLD, J) = ZMAT(KNEW, J);
    ZMAT(KNEW, J) = TEMP;
  end
  PTSID(KOLD) = PTSID(KNEW);
  PTSID(KNEW) = ZERO;
  W(NDIM + KNEW) = ZERO;
  NREM = NREM - 1;
  if (KNEW ~= KOPT)
    TEMP = VLAG(KOLD);
    VLAG(KOLD) = VLAG(KNEW);
    VLAG(KNEW) = TEMP;
    %
    %     Update the BMAT and ZMAT matrices so that the status of the KNEW-th
    %     interpolation point can be changed from provisional to original. The
    %     branch to label 350 occurs if all the original points are reinstated.
    %     The nonnegative values of W(NDIM+K) are required in the search below.
    %
    [BMAT, ZMAT, VLAG, W] = update (N, NPT, BMAT, ZMAT, VLAG, BETA, DENOM, ...
    KNEW, W);
    if (NREM == 0)
      return
    end
    for K = 1:NPT
      W(NDIM + K) = abs(W(NDIM + K));
    end
  end
  %
  %     Pick the index KNEW of an original interpolation point that has not
  %     yet replaced one of the provisional interpolation points, giving
  %     attention to the closeness to XOPT and to previous tries with KNEW.
  %
  flag = false;
  while (true)
    DSQMIN = ZERO;
    for K = 1:NPT
      if (W(NDIM + K) > ZERO)
        if (DSQMIN == ZERO || W(NDIM + K) < DSQMIN)
          KNEW = K;
          DSQMIN = W(NDIM + K);
        end
      end
    end
    if (DSQMIN == ZERO)
      flag = true;
      break
    end
    %
    %     Form the W-vector of the chosen original interpolation point.
    %
    for J = 1:N
      W(NPT + J) = XPT(KNEW, J);
    end
    for K = 1:NPT
      SUM = ZERO;
      if (K == KOPT)
        continue
      elseif (PTSID(K) == ZERO)
        for J = 1:N
          SUM = SUM + W(NPT + J) * XPT(K, J);
        end
      else
        IP = fix(PTSID(K));
        if (IP > 0)
          SUM = W(NPT + IP) * PTSAUX(1, IP);
        end
        IQ = fix(NP * PTSID(K) - (IP * NP));
        if (IQ > 0)
          IW = 1;
          if (IP == 0)
            IW = 2;
          end
          SUM = SUM + W(NPT + IQ) * PTSAUX(IW, IQ);
        end
      end
      W(K) = HALF * SUM * SUM;
    end
    %
    %     Calculate VLAG and BETA for the required updating of the H matrix if
    %     XPT(KNEW,.) is reinstated in the set of interpolation points.
    %
    for K = 1:NPT
      SUM = ZERO;
      for J = 1:N
        SUM = SUM + BMAT(K, J) * W(NPT + J);
      end
      VLAG(K) = SUM;
    end
    BETA = ZERO;
    for J = 1:NPTM
      SUM = ZERO;
      for K = 1:NPT
        SUM = SUM + ZMAT(K, J) * W(K);
      end
      BETA = BETA - SUM * SUM;
      for K = 1:NPT
        VLAG(K) = VLAG(K) + SUM * ZMAT(K, J);
      end
    end
    BSUM = ZERO;
    DISTSQ = ZERO;
    for J = 1:N
      SUM = ZERO;
      for K = 1:NPT
        SUM = SUM + BMAT(K, J) * W(K);
      end
      JP = J + NPT;
      BSUM = BSUM + SUM * W(JP);
      for IP = (NPT + 1):NDIM
        SUM = SUM + BMAT(IP, J) * W(IP);
      end
      BSUM = BSUM + SUM * W(JP);
      VLAG(JP) = SUM;
      DISTSQ = DISTSQ + XPT(KNEW, J) ^ 2;
    end
    BETA = HALF * DISTSQ * DISTSQ + BETA - BSUM;
    VLAG(KOPT) = VLAG(KOPT) + ONE;
    %
    %     KOLD is set to the index of the provisional interpolation point that is
    %     going to be deleted to make way for the KNEW-th original interpolation
    %     point. The choice of KOLD is governed by the avoidance of a small value
    %     of the denominator in the updating calculation of UPDATE.
    %
    DENOM = ZERO;
    VLMXSQ = ZERO;
    for K = 1:NPT
      if (PTSID(K) ~= ZERO)
        HDIAG = ZERO;
        for J = 1:NPTM
          HDIAG = HDIAG + ZMAT(K, J) ^ 2;
        end
        DEN = BETA * HDIAG + VLAG(K) ^ 2;
        if (DEN > DENOM)
          KOLD = K;
          DENOM = DEN;
        end
      end
      VLMXSQ = max(VLMXSQ, VLAG(K) ^ 2);
    end
    if (DENOM <= 1.0e-2 * VLMXSQ)
      W(NDIM + KNEW) = -W(NDIM + KNEW) - WINC;
    else
      break
    end
  end
  if flag
    break
  end
end
%
%     When label 260 is reached, all the final positions of the interpolation
%     points have been chosen although any changes have not been included yet
%     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
%     from the shift of XBASE, the updating of the quadratic model remains to
%     be done. The following cycle through the new interpolation points begins
%     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero,
%     except that a return occurs if MAXFUN prohibits another value of F.
%
for KPT = 1:NPT
  if (PTSID(KPT) == ZERO)
    continue
  end
  if (NF >= MAXFUN)
    NF = -1;
    return
  end
  IH = 0;
  for J = 1:N
    W(J) = XPT(KPT, J);
    XPT(KPT, J) = ZERO;
    TEMP = PQ(KPT) * W(J);
    for I = 1:J
      IH = IH + 1;
      HQ(IH) = HQ(IH) + TEMP * W(I);
    end
  end
  PQ(KPT) = ZERO;
  IP = fix(PTSID(KPT));
  IQ = fix(NP * PTSID(KPT) - (IP * NP));
  if (IP > 0)
    XP = PTSAUX(1, IP);
    XPT(KPT, IP) = XP;
  end
  if (IQ > 0)
    XQ = PTSAUX(1, IQ);
    if (IP == 0)
      XQ = PTSAUX(2, IQ);
    end
    XPT(KPT, IQ) = XQ;
  end
  %
  %     Set VQUAD to the value of the current model at the new point.
  %
  VQUAD = FBASE;
  if (IP > 0)
    IHP = (IP + IP * IP) / 2;
    VQUAD = VQUAD + XP * (GOPT(IP) + HALF * XP * HQ(IHP));
  end
  if (IQ > 0)
    IHQ = (IQ + IQ * IQ) / 2;
    VQUAD = VQUAD + XQ * (GOPT(IQ) + HALF * XQ * HQ(IHQ));
    if (IP > 0)
      IW = fix(max(IHP, IHQ)) - fix(abs(IP - IQ));
      VQUAD = VQUAD + XP * XQ * HQ(IW);
    end
  end
  for K = 1:NPT
    TEMP = ZERO;
    if (IP > 0)
      TEMP = TEMP + XP * XPT(K, IP);
    end
    if (IQ > 0)
      TEMP = TEMP + XQ * XPT(K, IQ);
    end
    VQUAD = VQUAD + HALF * PQ(K) * TEMP * TEMP;
  end
  %
  %     Calculate F at the new interpolation point, and set DIFF to the factor
  %     that is going to multiply the KPT-th Lagrange function when the model
  %     is updated to provide interpolation to the new function value.
  %
  for I = 1:N
    W(I) = min(max(XL(I), XBASE(I) + XPT(KPT, I)), XU(I));
    if (XPT(KPT, I) == SL(I))
      W(I) = XL(I);
    end
    if (XPT(KPT, I) == SU(I))
      W(I) = XU(I);
    end
  end
  NF = NF + 1;
  F = calfun(N, W);
  if (IPRINT == 3)
    disp(['Function number', num2str(NF), '    F =', num2str(F), ...
            '    The corresponding X is: [', num2str(X(1:N)), ']'])
  end
  FVAL(KPT) = F;
  if (F < FVAL(KOPT))
    KOPT = KPT;
  end
  DIFF = F - VQUAD;
  %
  %     Update the quadratic model. The return from the subroutine occurs when
  %     all the new interpolation points are included in the model.
  %
  for I = 1:N
    GOPT(I) = GOPT(I) + DIFF * BMAT(KPT, I);
  end
  for K = 1:NPT
    SUM = ZERO;
    for J = 1:NPTM
      SUM = SUM + ZMAT(K, J) * ZMAT(KPT, J);
    end
    TEMP = DIFF * SUM;
    if (PTSID(K) == ZERO)
      PQ(K) = PQ(K) + TEMP;
    else
      IP = fix(PTSID(K));
      IQ = fix((NP) * PTSID(K) - (IP * NP));
      IHQ = (IQ * IQ + IQ) / 2;
      if (IP == 0)
        HQ(IHQ) = HQ(IHQ) + TEMP * PTSAUX(2, IQ) ^ 2;
      else
        IHP = (IP * IP + IP) / 2;
        HQ(IHP) = HQ(IHP) + TEMP * PTSAUX(1, IP) ^ 2;
        if (IQ > 0)
          HQ(IHQ) = HQ(IHQ) + TEMP * PTSAUX(1, IQ) ^ 2;
          IW = fix(max(IHP, IHQ)) - fix(abs(IQ - IP));
          HQ(IW) = HQ(IW) + TEMP * PTSAUX(1, IP) * PTSAUX(1, IQ);
        end
      end
    end
  end
  PTSID(KPT) = ZERO;
end
return
end
