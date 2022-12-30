function [X, XBASE, XPT, FVAL, GOPT, HQ, PQ, BMAT, ZMAT, NF, KOPT] ...
    = prelim (N, NPT, X, XL, XU, RHOBEG, IPRINT, MAXFUN, XBASE, ...
    XPT, FVAL, GOPT, HQ, PQ, BMAT, ZMAT, SL, SU)
  %          IMPLICIT REAL*8 (A-H,O-Z)
  %          DIMENSION X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),GOPT(*),
  %      1     HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*)
  %
  %     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
  %       same as the corresponding arguments in SUBROUTINE BOBYQA.
  %     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
  %       are the same as the corresponding arguments in BOBYQB, the elements
  %       of SL and SU being set in BOBYQA.
  %     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
  %       it is set by PRELIM to the gradient of the quadratic model at XBASE.
  %       If XOPT is nonzero, BOBYQB will change it to its usual value later.
  %     NF is maintaned as the number of calls of CALFUN so far.
  %     KOPT will be such that the least calculated value of F so far is at
  %       the point XPT(KOPT,.)+XBASE in the space of the variables.
  %
  %     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
  %     BMAT and ZMAT for the first iteration, and it maintains the values of
  %     NF and KOPT. The vector X is also changed by PRELIM.
  %
  %     Set some constants.
  %
  HALF = 0.5e0;
  ONE = 1.0e0;
  TWO = 2.0e0;
  ZERO = 0.0e0;
  RHOSQ = RHOBEG * RHOBEG;
  RECIP = ONE / RHOSQ;
  NP = N + 1;
  %
  %     Set XBASE to the initial vector of variables, and set the initial
  %     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
  %     elements in matlab are all set to zeros() omit init
  %
  for J = 1:N
    XBASE(J) = X(J);
  end
  %
  %     Begin the initialization procedure. NF becomes one more than the number
  %     of function values so far. The coordinates of the displacement of the
  %     next initial interpolation point from XBASE are set in XPT(NF+1,.).
  %
  NF = 0;
  while (NF < NPT && NF < MAXFUN)
    NFM = NF;
    NFX = NF - N;
    NF = NF + 1;
    if (NFM <= 2 * N)
      if (NFM >= 1 && NFM <= N)
        STEPA = RHOBEG;
        if (SU(NFM) == ZERO)
          STEPA = -STEPA;
        end
        XPT(NF, NFM) = STEPA;
      elseif (NFM > N)
        STEPA = XPT(NF - N, NFX);
        STEPB = -RHOBEG;
        if (SL(NFX) == ZERO)
          STEPB = min(TWO * RHOBEG, SU(NFX));
        end
        if (SU(NFX) == ZERO)
          STEPB = max(-TWO * RHOBEG, SL(NFX));
        end
        XPT(NF, NFX) = STEPB;
      end
    else
      ITEMP = (NFM - NP) / N;
      JPT = NFM - ITEMP * N - N;
      IPT = JPT + ITEMP;
      if (IPT > N)
        ITEMP = JPT;
        JPT = IPT - N;
        IPT = ITEMP;
      end
      XPT(NF, IPT) = XPT(IPT + 1, IPT);
      XPT(NF, JPT) = XPT(JPT + 1, JPT);
    end
    %
    %     Calculate the next value of F. The least function value so far and
    %     its index are required.
    %
    for J = 1:N
      X(J) = min(max(XL(J), XBASE(J) + XPT(NF, J)), XU(J));
      if (XPT(NF, J) == SL(J))
        X(J) = XL(J);
      end
      if (XPT(NF, J) == SU(J))
        X(J) = XU(J);
      end
    end
    F = calfun (N, X);
    if (IPRINT == 3)
      disp(['Function number', num2str(NF), '    F =', num2str(F), ...
              ' The corresponding X is:[', num2str(X(1:N)), ']'])
    end
    FVAL(NF) = F;
    if (NF == 1)
      FBEG = F;
      KOPT = 1;
    elseif (F < FVAL(KOPT))
      KOPT = NF;
    end
    %
    %     Set the nonzero initial elements of BMAT and the quadratic model in the
    %     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
    %     of the NF-th and (NF-N)-th interpolation points may be switched, in
    %     order that the function value at the first of them contributes to the
    %     off-diagonal second derivative terms of the initial quadratic model.
    %
    if (NF <= 2 * N + 1)
      if (NF >= 2 && NF <= N + 1)
        GOPT(NFM) = (F - FBEG) / STEPA;
        if (NPT < NF + N)
          BMAT(1, NFM) = -ONE / STEPA;
          BMAT(NF, NFM) = ONE / STEPA;
          BMAT(NPT + NFM, NFM) = -HALF * RHOSQ;
        end
      elseif (NF >= N + 2)
        IH = (NFX * (NFX + 1)) / 2;
        TEMP = (F - FBEG) / STEPB;
        DIFF = STEPB - STEPA;
        HQ(IH) = TWO * (TEMP - GOPT(NFX)) / DIFF;
        GOPT(NFX) = (GOPT(NFX) * STEPB - TEMP * STEPA) / DIFF;
        if (STEPA * STEPB < ZERO)
          if (F < FVAL(NF - N))
            FVAL(NF) = FVAL(NF - N);
            FVAL(NF - N) = F;
            if (KOPT == NF)
              KOPT = NF - N;
            end
            XPT(NF - N, NFX) = STEPB;
            XPT(NF, NFX) = STEPA;
          end
        end
        BMAT(1, NFX) =- (STEPA + STEPB) / (STEPA * STEPB);
        BMAT(NF, NFX) = -HALF / XPT(NF - N, NFX);
        BMAT(NF - N, NFX) = -BMAT(1, NFX) - BMAT(NF, NFX);
        ZMAT(1, NFX) = sqrt(TWO) / (STEPA * STEPB);
        ZMAT(NF, NFX) = sqrt(HALF) / RHOSQ;
        ZMAT(NF - N, NFX) = -ZMAT(1, NFX) - ZMAT(NF, NFX);
      end
      %
      %     Set the off-diagonal second derivatives of the Lagrange functions and
      %     the initial quadratic model.
      %
    else
      IH = (IPT * (IPT - 1)) / 2 + JPT;
      ZMAT(1, NFX) = RECIP;
      ZMAT(NF, NFX) = RECIP;
      ZMAT(IPT + 1, NFX) = -RECIP;
      ZMAT(JPT + 1, NFX) = -RECIP;
      TEMP = XPT(NF, IPT) * XPT(NF, JPT);
      HQ(IH) = (FBEG - FVAL(IPT + 1) - FVAL(JPT + 1) + F) / TEMP;
    end
  end
end
