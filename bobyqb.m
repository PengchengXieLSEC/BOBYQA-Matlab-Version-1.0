function bobyqb(N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT, ...
    MAXFUN, SL, SU)
  %          IMPLICIT REAL*8 (A-H,O-Z)
  %          DIMENSION X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),
  %      1     XOPT(*),GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),
  %      2     SL(*),SU(*),XNEW(*),XALT(*),D(*),VLAG(*),W(*)
  %
  %     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN
  %       are identical to the corresponding arguments in SUBROUTINE BOBYQA.
  %     XBASE holds a shift of origin that should reduce the contributions
  %       from rounding errors to values of the model and Lagrange functions.
  %     XPT is a two-dimensional array that holds the coordinates of the
  %       interpolation points relative to XBASE.
  %     FVAL holds the values of F at the interpolation points.
  %     XOPT is set to the displacement from XBASE of the trust region centre.
  %     GOPT holds the gradient of the quadratic model at XBASE+XOPT.
  %     HQ holds the explicit second derivatives of the quadratic model.
  %     PQ contains the parameters of the implicit second derivatives of the
  %       quadratic model.
  %     BMAT holds the last N columns of H.
  %     ZMAT holds the factorization of the leading NPT by NPT submatrix of H,
  %       this factorization being ZMAT times ZMAT^T, which provides both the
  %       correct rank and positive semi-definiteness.
  %     NDIM is the first dimension of BMAT and has the value NPT+N.
  %     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively.
  %       All the components of every XOPT are going to satisfy the bounds
  %       SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when
  %       XOPT is on a constraint boundary.
  %     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the
  %       vector of variables for the next call of CALFUN. XNEW also satisfies
  %       the SL and SU constraints in the way that has just been mentioned.
  %     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW
  %       in order to increase the denominator in the updating of UPDATE.
  %     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.
  %     VLAG contains the values of the Lagrange functions at a new point X.
  %       They are part of a product that requires VLAG to be of length NDIM.
  %     W is a one-dimensional array thatNP is used for working space. Its length
  %       must be at least 3*NDIM = 3*(NPT+N).
  %
  %     Set some constants.
  %
  HALF = 0.5e0;
  ONE = 1.0e0;
  TEN = 10.0e0;
  TENTH = 0.1e0;
  TWO = 2.0e0;
  ZERO = 0.0e0;
  CAUCHY = ZERO;
  NP = N + 1;
  NDIM = NPT + N;
  NPTM = NPT - NP;
  NH = (N * NP) / 2;
  %
  % Init worspace
  %
  XBASE = zeros(1, N);
  XPT = zeros(NPT, N);
  FVAL = zeros(1, NPT);
  XOPT = zeros(1, N);
  GOPT = zeros(1, N);
  HQ = zeros(1, (N * NP) / 2);
  PQ = zeros(1, NPT);
  BMAT = zeros(NDIM, N);
  ZMAT = zeros(NPT, (NPT - NP));
  XNEW = zeros(1, N);
  XALT = zeros(1, N);
  D = zeros(1, N);
  VLAG = zeros(1, NDIM);
  W = zeros(1, max(5 * N, 3 * NDIM));
  %
  %     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
  %     BMAT and ZMAT for the first iteration, with the corresponding values of
  %     of NF and KOPT, which are the number of calls of CALFUN so far and the
  %     index of the interpolation point at the trust region centre. Then the
  %     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
  %     less than NPT. GOPT will be updated if KOPT is different from KBASE.
  %
  [X, XBASE, XPT, FVAL, GOPT, HQ, PQ, BMAT, ZMAT, NF, KOPT] ...
  = prelim(N, NPT, X, XL, XU, RHOBEG, IPRINT, MAXFUN, XBASE, XPT, ...
    FVAL, GOPT, HQ, PQ, BMAT, ZMAT, SL, SU);
  XOPTSQ = ZERO;
  for I = 1:N
    XOPT(I) = XPT(KOPT, I);
    XOPTSQ = XOPTSQ + XOPT(I) ^ 2;
  end
  FSAVE = FVAL(1);
  if (NF < NPT)
    if (IPRINT > 0)
      disp(['Return from BOBYQA because CALFUN has been', ...
            ' called MAXFUN times.'])
    end
    if (FVAL(KOPT) <= FSAVE)
      for I = 1:N
        X(I) = min(max(XL(I), XBASE(I) + XOPT(I)), XU(I));
        if (XOPT(I) == SL(I))
          X(I) = XL(I);
        end
        if (XOPT(I) == SU(I))
          X(I) = XU(I);
        end
      end
      F = FVAL(KOPT);
    end
    if (IPRINT >= 1)
      disp(['At the return from BOBYQA. Number of function values =', ...
              num2str(NF)])
      disp(['Least value of F =', num2str(F)])
      disp(['The corresponding X is:[', num2str(X(1:N)), ']'])
    end
    return
  end
  KBASE = 1;
  %
  %     Complete the settings that are required for the iterative procedure.
  %
  RHO = RHOBEG;
  DELTA = RHO;
  NRESC = NF;
  NTRITS = 0;
  DIFFA = ZERO;
  DIFFB = ZERO;
  ITEST = 0;
  NFSAV = NF;
  %
  %     Update GOPT if necessary before the first iteration and after each
  %     call of RESCUE that makes a call of CALFUN.
  %
  flag = 20;
  while (flag)
    switch flag
      case 20
        if (KOPT ~= KBASE)
          IH = 0;
          for J = 1:N
            for I = 1:J
              IH = IH + 1;
              if (I < J)
                GOPT(J) = GOPT(J) + HQ(IH) * XOPT(I);
              end
              GOPT(I) = GOPT(I) + HQ(IH) * XOPT(J);
            end
          end
          if (NF > NPT)
            for K = 1:NPT
              TEMP = ZERO;
              for J = 1:N
                TEMP = TEMP + XPT(K, J) * XOPT(J);
              end
              TEMP = PQ(K) * TEMP;
              for I = 1:N
                GOPT(I) = GOPT(I) + TEMP * XPT(K, I);
              end
            end
          end
        end
        flag = 60;
        %
        %     Generate the next point in the trust region that provides a small value
        %     of the quadratic model subject to the constraints on the variables.
        %     The integer NTRITS is set to the number "trust region" iterations that
        %     have occurred since the last "alternative" iteration. If the length
        %     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
        %     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.
        %
      case 60
        [XNEW, D, W(1:N), W(NP:2 * N), W(NP + N:3 * N), W(NP + 2 * N:4 * N), ...
           W(NP + 3 * N:5 * N), DSQ, CRVMIN] = trsbox(N, NPT, XPT, XOPT, GOPT, HQ, ...
          PQ, SL, SU, DELTA, XNEW, D, W(1:N), W(NP:2 * N), W(NP + N:3 * N), ...
          W(NP + 2 * N:4 * N), W(NP + 3 * N:5 * N));
        DNORM = min(DELTA, sqrt(DSQ));
        if (DNORM < HALF * RHO)
          while (true)
            NTRITS = -1;
            DISTSQ = (TEN * RHO) ^ 2;
            if (NF <= NFSAV + 2)
              flag = 650;
              break
            end
            %
            %     The following choice between labels 650 and 680 depends on whether or
            %     not our work with the current RHO seems to be complete. Either RHO is
            %     decreased or termination occurs if the errors in the quadratic model at
            %     the last three interpolation points compare favourably with predictions
            %     of likely improvements to the model within distance HALF*RHO of XOPT.
            %
            ERRBIG = max([DIFFA, DIFFB, DIFFC]);
            FRHOSQ = 0.125e0 * RHO * RHO;
            if (CRVMIN > ZERO && ERRBIG > FRHOSQ * CRVMIN)
              flag = 650;
              break
            end
            BDTOL = ERRBIG / RHO;
            for J = 1:N
              BDTEST = BDTOL;
              if (XNEW(J) == SL(J))
                BDTEST = W(J);
              end
              if (XNEW(J) == SU(J))
                BDTEST = -W(J);
              end
              if (BDTEST < BDTOL)
                CURV = HQ((J + J * J) / 2);
                for K = 1:NPT
                  CURV = CURV + PQ(K) * XPT(K, J) ^ 2;
                end
                BDTEST = BDTEST + HALF * CURV * RHO;
                if (BDTEST < BDTOL)
                  flag = 99999;
                  break
                end
              end
            end
            if (flag == 99999)
              flag = 650;
            else
              flag = 680;
            end
            break
          end
        else
          NTRITS = NTRITS + 1;
          flag = 90;
        end
        %
        %     Severe cancellation is likely to occur if XOPT is too far from XBASE.
        %     If the following test holds, then XBASE is shifted so that XOPT becomes
        %     zero. The appropriate changes are made to BMAT and to the second
        %     derivatives of the current model, beginning with the changes to BMAT
        %     that do not depend on ZMAT. VLAG is used temporarily for working space.
        %
      case 90
        if (DSQ <= 1.0e-3 * XOPTSQ)
          FRACSQ = 0.25e0 * XOPTSQ;
          SUMPQ = ZERO;
          for K = 1:NPT
            SUMPQ = SUMPQ + PQ(K);
            SUM = -HALF * XOPTSQ;
            for I = 1:N
              SUM = SUM + XPT(K, I) * XOPT(I);
            end
            W(NPT + K) = SUM;
            TEMP = FRACSQ - HALF * SUM;
            for I = 1:N
              W(I) = BMAT(K, I);
              VLAG(I) = SUM * XPT(K, I) + TEMP * XOPT(I);
              IP = NPT + I;
              for J = 1:I
                BMAT(IP, J) = BMAT(IP, J) + W(I) * VLAG(J) + VLAG(I) * W(J);
              end
            end
          end
          %
          %     Then the revisions of BMAT that depend on ZMAT are calculated.
          %
          for JJ = 1:NPTM
            SUMZ = ZERO;
            SUMW = ZERO;
            for K = 1:NPT
              SUMZ = SUMZ + ZMAT(K, JJ);
              VLAG(K) = W(NPT + K) * ZMAT(K, JJ);
              SUMW = SUMW + VLAG(K);
            end
            for J = 1:N
              SUM = (FRACSQ * SUMZ - HALF * SUMW) * XOPT(J);
              for K = 1:NPT
                SUM = SUM + VLAG(K) * XPT(K, J);
              end
              W(J) = SUM;
              for K = 1:NPT
                BMAT(K, J) = BMAT(K, J) + SUM * ZMAT(K, JJ);
              end
            end
            for I = 1:N
              IP = I + NPT;
              TEMP = W(I);
              for J = 1:I
                BMAT(IP, J) = BMAT(IP, J) + TEMP * W(J);
              end
            end
          end
          %
          %     The following instructions complete the shift, including the changes
          %     to the second derivative parameters of the quadratic model.
          %
          IH = 0;
          for J = 1:N
            W(J) = -HALF * SUMPQ * XOPT(J);
            for K = 1:NPT
              W(J) = W(J) + PQ(K) * XPT(K, J);
              XPT(K, J) = XPT(K, J) - XOPT(J);
            end
            for I = 1:J
              IH = IH + 1;
              HQ(IH) = HQ(IH) + W(I) * XOPT(J) + XOPT(I) * W(J);
              BMAT(NPT + I, J) = BMAT(NPT + J, I);
            end
          end

          for I = 1:N
            XBASE(I) = XBASE(I) + XOPT(I);
            XNEW(I) = XNEW(I) - XOPT(I);
            SL(I) = SL(I) - XOPT(I);
            SU(I) = SU(I) - XOPT(I);
            XOPT(I) = ZERO;
          end
          XOPTSQ = ZERO;
        end
        if (NTRITS == 0)
          flag = 210;
        else
          flag = 230;
        end

        %
        %     XBASE is also moved to XOPT by a call of RESCUE. This calculation is
        %     more expensive than the previous shift, because new matrices BMAT and
        %     ZMAT are generated from scratch, which may include the replacement of
        %     interpolation points whose positions seem to be causing near linear
        %     dependence in the interpolation conditions. Therefore RESCUE is called
        %     only if rounding errors have reduced by at least a factor of two the
        %     denominator of the formula for updating the H matrix. It provides a
        %     useful safeguard, but is not invoked in most applications of BOBYQA.
        %
      case 190
        NFSAV = NF;
        KBASE = KOPT;
        [XBASE, XPT, FVAL, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, SL, SU, NF, KOPT, ...
           VLAG, W(1:2 * N), W(NP + N:NDIM + N), W(NDIM + NP:3 * NDIM)] = ...
          rescue(N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT, ...
          GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL, SU, NF, DELTA, KOPT, VLAG, ...
          W(1:2 * N), W(NP + N:NDIM + N), W(NDIM + NP:3 * NDIM));
        %
        %     XOPT is updated now in case the branch below to label 720 is taken.
        %     Any updating of GOPT occurs after the branch below to label 20, which
        %     leads to a trust region iteration as does the branch to label 60.
        %
        XOPTSQ = ZERO;
        if (KOPT ~= KBASE)
          for I = 1:N
            XOPT(I) = XPT(KOPT, I);
            XOPTSQ = XOPTSQ + XOPT(I) ^ 2;
          end
        end
        if (NF < 0)
          NF = MAXFUN;
          if (IPRINT > 0)
            disp(['Return from BOBYQA because CALFUN has been', ...
                  ' called MAXFUN times.'])
          end
          break
        end
        NRESC = NF;
        if (NFSAV < NF)
          NFSAV = NF;
          flag = 20;
        elseif (NTRITS > 0)
          flag = 60;
        else
          flag = 210;
        end
        %
        %     Pick two alternative vectors of variables, relative to XBASE, that
        %     are suitable as new positions of the KNEW-th interpolation point.
        %     Firstly, XNEW is set to the point on a line through XOPT and another
        %     interpolation point that minimizes the predicted value of the next
        %     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL
        %     and SU bounds. Secondly, XALT is set to the best feasible point on
        %     a constrained version of the Cauchy step of the KNEW-th Lagrange
        %     function, the corresponding value of the square of this function
        %     being returned in CAUCHY. The choice between these alternatives is
        %     going to be made when the denominator is calculated.
        %
      case 210
        [XNEW, XALT, ALPHA, CAUCHY, W(1:N), W(NP:NDIM), W(NDIM + 1:NDIM + 2 * N)] ...
          = altmov (N, NPT, XPT, XOPT, BMAT, ZMAT, SL, SU, KOPT, KNEW, ...
          ADELT, XNEW, XALT, CAUCHY, W(1:N), W(NP:NDIM), W(NDIM + 1:NDIM + 2 * N));
        for I = 1:N
          D(I) = XNEW(I) - XOPT(I);
        end
        flag = 230;
        %
        %     Calculate VLAG and BETA for the current choice of D. The scalar
        %     product of D with XPT(K,.) is going to be held in W(NPT+K) for
        %     use when VQUAD is calculated.
        %
      case 230
        for K = 1:NPT
          SUMA = ZERO;
          SUMB = ZERO;
          SUM = ZERO;
          for J = 1:N
            SUMA = SUMA + XPT(K, J) * D(J);
            SUMB = SUMB + XPT(K, J) * XOPT(J);
            SUM = SUM + BMAT(K, J) * D(J);
          end
          W(K) = SUMA * (HALF * SUMA + SUMB);
          VLAG(K) = SUM;
          W(NPT + K) = SUMA;
        end
        BETA = ZERO;
        for JJ = 1:NPTM
          SUM = ZERO;
          for K = 1:NPT
            SUM = SUM + ZMAT(K, JJ) * W(K);
          end
          BETA = BETA - SUM * SUM;
          for K = 1:NPT
            VLAG(K) = VLAG(K) + SUM * ZMAT(K, JJ);
          end
        end
        DSQ = ZERO;
        BSUM = ZERO;
        DX = ZERO;
        for J = 1:N
          DSQ = DSQ + D(J) ^ 2;
          SUM = ZERO;
          for K = 1:NPT
            SUM = SUM + W(K) * BMAT(K, J);
          end
          BSUM = BSUM + SUM * D(J);
          JP = NPT + J;
          for I = 1:N
            SUM = SUM + BMAT(JP, I) * D(I);
          end
          VLAG(JP) = SUM;
          BSUM = BSUM + SUM * D(J);
          DX = DX + D(J) * XOPT(J);
        end
        BETA = DX * DX + DSQ * (XOPTSQ + DX + DX + HALF * DSQ) + BETA - BSUM;
        VLAG(KOPT) = VLAG(KOPT) + ONE;
        %
        %     If NTRITS is zero, the denominator may be increased by replacing
        %     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
        %     rounding errors have damaged the chosen denominator.
        %
        if (NTRITS == 0)
          DENOM = VLAG(KNEW) ^ 2 + ALPHA * BETA;
          if (DENOM < CAUCHY && CAUCHY > ZERO)
            for I = 1:N
              XNEW(I) = XALT(I);
              D(I) = XNEW(I) - XOPT(I);
            end
            CAUCHY = ZERO;
            flag = 230;
          elseif (DENOM <= HALF * VLAG(KNEW) ^ 2)
            if (NF > NRESC)
              flag = 190;
            else
              if (IPRINT > 0)
                disp(['Return from BOBYQA because of much', ...
                      ' cancellation in a denominator.'])
              end
              break
            end
          else
            flag = 360;
          end
          %
          %     Alternatively, if NTRITS is positive, then set KNEW to the index of
          %     the next interpolation point to be deleted to make room for a trust
          %     region step. Again RESCUE may be called if rounding errors have damaged
          %     the chosen denominator, which is the reason for attempting to select
          %     KNEW before calculating the next value of the objective function.
          %
        else
          DELSQ = DELTA * DELTA;
          SCADEN = ZERO;
          BIGLSQ = ZERO;
          KNEW = 0;
          for K = 1:NPT
            if (K == KOPT)
              continue
            end
            HDIAG = ZERO;
            for JJ = 1:NPTM
              HDIAG = HDIAG + ZMAT(K, JJ) ^ 2;
            end
            DEN = BETA * HDIAG + VLAG(K) ^ 2;
            DISTSQ = ZERO;
            for J = 1:N
              DISTSQ = DISTSQ + (XPT(K, J) - XOPT(J)) ^ 2;
            end
            TEMP = max(ONE, (DISTSQ / DELSQ) ^ 2);
            if (TEMP * DEN > SCADEN)
              SCADEN = TEMP * DEN;
              KNEW = K;
              DENOM = DEN;
            end
            BIGLSQ = max(BIGLSQ, TEMP * VLAG(K) ^ 2);
          end
          if (SCADEN <= HALF * BIGLSQ)
            if (NF > NRESC)
              flag = 190;
            else
              if (IPRINT > 0)
                disp('Return from BOBYQA because of much', ...
                ' cancellation in a denominator.')
              end
              break
            end
          else
            flag = 360;
          end
        end
        %
        %     Put the variables for the next calculation of the objective function
        %       in XNEW, with any adjustments for the bounds.
        %
        %
        %     Calculate the value of the objective function at XBASE+XNEW, unless
        %       the limit on the number of calculations of F has been reached.
        %
      case 360
        for I = 1:N
          X(I) = min(max(XL(I), XBASE(I) + XNEW(I)), XU(I));
          if (XNEW(I) == SL(I))
            X(I) = XL(I);
          end
          if (XNEW(I) == SU(I))
            X(I) = XU(I);
          end
        end
        if (NF >= MAXFUN)
          if (IPRINT > 0)
            disp(['Return from BOBYQA because CALFUN has been', ...
                  ' called MAXFUN times.'])
          end
          break
        end
        NF = NF + 1;
        F = calfun(N, X);
        if (IPRINT == 3)
          disp(['Function number', num2str(NF), '    F =', num2str(F), ...
                  '    The corresponding X is:[', num2str(X(1:N)), ']'])
        end
        if (NTRITS == -1)
          FSAVE = F;
          break
        end
        %
        %     Use the quadratic model to predict the change in F due to the step D,
        %       and set DIFF to the error of this prediction.
        %
        FOPT = FVAL(KOPT);
        VQUAD = ZERO;
        IH = 0;
        for J = 1:N
          VQUAD = VQUAD + D(J) * GOPT(J);
          for I = 1:J
            IH = IH + 1;
            TEMP = D(I) * D(J);
            if (I == J)
              TEMP = HALF * TEMP;
            end
            VQUAD = VQUAD + HQ(IH) * TEMP;
          end
        end
        for K = 1:NPT
          VQUAD = VQUAD + HALF * PQ(K) * W(NPT + K) ^ 2;
        end
        DIFF = F - FOPT - VQUAD;
        DIFFC = DIFFB;
        DIFFB = DIFFA;
        DIFFA = abs(DIFF);
        if (DNORM > RHO)
          NFSAV = NF;
        end
        %
        %     Pick the next value of DELTA after a trust region step.
        %
        if (NTRITS > 0)
          if (VQUAD >= ZERO)
            if (IPRINT > 0)
              disp(['Return from BOBYQA because a trust', ...
                    ' region step has failed to reduce Q.'])
            end
            break
          end
          RATIO = (F - FOPT) / VQUAD;
          if (RATIO <= TENTH)
            DELTA = min(HALF * DELTA, DNORM);
          elseif (RATIO <= 0.7e0)
            DELTA = max(HALF * DELTA, DNORM);
          else
            DELTA = max(HALF * DELTA, DNORM + DNORM);
          end
          if (DELTA <= 1.5e0 * RHO)
            DELTA = RHO;
          end
          %
          %     Recalculate KNEW and DENOM if the new F is less than FOPT.
          %
          if (F < FOPT)
            KSAV = KNEW;
            DENSAV = DENOM;
            DELSQ = DELTA * DELTA;
            SCADEN = ZERO;
            BIGLSQ = ZERO;
            KNEW = 0;
            for K = 1:NPT
              HDIAG = ZERO;
              for JJ = 1:NPTM
                HDIAG = HDIAG + ZMAT(K, JJ) ^ 2;
              end
              DEN = BETA * HDIAG + VLAG(K) ^ 2;
              DISTSQ = ZERO;
              for J = 1:N
                DISTSQ = DISTSQ + (XPT(K, J) - XNEW(J)) ^ 2;
              end
              TEMP = max(ONE, (DISTSQ / DELSQ) ^ 2);
              if (TEMP * DEN > SCADEN)
                SCADEN = TEMP * DEN;
                KNEW = K;
                DENOM = DEN;
              end
              BIGLSQ = max(BIGLSQ, TEMP * VLAG(K) ^ 2);
            end
            if (SCADEN <= HALF * BIGLSQ)
              KNEW = KSAV;
              DENOM = DENSAV;
            end
          end
        end
        %
        %     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
        %     moved. Also update the second derivative terms of the model.
        %
        [BMAT, ZMAT, VLAG, W] = update (N, NPT, BMAT, ZMAT, VLAG, BETA, ...
        DENOM, KNEW, W);
        IH = 0;
        PQOLD = PQ(KNEW);
        PQ(KNEW) = ZERO;
        for I = 1:N
          TEMP = PQOLD * XPT(KNEW, I);
          for J = 1:I
            IH = IH + 1;
            HQ(IH) = HQ(IH) + TEMP * XPT(KNEW, J);
          end
        end

        for JJ = 1:NPTM
          TEMP = DIFF * ZMAT(KNEW, JJ);
          for K = 1:NPT
            PQ(K) = PQ(K) + TEMP * ZMAT(K, JJ);
          end
        end
        %
        %     Include the new interpolation point, and make the changes to GOPT at
        %     the old XOPT that are caused by the updating of the quadratic model.
        %
        FVAL(KNEW) = F;
        for I = 1:N
          XPT(KNEW, I) = XNEW(I);
          W(I) = BMAT(KNEW, I);
        end

        for K = 1:NPT
          SUMA = ZERO;
          for JJ = 1:NPTM
            SUMA = SUMA + ZMAT(KNEW, JJ) * ZMAT(K, JJ);
          end
          SUMB = ZERO;
          for J = 1:N
            SUMB = SUMB + XPT(K, J) * XOPT(J);
          end
          TEMP = SUMA * SUMB;
          for I = 1:N
            W(I) = W(I) + TEMP * XPT(K, I);
          end
        end
        for I = 1:N
          GOPT(I) = GOPT(I) + DIFF * W(I);
        end
        %
        %     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
        %
        if (F < FOPT)
          KOPT = KNEW;
          XOPTSQ = ZERO;
          IH = 0;
          for J = 1:N
            XOPT(J) = XNEW(J);
            XOPTSQ = XOPTSQ + XOPT(J) ^ 2;
            for I = 1:J
              IH = IH + 1;
              if (I < J)
                GOPT(J) = GOPT(J) + HQ(IH) * D(I);
              end
              GOPT(I) = GOPT(I) + HQ(IH) * D(J);
            end
          end
          for K = 1:NPT
            TEMP = ZERO;
            for J = 1:N
              TEMP = TEMP + XPT(K, J) * D(J);
            end
            TEMP = PQ(K) * TEMP;
            for I = 1:N
              GOPT(I) = GOPT(I) + TEMP * XPT(K, I);
            end
          end
        end
        %
        %     Calculate the parameters of the least Frobenius norm interpolant to
        %     the current data, the gradient of this interpolant at XOPT being put
        %     into VLAG(NPT+I), I=1,2,...,N.
        %
        if (NTRITS > 0)
          for K = 1:NPT
            VLAG(K) = FVAL(K) - FVAL(KOPT);
            W(K) = ZERO;
          end
          for J = 1:NPTM
            SUM = ZERO;
            for K = 1:NPT
              SUM = SUM + ZMAT(K, J) * VLAG(K);
            end
            for K = 1:NPT
              W(K) = W(K) + SUM * ZMAT(K, J);
            end
          end
          for K = 1:NPT
            SUM = ZERO;
            for J = 1:N
              SUM = SUM + XPT(K, J) * XOPT(J);
            end
            W(K + NPT) = W(K);
            W(K) = SUM * W(K);
          end
          GQSQ = ZERO;
          GISQ = ZERO;
          for I = 1:N
            SUM = ZERO;
            for K = 1:NPT
              SUM = SUM + BMAT(K, I) * VLAG(K) + XPT(K, I) * W(K);
            end
            if (XOPT(I) == SL(I))
              GQSQ = GQSQ + min(ZERO, GOPT(I)) ^ 2;
              GISQ = GISQ + min(ZERO, SUM) ^ 2;
            elseif (XOPT(I) == SU(I))
              GQSQ = GQSQ + max(ZERO, GOPT(I)) ^ 2;
              GISQ = GISQ + max(ZERO, SUM) ^ 2;
            else
              GQSQ = GQSQ + GOPT(I) ^ 2;
              GISQ = GISQ + SUM * SUM;
            end
            VLAG(NPT + I) = SUM;
          end
          %
          %     Test whether to replace the new quadratic model by the least Frobenius
          %     norm interpolant, making the replacement if the test is satisfied.
          %
          ITEST = ITEST + 1;
          if (GQSQ < TEN * GISQ)
            ITEST = 0;
          end
          if (ITEST >= 3)
            for I = 1:fix(max(NPT, NH))
              if (I <= N)
                GOPT(I) = VLAG(NPT + I);
              end
              if (I <= NPT)
                PQ(I) = W(NPT + I);
              end
              if (I <= NH)
                HQ(I) = ZERO;
              end
              ITEST = 0;
            end
          end
        end
        %
        %     If a trust region step has provided a sufficient decrease in F, then
        %     branch for another trust region calculation. The case NTRITS=0 occurs
        %     when the new interpolation point was reached by an alternative step.
        %
        if (NTRITS == 0 || F <= FOPT + TENTH * VQUAD)
          flag = 60;
        else
          %
          %     Alternatively, find out if the interpolation points are close enough
          %       to the best point so far.
          %
          DISTSQ = max((TWO * DELTA) ^ 2, (TEN * RHO) ^ 2);
          flag = 650;
        end
      case 650
        KNEW = 0;
        for K = 1:NPT
          SUM = ZERO;
          for J = 1:N
            SUM = SUM + (XPT(K, J) - XOPT(J)) ^ 2;
          end
          if (SUM > DISTSQ)
            KNEW = K;
            DISTSQ = SUM;
          end
        end
        %
        %     If KNEW is positive, then ALTMOV finds alternative new positions for
        %     the KNEW-th interpolation point within distance ADELT of XOPT. It is
        %     reached via label 90. Otherwise, there is a branch to label 60 for
        %     another trust region iteration, unless the calculations with the
        %     current RHO are complete.
        %
        if (KNEW > 0)
          DIST = sqrt(DISTSQ);
          if (NTRITS == -1)
            DELTA = min(TENTH * DELTA, HALF * DIST);
            if (DELTA <= 1.5e0 * RHO)
              DELTA = RHO;
            end
          end
          NTRITS = 0;
          ADELT = max(min(TENTH * DIST, DELTA), RHO);
          DSQ = ADELT * ADELT;
          flag = 90;
        elseif (NTRITS == -1)
          flag = 680;
        elseif (RATIO > ZERO || max(DELTA, DNORM) > RHO)
          flag = 60;
        else
          flag = 680;
        end
        %
        %     The calculations with the current value of RHO are complete. Pick the
        %       next values of RHO and DELTA.
        %
      case 680
        if (RHO > RHOEND)
          DELTA = HALF * RHO;
          RATIO = RHO / RHOEND;
          if (RATIO <= 16.0e0)
            RHO = RHOEND;
          elseif (RATIO <= 250.0e0)
            RHO = sqrt(RATIO) * RHOEND;
          else
            RHO = TENTH * RHO;
          end
          DELTA = max(DELTA, RHO);
          if (IPRINT >= 2)
            disp(['New RHO =', num2str(RHO), ...
                    ' Number of function values =', num2str(NF)])
            disp(['Least value of F =', num2str(FVAL(KOPT))])
            disp(['The corresponding X is:[', num2str(XBASE(1:N) + XOPT(1:N)), ']'])
          end
          NTRITS = 0;
          NFSAV = NF;
          flag = 60;
          %
          %     Return from the calculation, after another Newton-Raphson step, if
          %       it is too short to have been tried before.
          %
        elseif (NTRITS == -1)
          flag = 360;
        else
          break
        end
    end
  end

  if (FVAL(KOPT) <= FSAVE)
    for I = 1:N
      X(I) = min(max(XL(I), XBASE(I) + XOPT(I)), XU(I));
      if (XOPT(I) == SL(I))
        X(I) = XL(I);
      end
      if (XOPT(I) == SU(I))
        X(I) = XU(I);
      end
    end
    F = FVAL(KOPT);
  end
  if (IPRINT >= 1)
    disp(['At the return from BOBYQA. Number of function values =', ...
            num2str(NF)])
    disp(['Least value of F =', num2str(F)])
    disp(['The corresponding X is:[', num2str(X(1:N)), ']'])
  end
end
