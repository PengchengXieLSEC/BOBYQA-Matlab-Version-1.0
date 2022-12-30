function [XNEW, D, GNEW, XBDI, S, HS, HRED, DSQ, CRVMIN] ...
    = trsbox (N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL, SU, DELTA, ...
    XNEW, D, GNEW, XBDI, S, HS, HRED)
  %          IMPLICIT REAL*8 (A-H,O-Z)
  %          DIMENSION XPT(NPT,*),XOPT(*),GOPT(*),HQ(*),PQ(*),SL(*),SU(*),
  %      1     XNEW(*),D(*),GNEW(*),XBDI(*),S(*),HS(*),HRED(*)
  %
  %     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same
  %       meanings as the corresponding arguments of BOBYQB.
  %     DELTA is the trust region radius for the present calculation, which
  %       seeks a small value of the quadratic model within distance DELTA of
  %       XOPT subject to the bounds on the variables.
  %     XNEW will be set to a new vector of variables that is approximately
  %       the one that minimizes the quadratic model within the trust region
  %       subject to the SL and SU constraints on the variables. It satisfies
  %       as equations the bounds that become active during the calculation.
  %     D is the calculated trial step from XOPT, generated iteratively from an
  %       initial value of zero. Thus XNEW is XOPT+D after the final iteration.
  %     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated
  %       when D is updated.
  %     XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is
  %       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the
  %       I-th variable has become fixed at a bound, the bound being SL(I) or
  %       SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This
  %       information is accumulated during the construction of XNEW.
  %     The arrays S, HS and HRED are also used for working space. They hold the
  %       current search direction, and the changes in the gradient of Q along S
  %       and the reduced D, respectively, where the reduced D is the same as D,
  %       except that the components of the fixed variables are zero.
  %     DSQ will be set to the square of the length of XNEW-XOPT.
  %     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise
  %       it is set to the least curvature of H that occurs in the conjugate
  %       gradient searches that are not restricted by any constraints. The
  %       value CRVMIN=-1.0e0 is set, however, if all of these searches are
  %       constrained.
  %
  %     A version of the truncated conjugate gradient is applied. If a line
  %     search is restricted by a constraint, then the procedure is restarted,
  %     the values of the variables that are at their bounds being fixed. If
  %     the trust region boundary is reached, then further changes may be made
  %     to D, each one being in the two dimensional space that is spanned
  %     by the current D and the gradient of Q at XOPT+D, staying on the trust
  %     region boundary. Termination occurs when the reduction in Q seems to
  %     be close to the greatest reduction that can be achieved.
  %
  %     Set some constants.
  %
  HALF = 0.5e0;
  ONE = 1.0e0;
  ONEMIN = -1.0e0;
  ZERO = 0.0e0;
  %
  %     The sign of GOPT(I) gives the sign of the change to the I-th variable
  %     that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether
  %     or not to fix the I-th variable at one of its bounds initially, with
  %     NACT being set to the number of fixed variables. D and GNEW are also
  %     set for the first iteration. DELSQ is the upper bound on the sum of
  %     squares of the free variables. QRED is the reduction in Q so far.
  %
  ITERC = 0;
  NACT = 0;
  for I = 1:N
    XBDI(I) = ZERO;
    if (XOPT(I) <= SL(I))
      if (GOPT(I) >= ZERO)
        XBDI(I) = ONEMIN;
      end
    elseif (XOPT(I) >= SU(I))
      if (GOPT(I) <= ZERO)
        XBDI(I) = ONE;
      end
    end
    if (XBDI(I) ~= ZERO)
      NACT = NACT + 1;
    end
    D(I) = ZERO;
    GNEW(I) = GOPT(I);
  end
  DELSQ = DELTA * DELTA;
  QRED = ZERO;
  CRVMIN = ONEMIN;
  %
  %     Set the next search direction of the conjugate gradient method. It is
  %     the steepest descent direction initially and when the iterations are
  %     restarted because a variable has just been fixed by a bound, and of
  %     course the components of the fixed variables are zero. ITERMAX is an
  %     upper bound on the indices of the conjugate gradient iterations.
  %
  flag = 20;
  while (true)
    switch flag
      case 20
        BETA = ZERO;
        flag = 30;
      case 30
        STEPSQ = ZERO;
        for I = 1:N
          if (XBDI(I) ~= ZERO)
            S(I) = ZERO;
          elseif (BETA == ZERO)
            S(I) = -GNEW(I);
          else
            S(I) = BETA * S(I) - GNEW(I);
          end
          STEPSQ = STEPSQ + S(I) ^ 2;
        end
        if (STEPSQ == ZERO)
          break
        end
        if (BETA == ZERO)
          GREDSQ = STEPSQ;
          ITERMAX = ITERC + N - NACT;
        end
        if (GREDSQ * DELSQ <= 1.0e-4 * QRED * QRED)
          break
        end
        %
        %     Multiply the search direction by the second derivative matrix of Q and
        %     calculate some scalars for the choice of steplength. Then set BLEN to
        %     the length of the the step to the trust region boundary and STPLEN to
        %     the steplength, ignoring the simple bounds.
        %
        flag = 210;
      case 50
        RESID = DELSQ;
        DS = ZERO;
        SHS = ZERO;
        for I = 1:N
          if (XBDI(I) == ZERO)
            RESID = RESID - D(I) ^ 2;
            DS = DS + S(I) * D(I);
            SHS = SHS + S(I) * HS(I);
          end
        end
        if (RESID <= ZERO)
          flag = 90;
        else
          TEMP = sqrt(STEPSQ * RESID + DS * DS);
          if (DS < ZERO)
            BLEN = (TEMP - DS) / STEPSQ;
          else
            BLEN = RESID / (TEMP + DS);
          end
          STPLEN = BLEN;
          if (SHS > ZERO)
            STPLEN = min(BLEN, GREDSQ / SHS);
          end
          %
          %     Reduce STPLEN if necessary in order to preserve the simple bounds,
          %     letting IACT be the index of the new constrained variable.
          %
          IACT = 0;
          for I = 1:N
            if (S(I) ~= ZERO)
              XSUM = XOPT(I) + D(I);
              if (S(I) > ZERO)
                TEMP = (SU(I) - XSUM) / S(I);
              else
                TEMP = (SL(I) - XSUM) / S(I);
              end
              if (TEMP < STPLEN)
                STPLEN = TEMP;
                IACT = I;
              end
            end
          end
          %
          %     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.
          %
          SDEC = ZERO;
          if (STPLEN > ZERO)
            ITERC = ITERC + 1;
            TEMP = SHS / STEPSQ;
            if (IACT == 0 && TEMP > ZERO)
              CRVMIN = min(CRVMIN, TEMP);
              if (CRVMIN == ONEMIN)
                CRVMIN = TEMP;
              end
            end
            GGSAV = GREDSQ;
            GREDSQ = ZERO;
            for I = 1:N
              GNEW(I) = GNEW(I) + STPLEN * HS(I);
              if (XBDI(I) == ZERO)
                GREDSQ = GREDSQ + GNEW(I) ^ 2;
              end
              D(I) = D(I) + STPLEN * S(I);
            end
            SDEC = max(STPLEN * (GGSAV - HALF * STPLEN * SHS), ZERO);
            QRED = QRED + SDEC;
          end
          %
          %     Restart the conjugate gradient method if it has hit a new bound.
          %
          if (IACT > 0)
            NACT = NACT + 1;
            XBDI(IACT) = ONE;
            if (S(IACT) < ZERO)
              XBDI(IACT) = ONEMIN;
            end
            DELSQ = DELSQ - D(IACT) ^ 2;
            if (DELSQ <= ZERO)
              flag = 90;
            else
              flag = 20;
            end
          else
            %
            %     If STPLEN is less than BLEN, then either apply another conjugate
            %     gradient iteration or RETURN.
            %
            if (STPLEN < BLEN)
              if (ITERC == ITERMAX || SDEC <= 0.01e0 * QRED)
                break
              end
              BETA = GREDSQ / GGSAV;
              flag = 30;
            else
              flag = 90;
            end
          end
        end
      case 90
        CRVMIN = ZERO;
        flag = 100;
        %
        %     Prepare for the alternative iteration by calculating some scalars and
        %     by multiplying the reduced D by the second derivative matrix of Q.
        %
      case 100
        if (NACT >= N - 1)
          break
        end
        DREDSQ = ZERO;
        DREDG = ZERO;
        GREDSQ = ZERO;
        for I = 1:N
          if (XBDI(I) == ZERO)
            DREDSQ = DREDSQ + D(I) ^ 2;
            DREDG = DREDG + D(I) * GNEW(I);
            GREDSQ = GREDSQ + GNEW(I) ^ 2;
            S(I) = D(I);
          else
            S(I) = ZERO;
          end
        end
        ITCSAV = ITERC;
        flag = 210;
        %
        %     The following instructions multiply the current S-vector by the second
        %     derivative matrix of the quadratic model, putting the product in HS.
        %     They are reached from three different parts of the software above and
        %     they can be regarded as an external subroutine.
        %
      case 210
        IH = 0;
        for J = 1:N
          HS(J) = ZERO;
          for I = 1:J
            IH = IH + 1;
            if (I < J)
              HS(J) = HS(J) + HQ(IH) * S(I);
            end
            HS(I) = HS(I) + HQ(IH) * S(J);
          end
        end
        for K = 1:NPT
          if (PQ(K) ~= ZERO)
            TEMP = ZERO;
            for J = 1:N
              TEMP = TEMP + XPT(K, J) * S(J);
            end
            TEMP = TEMP * PQ(K);
            for I = 1:N
              HS(I) = HS(I) + TEMP * XPT(K, I);
            end
          end
        end
        if (CRVMIN ~= ZERO)
          flag = 50;
        elseif (ITERC > ITCSAV)
          flag = 150;
        else
          for I = 1:N
            HRED(I) = HS(I);
          end
          flag = 120;
        end

        %
        %     Let the search direction S be a linear combination of the reduced D
        %     and the reduced G that is orthogonal to the reduced D.
        %
      case 120
        ITERC = ITERC + 1;
        TEMP = GREDSQ * DREDSQ - DREDG * DREDG;
        if (TEMP <= 1.0e-4 * QRED * QRED)
          break
        end
        TEMP = sqrt(TEMP);
        for I = 1:N
          if (XBDI(I) == ZERO)
            S(I) = (DREDG * D(I) - DREDSQ * GNEW(I)) / TEMP;
          else
            S(I) = ZERO;
          end
        end
        SREDG = -TEMP;
        %
        %     By considering the simple bounds on the variables, calculate an upper
        %     bound on the tangent of half the angle of the alternative iteration,
        %     namely ANGBD, except that, if already a free variable has reached a
        %     bound, there is a branch back to label 100 after fixing that variable.
        %
        ANGBD = ONE;
        IACT = 0;
        for I = 1:N
          if (XBDI(I) == ZERO)
            TEMPA = XOPT(I) + D(I) - SL(I);
            TEMPB = SU(I) - XOPT(I) - D(I);
            if (TEMPA <= ZERO)
              NACT = NACT + 1;
              XBDI(I) = ONEMIN;
              flag = 100;
              break
            elseif (TEMPB <= ZERO)
              NACT = NACT + 1;
              XBDI(I) = ONE;
              flag = 100;
              break
            end
            SSQ = D(I) ^ 2 + S(I) ^ 2;
            TEMP = SSQ - (XOPT(I) - SL(I)) ^ 2;
            if (TEMP > ZERO)
              TEMP = sqrt(TEMP) - S(I);
              if (ANGBD * TEMP > TEMPA)
                ANGBD = TEMPA / TEMP;
                IACT = I;
                XSAV = ONEMIN;
              end
            end
            TEMP = SSQ - (SU(I) - XOPT(I)) ^ 2;
            if (TEMP > ZERO)
              TEMP = sqrt(TEMP) + S(I);
              if (ANGBD * TEMP > TEMPB)
                ANGBD = TEMPB / TEMP;
                IACT = I;
                XSAV = ONE;
              end
            end
          end
          flag = 210;
        end
        %
        %     Calculate HHD and some curvatures for the alternative iteration.
        %
      case 150
        SHS = ZERO;
        DHS = ZERO;
        DHD = ZERO;
        for I = 1:N
          if (XBDI(I) == ZERO)
            SHS = SHS + S(I) * HS(I);
            DHS = DHS + D(I) * HS(I);
            DHD = DHD + D(I) * HRED(I);
          end
        end
        %
        %     Seek the greatest reduction in Q for a range of equally spaced values
        %     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
        %     the alternative iteration.
        %
        REDMAX = ZERO;
        ISAV = 0;
        REDSAV = ZERO;
        IU = 17.0e0 * ANGBD +3.1e0;
        for I = 1:IU
          ANGT = ANGBD * I / IU;
          STH = (ANGT + ANGT) / (ONE + ANGT * ANGT);
          TEMP = SHS + ANGT * (ANGT * DHD - DHS - DHS);
          REDNEW = STH * (ANGT * DREDG - SREDG - HALF * STH * TEMP);
          if (REDNEW > REDMAX)
            REDMAX = REDNEW;
            ISAV = I;
            RDPREV = REDSAV;
          elseif (I == ISAV + 1)
            RDNEXT = REDNEW;
          end
          REDSAV = REDNEW;
        end
        %
        %     Return if the reduction is zero. Otherwise, set the sine and cosine
        %     of the angle of the alternative iteration, and calculate SDEC.
        %
        if (ISAV == 0)
          break
        end
        if (ISAV < IU)
          TEMP = (RDNEXT - RDPREV) / (REDMAX + REDMAX - RDPREV - RDNEXT);
          ANGT = ANGBD * (ISAV + HALF * TEMP) / IU;
        end
        CTH = (ONE - ANGT * ANGT) / (ONE + ANGT * ANGT);
        STH = (ANGT + ANGT) / (ONE + ANGT * ANGT);
        TEMP = SHS + ANGT * (ANGT * DHD - DHS - DHS);
        SDEC = STH * (ANGT * DREDG - SREDG - HALF * STH * TEMP);
        if (SDEC <= ZERO)
          break
        end
        %
        %     Update GNEW, D and HRED. If the angle of the alternative iteration
        %     is restricted by a bound on a free variable, that variable is fixed
        %     at the bound.
        %
        DREDG = ZERO;
        GREDSQ = ZERO;
        for I = 1:N
          GNEW(I) = GNEW(I) + (CTH - ONE) * HRED(I) + STH * HS(I);
          if (XBDI(I) == ZERO)
            D(I) = CTH * D(I) + STH * S(I);
            DREDG = DREDG + D(I) * GNEW(I);
            GREDSQ = GREDSQ + GNEW(I) ^ 2;
          end
          HRED(I) = CTH * HRED(I) + STH * HS(I);
        end
        QRED = QRED + SDEC;
        if (IACT > 0 && ISAV == IU)
          NACT = NACT + 1;
          XBDI(IACT) = XSAV;
          flag = 100;
          %
          %     If SDEC is sufficiently small, then RETURN after setting XNEW to
          %     XOPT+D, giving careful attention to the bounds.
          %
        elseif (SDEC > 0.01e0 * QRED)
          flag = 120;
        else
          break;
        end
    end
  end
  DSQ = ZERO;
  for I = 1:N
    XNEW(I) = max(min(XOPT(I) + D(I), SU(I)), SL(I));
    if (XBDI(I) == ONEMIN)
      XNEW(I) = SL(I);
    end
    if (XBDI(I) == ONE)
      XNEW(I) = SU(I);
    end
    D(I) = XNEW(I) - XOPT(I);
    DSQ = DSQ + D(I) ^ 2;
  end
end
