!*************************************************
!This part is the complex LU decomposition
!And the linear matrix solver (complex number) based on LU decomposition
!The original code is from J-P Moreau (see below) which is only for real numbers
!The complex subroutines are based on Evert's code that are modified to be complex*16
!z~ -> complex LU; d~ -> real LU


!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!*                        (www.jpmoreau.fr)            *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     *
!*******************************************************

!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************


!----------------------------------------
!LU decomposition part

    subroutine zLUDCMP(A,N,INDX,D,CODE)

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        COMPLEX(KIND=KIND(1.0D0)), INTENT(INOUT) :: A(N,N)
        INTEGER, INTENT(OUT) :: CODE, D, INDX(N)
        COMPLEX(KIND=KIND(1.0D0)) :: VV(N)     !complex vector (n)
        INTEGER :: I, J, K, IMAX
        DOUBLE PRECISION, PARAMETER :: dTINY=2.2D-16
        COMPLEX(KIND=KIND(1.0D0)) :: AMAX,DUM,zSUM

        D=1
        CODE=0

        DO I=1,N
            AMAX=DCMPLX(0.0d0,0.0d0)
            DO J=1,N
                IF (CDABS(A(I,J)).GT.CDABS(AMAX)) AMAX=A(I,J)
            END DO ! j loop
            IF(CDABS(AMAX).LT.dTINY) THEN
                CODE = 1
                RETURN
            END IF
            VV(I) = 1.0d0 / AMAX
        END DO

        DO J=1,N
            DO I=1,J-1
                zSUM = A(I,J)
                DO K=1,I-1
                    zSUM = zSUM - A(I,K)*A(K,J)
                END DO ! k loop
                A(I,J) = zSUM
            END DO ! i loop
            AMAX = DCMPLX(0.0d0,0.0d0)
            DO I=J,N
                zSUM = A(I,J)
                DO K=1,J-1
                    zSUM = zSUM - A(I,K)*A(K,J)
                END DO ! k loop
                A(I,J) = zSUM
                DUM = VV(I)*CDABS(zSUM)
                IF(CDABS(DUM).GE.CDABS(AMAX)) THEN
                    IMAX = I
                    AMAX = DUM
                END IF
            END DO ! i loop

            IF(J.NE.IMAX) THEN
                DO K=1,N
                    DUM = A(IMAX,K)
                    A(IMAX,K) = A(J,K)
                    A(J,K) = DUM
                END DO ! k loop
                D = -D
                VV(IMAX) = VV(J)
            END IF

            INDX(J) = IMAX
            IF(CDABS(A(J,J)) < dTINY)  A(J,J) = DCMPLX(dTINY,0.0d0)

            IF(J.NE.N) THEN
                DUM = 1.0d0 / A(J,J)
                DO I=J+1,N
                    A(I,J) = A(I,J)*DUM
                END DO ! i loop
            END IF
        END DO ! j loop

        RETURN
    end Subroutine


!********************************************************************
!* Solves the set of N complex linear equations A . X = B.  Here A  *
!* is input, not as the matrix A but rather as its LU decomposition,*
!* determined by the routine CLUDCMP. INDX is input as the permuta- *
!* tion vector returned by CLUDCMP. B is input as the right-hand    *
!* side complex vector B, and returns with the solution vector X. A,*
!* N and INDX are not modified by this routine and can be used for  *
!* successive calls with different right-hand sides. This routine is*
!* also efficient for plain complex matrix inversion.               *
!********************************************************************
    subroutine zLUBKSB(A,N,INDX,B)

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        COMPLEX(KIND=KIND(1.0D0)), INTENT(IN) :: A(N,N)
        COMPLEX(KIND=KIND(1.0D0)), INTENT(OUT) :: B(N)
        INTEGER, INTENT (IN) :: INDX(N)
        COMPLEX(KIND=KIND(1.0D0)) :: zSUM
        INTEGER :: II, I, LL, J

        II = 0

        DO I=1,N
            LL = INDX(I)
            zSUM = B(LL)
            B(LL) = B(I)
            IF(II.NE.0) THEN
                DO J=II,I-1
                    zSUM = zSUM - A(I,J)*B(J)
                END DO ! j loop
            ELSE IF(CDABS(zSUM).NE.0.) THEN
                II = I
            END IF
            B(I) = zSUM
        END DO ! i loop

        DO I=N,1,-1
            zSUM = B(I)
            IF(I < N) THEN
                DO J=I+1,N
                    zSUM = zSUM - A(I,J)*B(J)
                END DO ! j loop
            END IF
            B(I) = zSUM / A(I,I)
        END DO ! i loop

        RETURN
    end Subroutine

    subroutine dLUDCMP(A,N,INDX,D,CODE)

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        DOUBLE PRECISION, INTENT(INOUT) :: A(N,N)
        INTEGER, INTENT(OUT) :: CODE, D, INDX(N)
        DOUBLE PRECISION :: VV(N)     !complex vector (n)
        INTEGER :: I, J, K, IMAX
        DOUBLE PRECISION, PARAMETER :: dTINY=2.2D-16
        DOUBLE PRECISION :: AMAX,DUM,dSUM


        D=1
        CODE=0

        DO I=1,N
            AMAX=0.d0
            DO J=1,N
                IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
            END DO ! j loop
            IF(AMAX.LT.dTINY) THEN
                CODE = 1
                RETURN
            END IF
            VV(I) = 1.d0 / AMAX
        END DO ! i loop

        DO J=1,N
            DO I=1,J-1
                dSUM = A(I,J)
                DO K=1,I-1
                    dSUM = dSUM - A(I,K)*A(K,J)
                END DO ! k loop
                A(I,J) = dSUM
            END DO ! i loop
            AMAX = 0.d0
            DO I=J,N
                dSUM = A(I,J)
                DO K=1,J-1
                    dSUM = dSUM - A(I,K)*A(K,J)
                END DO ! k loop
                A(I,J) = dSUM
                DUM = VV(I)*DABS(dSUM)
                IF(DUM.GE.AMAX) THEN
                    IMAX = I
                    AMAX = DUM
                END IF
            END DO ! i loop

            IF(J.NE.IMAX) THEN
                DO K=1,N
                    DUM = A(IMAX,K)
                    A(IMAX,K) = A(J,K)
                    A(J,K) = DUM
                END DO ! k loop
                D = -D
                VV(IMAX) = VV(J)
            END IF

            INDX(J) = IMAX
            IF(DABS(A(J,J)) < dTINY) A(J,J) = dTINY

            IF(J.NE.N) THEN
                DUM = 1.d0 / A(J,J)
                DO I=J+1,N
                    A(I,J) = A(I,J)*DUM
                END DO ! i loop
            END IF
        END DO ! j loop

        RETURN
    END subroutine


!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
    subroutine dLUBKSB(A,N,INDX,B)

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        DOUBLE PRECISION, INTENT(IN) :: A(N,N)
        DOUBLE PRECISION, INTENT(OUT) :: B(N)
        INTEGER, INTENT (IN) :: INDX(N)
        DOUBLE PRECISION :: dSUM
        INTEGER :: II, I, LL, J
        II = 0

        DO I=1,N
            LL = INDX(I)
            dSUM = B(LL)
            B(LL) = B(I)
            IF(II.NE.0) THEN
                DO J=II,I-1
                    dSUM = dSUM - A(I,J)*B(J)
                END DO ! j loop
            ELSE IF(dSUM.NE.0.d0) THEN
                II = I
            END IF
            B(I) = dSUM
        END DO ! i loop

        DO I=N,1,-1
            dSUM = B(I)
            IF(I < N) THEN
                DO J=I+1,N
                    dSUM = dSUM - A(I,J)*B(J)
                END DO ! j loop
            END IF
            B(I) = dSUM / A(I,I)
        END DO ! i loop

        RETURN
    END subroutine

!LU decompositions end
!*************************************************


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine dbicgstab2 (okprint,l,n,x,nonzero_x,rhs,matvec,precond,toler,typestop, &
                        & mxmatvec,info)
!
! Simple BiCGstab(\ell) iterative method, \ell <= 2
! By M.A.Botchev, Jan.'98
! report bugs to botchev@cwi.nl or botchev@excite.com
!
! Copyright (c) 1998 by M.A.Botchev
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.
!
! This is the "vanilla" version of BiCGstab(\ell) as described
! in PhD thesis of D.R.Fokkema, Chapter 3.  It includes two enhancements
! to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
! 1) G.Sleijpen and H.van der Vorst "Maintaining convergence 
!    properties of BiCGstab methods in finite precision arithmetic",
!    Numerical Algorithms, 10, 1995, pp.203-223
! 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
!    hybrid BiCG methods", Computing, 56, 1996, pp.141-163
!
! {{ This code based on:
! subroutine bistbl v1.0 1995
!
! Copyright (c) 1995 by D.R. Fokkema.
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.  }}
!
! okprint  == (input) LOGICAL. If okprint=.true. residual norm
!            will be printed to *
! l        == (input) INTEGER BiCGstab's dimension <= 2
!            Set l=2 for highly nonsymmetric problems
! n        == (input) INTEGER size of the system to solve
! x        == (input/output) DOUBLE PRECISION array dimension n
!            initial guess on input, solution on output
! rhs      == (input) DOUBLE PRECISION array dimension n
!            right-hand side (rhs) vector
! matvec   == (input) EXTERNAL name of matrix vector subroutine
!            to deliver y:=A*x by CALL matvec(n,x,y)
! nonzero_x== (input) LOGICAL tells
!            BiCGstab(\ell) if initial guess in x is zero or not.
!            If nonzero_x is .FALSE., one MATVEC call is saved.
! toler    == (input/output) DOUBLE PRECISION tolerance for all possible
!            stopping criteria (see the 'typestop' parameter)
!            On output, if info=0 or 1, toler is actually achieved
!            residual reduction or whatever (see the 'typestop' parameter)
! typestop == (input) CHARACTER*3 stopping criterion (||.|| denotes
!            the 2-norm):
!            typestop='rel' -- relative stopping crit.: ||res|| < toler*||res0||
!            typestop='abs' -- absolute stopping crit.: ||res||<toler
!            typestop='max' -- maximum  stopping crit.: max(abs(res))<toler
! NOTE(for typestop='rel' and 'abs'): To save computational work, the value of
!            residual norm used to check the convergence inside the main iterative
!            loop is computed from
!            projections, i.e. it can be smaller than the true residual norm
!            (it may happen when e.g. the 'matrix-free' approach is used).
!            Thus, it is possible that the true residual does NOT satisfy
!            the stopping criterion ('rel' or 'abs').
!            The true residual norm (or residual reduction) is reported on
!            output in parameter TOL -- this can be changed to save 1 MATVEC
!            (see comments at the end of the subroutine)
! mxmatvec==  (input/output) INTEGER.  On input: maximum number of matrix
!            vector multiplications allowed to be done.  On output:
!            actual number of matrix vector multiplications done
! work    ==  (workspace) DOUBLE PRECISION array dimension (n,2*l+5))
! info    ==  (output) INTEGER.  info = 0 in case of normal computations
!            and
!            info = -m (<0) - means paramater number m has an illegal value
!            info = 1 - means no convergence achieved (stopping criterion
!            is not fulfilled)
!            info = 2 - means breakdown of the algorithm (try to enlarge
!            parameter l=\ell to get rid of this)
! ----------------------------------------------------------
        implicit none
        external matvec, precond
        integer l, n, mxmatvec, info
        DOUBLE PRECISION  x(n), rhs(n), toler
        logical   okprint,nonzero_x
        character*3   typestop
!     -----------------------------------------
        DOUBLE PRECISION  work(n,5+2*l)
        integer   lmax
        parameter(lmax=2)
        DOUBLE PRECISION rwork(lmax+1,3+2*(lmax+1))

        logical GoOn, rcmp, xpdt
        integer ii, i1, jj, kk, nmatvec
        DOUBLE PRECISION alpha,beta,hatgamma,kappa0, kappal,maxval1, &
                        &mxnrmr,mxnrmx,omega,rho0,rho1,rnrm0,rnrm,rnrmMax, &
                        &sigma,sum1,varrho
        integer z, zz, y0, yl, y
        integer rr, r, u, xp, bp

        DOUBLE PRECISION    zero, one, delta
        parameter(zero=0d0,one=1d0,delta=1d-2)

        DOUBLE PRECISION    ddot, dnorm2

        logical filexist

        info = 0

        if (l.gt.lmax .or. l.lt.1) info = -2
        if (l.lt.1)        info = -2
        if (toler.le.zero) info = -8
        if (mxmatvec.lt.0) info = -10

        rr = 1
        r = rr+1
        u = r+(l+1)
        xp = u+(l+1)
        bp = xp+1

        z = 1
        zz = z+(l+1)
        y0 = zz+(l+1)
        yl = y0+1
        y = yl+1

        if (info.ne.0) return

!      Initialize first residual
        if (nonzero_x) then
            call matvec (n, x, work(1:n,r) )

!$OMP PARALLEL
!$OMP DO
            do ii=1,n
                work(ii,r) = rhs(ii) - work(ii,r)
            enddo
!$OMP END DO
!$OMP END PARALLEL

            nmatvec = 1
        else

!$OMP PARALLEL
!$OMP DO
            do ii=1,n
                work(ii,r) = rhs(ii)
            enddo
!$OMP END DO
!$OMP END PARALLEL

            nmatvec = 0
        endif
        call precond (n, work(1:n,r))

!     Initialize iteration loop

!$OMP PARALLEL
!$OMP DO
        do ii=1,n
            work(ii,rr) = work(ii,r)
            work(ii,bp) = work(ii,r)
            work(ii,xp) = x(ii)
            x(ii) = zero

        enddo
!$OMP END DO
!$OMP END PARALLEL

        rnrm0 = dnorm2(n, work(:,r))
        rnrm = rnrm0

        mxnrmx = rnrm0
        mxnrmr = rnrm0
        rcmp = .false.
        xpdt = .false.

        alpha = zero
        omega = one
        sigma = one
        rho0 =  one

        if (okprint) then
            inquire (file="Resrcd_dbicg.dat", exist=filexist)
            if (filexist) then
                open (111, file="Resrcd_dbicg.dat", status="OLD", &
                    & position="APPEND", action="WRITE")
            else
                open (111, file="Resrcd_dbicg.dat", status="NEW", &
                    & action="WRITE")
            end if
            write (111, *)
            write (111, *) '========================================='
            write (111, *) 'start bicgstab2!'
            write (111, *) typestop
            write (111, *) 'toler: ', toler
            write (111, *)
            close (111)
        end if

!     Iterate
        if(typestop.eq.'rel')then
            GoOn = rnrm.ge.toler*rnrm0 .and. nmatvec.lt.mxmatvec
            if (okprint) then
                inquire (file="Resrcd_dbicg.dat", exist=filexist)
                if (filexist) then
                    open (111, file="Resrcd_dbicg.dat", status="OLD", &
                        & position="APPEND", action="WRITE")
                else
                    open (111, file="Resrcd_dbicg.dat", status="NEW", &
                        & action="WRITE")
                end if
                write (111, *) nmatvec, rnrm/rnrm0
                close (111)
            end if
        else if(typestop.eq.'abs')then
            GoOn = rnrm.ge.toler       .and. nmatvec.lt.mxmatvec
            if (okprint) then
                inquire (file="Resrcd_dbicg.dat", exist=filexist)
                if (filexist) then
                    open (111, file="Resrcd_dbicg.dat", status="OLD", &
                        & position="APPEND", action="WRITE")
                else
                    open (111, file="Resrcd_dbicg.dat", status="NEW", &
                        & action="WRITE")
                end if
                write (111, *) nmatvec, rnrm
                close (111)
            end if
        else
            info = -9
            return
        end if

        do while (GoOn)
!     =====================
!     --- The BiCG part ---
!     =====================
            rho0 = -omega*rho0
            do kk=1,l

                sum1 = ddot(n, work(:,rr), work(:,r+kk-1))
                rho1 = sum1

                if (rho0.eq.zero) then
                    info = 2
                    return
                endif
                beta = alpha*(rho1/rho0)
                rho0 = rho1
                do jj=0,kk-1

!$OMP PARALLEL
!$OMP DO
                    do ii=1,n
                        work(ii,u+jj) = work(ii,r+jj) - beta*work(ii,u+jj)
                    enddo
!$OMP END DO
!$OMP END PARALLEL

                enddo

                call matvec(n, work(1:n,u+kk-1), work(1:n,u+kk))
                call precond (n, work(1:n,u+kk))
                nmatvec = nmatvec+1

                sum1 = ddot(n, work(:,rr), work(:,u+kk))
                sigma = sum1

                if (sigma.eq.zero) then
                    info = 2
                    return
                endif

                alpha = rho1/sigma

!$OMP PARALLEL
!$OMP DO
                do ii=1,n
                    x(ii) = alpha*work(ii,u) + x(ii)
                enddo
!$OMP END DO
!$OMP END PARALLEL

                do jj=0,kk-1

!$OMP PARALLEL
!$OMP DO
                    do ii=1,n
                        work(ii,r+jj) = -alpha*work(ii,u+jj+1) + work(ii,r+jj)
                    enddo
!$OMP END DO
!$OMP END PARALLEL

                enddo

                call matvec (n, work(1:n,r+kk-1), work(1:n,r+kk))
                call precond (n, work(1:n,r+kk))
                nmatvec = nmatvec+1

                rnrm = dnorm2(n, work(:,r))

                mxnrmx = max (mxnrmx, rnrm)
                mxnrmr = max (mxnrmr, rnrm)
            enddo

!     ==================================
!     --- The convex polynomial part ---
!     ==================================

!        Z = R'R
            do i1=1,l+1
                do jj=i1-1,l
                    sum1 = ddot(n, work(:,r+jj), work(:,r+i1-1))
                    rwork(jj+1,z+i1-1) = sum1
                    rwork(z+i1-1,jj+1) = rwork(jj+1,z+i1-1)
                enddo
            enddo

            do i1=zz,zz+l
                do ii=1,l+1
                    rwork(ii,i1)   = rwork(ii,i1+(z-zz))
                enddo
            enddo
!        tilde r0 and tilde rl (small vectors)

            rwork(1,y0) = -one
            rwork(2,y0) = rwork(2,z) / rwork(2,zz+1)
            rwork(l+1,y0) = zero

            rwork(1,yl) = zero
            rwork(2,yl) = rwork(2,z+l) / rwork(2,zz+1)
            rwork(l+1,yl) = -one

!        Convex combination
            do ii=1,l+1
                rwork(ii,y) = zero
            enddo
            do jj=1,l+1
                do ii=1,l+1
                    rwork(ii,y) = rwork(ii,y) + rwork(jj,yl)*rwork(ii,z+jj-1)
                enddo
            enddo
            sum1 = zero
            do ii=1,l+1
                sum1=sum1+ rwork(ii,yl)*rwork(ii,y)
            enddo
            kappal = dsqrt( sum1 )

            do ii=1,l+1
                rwork(ii,y) = zero
            enddo
            do jj=1,l+1
                do ii=1,l+1
                    rwork(ii,y) = rwork(ii,y) + rwork(jj,y0)*rwork(ii,z+jj-1)
                enddo
            enddo
            sum1 = zero
                do ii=1,l+1
                    sum1=sum1+ rwork(ii,y0)*rwork(ii,y)
                enddo
            kappa0 = dsqrt( sum1 )

            sum1 = zero
            do ii=1,l+1
                sum1=sum1+ rwork(ii,yl)*rwork(ii,y)
            enddo
            varrho = sum1
            varrho = varrho / (kappa0*kappal)

            hatgamma = sign(1d0,varrho)*max(abs(varrho),5d-1)*(kappa0/kappal)

            do ii=1,l+1
                rwork(ii,y0) = -hatgamma*rwork(ii,yl) +      rwork(ii,y0)
            enddo

!        Update
            omega = rwork(l+1,y0)
            do jj=1,l

!$OMP PARALLEL
!$OMP DO
                do ii=1,n
                    work(ii,u) = work(ii,u) - rwork(1+jj,y0)*work(ii,u+jj)
                    x(ii)      = x(ii)      + rwork(1+jj,y0)*work(ii,r+jj-1)
                    work(ii,r) = work(ii,r) - rwork(1+jj,y0)*work(ii,r+jj)
                enddo
!$OMP END DO
!$OMP END PARALLEL

            enddo

            do ii=1,l+1
                rwork(ii,y) = zero
            enddo
            do jj=1,l+1
                do ii=1,l+1
                    rwork(ii,y) = rwork(ii,y) + rwork(jj,y0)*rwork(ii,z+jj-1)
                enddo
            enddo

            sum1 = zero
            do ii=1,l+1
                sum1=sum1+ rwork(ii,y0)*rwork(ii,y)
            enddo
            rnrm = dsqrt( sum1 )

!     ================================
!     --- The reliable update part ---
!     ================================
            mxnrmx = max (mxnrmx, rnrm)
            mxnrmr = max (mxnrmr, rnrm)
            xpdt = (rnrm.lt.delta*rnrm0.and.rnrm0.lt.mxnrmx)
            rcmp = ((rnrm.lt.delta*mxnrmr.and.rnrm0.lt.mxnrmr) .or.xpdt)
            if (rcmp) then
                call matvec (n, x, work(1:n,r) )
                call precond (n, work(1:n,r))
                nmatvec = nmatvec + 1

!$OMP PARALLEL
!$OMP DO
                do ii=1,n
                    work(ii,r) = work(ii,bp) - work(ii,r)
                enddo
!$OMP END DO
!$OMP END PARALLEL

                mxnrmr = rnrm
                if (xpdt) then

!$OMP PARALLEL
!$OMP DO
                    do ii=1,n
                        work(ii,xp) = x(ii) + work(ii,xp)
                        x(ii) = zero
                        work(ii,bp) = work(ii,r)
                    enddo
!$OMP END DO
!$OMP END PARALLEL

                    mxnrmx = rnrm
                endif
            endif

            if(typestop.eq.'rel')then
                GoOn = rnrm.ge.toler*rnrm0 .and. nmatvec.lt.mxmatvec
                if (okprint) then
                    inquire (file="Resrcd_dbicg.dat", exist=filexist)
                    if (filexist) then
                        open (111, file="Resrcd_dbicg.dat", status="OLD", &
                            & position="APPEND", action="WRITE")
                    else
                        open (111, file="Resrcd_dbicg.dat", status="NEW", &
                            & action="WRITE")
                    end if
                    write (111, *) nmatvec, rnrm/rnrm0
                    close (111)
                end if
            else if(typestop.eq.'abs')then
                GoOn = rnrm.ge.toler       .and. nmatvec.lt.mxmatvec
                if (okprint) then
                    inquire (file="Resrcd_dbicg.dat", exist=filexist)
                    if (filexist) then
                        open (111, file="Resrcd_dbicg.dat", status="OLD", &
                            & position="APPEND", action="WRITE")
                    else
                        open (111, file="Resrcd_dbicg.dat", status="NEW", &
                            & action="WRITE")
                    end if
                    write (111, *) nmatvec, rnrm
                    close (111)
                end if
            end if

        enddo

!     =========================
!     --- End of iterations ---
!     =========================

!$OMP PARALLEL
!$OMP DO
        do ii=1,n
            x(ii) = work(ii,xp) + x(ii)
        enddo
!$OMP END DO
!$OMP END PARALLEL

! --------------------- One matvec can be saved by commenting out this:
!                       (this is to compute the true residual)
!       call matvec (n, x, work(1:n,r) )
!       do ii=1,n
!           work(ii,r) = rhs(ii) - work(ii,r)
!       enddo
!       call precond (n, work(1:n,r))
!
!       sum1 = zero
!       do ii=1,n
!           sum1=sum1+ work(ii,r)**2
!       enddo
!       rnrm = sqrt( sum1 )
! --------------------- One matvec can be saved by commenting out this^
!
!       if(typestop.eq.'rel')then
!           if (rnrm.gt.toler*rnrm0) info = 1
!           toler = rnrm/rnrm0
!       else if(typestop.eq.'abs')then
!           if (rnrm.gt.toler) info = 1
!           toler = rnrm
!       end if
!
!       mxmatvec = nmatvec

        return

    end


    subroutine zbicgstab2 (okprint,l,n,x,nonzero_x,rhs,matvec,precond,toler,typestop, &
                         & mxmatvec,info)

! subroutine zbcg2 (okprint,l,n,x,nonzero_x,rhs,matvec,precond,toler, &
!                   mxmatvec,work,info)
!
! Improved "vanilla" BiCGStab(2) iterative method
!
! Copyright (c) 2001 by M.A.Botchev (http://www.math.utwente.nl/~botchev/),
!                       University of Twente
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.
!
! This is the "vanilla" version of BiCGstab(\ell) as described
! in PhD thesis of D.R.Fokkema, Chapter 3 (also available as
! Preprint 976, Dept. of Mathematics, Utrecht University, URL
! http://www.math.uu.nl/publications/).  It includes two enhancements
! to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
! 1) G.Sleijpen and H.van der Vorst "Maintaining convergence 
!    properties of BiCGstab methods in finite precision arithmetic",
!    Numerical Algorithms, 10, 1995, pp.203-223
! 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
!    hybrid BiCG methods", Computing, 56, 1996, pp.141-163
!
! {{ This code based on original work of D.R.Fokkema:
!
! subroutine zbistbl v1.1 1998
! Copyright (c) 1995-1998 by D.R. Fokkema.
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.
!
! }}
!
! Your bug reports, comments, etc. are welcome:
! m.a.botchev@math.utwente.nl
!
! ------------------------------
! Description of the parameters:
! ------------------------------
!
! okprint    (input) LOGICAL. If okprint=.true. the number of
!            matrix-vector multiplications done so far and residual norm will
!            be printed to the standard output each iteration
!
! l          (input) INTEGER the dimension \ell of BiCGstab(\ell)
!            in this simple version it is required that l <= 2
!            l=2 is often useful for systems with nonsymmetric matrices
!
! n          (input) INTEGER size of the linear system to solve
!
! x          (input/output) COMPLEX(KIND=KIND(1.0D0)) array dimension n
!            initial guess on input, solution on output
!
! rhs        (input) COMPLEX(KIND=KIND(1.0D0)) array dimension n
!            the right-hand side (r.h.s.) vector
!
! matvec     (input) EXTERNAL name of matrix vector subroutine
!            to deliver y:=A*x by CALL matvec(n,x,y)
!
! nonzero_x  (input) LOGICAL tells
!            BiCGstab(\ell) if the initial guess x is zero or not.
!            If nonzero_x is .FALSE., initial residual equals r.h.s. vector
!            and one MATVEC call is avoided
!
! toler      (input/output) DOUBLE PRECISION tolerance: the iterations are
!            stopped as soon as || residual ||/|| initial residual|| <= toler,
!            the norm is Euclidean.  On output, if info>=0, the value of
!            toler is set to the actually achieved residual reduction
!
! mxmatvec   (input/output) INTEGER.  On input: maximum number of matrix
!            vector multiplications allowed to be done.  On output:
!            if info>=0, mxmatvec is set to the actual number of matrix
!            vector multiplications done
!
! work       (workspace) COMPLEX(KIND=KIND(1.0D0)) array of dimension (n,2*l+5)
!
! info       (output) INTEGER.  info = 0 in case of succesful computations
!            and
!            info = -m (<0) - means paramater number m has an illegal value
!            info = 1 - means no convergence achieved (stopping criterion
!            is not fulfilled)
!            info = 2 - means breakdown of the algorithm (taking a larger
!            value of parameter l usually helps)
!
! WARNING: If the iterations are ended normally (info=0 or info=1),
! the true residual norm is computed and returned as an output value
! of the parameter toler.  The true residual norm can be slightly larger
! than the projected residual norm used by the algorithm to stop the
! iterations.  It may thus happen that on output info=0 but the value
! of toler is (slightly) larger than tolerance prescribed on input.
! ----------------------------------------------------------
        implicit none

! Parameters:
        logical,    intent(in)   :: okprint,nonzero_x
        integer,    intent(in)   :: l, n
        integer,    intent(in)   :: mxmatvec
        integer,    intent(out)  :: info
        COMPLEX(KIND=KIND(1.0D0)), intent(inout):: x(n)
        COMPLEX(KIND=KIND(1.0D0)), intent(in)   :: rhs(n)
        DOUBLE PRECISION,     intent(in)   :: toler
        character*3,intent(in)   :: typestop
        external                    matvec,precond

! Local variables:
        COMPLEX(KIND=KIND(1.0D0)) :: work(n,3+2*(l+1))
        COMPLEX(KIND=KIND(1.0D0)) :: matrix_z(l+1,l+1),y0(l+1),yl(l+1),zy0(l+1),zyl(l+1)
        logical    :: GoOn, rcmp, xpdt
        integer    :: i, j, k, ii, jj, kk, nmatvec
        COMPLEX(KIND=KIND(1.0D0)) :: alpha, beta, omega, rho0, rho1, sigma
        COMPLEX(KIND=KIND(1.0D0)) :: varrho, hatgamma
        DOUBLE PRECISION     :: rnrm0, rnrm, rnrmMax
        DOUBLE PRECISION     :: mxnrmx, mxnrmr, maxval1
        COMPLEX(KIND=KIND(1.0D0)) :: kappa0, kappal
        logical    :: filexist

! Aliases for the parts of the work array:
        integer              :: rr, r, u, xp, bp

! Constants:
        DOUBLE PRECISION,    parameter :: zero = 0.0d0, one = 1.0d0, delta = 1d-2
        COMPLEX(KIND=KIND(1.0D0)),parameter :: zzero = (0d0,0d0), zone = (1d0,0d0)

! Functions:
        DOUBLE PRECISION     :: zdnorm2
        COMPLEX(KIND=KIND(1.0D0)) :: zdot


        info = 0

        if (l<1 .or. l>2) info = -2
        if (n<1) info = -3
        if (toler<=0d0) info = -9
        if (mxmatvec<0) info = -10

        rr = 1
        r = rr+1
        u = r+(l+1)
        xp = u+(l+1)
        bp = xp+1

        if (info/=0) return

! Initialize first residual

        if (nonzero_x) then

            call matvec (n, x, work(1:n,r))

!$OMP PARALLEL
!$OMP DO
            do ii = 1, n
                work(ii,r) = rhs(ii) - work(ii,r)
            end do
!$OMP END DO
!$OMP END PARALLEL

            nmatvec = 1

        else

!$OMP PARALLEL
!$OMP DO
            do ii = 1, n
                work(ii,r) = rhs(ii)
            end do
!$OMP END DO
!$OMP END PARALLEL

            nmatvec = 0

        end if

        call precond (n, work(1:n,r))

! Initialize iteration loop

!$OMP PARALLEL
!$OMP DO
        do ii = 1, n
            work(ii,rr) = work(ii,r)
            work(ii,bp) = work(ii,r)
            work(ii,xp) = x(ii)
            x(ii) = zzero
        end do
!$OMP END DO
!$OMP END PARALLEL

        rnrm0 = zdnorm2 (n, work(1:n,r))
        rnrm = rnrm0

        mxnrmx = rnrm0
        mxnrmr = rnrm0
        rcmp = .false.
        xpdt = .false.

        alpha = zzero
        omega = zone
        sigma = zone
        rho0  = zone

        if (okprint) then
            inquire (file="Resrcd_dbicg.dat", exist=filexist)
            if (filexist) then
                open (111, file="Resrcd_dbicg.dat", status="OLD", &
                    & position="APPEND", action="WRITE")
            else
                open (111, file="Resrcd_dbicg.dat", status="NEW", &
                    & action="WRITE")
            end if
            write (111, *)
            write (111, *) '========================================='
            write (111, *) 'start bicgstab2!'
            write (111, *) typestop
            write (111, *) 'toler: ', toler
            write (111, *)
            close (111)
        end if

! Iterate

        if(typestop.eq.'rel')then
            GoOn = rnrm.ge.toler*rnrm0 .and. nmatvec.lt.mxmatvec
            if (okprint) then
                inquire (file="Resrcd_dbicg.dat", exist=filexist)
                if (filexist) then
                    open (111, file="Resrcd_dbicg.dat", status="OLD", &
                        & position="APPEND", action="WRITE")
                else
                    open (111, file="Resrcd_dbicg.dat", status="NEW", &
                        & action="WRITE")
                end if
                write (111, *) nmatvec, rnrm/rnrm0
                close (111)
            end if
        else if(typestop.eq.'abs')then
            GoOn = rnrm.ge.toler       .and. nmatvec.lt.mxmatvec
            if (okprint) then
                inquire (file="Resrcd_dbicg.dat", exist=filexist)
                if (filexist) then
                    open (111, file="Resrcd_dbicg.dat", status="OLD", &
                        & position="APPEND", action="WRITE")
                else
                    open (111, file="Resrcd_dbicg.dat", status="NEW", &
                        & action="WRITE")
                end if
                write (111, *) nmatvec, rnrm
                close (111)
            end if
        else
            info = -9
            return
        end if

        do while (GoOn)

! =====================
! The BiCG part ---
! =====================

            rho0 = -omega*rho0
            do k=1,l
                rho1 = zdot (n, work(1:n,rr), work(1:n,r+k-1))
                if (rho0.eq.zzero) then
                    info = 2
!                   toler = rnrm/rnrm0
!                   mxmatvec = nmatvec
                    return
                endif
                beta = alpha*(rho1/rho0)
                rho0 = rho1
                do j=0,k-1

!$OMP PARALLEL
!$OMP DO
                    do ii = 1, n
                        work(ii,u+j) = work(ii,r+j) - beta*work(ii,u+j)
                    end do
!$OMP END DO
!$OMP END PARALLEL

                enddo
                call matvec (n, work(1:n,u+k-1), work(1:n,u+k))
                call precond (n, work(1:n,u+k))
                nmatvec = nmatvec+1

                sigma = zdot (n, work(1:n,rr), work(1:n,u+k))
                if (sigma.eq.zzero) then
                    info = 2
!                   toler = rnrm/rnrm0
!                   mxmatvec = nmatvec
                    return
                endif
                alpha = rho1/sigma

!$OMP PARALLEL
!$OMP DO
                do ii = 1, n
                    x(ii) = alpha*work(ii,u) + x(ii)
                end do
!$OMP END DO
!$OMP END PARALLEL

                do j=0,k-1

!$OMP PARALLEL
!$OMP DO
                    do ii = 1, n
                        work(ii,r+j) = -alpha*work(ii,u+j+1) + work(ii,r+j)
                    end do
!$OMP END DO
!$OMP END PARALLEL

                enddo
                call matvec (n, work(1:n,r+k-1), work(1:n,r+k))
                call precond (n, work(1:n,r+k))
                nmatvec = nmatvec+1
                rnrm = zdnorm2 (n, work(1,r))
                mxnrmx = max (mxnrmx, rnrm)
                mxnrmr = max (mxnrmr, rnrm)
            enddo

! ==================================
! The convex polynomial part ---
! ==================================

!  --- Z = R'R

            do i=1,l+1
                do j=1,i
                    matrix_z(i,j) = dconjg(zdot(n, work(1:n,r+j-1), work(1:n,r+i-1)))
                end do
            end do

!  lower triangular part of Z is computed; compute the rest knowing that Z^H=Z
            do j=2,l+1
                matrix_z(1:j-1,j) = dconjg( matrix_z(j,1:j-1) )
            end do

!  small vectors y0 and yl

            y0(1) = -zone
            y0(2) =      ( matrix_z(2,1) / matrix_z(2,2) )   ! works only for l=2
            y0(l+1) = zzero

            yl(1) = zzero
            yl(2) =      ( matrix_z(2,3) / matrix_z(2,2) )   ! works only for l=2
            yl(l+1) = -zone

!  --- Convex combination

! compute Z*y0 and Z*yl
            zy0 = zzero
            zyl = zzero
            do j=1,l+1
                zy0 = zy0 + matrix_z(:,j)*y0(j)
                zyl = zyl + matrix_z(:,j)*yl(j)
            end do

            kappa0 = dsqrt( abs(zdot(l+1,y0,zy0)) )
            kappal = dsqrt( abs(zdot(l+1,yl,zyl)) )

            varrho = zdot(l+1,yl,zy0) / (kappa0*kappal)

            hatgamma = varrho/abs(varrho) * max( abs(varrho),5d-1 )

            y0 = y0 - (hatgamma*kappa0/kappal)*yl


!  --- Update

            omega = y0(l+1)

            do j=1,l

!$OMP PARALLEL
!$OMP DO
                do ii = 1, n
                    work(ii,u) = work(ii,u) - y0(j+1)*work(ii,u+j)
                    x(ii)      = x(ii)      + y0(j+1)*work(ii,r+j-1)
                    work(ii,r) = work(ii,r) - y0(j+1)*work(ii,r+j)
                end do
!$OMP END DO
!$OMP END PARALLEL

            enddo

! y0 has changed; compute Z*y0 once more
            zy0 = zzero
            do j=1,l+1
                zy0 = zy0 + matrix_z(:,j)*y0(j)
            end do

            rnrm = dsqrt( abs(zdot(l+1,y0,zy0)) )

! ================================
! The reliable update part ---
! ================================

            mxnrmx = max (mxnrmx, rnrm)
            mxnrmr = max (mxnrmr, rnrm)
            xpdt =  (rnrm < delta*rnrm0  .and. rnrm0 < mxnrmx)
            rcmp = ((rnrm < delta*mxnrmr .and. rnrm0 < mxnrmr) .or. xpdt)
            if (rcmp) then
                call matvec (n, x, work(1:n,r))
                call precond (n, work(1:n,r))
                nmatvec = nmatvec + 1

!$OMP PARALLEL
!$OMP DO
                do ii = 1, n
                    work(ii,r) =  work(ii,bp) - work(ii,r)
                end do
!$OMP END DO
!$OMP END PARALLEL

                mxnrmr = rnrm
                if (xpdt) then

!$OMP PARALLEL
!$OMP DO
                    do ii = 1, n
                        work(ii,xp) = x(ii) + work(ii,xp)
                        x(ii) = zzero
                        work(ii,bp) = work(ii,r)
                    end do
!$OMP END DO
!$OMP END PARALLEL

                    mxnrmx = rnrm
                endif
            endif

            if(typestop.eq.'rel')then
                GoOn = rnrm.ge.toler*rnrm0 .and. nmatvec.lt.mxmatvec
                if (okprint) then
                    inquire (file="Resrcd_dbicg.dat", exist=filexist)
                    if (filexist) then
                        open (111, file="Resrcd_dbicg.dat", status="OLD", &
                            & position="APPEND", action="WRITE")
                    else
                        open (111, file="Resrcd_dbicg.dat", status="NEW", &
                            & action="WRITE")
                    end if
                    write (111, *) nmatvec, rnrm/rnrm0
                    close (111)
                end if
            else if(typestop.eq.'abs')then
                GoOn = rnrm.ge.toler       .and. nmatvec.lt.mxmatvec
                if (okprint) then
                    inquire (file="Resrcd_dbicg.dat", exist=filexist)
                    if (filexist) then
                        open (111, file="Resrcd_dbicg.dat", status="OLD", &
                            & position="APPEND", action="WRITE")
                    else
                        open (111, file="Resrcd_dbicg.dat", status="NEW", &
                            & action="WRITE")
                    end if
                    write (111, *) nmatvec, rnrm
                    close (111)
                end if
            end if

        enddo

! =========================
! End of iterations ---
! =========================

!$OMP PARALLEL
!$OMP DO
        do ii = 1, n
            x(ii) = x(ii) + work(ii,xp)
        end do
!$OMP END DO
!$OMP END PARALLEL

! compute the true residual:

! --------------------- One matvec can be saved by commenting out this:
!       call matvec (n, x, work(1:n,r) )
!       work(1:n,r) = rhs(1:n) - work(1:n,r)
!       call precond (n, work(1:n,r))
!       rnrm = zdnorm2(n,work(1:n,r))
!       nmatvec = nmatvec+1
! --------------------- One matvec can be saved by commenting out this^
!
!       if(typestop.eq.'rel')then
!           if (rnrm.gt.toler*rnrm0) info = 1
!           toler = rnrm/rnrm0
!       else if(typestop.eq.'abs')then
!           if (rnrm.gt.toler) info = 1
!           toler = rnrm
!       end if
!
!       mxmatvec = nmatvec

        return

    end subroutine

    DOUBLE PRECISION function ddot(n,zx,zy)

        implicit none
        integer,       intent(in):: n
        DOUBLE PRECISION,intent(in):: zx(n),zy(n)
        DOUBLE PRECISION :: gmre_ctpsum
        integer :: i

        ddot = 0.0d0
!$OMP PARALLEL PRIVATE(gmre_ctpsum)
        gmre_ctpsum = 0.0d0
!$OMP DO
        do i = 1, n
            gmre_ctpsum = gmre_ctpsum + zx(i) * zy(i)
        end do
!$OMP END DO
!$OMP ATOMIC
        ddot = ddot + gmre_ctpsum
!$OMP END PARALLEL

    end function ddot

    DOUBLE PRECISION function dnorm2(n,zx)

        implicit none
        integer :: n
        DOUBLE PRECISION, intent(in) :: zx(n)
        DOUBLE PRECISION, external :: ddot

        dnorm2 = dsqrt(ddot(n,zx,zx))

    end function

    COMPLEX(KIND=KIND(1.0D0)) function zdot(n,zx,zy)

! complex inner product function

        implicit none
        integer,       intent(in):: n
        COMPLEX(KIND=KIND(1.0D0)),intent(in):: zx(n),zy(n)
        COMPLEX(KIND=KIND(1.0D0)) :: gmre_ctpsum
        integer :: i

        !zdot = sum( conjg(zx) * zy )

        zdot = dcmplx(0.0d0,0.0d0)
!$OMP PARALLEL PRIVATE(gmre_ctpsum)
        gmre_ctpsum = dcmplx(0.0d0,0.0d0)
!$OMP DO
        do i = 1, n
            gmre_ctpsum = gmre_ctpsum + dconjg(zx(i)) * zy(i)
        end do
!$OMP END DO
!$OMP ATOMIC
        zdot = zdot+gmre_ctpsum
!$OMP END PARALLEL

    end function zdot


    DOUBLE PRECISION function zdnorm2(n,zx)

! l2 norm function

        implicit none
        integer,       intent(in):: n
        COMPLEX(KIND=KIND(1.0D0)),intent(in):: zx(n)
        COMPLEX(KIND=KIND(1.0D0)),external  :: zdot

        zdnorm2 = dsqrt( abs( zdot(n, zx, zx) ) )

    end function zdnorm2


!   zGETdet and dGETdet use to calculate determinant of complex and real matrix
!   The sample code is from the following link
!   https://dualm.wordpress.com/2012/01/06/computing-determinant-in-fortran/
!   Modified by SQ on 04-Oct-2018

    SUBROUTINE zGETdet (N, MAT, zdet)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: N

        COMPLEX(KIND=KIND(1.0D0)), INTENT(IN) :: MAT(N,N)

        COMPLEX(KIND=KIND(1.0D0)), INTENT(OUT) :: zdet

        INTEGER :: i, info

        INTEGER, ALLOCATABLE :: ipiv(:)

        DOUBLE PRECISION :: sgn

        ALLOCATE(ipiv(N))

        ipiv = 0

        CALL zgetrf(N, N, MAT, N, ipiv, info)

        zdet = DCMPLX(1.0d0, 0.0d0)

        DO i = 1, N

            zdet = zdet*MAT(i, i)

        END DO

        sgn = 1.0d0

        DO i = 1, N

            IF (ipiv(i) /= i) THEN

                sgn = -sgn

            END IF

        END DO

        zdet = sgn*zdet

    END SUBROUTINE


    SUBROUTINE dGETdet (N, MAT, ddet)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: N

        DOUBLE PRECISION, INTENT(IN) :: MAT(N,N)

        DOUBLE PRECISION, INTENT(OUT) :: ddet

        INTEGER :: i, info

        INTEGER, ALLOCATABLE :: ipiv(:)

        DOUBLE PRECISION :: sgn

        ALLOCATE(ipiv(N))

        ipiv = 0

        CALL dgetrf(N, N, MAT, N, ipiv, info)

        ddet = 1.0d0

        DO i = 1, N

            ddet = ddet*MAT(i, i)

        END DO

        sgn = 1.0d0

        DO i = 1, N

            IF (ipiv(i) /= i) THEN

                sgn = -sgn

            END IF

        END DO

        ddet = sgn*ddet

    END SUBROUTINE

