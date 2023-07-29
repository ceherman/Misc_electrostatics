
MODULE DH_SurfCal_Solver

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData

    USE BRIEFGHReal

    USE DH_SurfCal_GlobalData

    IMPLICIT NONE

    INTEGER :: rowid,colid,rowscl,colscl

    DOUBLE PRECISION,PRIVATE, ALLOCATABLE, DIMENSION (:,:)   :: &
    &   bimGGEX, bimHHEX, bimGGIN, bimHHIN
    DOUBLE PRECISION,PRIVATE, ALLOCATABLE, DIMENSION (:,:)   :: &
    &   bimBMGGEX,bimBMHHEX,bimBMaxGGEX,bimax0GGEX,bimax0HHEX

    DOUBLE PRECISION,PRIVATE, ALLOCATABLE, DIMENSION (:,:)   :: bimAA,tpPreCD,tpPrercd
    DOUBLE PRECISION,PRIVATE, ALLOCATABLE, DIMENSION (:)     :: bimX,bimBB,tpPreBB,tpPreX
    DOUBLE PRECISION,PRIVATE, ALLOCATABLE, DIMENSION (:,:,:) :: bimPreCD

    CONTAINS

    SUBROUTINE SlvPrblm_DH

        CHARACTER (LEN=1) :: ItrtnorElmntn = 'L'
        INTEGER, ALLOCATABLE, DIMENSION (:) :: slv_ipiv
        INTEGER :: slv_info, slv_D

        INTEGER :: ithprtl, jthprtl, kthprtl, i, j, k, ii, jj, kk, &
        &          icmpt, id_tp, slf_i, slf_j, slf_k, GLQi, icnt

        INTEGER :: NdA, NdB, NdC, NdD, NdE, NdF, ndst, nded, elst, eled

        INTEGER :: Ttlbim, MatSize
        INTEGER, ALLOCATABLE, DIMENSION (:) :: mtidst

        DOUBLE PRECISION :: mdm_oi, mdm_io
        DOUBLE PRECISION :: tp, tp1, tp2, tp3, tp4, tp5, tp6

        CHARACTER (LEN=99) :: filename

        Ttlbim = ttlnmbrnd
        ALLOCATE (mtidst(ttlnmbrnd))

        MatSize = 0
        DO ithprtl = 1, nmbrprtl
            IF (BCType_DH(ithprtl) /= '2SD') THEN
                id_tp = 1
                DO i = ndstaID(ithprtl), ndendID(ithprtl)
                    mtidst(i) = MatSize + id_tp
                    id_tp = id_tp + 1
                END DO
                MatSize = MatSize + nmbrnd(ithprtl)
            END IF
            IF (BCType_DH(ithprtl) == '2SD') THEN
                id_tp = 1
                DO i = ndstaID(ithprtl), ndendID(ithprtl)
                    mtidst(i) = MatSize + id_tp
                    id_tp = id_tp + 2
                END DO
                MatSize = MatSize + 2*nmbrnd(ithprtl)
            END IF
        END DO

        ALLOCATE(bimAA(MatSize,MatSize))
        ALLOCATE(bimBB(MatSize))
        ALLOCATE(bimX(MatSize))

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
        DO i = 1, MatSize
            DO j = 1, MatSize
                bimAA(i,j) = 0.0d0
            END DO
            bimBB(i) = 0.0d0
            bimX(i) = 0.0d0
        END DO
!$OMP END DO
!$OMP END PARALLEL

        CALL Get_BRIEFGGHHEX_DH
        CALL Get_BRIEFGGHHIN_DH


!$OMP PARALLEL PRIVATE (i,j,ii,jj,kk,NdA) &
!$OMP & PRIVATE (ithprtl,jthprtl,id_tp) &
!$OMP & PRIVATE (tp,tp1,tp2,tp3,tp4,tp5,tp6,mdm_oi,mdm_io)
!$OMP DO
        DO i = 1, ttlnmbrnd

            ithprtl = 1
            IF (nmbrprtl > 1 .AND. i > ndendID(1)) THEN
                DO id_tp = 2, nmbrprtl
                    IF (i >= ndstaID(id_tp) .AND. i <= ndendID(id_tp)) THEN
                        ithprtl = id_tp
                        EXIT
                    END IF
                END DO
            END IF

!******
!ithparticle with double-sided Condition:

            IF (BCType_DH(ithprtl) == '2SD') THEN

                IF (srcendID_DH(corelnkshell(ithprtl)) > 0) THEN
                    id_tp = corelnkshell(ithprtl)
                    tp = 0.0d0
                    ii = srcstaID_DH(corelnkshell(ithprtl))
                    jj = srcendID_DH(corelnkshell(ithprtl))
                    DO kk = ii, jj
                        tp1 = xnd(i) - xsrc_DH(kk)
                        tp2 = ynd(i) - ysrc_DH(kk)
                        tp3 = znd(i) - zsrc_DH(kk)
                        tp4 = DSQRT(tp1*tp1+tp2*tp2+tp3*tp3)
                        IF (id_tp == 0) THEN
                            tp5 = DEXP(-exK_DH*tp4)/tp4
                        END IF
                        IF (id_tp > 0) THEN
                            tp5 = DEXP(-inK_DH(id_tp)*tp4)/tp4
                        END IF
                        tp = tp + srcStrength_DH(kk)*tp5
                    END DO
                    IF (id_tp == 0) THEN
                        tp = tp/exeps_DH
                    ELSE
                        tp = tp/ineps_DH(id_tp)
                    END IF
                    IF (NrmlInOut(ithprtl) == 1) THEN
                        bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) + tp
                    END IF
                    IF (NrmlInOut(ithprtl) ==-1) THEN
                        bimBB(mtidst(i)  ) = bimBB(mtidst(i)  ) - tp
                    END IF
                END IF

                IF (srcendID_DH(ithprtl) > 0) THEN
                    tp = 0.0d0
                    ii = srcstaID_DH(ithprtl)
                    jj = srcendID_DH(ithprtl)
                    DO kk = ii, jj
                        tp1 = xnd(i) - xsrc_DH(kk)
                        tp2 = ynd(i) - ysrc_DH(kk)
                        tp3 = znd(i) - zsrc_DH(kk)
                        tp4 = DSQRT(tp1*tp1+tp2*tp2+tp3*tp3)
                        tp5 = DEXP(-inK_DH(ithprtl)*tp4)/tp4
                        tp = tp + srcStrength_DH(kk)*tp5
                    END DO
                    tp = tp/ineps_DH(ithprtl)
                    IF (NrmlInOut(ithprtl) == 1) THEN
                        bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) - tp
                    END IF
                    IF (NrmlInOut(ithprtl) ==-1) THEN
                        bimBB(mtidst(i)+1) = bimBB(mtidst(i)+1) + tp
                    END IF
                END IF

!***
!External field

!---
!DH(exPhi2_DH) = 0
!ROW 1 with GGEX and HHEX
!Part 1: surfaces share the same 'corelnkshell'

                DO jthprtl = 1, nmbrprtl

                    id_tp = corelnkshell(jthprtl)
                    IF (id_tp == 0) THEN
                        mdm_oi = exeps_DH/ineps_DH(jthprtl)
                    ELSE
                        mdm_oi = ineps_DH(id_tp)/ineps_DH(jthprtl)
                    END IF

                    IF (corelnkshell(jthprtl) == corelnkshell(ithprtl)) THEN

                        DO j = ndstaID(jthprtl), ndendID(jthprtl)

                            IF (BCType_DH(jthprtl) == '2SD') THEN

                                bimAA(mtidst(i)  ,mtidst(j)  ) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)  ) + bimHHEX(i,j)

                                bimAA(mtidst(i)  ,mtidst(j)+1) &
                            &=  bimAA(mtidst(i)  ,mtidst(j)+1) - bimGGEX(i,j)

                            END IF

                        END DO

                    END IF

                END DO

!---

!***
!Internal field
!DH(inPhi) = 0 at Eq. (2)
!ROW 2 with GGIN and HHIN

                DO jthprtl = 1, nmbrprtl

                    id_tp = corelnkshell(jthprtl)
                    IF (id_tp == 0) THEN
                        mdm_oi = exeps_DH/ineps_DH(jthprtl)
                    ELSE
                        mdm_oi = ineps_DH(id_tp)/ineps_DH(jthprtl)
                    END IF

!---
!Part1: on the surface where the host node x0 locates
!!NOTE: 'phi' is 'inPhi'

                    IF (jthprtl == ithprtl) THEN

                        DO j = ndstaID(jthprtl), ndendID(jthprtl)

                            bimAA(mtidst(i)+1,mtidst(j)  ) &
                        &=  bimAA(mtidst(i)+1,mtidst(j)  ) + bimHHIN(i,j)

                            bimAA(mtidst(i)+1,mtidst(j)+1) &
                        &=  bimAA(mtidst(i)+1,mtidst(j)+1) - bimGGIN(i,j) * mdm_oi

                            IF (NrmlInOut(jthprtl) == 1) THEN

                                tp1 = exPhi1_DH(j) - inPhi1_DH(j)
                                tp2 =-ndBCKnown_DH(j)/ineps_DH(jthprtl) &
                                &    +(mdm_oi*exPhi1dn_DH(j) - inPhi1dn_DH(j))

                                    bimBB(mtidst(i)+1) &
                                &=  bimBB(mtidst(i)+1) - bimHHIN(i,j) * tp1 &
                                &                      + bimGGIN(i,j) * tp2

                            END IF

                            IF (NrmlInOut(jthprtl) ==-1) THEN

                                tp1 = exPhi1_DH(j) - inPhi1_DH(j)
                                tp2 = ndBCKnown_DH(j)/ineps_DH(jthprtl) &
                                &    +(mdm_oi*exPhi1dn_DH(j) - inPhi1dn_DH(j))

                                    bimBB(mtidst(i)+1) &
                                &=  bimBB(mtidst(i)+1) - bimHHIN(i,j) * tp1 &
                                &                      + bimGGIN(i,j) * tp2

                            END IF

                        END DO

                    END IF
                    
                END DO
                
            END IF

!---
!End Internal field
!***

!end with ithparticle with double-sided Condition
!******

        END DO
!$OMP END DO
!$OMP END PARALLEL

        DEALLOCATE (bimGGEX, bimHHEX, bimGGIN, bimHHIN)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!second argument = 0 -> LU decomposition
        CALL Slv_LinearMatrixREAL_LU(MatSize, 0)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!reload the results

!$OMP PARALLEL PRIVATE (i,NdA) &
!$OMP & PRIVATE (ithprtl,id_tp) &
!$OMP & PRIVATE (tp,tp1,tp2,tp3,tp4,tp5,tp6,mdm_oi,mdm_io)
!$OMP DO
        DO i = 1, ttlnmbrnd

            ithprtl = 1
            IF (nmbrprtl > 1 .AND. i > ndendID(1)) THEN
                DO id_tp = 2, nmbrprtl
                    IF (i >= ndstaID(id_tp) .AND. i <= ndendID(id_tp)) THEN
                        ithprtl = id_tp
                        EXIT
                    END IF
                END DO
            END IF

            id_tp = corelnkshell(ithprtl)
            IF (id_tp == 0) THEN
                mdm_oi = exeps_DH/ineps_DH(ithprtl)
            ELSE
                mdm_oi = ineps_DH(id_tp)/ineps_DH(ithprtl)
            END IF

            IF (BCType_DH(ithprtl) == '2SD') THEN

                exPhi2_DH(i) = bimX(mtidst(i)  )
                exPhi3_DH(i) = exPhi1_DH(i) + exPhi2_DH(i)

                exPhi2dn_DH(i) = bimX(mtidst(i)+1)
                exPhi3dn_DH(i) = exPhi1dn_DH(i) + exPhi2dn_DH(i)

                inPhi3_DH(i) = exPhi3_DH(i)
                inPhi2_DH(i) = inPhi3_DH(i) - inPhi1_DH(i)

                IF (NrmlInOut(ithprtl) == 1) THEN
                    tp2 =-ndBCKnown_DH(i)/ineps_DH(ithprtl) &
                    &    +(mdm_oi*exPhi3dn_DH(i) - inPhi1dn_DH(i))
                    inPhi2dn_DH(i) = tp2
                END IF
                IF (NrmlInOut(ithprtl) ==-1) THEN
                    tp2 = ndBCKnown_DH(i)/ineps_DH(ithprtl) &
                    &    +(mdm_oi*exPhi3dn_DH(i) - inPhi1dn_DH(i))
                    inPhi2dn_DH(i) = tp2
                END IF
                inPhi3dn_DH(i) = inPhi1dn_DH(i) + inPhi2dn_DH(i)

                IF (NrmlInOut(ithprtl) == 1) THEN
                    IF (corelnkshell(ithprtl) == 0) &
                    &   sigma_DH(i) = exPhi3dn_DH(i)*exeps_DH &
                    &                -inPhi3dn_DH(i)*ineps_DH(ithprtl)
                    IF (corelnkshell(ithprtl) /= 0) &
                    &   sigma_DH(i) = exPhi3dn_DH(i)*ineps_DH(corelnkshell(ithprtl)) &
                    &                -inPhi3dn_DH(i)*ineps_DH(ithprtl)
                END IF

                IF (NrmlInOut(ithprtl) ==-1) THEN
                    IF (corelnkshell(ithprtl) == 0) &
                    &   sigma_DH(i) =-exPhi3dn_DH(i)*exeps_DH &
                    &                +inPhi3dn_DH(i)*ineps_DH(ithprtl)
                    IF (corelnkshell(ithprtl) /= 0) &
                    &   sigma_DH(i) =-exPhi3dn_DH(i)*ineps_DH(corelnkshell(ithprtl)) &
                    &                +inPhi3dn_DH(i)*ineps_DH(ithprtl)
                END IF

            END IF

        END DO
!$OMP END DO
!$OMP END PARALLEL

        DEALLOCATE (mtidst, bimX)

        DO ithprtl = 1, nmbrprtl

            surfQ_DH(ithprtl) = 0.0d0

            IF (MeshType == 'L') THEN

                DO k = elstaID(ithprtl), elendID(ithprtl)

                    tp1 = sigma_DH(elmntlnknd(k,1))
                    tp2 = sigma_DH(elmntlnknd(k,2))
                    tp3 = sigma_DH(elmntlnknd(k,3))

                    DO GLQi = 1, n_glqtr2d

                        icnt = n_glqtr2d*(k-1)+GLQi

                        surfQ_DH(ithprtl) = surfQ_DH(ithprtl) &
                        &                  +srcfmm_wtnd(1,icnt) * tp1 &
                        &                  +srcfmm_wtnd(2,icnt) * tp2 &
                        &                  +srcfmm_wtnd(3,icnt) * tp3

                    END DO

                END DO

            END IF

            IF (MeshType == 'Q') THEN

                DO k = elstaID(ithprtl), elendID(ithprtl)

                    tp1 = sigma_DH(elmntlnknd(k,1))
                    tp2 = sigma_DH(elmntlnknd(k,2))
                    tp3 = sigma_DH(elmntlnknd(k,3))
                    tp4 = sigma_DH(elmntlnknd(k,4))
                    tp5 = sigma_DH(elmntlnknd(k,5))
                    tp6 = sigma_DH(elmntlnknd(k,6))

                    DO GLQi = 1, n_glqtr2d

                        icnt = n_glqtr2d*(k-1)+GLQi

                        surfQ_DH(ithprtl) = surfQ_DH(ithprtl) &
                        &                  +srcfmm_wtnd(1,icnt) * tp1 &
                        &                  +srcfmm_wtnd(2,icnt) * tp2 &
                        &                  +srcfmm_wtnd(3,icnt) * tp3 &
                        &                  +srcfmm_wtnd(4,icnt) * tp4 &
                        &                  +srcfmm_wtnd(5,icnt) * tp5 &
                        &                  +srcfmm_wtnd(6,icnt) * tp6

                    END DO


                END DO

            END IF

        END DO

    END SUBROUTINE

    SUBROUTINE Slv_LinearMatrixREAL_LU(TtlbimA,CalMtrxCNDet_On)

        INTEGER, INTENT(IN) :: TtlbimA,CalMtrxCNDet_On

        CHARACTER (LEN=1) :: ItrtnorElmntn = 'I'

        INTEGER :: i, j, k, num_PreCD

        INTEGER, ALLOCATABLE, DIMENSION (:) :: slv_ipiv
        INTEGER :: slv_info, slv_D
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: &
        &   bimMtrxInvCNDet,bimMtrxVL,bimMtrxVR
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: bimMtrxwork, bimMtrxW
        DOUBLE PRECISION :: bimMtrxRCN, bimMtrxanorm
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: bimMtrxrwork

        DOUBLE PRECISION :: tp, tp1, tp2, tp3, tp4, tp5, tp6

        CHARACTER (LEN=99) :: filename

!CalMtrxCNDet_On = 0 -> LU decomposition
!CalMtrxCNDet_On = 1 -> LU decomposition & determinant
!CalMtrxCNDet_On = 2 -> LU decomposition & eigenvalues and eigenvectors

!
!solve the matrix by LU decomposition method
!
        IF (ItrtnorElmntn == "L") THEN
            ALLOCATE (slv_ipiv(TtlbimA))
            CALL dLUDCMP(bimAA,TtlbimA,slv_ipiv,slv_D,slv_info)
            CALL dLUBKSB(bimAA,TtlbimA,slv_ipiv,bimBB)
!$OMP PARALLEL
!$OMP DO
            DO i = 1, TtlbimA
                bimX(i) = bimBB(i)
            END DO
!$OMP END DO
!$OMP END PARALLEL
            DEALLOCATE (slv_ipiv)
            PRINT *, slv_info
        END IF

!
!solve the matrix by Lapack LU decomposition
!

        IF (ItrtnorElmntn == "I") THEN
            ALLOCATE (slv_ipiv(TtlbimA))
            !Lapack
            CALL dgesv (TtlbimA, 1, bimAA, TtlbimA, slv_ipiv, bimBB, TtlbimA, slv_info)
            IF (slv_info == 0) THEN
!$OMP PARALLEL
!$OMP DO
                DO i = 1, TtlbimA
                    bimX(i) = bimBB(i)
                END DO
!$OMP END DO
!$OMP END PARALLEL
            END IF
            PRINT *, slv_info
            DEALLOCATE (slv_ipiv)
        END IF

        DEALLOCATE (bimAA,bimBB)

    END SUBROUTINE

    SUBROUTINE Get_BRIEFGGHHIN_DH

        INTEGER :: ithprtl, jthprtl, kthprtl, i, j, k
        INTEGER :: NdA, NdB, NdC, NdD, NdE, NdF, id_tp

        INTEGER :: Ttlbim

        DOUBLE PRECISION :: pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz

        DOUBLE PRECISION :: g1, g2, g3, g4, g5, g6
        DOUBLE PRECISION :: h1, h2, h3, h4, h5, h6
        DOUBLE PRECISION :: gnsbim, hnsbim

        DOUBLE PRECISION :: bimGGINNS,bimHHINNS, tp

!*********

        Ttlbim = ttlnmbrnd

!*********
!BRIEF - G & H

        ALLOCATE (bimGGIN(Ttlbim,Ttlbim))
        ALLOCATE (bimHHIN(Ttlbim,Ttlbim))

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
        DO i = 1, Ttlbim
            DO j = 1, Ttlbim
                bimGGIN(i,j) = 0.0d0
                bimHHIN(i,j) = 0.0d0
            END DO
        END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE (i,k,ithprtl,jthprtl,id_tp) &
!$OMP & PRIVATE (NdA,NdB,NdC,NdD,NdE,NdF) &
!$OMP & PRIVATE (pr0x,pr0y,pr0z,pr0nnx,pr0nny,pr0nnz) &
!$OMP & PRIVATE (bimGGINNS,bimHHINNS,tp) &
!$OMP & PRIVATE (g1,g2,g3,g4,g5,g6,h1,h2,h3,h4,h5,h6,gnsbim,hnsbim)
!$OMP DO
        DO i = 1, ttlnmbrnd

            ithprtl = 1
            IF (nmbrprtl > 1 .AND. i > ndendID(1)) THEN
                DO id_tp = 2, nmbrprtl
                    IF (i >= ndstaID(id_tp) .AND. i <= ndendID(id_tp)) THEN
                        ithprtl = id_tp
                        EXIT
                    END IF
                END DO
            END IF

            pr0x = xnd(i)
            pr0y = ynd(i)
            pr0z = znd(i)
            pr0nnx = nnx(i)
            pr0nny = nny(i)
            pr0nnz = nnz(i)

!***
!Internal GG and HH
!   The surface where host node x0 locates
! + The surfaces that are enclosed by the surface where x0 locates
!   (just one level down needed)

            bimGGINNS = 0.0d0
            bimHHINNS = 0.0d0

            DO k = 1, ttlnmbrelmnt

                jthprtl = 1
                IF (nmbrprtl > 1 .AND. k > elendID(1)) THEN
                    DO id_tp = 2, nmbrprtl
                        IF (k >= elstaID(id_tp) .AND. k <= elendID(id_tp)) THEN
                            jthprtl = id_tp
                            EXIT
                        END IF
                    END DO
                END IF

                IF (corelnkshell(jthprtl) == ithprtl .OR. jthprtl == ithprtl) THEN

                    IF (MeshType == "L") THEN

                        NdA = elmntlnknd(k,1)
                        NdB = elmntlnknd(k,2)
                        NdC = elmntlnknd(k,3)

                        CALL CalGHLnrBRIEFLnrREAL( ink_DH(ithprtl), k, &
                        &       pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, &
                        &       g1, g2, g3, h1, h2, h3, gnsbim, hnsbim)

                        bimGGIN(i,NdA) = bimGGIN(i,NdA) + g1
                        bimGGIN(i,NdB) = bimGGIN(i,NdB) + g2
                        bimGGIN(i,NdC) = bimGGIN(i,NdC) + g3

                        bimHHIN(i,NdA) = bimHHIN(i,NdA) + h1
                        bimHHIN(i,NdB) = bimHHIN(i,NdB) + h2
                        bimHHIN(i,NdC) = bimHHIN(i,NdC) + h3

                        bimGGINNS = bimGGINNS + gnsbim
                        bimHHINNS = bimHHINNS + hnsbim

                   END IF


                    IF (MeshType == "Q") THEN

                        NdA = elmntlnknd(k,1)
                        NdB = elmntlnknd(k,2)
                        NdC = elmntlnknd(k,3)
                        NdD = elmntlnknd(k,4)
                        NdE = elmntlnknd(k,5)
                        NdF = elmntlnknd(k,6)

                        CALL CalGHQdrBRIEFLnrREAL( ink_DH(ithprtl), k, &
                        &       pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, &
                        &       g1, g2, g3, g4, g5, g6, &
                        &       h1, h2, h3, h4, h5, h6, gnsbim, hnsbim)

                        bimGGIN(i,NdA) = bimGGIN(i,NdA) + g1
                        bimGGIN(i,NdB) = bimGGIN(i,NdB) + g2
                        bimGGIN(i,NdC) = bimGGIN(i,NdC) + g3
                        bimGGIN(i,NdD) = bimGGIN(i,NdD) + g4
                        bimGGIN(i,NdE) = bimGGIN(i,NdE) + g5
                        bimGGIN(i,NdF) = bimGGIN(i,NdF) + g6

                        bimHHIN(i,NdA) = bimHHIN(i,NdA) + h1
                        bimHHIN(i,NdB) = bimHHIN(i,NdB) + h2
                        bimHHIN(i,NdC) = bimHHIN(i,NdC) + h3
                        bimHHIN(i,NdD) = bimHHIN(i,NdD) + h4
                        bimHHIN(i,NdE) = bimHHIN(i,NdE) + h5
                        bimHHIN(i,NdF) = bimHHIN(i,NdF) + h6

                        bimGGINNS = bimGGINNS + gnsbim
                        bimHHINNS = bimHHINNS + hnsbim

                    END IF

                END IF

            END DO

            bimGGIN(i,i) = bimGGIN(i,i) + bimGGINNS
            bimHHIN(i,i) = bimHHIN(i,i) + bimHHINNS

!End internal GG and HH
!***
        END DO
!$OMP END DO
!$OMP END PARALLEL

    END SUBROUTINE

    SUBROUTINE Get_BRIEFGGHHEX_DH

        INTEGER :: ithprtl, jthprtl, kthprtl, i, j, k
        INTEGER :: NdA, NdB, NdC, NdD, NdE, NdF, id_tp

        INTEGER :: Ttlbim

        DOUBLE PRECISION :: pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz

        DOUBLE PRECISION :: g1, g2, g3, g4, g5, g6
        DOUBLE PRECISION :: h1, h2, h3, h4, h5, h6
        DOUBLE PRECISION :: gnsbim, hnsbim

        DOUBLE PRECISION :: bimGGEXNS,bimHHEXNS, tp

!*********

        Ttlbim = ttlnmbrnd

!*********
!BRIEF - G & H

        ALLOCATE (bimGGEX(Ttlbim,Ttlbim))
        ALLOCATE (bimHHEX(Ttlbim,Ttlbim))

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
        DO i = 1, Ttlbim
            DO j = 1, Ttlbim
                bimGGEX(i,j) = 0.0d0
                bimHHEX(i,j) = 0.0d0
            END DO
        END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE (i,k,ithprtl,jthprtl,id_tp) &
!$OMP & PRIVATE (NdA,NdB,NdC,NdD,NdE,NdF) &
!$OMP & PRIVATE (pr0x,pr0y,pr0z,pr0nnx,pr0nny,pr0nnz) &
!$OMP & PRIVATE (bimGGEXNS,bimHHEXNS,tp) &
!$OMP & PRIVATE (g1,g2,g3,g4,g5,g6,h1,h2,h3,h4,h5,h6,gnsbim,hnsbim)
!$OMP DO
        DO i = 1, ttlnmbrnd

            ithprtl = 1
            IF (nmbrprtl > 1 .AND. i > ndendID(1)) THEN
                DO id_tp = 2, nmbrprtl
                    IF (i >= ndstaID(id_tp) .AND. i <= ndendID(id_tp)) THEN
                        ithprtl = id_tp
                        EXIT
                    END IF
                END DO
            END IF

            pr0x = xnd(i)
            pr0y = ynd(i)
            pr0z = znd(i)
            pr0nnx = nnx(i)
            pr0nny = nny(i)
            pr0nnz = nnz(i)

!***
!External GG and HH

!---
!There are two cases:
!Case 1: the very external domain with surface at infinity
!Case 2: an domain with external bounded surface: corelnkshell(ithprtl)
!For both cases, the surfaces with the same 'corelnkshell(ithprtl)' should be involved

            IF (corelnkshell(ithprtl) == 0) THEN         !Case 1

                IF (NrmlInOut(ithprtl) == 1) THEN
                    bimGGEXNS = 0.0d0
                    bimHHEXNS = 4.0d0*pai     !contribution from infinity
                END IF

                IF (NrmlInOut(ithprtl) ==-1) THEN
                    bimGGEXNS = 0.0d0
                    bimHHEXNS =-4.0d0*pai     !contribution from infinity
                END IF

                tp = exk_DH

            END IF

            IF (corelnkshell(ithprtl) > 0) THEN          !Case 2

                bimGGEXNS = 0.0d0
                bimHHEXNS = 0.0d0

                tp = ink_DH(corelnkshell(ithprtl))

            END IF

            DO k = 1, ttlnmbrelmnt

                jthprtl = 1
                IF (nmbrprtl > 1 .AND. k > elendID(1)) THEN
                    DO id_tp = 2, nmbrprtl
                        IF (k >= elstaID(id_tp) .AND. k <= elendID(id_tp)) THEN
                            jthprtl = id_tp
                            EXIT
                        END IF
                    END DO
                END IF

                IF (  corelnkshell(jthprtl) == corelnkshell(ithprtl) &
                &.OR. jthprtl == corelnkshell(ithprtl)) THEN

                    IF (MeshType == "L") THEN

                        NdA = elmntlnknd(k,1)
                        NdB = elmntlnknd(k,2)
                        NdC = elmntlnknd(k,3)

                        CALL CalGHLnrBRIEFLnrREAL( tp, k, &
                        &       pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, &
                        &       g1, g2, g3, h1, h2, h3, gnsbim, hnsbim)

                        bimGGEX(i,NdA)=bimGGEX(i,NdA)+g1
                        bimGGEX(i,NdB)=bimGGEX(i,NdB)+g2
                        bimGGEX(i,NdC)=bimGGEX(i,NdC)+g3

                        bimHHEX(i,NdA)=bimHHEX(i,NdA)+h1
                        bimHHEX(i,NdB)=bimHHEX(i,NdB)+h2
                        bimHHEX(i,NdC)=bimHHEX(i,NdC)+h3

                        bimGGEXNS = bimGGEXNS + gnsbim
                        bimHHEXNS = bimHHEXNS + hnsbim

                    END IF


                    IF (MeshType == "Q") THEN

                        NdA = elmntlnknd(k,1)
                        NdB = elmntlnknd(k,2)
                        NdC = elmntlnknd(k,3)
                        NdD = elmntlnknd(k,4)
                        NdE = elmntlnknd(k,5)
                        NdF = elmntlnknd(k,6)

                        CALL CalGHQdrBRIEFLnrREAL( tp, k, &
                        &       pr0x, pr0y, pr0z, pr0nnx, pr0nny, pr0nnz, &
                        &       g1, g2, g3, g4, g5, g6, &
                        &       h1, h2, h3, h4, h5, h6, gnsbim, hnsbim)

                        bimGGEX(i,NdA)=bimGGEX(i,NdA)+g1
                        bimGGEX(i,NdB)=bimGGEX(i,NdB)+g2
                        bimGGEX(i,NdC)=bimGGEX(i,NdC)+g3
                        bimGGEX(i,NdD)=bimGGEX(i,NdD)+g4
                        bimGGEX(i,NdE)=bimGGEX(i,NdE)+g5
                        bimGGEX(i,NdF)=bimGGEX(i,NdF)+g6

                        bimHHEX(i,NdA)=bimHHEX(i,NdA)+h1
                        bimHHEX(i,NdB)=bimHHEX(i,NdB)+h2
                        bimHHEX(i,NdC)=bimHHEX(i,NdC)+h3
                        bimHHEX(i,NdD)=bimHHEX(i,NdD)+h4
                        bimHHEX(i,NdE)=bimHHEX(i,NdE)+h5
                        bimHHEX(i,NdF)=bimHHEX(i,NdF)+h6

                        bimGGEXNS = bimGGEXNS + gnsbim
                        bimHHEXNS = bimHHEXNS + hnsbim


                    END IF

                END IF

            END DO

            bimGGEX(i,i) = bimGGEX(i,i) + bimGGEXNS
            bimHHEX(i,i) = bimHHEX(i,i) + bimHHEXNS

!---
!End calculating external GG and HH
!***

        END DO
!$OMP END DO
!$OMP END PARALLEL

    END SUBROUTINE

END MODULE DH_SurfCal_Solver
