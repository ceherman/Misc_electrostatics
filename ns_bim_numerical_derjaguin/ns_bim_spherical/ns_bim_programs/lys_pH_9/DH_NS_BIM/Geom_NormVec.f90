!   (updated on 11-Sep-2019)
!
!   The normal vector calculation at nodes is changed
!
!   The d/dt weights are changes
!
!   (updated on 11-Sep-2019)


MODULE Geom_NormVec

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData

    IMPLICIT NONE

    INTEGER, PRIVATE :: dt_option = 2

    CONTAINS

    SUBROUTINE Getndnrml

        INTEGER :: ithprtl, i, j, k, ii, ij, ik, jj, kk, itmp, elmntA, ndA, ndB, ndC
        INTEGER :: ndst,nded,elst,eled,id_tp

        DOUBLE PRECISION :: tpvctACx, tpvctACy, tpvctACz, &
                         &  tpvctBCx, tpvctBCy, tpvctBCz, &
                         &  tpelmntnrmlx, tpelmntnrmly, tpelmntnrmlz, &
                         &  tpelmntareanrml

        DOUBLE PRECISION :: tpsumndnnx, tpsumndnny, tpsumndnnz, tpsumndnnmdl
        DOUBLE PRECISION :: tpnx, tpny, tpnz, tparea

        DOUBLE PRECISION ::  dmnx, dmny, dmnz, tplmn1, tplmn2, tpll1, tpll2, tpll3, &
                            &tpmm1, tpmm2, tpmm3, tpnn1, tpnn2, tpnn3

        DOUBLE PRECISION :: dtetax,dtetay,dtetaz,dleta,dt_xix,dt_xiy,dt_xiz,dt_xi,&
                           &dteta1,dteta2,dt_xi1,dt_xi2,dfac,dl_xi
        INTEGER :: slfNdA,slfNdB,slfNdC

        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: &
        &   cbcmtrxA,cbcmtrxAA,cbcmtrxAArcd,cbcmtrxAA_inv,cbcmtrxEnd,cbcmtrxA_tr
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: &
        &   cbcmtrxB,cbcmtrxW,cbcmtrxBB,cbcmtrxXX,cbcmtrxS
        INTEGER :: slv_ipiv(6),slv_D,slv_info
        DOUBLE PRECISION :: tp, tpcbcx, tpcbcy, tpcbcz, lcbcx, lcbcy, lcbcz

        INTEGER ::  GLQi, icnt
        DOUBLE PRECISION :: tp1,tp2,tp3,tp4,tp5,tp6,tp_a,tp_b,tp_g,tpxi,tpet
        DOUBLE PRECISION :: tprx,tpry,tprz,tpEnx,tpEny,tpEnz,JcbDtmn
        DOUBLE PRECISION :: tpendx1,tpendx2,tpendx3,tpendx4,tpendx5,tpendx6
        DOUBLE PRECISION :: tpendy1,tpendy2,tpendy3,tpendy4,tpendy5,tpendy6
        DOUBLE PRECISION :: tpendz1,tpendz2,tpendz3,tpendz4,tpendz5,tpendz6
        DOUBLE PRECISION :: drx_deps,drx_dyet,dry_deps,dry_dyet,drz_deps,drz_dyet
        DOUBLE PRECISION :: dphi_deps1,dphi_deps2,dphi_deps3,&
        &                   dphi_deps4,dphi_deps5,dphi_deps6
        DOUBLE PRECISION :: dphi_dyet1,dphi_dyet2,dphi_dyet3,&
        &                   dphi_dyet4,dphi_dyet5,dphi_dyet6
        DOUBLE PRECISION :: h_eps,h_yet,tpkappa,t_teps,t_tyet

        DOUBLE PRECISION :: tpAx, tpAy, tpAz, tpBx, tpBy, tpBz, tpCx, tpCy, tpCz

!==============
!   Element area calculation

!$OMP PARALLEL PRIVATE(k,ndA,ndB,ndC) &
!$OMP & PRIVATE(tpvctACx,tpvctACy,tpvctACz) &
!$OMP & PRIVATE(tpvctBCx,tpvctBCy,tpvctBCz) &
!$OMP & PRIVATE(tpelmntnrmlx,tpelmntnrmly,tpelmntnrmlz) &
!$OMP & PRIVATE(tpelmntareanrml)
!$OMP DO
        DO k = 1, ttlnmbrelmnt

            ndA = elmntlnknd(k,1)
            ndB = elmntlnknd(k,2)
            ndC = elmntlnknd(k,3)

            tpvctACx = xnd(ndA) - xnd(ndC)
            tpvctACy = ynd(ndA) - ynd(ndC)
            tpvctACz = znd(ndA) - znd(ndC)

            tpvctBCx = xnd(ndB) - xnd(ndC)
            tpvctBCy = ynd(ndB) - ynd(ndC)
            tpvctBCz = znd(ndB) - znd(ndC)

            tpelmntnrmlx = tpvctACy*tpvctBCz - tpvctACz*tpvctBCy
            tpelmntnrmly = tpvctACz*tpvctBCx - tpvctACx*tpvctBCz
            tpelmntnrmlz = tpvctACx*tpvctBCy - tpvctACy*tpvctBCx
            tpelmntareanrml = DSQRT(  tpelmntnrmlx**2 &
            &                       + tpelmntnrmly**2 &
            &                       + tpelmntnrmlz**2   )

            elmntarea(k) = 0.5d0 * tpelmntareanrml

            nnxelmnt(k) = tpelmntnrmlx / tpelmntareanrml
            nnyelmnt(k) = tpelmntnrmly / tpelmntareanrml
            nnzelmnt(k) = tpelmntnrmlz / tpelmntareanrml

        END DO
!$OMP END DO
!$OMP END PARALLEL

!==============

!==============
!   Gauss point weight, nodal position, nodal normal vector on each element

!$OMP PARALLEL PRIVATE(k,GLQi,icnt) &
!$OMP & PRIVATE(tp1,tp2,tp3,tpxi,tpet)&
!$OMP & PRIVATE(tprx,tpry,tprz,tpEnx,tpEny,tpEnz,JcbDtmn) &
!$OMP & PRIVATE(tpendx1,tpendx2,tpendx3) &
!$OMP & PRIVATE(tpendy1,tpendy2,tpendy3) &
!$OMP & PRIVATE(tpendz1,tpendz2,tpendz3)
!$OMP DO
        DO k = 1, ttlnmbrelmnt

            tpendx1 = xnd(elmntlnknd(k,1))
            tpendx2 = xnd(elmntlnknd(k,2))
            tpendx3 = xnd(elmntlnknd(k,3))

            tpendy1 = ynd(elmntlnknd(k,1))
            tpendy2 = ynd(elmntlnknd(k,2))
            tpendy3 = ynd(elmntlnknd(k,3))

            tpendz1 = znd(elmntlnknd(k,1))
            tpendz2 = znd(elmntlnknd(k,2))
            tpendz3 = znd(elmntlnknd(k,3))

            DO GLQi = 1, n_glqtr2d

                icnt = n_glqtr2d*(k-1)+GLQi

                tpxi = xg_glqtr2d(GLQi)
                tpet = yg_glqtr2d(GLQi)

                JcbDtmn = 2.0d0*elmntarea(k)

                tpEnx = nnxelmnt(k)
                tpEny = nnyelmnt(k)
                tpEnz = nnzelmnt(k)

                tp2 = tpxi
                tp3 = tpet
                tp1 = 1.0d0 - tp2 - tp3

                tprx = tp1*tpendx1+tp2*tpendx2+tp3*tpendx3
                tpry = tp1*tpendy1+tp2*tpendy2+tp3*tpendy3
                tprz = tp1*tpendz1+tp2*tpendz2+tp3*tpendz3

                srcfmm_vec(1, icnt) = tprx
                srcfmm_vec(2, icnt) = tpry
                srcfmm_vec(3, icnt) = tprz

                srcfmm_nrm(1, icnt) = tpEnx
                srcfmm_nrm(2, icnt) = tpEny
                srcfmm_nrm(3, icnt) = tpEnz

                srcfmm_wght(icnt) = wg_glqtr2d(GLQi)*JcbDtmn

                srcfmm_wtnd(1,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp1
                srcfmm_wtnd(2,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp2
                srcfmm_wtnd(3,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp3

            END DO

        END DO
!$OMP END DO
!$OMP END PARALLEL

!==============

!==============
!   Normal, tangential vector calculation
!   Normal is based on the facet elements connecting to the node of interest
!   Weights are chosen as the combination of the element area
!   and how far the element is relative to the node
!   combine "MAX, N. Journal of Graphics Tools, Vol. 4, No. 2." &
!   "Chen, Wu Computer Aided Geometric Design 21 (2004) 447–458"
!   as: Area_j / (|g_j - v|^2) where g_j is the center of the triangle

!$OMP PARALLEL PRIVATE(i,k,elmntA,ndA,ndB,ndC) &
!$OMP & PRIVATE(tpsumndnnx,tpsumndnny,tpsumndnnz)&
!$OMP & PRIVATE(tpvctACx,tpvctACy,tpvctACz) &
!$OMP & PRIVATE(tpvctBCx,tpvctBCy,tpvctBCz) &
!$OMP & PRIVATE(tpelmntnrmlx,tpelmntnrmly,tpelmntnrmlz) &
!$OMP & PRIVATE(tpelmntareanrml) &
!$OMP & PRIVATE(tpnx,tpny,tpnz,tparea,tpsumndnnmdl) &
!$OMP & PRIVATE(tpll1,tpmm1,tpnn1,tplmn1,tpll2,tpmm2,tpnn2,tplmn2)
!$OMP DO
        DO i = 1, ttlnmbrnd

            tpsumndnnx = 0.0d0
            tpsumndnny = 0.0d0
            tpsumndnnz = 0.0d0

            DO k = 1, mxnmbrndlnkelmnt

                IF (ndlnkelmnt(i,k) /= 0) THEN

                    elmntA = ndlnkelmnt(i,k)

                    IF (i == elmntlnknd(elmntA,1)) THEN
                        ndC = elmntlnknd(elmntA,1)
                        ndA = elmntlnknd(elmntA,2)
                        ndB = elmntlnknd(elmntA,3)
                    END IF

                    IF (i == elmntlnknd(elmntA,2)) THEN
                        ndC = elmntlnknd(elmntA,2)
                        ndA = elmntlnknd(elmntA,3)
                        ndB = elmntlnknd(elmntA,1)
                    END IF

                    IF (i == elmntlnknd(elmntA,3)) THEN
                        ndC = elmntlnknd(elmntA,3)
                        ndA = elmntlnknd(elmntA,1)
                        ndB = elmntlnknd(elmntA,2)
                    END IF

                    tpvctACx = xnd(ndA) - xnd(ndC)
                    tpvctACy = ynd(ndA) - ynd(ndC)
                    tpvctACz = znd(ndA) - znd(ndC)

                    tpvctBCx = xnd(ndB) - xnd(ndC)
                    tpvctBCy = ynd(ndB) - ynd(ndC)
                    tpvctBCz = znd(ndB) - znd(ndC)

                    tpelmntnrmlx = tpvctACy*tpvctBCz - tpvctACz*tpvctBCy
                    tpelmntnrmly = tpvctACz*tpvctBCx - tpvctACx*tpvctBCz
                    tpelmntnrmlz = tpvctACx*tpvctBCy - tpvctACy*tpvctBCx

                    tpelmntareanrml=DSQRT(tpelmntnrmlx**2+tpelmntnrmly**2+tpelmntnrmlz**2)

                    tpnx = tpelmntnrmlx / tpelmntareanrml
                    tpny = tpelmntnrmly / tpelmntareanrml
                    tpnz = tpelmntnrmlz / tpelmntareanrml

                    tpelmntnrmlx = (xnd(ndA) + xnd(ndB) + xnd(ndC))/3.0d0
                    tpelmntnrmly = (ynd(ndA) + ynd(ndB) + ynd(ndC))/3.0d0
                    tpelmntnrmlz = (znd(ndA) + znd(ndB) + znd(ndC))/3.0d0

                    tparea = (tpelmntnrmlx-xnd(ndC))**2 &
                    &       +(tpelmntnrmly-ynd(ndC))**2 &
                    &       +(tpelmntnrmlz-znd(ndC))**2

                    tpsumndnnx = tpsumndnnx + tpnx*tpelmntareanrml/tparea
                    tpsumndnny = tpsumndnny + tpny*tpelmntareanrml/tparea
                    tpsumndnnz = tpsumndnnz + tpnz*tpelmntareanrml/tparea

                END IF

            END DO

            tpsumndnnmdl = DSQRT(tpsumndnnx**2 + tpsumndnny**2 + tpsumndnnz**2)
            nnx(i) = tpsumndnnx/tpsumndnnmdl
            nny(i) = tpsumndnny/tpsumndnnmdl
            nnz(i) = tpsumndnnz/tpsumndnnmdl

            !t1 dot n = 0
            IF (DABS(nnx(i)).gt.0.2d0) THEN
                tpll1 = nny(i)
                tpmm1 =-nnx(i)
                tpnn1 = 0.0d0
            ELSE
                IF (DABS(nny(i)).gt.0.2d0) THEN
                    tpll1 = 0.0d0
                    tpmm1 =-nnz(i)
                    tpnn1 = nny(i)
                ELSE
                    tpll1 = nnz(i)
                    tpmm1 = 0.0d0
                    tpnn1 =-nnx(i)
                END IF
            END IF
            tplmn1 = 1.0d0/DSQRT(tpll1**2 + tpmm1**2 + tpnn1**2)
            t1x(i) = tpll1*tplmn1
            t1y(i) = tpmm1*tplmn1
            t1z(i) = tpnn1*tplmn1

            !t2 = n \cross t1
            tpll2 = nny(i)*t1z(i)&
            &      -nnz(i)*t1y(i)
            tpmm2 = nnz(i)*t1x(i)&
            &      -nnx(i)*t1z(i)
            tpnn2 = nnx(i)*t1y(i)&
            &      -nny(i)*t1x(i)
            tplmn2 = 1.0d0/DSQRT(tpll2**2 + tpmm2**2 + tpnn2**2)
            t2x(i) = tpll2*tplmn2
            t2y(i) = tpmm2*tplmn2
            t2z(i) = tpnn2*tplmn2

        END DO
!$OMP END DO
!$OMP END PARALLEL


!==============

        DO ithprtl = 1, nmbrprtl
            tp = 0.0d0
!$OMP PARALLEL PRIVATE(tp2)
            tp2 = 0.0d0
!$OMP DO SCHEDULE(GUIDED,4) PRIVATE(k,GLQi,icnt)
            DO k = elstaID(ithprtl), elendID(ithprtl)
                DO GLQi = 1, n_glqtr2d
                    icnt = n_glqtr2d*(k-1)+GLQi
                    tp2 = tp2 + srcfmm_wght(icnt)
                END DO
            END DO
!$OMP END DO
!$OMP ATOMIC
            tp = tp + tp2
!$OMP END PARALLEL
            surfarea(ithprtl) = tp
        END DO


        DO ithprtl = 1, nmbrprtl
            tp = 0.0d0
!$OMP PARALLEL PRIVATE(tp2)
            tp2 = 0.0d0
!$OMP DO SCHEDULE(GUIDED,4) PRIVATE(tpAx,tpAy,tpAz,tpBx,tpBy,tpBz,tpCx,tpCy,tpCz,tp1)
            DO k = elstaID(ithprtl), elendID(ithprtl)
                tpAx = xnd(elmntlnknd(k,1))
                tpAy = ynd(elmntlnknd(k,1))
                tpAz = znd(elmntlnknd(k,1))
                tpBx = xnd(elmntlnknd(k,2))
                tpBy = ynd(elmntlnknd(k,2))
                tpBz = znd(elmntlnknd(k,2))
                tpCx = xnd(elmntlnknd(k,3))
                tpCy = ynd(elmntlnknd(k,3))
                tpCz = znd(elmntlnknd(k,3))
                tp1 = - tpAx*tpBy*tpCz + tpAx*tpCy*tpBz + tpBx*tpAy*tpCz &
                    & - tpBx*tpCy*tpAz - tpCx*tpAy*tpBz + tpCx*tpBy*tpAz
                tp2 = tp2 + tp1
            END DO
!$OMP END DO
!$OMP ATOMIC
            tp = tp + tp2
!$OMP END PARALLEL
            IF (NrmlInOut(ithprtl) == 1) volume(ithprtl) = tp/6.0d0
            IF (NrmlInOut(ithprtl) ==-1) volume(ithprtl) =-tp/6.0d0
        END DO

    END SUBROUTINE



    SUBROUTINE GetndnrmlQdrtcLnr

        INTEGER :: ithprtl, i, j, k, ii, ij, ik, jj, kk, itmp, elmntA, ndA, ndB, ndC
        INTEGER :: ndst,nded,elst,eled,id_tp
        INTEGER, DIMENSION (mxnmbrndlnknd2nd,2) :: Nd2ndPairRcd

        DOUBLE PRECISION :: tpvctACx, tpvctACy, tpvctACz, &
                         &  tpvctBCx, tpvctBCy, tpvctBCz, &
                         &  tpelmntnrmlx, tpelmntnrmly, tpelmntnrmlz, &
                         &  tpelmntareanrml

        DOUBLE PRECISION :: tpsumndnnx, tpsumndnny, tpsumndnnz, tpsumndnnmdl
        DOUBLE PRECISION :: tpnx, tpny, tpnz, tparea

        DOUBLE PRECISION ::  dmnx, dmny, dmnz, tplmn1, tplmn2, tpll1, tpll2, tpll3, &
                            &tpmm1, tpmm2, tpmm3, tpnn1, tpnn2, tpnn3

        DOUBLE PRECISION :: dtetax,dtetay,dtetaz,dleta,dt_xix,dt_xiy,dt_xiz,dt_xi,&
                           &dteta1,dteta2,dt_xi1,dt_xi2,dfac,dl_xi
        INTEGER :: slfNdA,slfNdB,slfNdC,slfNdD,slfNdE,slfNdF

        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: &
        &   cbcmtrxA,cbcmtrxAA,cbcmtrxAArcd,cbcmtrxAA_inv,cbcmtrxEnd,cbcmtrxA_tr
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: &
        &   cbcmtrxB,cbcmtrxW,cbcmtrxBB,cbcmtrxXX,cbcmtrxS
        INTEGER :: slv_ipiv(6),slv_D,slv_info
        DOUBLE PRECISION :: tp, tpcbcx, tpcbcy, tpcbcz, lcbcx, lcbcy, lcbcz

        INTEGER ::  GLQi, icnt
        DOUBLE PRECISION :: tp1,tp2,tp3,tp4,tp5,tp6,tp_a,tp_b,tp_g,tpxi,tpet
        DOUBLE PRECISION :: tprx,tpry,tprz,tpEnx,tpEny,tpEnz,JcbDtmn
        DOUBLE PRECISION :: tpendx1,tpendx2,tpendx3,tpendx4,tpendx5,tpendx6
        DOUBLE PRECISION :: tpendy1,tpendy2,tpendy3,tpendy4,tpendy5,tpendy6
        DOUBLE PRECISION :: tpendz1,tpendz2,tpendz3,tpendz4,tpendz5,tpendz6
        DOUBLE PRECISION :: drx_deps,drx_dyet,dry_deps,dry_dyet,drz_deps,drz_dyet
        DOUBLE PRECISION :: dphi_deps1,dphi_deps2,dphi_deps3,&
        &                   dphi_deps4,dphi_deps5,dphi_deps6
        DOUBLE PRECISION :: dphi_dyet1,dphi_dyet2,dphi_dyet3,&
        &                   dphi_dyet4,dphi_dyet5,dphi_dyet6
        DOUBLE PRECISION :: h_eps,h_yet,tpkappa,t_teps,t_tyet

        DOUBLE PRECISION :: tpAx, tpAy, tpAz, tpBx, tpBy, tpBz, tpCx, tpCy, tpCz

        CHARACTER (LEN=1) :: ndnrml_mth = 'L'


!==============
!   Gauss point weight, nodal position, nodal normal vector on each element

!$OMP PARALLEL PRIVATE(k,GLQi,icnt) &
!$OMP & PRIVATE(tp1,tp2,tp3,tp4,tp5,tp6,tp_a,tp_b,tp_g,tpxi,tpet)&
!$OMP & PRIVATE(tprx,tpry,tprz,tpEnx,tpEny,tpEnz,JcbDtmn) &
!$OMP & PRIVATE(tpendx1,tpendx2,tpendx3,tpendx4,tpendx5,tpendx6) &
!$OMP & PRIVATE(tpendy1,tpendy2,tpendy3,tpendy4,tpendy5,tpendy6) &
!$OMP & PRIVATE(tpendz1,tpendz2,tpendz3,tpendz4,tpendz5,tpendz6) &
!$OMP & PRIVATE(drx_deps,drx_dyet,dry_deps,dry_dyet,drz_deps,drz_dyet)
!$OMP DO
        DO k = 1, ttlnmbrelmnt

            tpendx1 = xnd(elmntlnknd(k,1))
            tpendx2 = xnd(elmntlnknd(k,2))
            tpendx3 = xnd(elmntlnknd(k,3))
            tpendx4 = xnd(elmntlnknd(k,4))
            tpendx5 = xnd(elmntlnknd(k,5))
            tpendx6 = xnd(elmntlnknd(k,6))
            tpendy1 = ynd(elmntlnknd(k,1))
            tpendy2 = ynd(elmntlnknd(k,2))
            tpendy3 = ynd(elmntlnknd(k,3))
            tpendy4 = ynd(elmntlnknd(k,4))
            tpendy5 = ynd(elmntlnknd(k,5))
            tpendy6 = ynd(elmntlnknd(k,6))
            tpendz1 = znd(elmntlnknd(k,1))
            tpendz2 = znd(elmntlnknd(k,2))
            tpendz3 = znd(elmntlnknd(k,3))
            tpendz4 = znd(elmntlnknd(k,4))
            tpendz5 = znd(elmntlnknd(k,5))
            tpendz6 = znd(elmntlnknd(k,6))

            tp1 = DSQRT( (tpendx4-tpendx2)**2&
            &           +(tpendy4-tpendy2)**2&
            &           +(tpendz4-tpendz2)**2 )
            tp2 = DSQRT( (tpendx4-tpendx1)**2&
            &           +(tpendy4-tpendy1)**2&
            &           +(tpendz4-tpendz1)**2 )
            tp_a = 1.0d0/(1.0d0 + tp1/tp2)
            tp1 = DSQRT( (tpendx6-tpendx3)**2&
            &           +(tpendy6-tpendy3)**2&
            &           +(tpendz6-tpendz3)**2 )
            tp2 = DSQRT( (tpendx6-tpendx1)**2&
            &           +(tpendy6-tpendy1)**2&
            &           +(tpendz6-tpendz1)**2 )
            tp_b = 1.0d0/(1.0d0 + tp1/tp2)
            tp1 = DSQRT( (tpendx5-tpendx2)**2&
            &           +(tpendy5-tpendy2)**2&
            &           +(tpendz5-tpendz2)**2 )
            tp2 = DSQRT( (tpendx5-tpendx3)**2&
            &           +(tpendy5-tpendy3)**2&
            &           +(tpendz5-tpendz3)**2 )
            tp_g = 1.0d0/(1.0d0 + tp1/tp2)

            DO GLQi = 1, n_glqtr2d

                icnt = n_glqtr2d*(k-1)+GLQi

                tpxi = xg_glqtr2d(GLQi)
                tpet = yg_glqtr2d(GLQi)

                tp2= 1.0d0/(1.0d0-tp_a)*(tpxi-tp_a+(tp_a-tp_g)/(1.0d0-tp_g)*tpet+tpxi)
                tp3= 1.0d0/(1.0d0-tp_b)*tpet*(tp_b+tp_g-1.0d0)/(tp_g)
                tp4= 1.0d0/(tp_a*(1.0d0-tp_a))*(1.0d0-tpxi-tpet-tpxi)
                tp5= 1.0d0/(tp_g*(1.0d0-tp_g))*tpet
                tp6=-1.0d0/(tp_b*(1.0d0-tp_b))*tpet
                tp1=-tp2-tp3-tp4-tp5-tp6

                drx_deps = tp1*tpendx1+tp2*tpendx2+tp3*tpendx3 &
                &         +tp4*tpendx4+tp5*tpendx5+tp6*tpendx6
                dry_deps = tp1*tpendy1+tp2*tpendy2+tp3*tpendy3 &
                &         +tp4*tpendy4+tp5*tpendy5+tp6*tpendy6
                drz_deps = tp1*tpendz1+tp2*tpendz2+tp3*tpendz3 &
                &         +tp4*tpendz4+tp5*tpendz5+tp6*tpendz6

                tp2= 1.0d0/(1.0d0-tp_a)*tpxi*(tp_a-tp_g)/(1.0d0-tp_g)
                tp3= 1.0d0/(1.0d0-tp_b)*(tpet-tp_b+(tp_b+tp_g-1.0d0)/(tp_g)*tpxi+tpet)
                tp4=-1.0d0/(tp_a*(1.0d0-tp_a))*tpxi
                tp5= 1.0d0/(tp_g*(1.0d0-tp_g))*tpxi
                tp6= 1.0d0/(tp_b*(1.0d0-tp_b))*(1.0d0-tpxi-tpet-tpet)
                tp1=-tp2-tp3-tp4-tp5-tp6

                drx_dyet = tp1*tpendx1+tp2*tpendx2+tp3*tpendx3 &
                &         +tp4*tpendx4+tp5*tpendx5+tp6*tpendx6
                dry_dyet = tp1*tpendy1+tp2*tpendy2+tp3*tpendy3 &
                &         +tp4*tpendy4+tp5*tpendy5+tp6*tpendy6
                drz_dyet = tp1*tpendz1+tp2*tpendz2+tp3*tpendz3 &
                &         +tp4*tpendz4+tp5*tpendz5+tp6*tpendz6

                tpEnx = dry_deps*drz_dyet - drz_deps*dry_dyet
                tpEny = drz_deps*drx_dyet - drx_deps*drz_dyet
                tpEnz = drx_deps*dry_dyet - dry_deps*drx_dyet

                JcbDtmn = DSQRT(tpEnx**2+tpEny**2+tpEnz**2)

                tpEnx = tpEnx/JcbDtmn
                tpEny = tpEny/JcbDtmn
                tpEnz = tpEnz/JcbDtmn

                tp2 = 1.0d0/(1.0d0-tp_a)*tpxi*(tpxi-tp_a+(tp_a-tp_g)/(1.0d0-tp_g)*tpet)
                tp3 = 1.0d0/(1.0d0-tp_b)*tpet*(tpet-tp_b+(tp_b+tp_g-1.0d0)/(tp_g)*tpxi)
                tp4 = 1.0d0/(tp_a*(1.0d0-tp_a))*tpxi*(1.0d0-tpxi-tpet)
                tp5 = 1.0d0/(tp_g*(1.0d0-tp_g))*tpxi*tpet
                tp6 = 1.0d0/(tp_b*(1.0d0-tp_b))*tpet*(1.0d0-tpxi-tpet)
                tp1 = 1.0d0-tp2-tp3-tp4-tp5-tp6

                tprx = tp1*tpendx1+tp2*tpendx2+tp3*tpendx3 &
                &     +tp4*tpendx4+tp5*tpendx5+tp6*tpendx6
                tpry = tp1*tpendy1+tp2*tpendy2+tp3*tpendy3 &
                &     +tp4*tpendy4+tp5*tpendy5+tp6*tpendy6
                tprz = tp1*tpendz1+tp2*tpendz2+tp3*tpendz3 &
                &     +tp4*tpendz4+tp5*tpendz5+tp6*tpendz6

                srcfmm_vec(1, icnt) = tprx
                srcfmm_vec(2, icnt) = tpry
                srcfmm_vec(3, icnt) = tprz

                srcfmm_nrm(1, icnt) = tpEnx
                srcfmm_nrm(2, icnt) = tpEny
                srcfmm_nrm(3, icnt) = tpEnz

                srcfmm_wght(icnt) = wg_glqtr2d(GLQi)*JcbDtmn

                srcfmm_wtnd(1,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp1
                srcfmm_wtnd(2,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp2
                srcfmm_wtnd(3,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp3
                srcfmm_wtnd(4,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp4
                srcfmm_wtnd(5,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp5
                srcfmm_wtnd(6,icnt) = wg_glqtr2d(GLQi)*JcbDtmn*tp6

            END DO

        END DO
!$OMP END DO
!$OMP END PARALLEL
!==============

!==============
!   Normal, tangential vector calculation
!   Normal, tangential vector calculation
!   Normal is based on the facet elements connecting to the node of interest
!   Weights are chosen as the combination of the element area
!   and how far the element is relative to the node
!   combine "MAX, N. Journal of Graphics Tools, Vol. 4, No. 2." &
!   "Chen, Wu Computer Aided Geometric Design 21 (2004) 447–458"
!   as: Area_j / (|g_j - v|^2) where g_j is the center of the triangle

        IF (ndnrml_mth == 'L') THEN

!$OMP PARALLEL PRIVATE(i,k,elmntA,ndA,ndB,ndC) &
!$OMP & PRIVATE(tpsumndnnx,tpsumndnny,tpsumndnnz)&
!$OMP & PRIVATE(tpvctACx,tpvctACy,tpvctACz) &
!$OMP & PRIVATE(tpvctBCx,tpvctBCy,tpvctBCz) &
!$OMP & PRIVATE(tpelmntnrmlx,tpelmntnrmly,tpelmntnrmlz) &
!$OMP & PRIVATE(tpelmntareanrml) &
!$OMP & PRIVATE(tpnx,tpny,tpnz,tparea,tpsumndnnmdl) &
!$OMP & PRIVATE(tpll1,tpmm1,tpnn1,tplmn1,tpll2,tpmm2,tpnn2,tplmn2)
!$OMP DO
            DO i = 1, ttlnmbrnd

                tpsumndnnx = 0.0d0
                tpsumndnny = 0.0d0
                tpsumndnnz = 0.0d0

                DO k = 1, mxnmbrndlnkelmntlnr

                    IF (ndlnkelmntlnr(i,k) /= 0) THEN

                        elmntA = ndlnkelmntlnr(i,k)

                        IF (i == elmntlnkndlnr(elmntA,1)) THEN
                            ndC = elmntlnkndlnr(elmntA,1)
                            ndA = elmntlnkndlnr(elmntA,2)
                            ndB = elmntlnkndlnr(elmntA,3)
                        END IF

                        IF (i == elmntlnkndlnr(elmntA,2)) THEN
                            ndC = elmntlnkndlnr(elmntA,2)
                            ndA = elmntlnkndlnr(elmntA,3)
                            ndB = elmntlnkndlnr(elmntA,1)
                        END IF

                        IF (i == elmntlnkndlnr(elmntA,3)) THEN
                            ndC = elmntlnkndlnr(elmntA,3)
                            ndA = elmntlnkndlnr(elmntA,1)
                            ndB = elmntlnkndlnr(elmntA,2)
                        END IF

                        tpvctACx = xnd(ndA) - xnd(ndC)
                        tpvctACy = ynd(ndA) - ynd(ndC)
                        tpvctACz = znd(ndA) - znd(ndC)

                        tpvctBCx = xnd(ndB) - xnd(ndC)
                        tpvctBCy = ynd(ndB) - ynd(ndC)
                        tpvctBCz = znd(ndB) - znd(ndC)

                        tpelmntnrmlx = tpvctACy*tpvctBCz - tpvctACz*tpvctBCy
                        tpelmntnrmly = tpvctACz*tpvctBCx - tpvctACx*tpvctBCz
                        tpelmntnrmlz = tpvctACx*tpvctBCy - tpvctACy*tpvctBCx

                        tpelmntareanrml=DSQRT( tpelmntnrmlx**2 &
                        &                     +tpelmntnrmly**2 &
                        &                     +tpelmntnrmlz**2  )

                        tpnx = tpelmntnrmlx / tpelmntareanrml
                        tpny = tpelmntnrmly / tpelmntareanrml
                        tpnz = tpelmntnrmlz / tpelmntareanrml

                        tpelmntnrmlx = (xnd(ndA) + xnd(ndB) + xnd(ndC))/3.0d0
                        tpelmntnrmly = (ynd(ndA) + ynd(ndB) + ynd(ndC))/3.0d0
                        tpelmntnrmlz = (znd(ndA) + znd(ndB) + znd(ndC))/3.0d0

                        tparea = (tpelmntnrmlx-xnd(ndC))**2 &
                        &       +(tpelmntnrmly-ynd(ndC))**2 &
                        &       +(tpelmntnrmlz-znd(ndC))**2

                        tpsumndnnx = tpsumndnnx + tpnx*tpelmntareanrml/tparea
                        tpsumndnny = tpsumndnny + tpny*tpelmntareanrml/tparea
                        tpsumndnnz = tpsumndnnz + tpnz*tpelmntareanrml/tparea

                    END IF

                END DO

                tpsumndnnmdl = DSQRT(tpsumndnnx**2 + tpsumndnny**2 + tpsumndnnz**2)
                nnx(i) = tpsumndnnx/tpsumndnnmdl
                nny(i) = tpsumndnny/tpsumndnnmdl
                nnz(i) = tpsumndnnz/tpsumndnnmdl

                !t1 dot n = 0
                IF (DABS(nnx(i)).gt.0.2d0) THEN
                    tpll1 = nny(i)
                    tpmm1 =-nnx(i)
                    tpnn1 = 0.0d0
                ELSE
                    IF (DABS(nny(i)).gt.0.2d0) THEN
                        tpll1 = 0.0d0
                        tpmm1 =-nnz(i)
                        tpnn1 = nny(i)
                    ELSE
                        tpll1 = nnz(i)
                        tpmm1 = 0.0d0
                        tpnn1 =-nnx(i)
                    END IF
                END IF
                tplmn1 = 1.0d0/DSQRT(tpll1**2 + tpmm1**2 + tpnn1**2)
                t1x(i) = tpll1*tplmn1
                t1y(i) = tpmm1*tplmn1
                t1z(i) = tpnn1*tplmn1

                !t2 = n \cross t1
                tpll2 = nny(i)*t1z(i)&
                    &  -nnz(i)*t1y(i)
                tpmm2 = nnz(i)*t1x(i)&
                    &  -nnx(i)*t1z(i)
                tpnn2 = nnx(i)*t1y(i)&
                    &  -nny(i)*t1x(i)
                tplmn2 = 1.0d0/DSQRT(tpll2**2 + tpmm2**2 + tpnn2**2)
                t2x(i) = tpll2*tplmn2
                t2y(i) = tpmm2*tplmn2
                t2z(i) = tpnn2*tplmn2

            END DO
!$OMP END DO
!$OMP END PARALLEL

        END IF


!==============

        DO ithprtl = 1, nmbrprtl
            tp = 0.0d0
!$OMP PARALLEL PRIVATE(tp2)
            tp2 = 0.0d0
!$OMP DO SCHEDULE(GUIDED,4) PRIVATE(k,GLQi,icnt)
            DO k = elstaID(ithprtl), elendID(ithprtl)
                DO GLQi = 1, n_glqtr2d
                    icnt = n_glqtr2d*(k-1)+GLQi
                    tp2 = tp2 + srcfmm_wght(icnt)
                END DO
            END DO
!$OMP END DO
!$OMP ATOMIC
            tp = tp + tp2
!$OMP END PARALLEL
            surfarea(ithprtl) = tp
        END DO

        DO ithprtl = 1, nmbrprtl
            tp = 0.0d0
!$OMP PARALLEL PRIVATE(tp2)
            tp2 = 0.0d0
!$OMP DO SCHEDULE(GUIDED,4) PRIVATE(tpAx,tpAy,tpAz,tpBx,tpBy,tpBz,tpCx,tpCy,tpCz,tp1)
            DO k = elstaID(ithprtl), elendID(ithprtl)
                tpAx = xnd(elmntlnknd(k,1))
                tpAy = ynd(elmntlnknd(k,1))
                tpAz = znd(elmntlnknd(k,1))
                tpBx = xnd(elmntlnknd(k,4))
                tpBy = ynd(elmntlnknd(k,4))
                tpBz = znd(elmntlnknd(k,4))
                tpCx = xnd(elmntlnknd(k,6))
                tpCy = ynd(elmntlnknd(k,6))
                tpCz = znd(elmntlnknd(k,6))
                tp1 = - tpAx*tpBy*tpCz + tpAx*tpCy*tpBz + tpBx*tpAy*tpCz &
                    & - tpBx*tpCy*tpAz - tpCx*tpAy*tpBz + tpCx*tpBy*tpAz
                tp2 = tp2 + tp1
                tpAx = xnd(elmntlnknd(k,4))
                tpAy = ynd(elmntlnknd(k,4))
                tpAz = znd(elmntlnknd(k,4))
                tpBx = xnd(elmntlnknd(k,2))
                tpBy = ynd(elmntlnknd(k,2))
                tpBz = znd(elmntlnknd(k,2))
                tpCx = xnd(elmntlnknd(k,5))
                tpCy = ynd(elmntlnknd(k,5))
                tpCz = znd(elmntlnknd(k,5))
                tp1 = - tpAx*tpBy*tpCz + tpAx*tpCy*tpBz + tpBx*tpAy*tpCz &
                    & - tpBx*tpCy*tpAz - tpCx*tpAy*tpBz + tpCx*tpBy*tpAz
                tp2 = tp2 + tp1
                tpAx = xnd(elmntlnknd(k,6))
                tpAy = ynd(elmntlnknd(k,6))
                tpAz = znd(elmntlnknd(k,6))
                tpBx = xnd(elmntlnknd(k,5))
                tpBy = ynd(elmntlnknd(k,5))
                tpBz = znd(elmntlnknd(k,5))
                tpCx = xnd(elmntlnknd(k,3))
                tpCy = ynd(elmntlnknd(k,3))
                tpCz = znd(elmntlnknd(k,3))
                tp1 = - tpAx*tpBy*tpCz + tpAx*tpCy*tpBz + tpBx*tpAy*tpCz &
                    & - tpBx*tpCy*tpAz - tpCx*tpAy*tpBz + tpCx*tpBy*tpAz
                tp2 = tp2 + tp1
                tpAx = xnd(elmntlnknd(k,4))
                tpAy = ynd(elmntlnknd(k,4))
                tpAz = znd(elmntlnknd(k,4))
                tpBx = xnd(elmntlnknd(k,5))
                tpBy = ynd(elmntlnknd(k,5))
                tpBz = znd(elmntlnknd(k,5))
                tpCx = xnd(elmntlnknd(k,6))
                tpCy = ynd(elmntlnknd(k,6))
                tpCz = znd(elmntlnknd(k,6))
                tp1 = - tpAx*tpBy*tpCz + tpAx*tpCy*tpBz + tpBx*tpAy*tpCz &
                    & - tpBx*tpCy*tpAz - tpCx*tpAy*tpBz + tpCx*tpBy*tpAz
                tp2 = tp2 + tp1
            END DO
!$OMP END DO
!$OMP ATOMIC
            tp = tp + tp2
!$OMP END PARALLEL
            IF (NrmlInOut(ithprtl) == 1) volume(ithprtl) = tp/6.0d0
            IF (NrmlInOut(ithprtl) ==-1) volume(ithprtl) =-tp/6.0d0
        END DO

    END SUBROUTINE

END MODULE

