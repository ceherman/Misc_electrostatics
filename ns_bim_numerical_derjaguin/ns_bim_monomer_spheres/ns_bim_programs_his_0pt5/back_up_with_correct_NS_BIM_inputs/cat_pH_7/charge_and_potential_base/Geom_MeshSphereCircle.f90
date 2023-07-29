
MODULE Geom_MeshSphereCircle

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE GetMeshNdElmntIcshdrlPaper(dmnlvlprtl, dmRlvlStpScl, dmMeshType)

        IMPLICIT NONE

        INTEGER, INTENT (IN) :: dmnlvlprtl

        DOUBLE PRECISION, INTENT (IN) :: dmRlvlStpScl

        CHARACTER (LEN=1), INTENT (IN) :: dmMeshType

        DOUBLE PRECISION :: SphereMesh_pi = DATAN(1.0d0)*4.0d0

        DOUBLE PRECISION, DIMENSION (12) :: xcshd, ycshd, zcshd

        DOUBLE PRECISION :: tp, tp1, cshdxy, cshdz, &
        &                   tpAx, tpAy, tpAz, tpBx, tpBy, tpBz, &
        &                   tpi_n, tpxcdnt, tpycdnt, tpzcdnt

        INTEGER :: i, j, k, NdCnt, nlvl, ElmtCnt, lvlmk, lvlnd, lvli, pstni, pstnj

        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: xsphr, ysphr ,zsphr

        INTEGER, ALLOCATABLE, DIMENSION (:, :, :) :: NdIDRgbsCll

        INTEGER, ALLOCATABLE, DIMENSION (:, :) :: NdIDEg

        INTEGER, ALLOCATABLE, DIMENSION (:) :: ElmtLnk1, ElmtLnk2, ElmtLnk3

!=====================================================================================
! icosahedron vertrices coordinates

        tp = 2.0d0*DATAN( 2.0d0/( 1.0d0+DSQRT(5.0d0) ) )   !golden ratio with rotation
        cshdz  = DCOS(tp)
        cshdxy = DSIN(tp)

        tp = 0.4d0*SphereMesh_pi    !(2*pi/5)

        xcshd(1) = 0.0d0
        ycshd(1) = 0.0d0
        zcshd(1) = 1.0d0

        DO i = 2, 6

            xcshd(i) = cshdxy*DCOS( tp * DBLE(i-2)  )
            ycshd(i) = cshdxy*DSIN( tp * DBLE(i-2)  )
            zcshd(i) = cshdz

        END DO

        xcshd(7) =  0.0d0
        ycshd(7) =  0.0d0
        zcshd(7) = -1.0d0

        DO i = 8, 12

            xcshd(i) = cshdxy*DCOS( SphereMesh_pi + tp * DBLE(i-8)  )
            ycshd(i) = cshdxy*DSIN( SphereMesh_pi + tp * DBLE(i-8)  )
            zcshd(i) =-cshdz

        END DO

!=====================================================================================

!=====================================================================================
!sphere mesh

        nlvl = 2*dmnlvlprtl

        ALLOCATE (xsphr(10*nlvl*nlvl+2))
        ALLOCATE (ysphr(10*nlvl*nlvl+2))
        ALLOCATE (zsphr(10*nlvl*nlvl+2))

        ALLOCATE (ElmtLnk1(10*2*nlvl*nlvl))
        ALLOCATE (ElmtLnk2(10*2*nlvl*nlvl))
        ALLOCATE (ElmtLnk3(10*2*nlvl*nlvl))

        ALLOCATE (NdIDRgbsCll(10, nlvl+1, nlvl+1))

        ALLOCATE (NdIDEg(20, nlvl+1))

        DO i = 1, 12
            xsphr(i) = xcshd(i)
            ysphr(i) = ycshd(i)
            zsphr(i) = zcshd(i)
        END DO

        NdCnt = 12

    !NdIDEg
        DO k = 1, 20

            IF (k == 1) THEN ! 1_2
                NdIDEg(k,1) = 1
                NdIDEg(k,nlvl+1) = 2
            END IF
            IF (k == 2) THEN ! 1_3
                NdIDEg(k,1) = 1
                NdIDEg(k,nlvl+1) = 3
            END IF
            IF (k == 3) THEN ! 1_4
                NdIDEg(k,1) = 1
                NdIDEg(k,nlvl+1) = 4
            END IF
            IF (k == 4) THEN ! 1_5
                NdIDEg(k,1) = 1
                NdIDEg(k,nlvl+1) = 5
            END IF
            IF (k == 5) THEN ! 1_6
                NdIDEg(k,1) = 1
                NdIDEg(k,nlvl+1) = 6
            END IF

            IF (k == 6) THEN ! 2_11
                NdIDEg(k,1) = 2
                NdIDEg(k,nlvl+1) = 11
            END IF
            IF (k == 7) THEN ! 2_10
                NdIDEg(k,1) = 2
                NdIDEg(k,nlvl+1) = 10
            END IF
            IF (k == 8) THEN ! 3_11
                NdIDEg(k,1) = 3
                NdIDEg(k,nlvl+1) = 11
            END IF
            IF (k == 9) THEN ! 3_12
                NdIDEg(k,1) = 3
                NdIDEg(k,nlvl+1) = 12
            END IF
            IF (k == 10) THEN ! 4_12
                NdIDEg(k,1) = 4
                NdIDEg(k,nlvl+1) = 12
            END IF

            IF (k == 11) THEN ! 4_8
                NdIDEg(k,1) = 4
                NdIDEg(k,nlvl+1) = 8
            END IF
            IF (k == 12) THEN ! 5_8
                NdIDEg(k,1) = 5
                NdIDEg(k,nlvl+1) = 8
            END IF
            IF (k == 13) THEN ! 5_9
                NdIDEg(k,1) = 5
                NdIDEg(k,nlvl+1) = 9
            END IF
            IF (k == 14) THEN ! 6_9
                NdIDEg(k,1) = 6
                NdIDEg(k,nlvl+1) = 9
            END IF
            IF (k == 15) THEN ! 6_10
                NdIDEg(k,1) = 6
                NdIDEg(k,nlvl+1) = 10
            END IF

            IF (k == 16) THEN ! 8_7
                NdIDEg(k,1) = 8
                NdIDEg(k,nlvl+1) = 7
            END IF
            IF (k == 17) THEN ! 9_7
                NdIDEg(k,1) = 9
                NdIDEg(k,nlvl+1) = 7
            END IF
            IF (k == 18) THEN ! 10_7
                NdIDEg(k,1) = 10
                NdIDEg(k,nlvl+1) = 7
            END IF
            IF (k == 19) THEN ! 11_7
                NdIDEg(k,1) = 11
                NdIDEg(k,nlvl+1) = 7
            END IF
            IF (k == 20) THEN ! 12_7
                NdIDEg(k,1) = 12
                NdIDEg(k,nlvl+1) = 7
            END IF

            tpAx = xsphr(NdIDEg(k,1))
            tpAy = ysphr(NdIDEg(k,1))
            tpAz = zsphr(NdIDEg(k,1))
            tpBx = xsphr(NdIDEg(k,nlvl+1))
            tpBy = ysphr(NdIDEg(k,nlvl+1))
            tpBz = zsphr(NdIDEg(k,nlvl+1))

            DO i = 2, nlvl

                NdCnt = NdCnt + 1
                NdIDEg(k,i) = NdCnt


                tpi_n = DBLE(i-1)/DBLE(nlvl)
                CALL NdCodntssphr(tpAx, tpAy, tpAz, tpBx, tpBy, tpBz, &
                &                 tpi_n, tpxcdnt, tpycdnt, tpzcdnt)
                xsphr(NdIDEg(k,i)) = tpxcdnt
                ysphr(NdIDEg(k,i)) = tpycdnt
                zsphr(NdIDEg(k,i)) = tpzcdnt

            END DO

        END DO

        DO k = 1, 10

            IF (k == 1) THEN !1_2, 1_3, 2_11, 3_11
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(1,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(2,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(6,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(8,i)
                END DO
            END IF
            IF (k == 2) THEN !1_3, 1_4, 3_12, 4_12
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(2,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(3,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(9,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(10,i)
                END DO
            END IF
            IF (k == 3) THEN !1_4, 1_5, 4_8, 5_8
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(3,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(4,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(11,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(12,i)
                END DO
            END IF
            IF (k == 4) THEN !1_5, 1_6, 5_9, 6_9
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(4,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(5,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(13,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(14,i)
                END DO
            END IF
            IF (k == 5) THEN !1_6, 1_2, 6_10, 2_10
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(5,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(1,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(15,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(7,i)
                END DO
            END IF
            IF (k == 6) THEN !2_11, 2_10, 11_7, 10_7
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(6,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(7,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(19,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(18,i)
                END DO
            END IF
            IF (k == 7) THEN !3_11, 3_12, 11_7, 12_7
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(8,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(9,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(19,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(20,i)
                END DO
            END IF
            IF (k == 8) THEN !4_12, 4_8, 12_7, 8_7
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(10,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(11,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(20,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(16,i)
                END DO
            END IF
            IF (k == 9) THEN !5_8, 5_9, 8_7, 9_7
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(12,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(13,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(16,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(17,i)
                END DO
            END IF
            IF (k == 10) THEN !6_9, 6_10, 9_7, 10_7
                DO i = 1, nlvl+1
                    NdIDRgbsCll(k,i,1) = NdIDEg(14,i)
                    NdIDRgbsCll(k,1,i) = NdIDEg(15,i)
                    NdIDRgbsCll(k,nlvl+1,i) = NdIDEg(17,i)
                    NdIDRgbsCll(k,i,nlvl+1) = NdIDEg(18,i)
                END DO
            END IF

            DO lvlmk = 2, nlvl

                tpAx = xsphr(NdIDRgbsCll(k,lvlmk+1,1))
                tpAy = ysphr(NdIDRgbsCll(k,lvlmk+1,1))
                tpAz = zsphr(NdIDRgbsCll(k,lvlmk+1,1))
                tpBx = xsphr(NdIDRgbsCll(k,1,lvlmk+1))
                tpBy = ysphr(NdIDRgbsCll(k,1,lvlmk+1))
                tpBz = zsphr(NdIDRgbsCll(k,1,lvlmk+1))

                lvli = 0

                DO lvlnd = 1, lvlmk-1

                    pstni = lvlmk-(lvlnd-1)
                    pstnj = (lvlmk+2)-pstni

                    NdCnt = NdCnt + 1
                    NdIDRgbsCll(k,pstni,pstnj) = NdCnt

                    lvli = lvli + 1

                    tpi_n = DBLE(lvli)/DBLE(lvlmk)
                    CALL NdCodntssphr(tpAx, tpAy, tpAz, tpBx, tpBy, tpBz, &
                    &                 tpi_n, tpxcdnt, tpycdnt, tpzcdnt)
                    xsphr(NdIDRgbsCll(k,pstni,pstnj)) = tpxcdnt
                    ysphr(NdIDRgbsCll(k,pstni,pstnj)) = tpycdnt
                    zsphr(NdIDRgbsCll(k,pstni,pstnj)) = tpzcdnt

                END DO

            END DO
            DO lvlmk = nlvl-1, 2, -1

                tpAx = xsphr(NdIDRgbsCll(k,nlvl+1,nlvl+1-lvlmk))
                tpAy = ysphr(NdIDRgbsCll(k,nlvl+1,nlvl+1-lvlmk))
                tpAz = zsphr(NdIDRgbsCll(k,nlvl+1,nlvl+1-lvlmk))
                tpBx = xsphr(NdIDRgbsCll(k,nlvl+1-lvlmk,nlvl+1))
                tpBy = ysphr(NdIDRgbsCll(k,nlvl+1-lvlmk,nlvl+1))
                tpBz = zsphr(NdIDRgbsCll(k,nlvl+1-lvlmk,nlvl+1))

                lvli = 0

                DO lvlnd = 1, lvlmk-1

                    pstnj = (nlvl+1)-(lvlmk-lvlnd)
                    pstni = 2*(nlvl+1)-lvlmk-pstnj

                    NdCnt = NdCnt + 1
                    NdIDRgbsCll(k,pstni,pstnj) = NdCnt

                    lvli = lvli + 1

                    tpi_n = DBLE(lvli)/DBLE(lvlmk)
                    CALL NdCodntssphr(tpAx, tpAy, tpAz, tpBx, tpBy, tpBz, &
                    &                 tpi_n, tpxcdnt, tpycdnt, tpzcdnt)
                    xsphr(NdIDRgbsCll(k,pstni,pstnj)) = tpxcdnt
                    ysphr(NdIDRgbsCll(k,pstni,pstnj)) = tpycdnt
                    zsphr(NdIDRgbsCll(k,pstni,pstnj)) = tpzcdnt

                END DO

            END DO

        END DO

        OPEN (111, FILE = "Prtl_Orgnl.vrt", STATUS = "REPLACE")

        DO i = 1, 10*nlvl*nlvl+2
            WRITE (111, *) i, xsphr(i), ysphr(i), zsphr(i)
        END DO

        CLOSE (111)

        OPEN (121, FILE = "Prtl_Orgnl.cel", STATUS = "REPLACE")

        IF (dmMeshType == "L" .OR. dmMeshType == "C") THEN

            ElmtCnt = 0
            DO k =1, 10
                DO i = 1, nlvl
                    DO j = 1, nlvl
                        ElmtCnt = ElmtCnt + 1
                        WRITE (121, *) ElmtCnt, NdIDRgbsCll(k,i,j), &
                        &              NdIDRgbsCll(k,i,j+1), NdIDRgbsCll(k,i+1,j)
                        ElmtCnt = ElmtCnt + 1
                        WRITE (121, *) ElmtCnt, NdIDRgbsCll(k,i,j+1), &
                        &              NdIDRgbsCll(k,i+1,j+1), NdIDRgbsCll(k,i+1,j)
                    END DO
                END DO
            END DO

        END IF

        IF (dmMeshType == "Q") THEN

            ElmtCnt = 0
            DO k =1, 10
                DO i = 1, nlvl+1-2,2
                    DO j = 1, nlvl+1-2,2
                        ElmtCnt = ElmtCnt + 1
                        WRITE (121, *) ElmtCnt, &
                        &   NdIDRgbsCll(k,i,j), &
                        &   NdIDRgbsCll(k,i,j+2), &
                        &   NdIDRgbsCll(k,i+2,j), &
                        &   NdIDRgbsCll(k,i,j+1), &
                        &   NdIDRgbsCll(k,i+1,j+1), &
                        &   NdIDRgbsCll(k,i+1,j)
                        ElmtCnt = ElmtCnt + 1
                        WRITE (121, *) ElmtCnt, &
                        &   NdIDRgbsCll(k,i,j+2), &
                        &   NdIDRgbsCll(k,i+2,j+2), &
                        &   NdIDRgbsCll(k,i+2,j), &
                        &   NdIDRgbsCll(k,i+1,j+2), &
                        &   NdIDRgbsCll(k,i+2,j+1),&
                        &   NdIDRgbsCll(k,i+1,j+1)
                    END DO
                END DO
            END DO

        END IF

        CLOSE (121)

        OPEN (101, FILE = "Prtl_Orgnl.inp", STATUS = "REPLACE")

        WRITE (101, *)
        WRITE (101, *)
        WRITE (101, *) 10*nlvl*nlvl+2
        WRITE (101, *) ElmtCnt
        WRITE (101, *) 0.0d0, 0.0d0, 0.0d0

        CLOSE (101)

        DEALLOCATE (xsphr)
        DEALLOCATE (ysphr)
        DEALLOCATE (zsphr)
        DEALLOCATE (NdIDRgbsCll)
        DEALLOCATE (NdIDEg)
        DEALLOCATE (ElmtLnk1)
        DEALLOCATE (ElmtLnk2)
        DEALLOCATE (ElmtLnk3)

    END SUBROUTINE

    SUBROUTINE NdCodntssphr(dmtpAx, dmtpAy, dmtpAz, dmtpBx, dmtpBy, dmtpBz, &
                            &dmtpi_nlvl, dmtpxcdnt, dmtpycdnt, dmtpzcdnt)

        DOUBLE PRECISION, INTENT(IN) :: dmtpAx,dmtpAy,dmtpAz,dmtpBx,dmtpBy,dmtpBz, &
        &   dmtpi_nlvl
        DOUBLE PRECISION, INTENT(OUT) :: dmtpxcdnt, dmtpycdnt, dmtpzcdnt

        DOUBLE PRECISION, DIMENSION (3,4):: GEclM
        DOUBLE PRECISION, DIMENSION (3) :: GErslt
        DOUBLE PRECISION, DIMENSION (4) :: GETp

        DOUBLE PRECISION :: ArcAB, tp

        INTEGER :: i, j, k, p, remp

        ArcAB = DACOS(DABS(dmtpAx*dmtpBx+dmtpAy*dmtpBy+dmtpAz*dmtpBz))

        GEclM(1,4) = DCOS(dmtpi_nlvl*ArcAB)
        GEclM(2,4) = DCOS((1.0d0-dmtpi_nlvl)*ArcAB)
        GEclM(3,4) = 0.0d0

        GEclM(1,1) = dmtpAx
        GEclM(1,2) = dmtpAy
        GEclM(1,3) = dmtpAz

        GEclM(2,1) = dmtpBx
        GEclM(2,2) = dmtpBy
        GEclM(2,3) = dmtpBz

        GEclM(3,1) = (dmtpAy*dmtpBz - dmtpAz*dmtpBy)
        GEclM(3,2) = (dmtpAz*dmtpBx - dmtpAx*dmtpBz)
        GEclM(3,3) = (dmtpAx*dmtpBy - dmtpAy*dmtpBx)

        DO j=1,3-1
            tp=DABS(GEclM(j,j))
            remp = j
            DO k=j+1,3
                IF (tp<DABS(GEclM(k,j))) THEN
                    tp=DABS(GEclM(k,j))
                    remp=k
                END IF
            END DO
            IF (remp>j) THEN
                GETp(:)=GEclM(remp,:)
                GEclM(remp,:)=GEclM(j,:)
                GEclM(j,:)=GETp(:)
            END IF
            DO i=j+1,3
                GEclM(i,:)=GEclM(i,:)-GEclM(j,:)*GEclM(i,j)/GEclM(j,j)
            END DO
        END DO

        GErslt(3)=GEclM(3,3+1)/GEclM(3,3)
        DO i=3-1, 1, -1
            tp=0.0d0
            DO j=i+1,3
                tp=tp+GEclM(i,j)*GErslt(j)
            END DO
            GErslt(i)=(GEclM(i,3+1)-tp)/GEclM(i,i)
        END DO

        tp = DSQRT(GErslt(1)**2+GErslt(2)**2+GErslt(3)**2)
        dmtpxcdnt = GErslt(1)/tp
        dmtpycdnt = GErslt(2)/tp
        dmtpzcdnt = GErslt(3)/tp

    END SUBROUTINE

    

END MODULE
