
MODULE DH_SurfCal_VariableInt

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData

    USE DH_SurfCal_GlobalData

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE GetVariableInt_DH

        IMPLICIT NONE

        INTEGER :: i, j, ithprtl, IOS
        DOUBLE PRECISION :: tp

        IF (ttlnmbrsrc_DH > 0) THEN
            ALLOCATE(srcType_DH(ttlnmbrsrc_DH))
            ALLOCATE(srcStrength_DH(ttlnmbrsrc_DH))
            ALLOCATE(xsrc_DH(ttlnmbrsrc_DH))
            ALLOCATE(ysrc_DH(ttlnmbrsrc_DH))
            ALLOCATE(zsrc_DH(ttlnmbrsrc_DH))
            DO i = 1, ttlnmbrsrc_DH
                srcStrength_DH(i) = 0.0d0
                xsrc_DH(i) = 0.0d0
                ysrc_DH(i) = 0.0d0
                zsrc_DH(i) = 0.0d0
            END DO
        ELSE
            ALLOCATE(srcType_DH(1))
            ALLOCATE(srcStrength_DH(1))
            ALLOCATE(xsrc_DH(1))
            ALLOCATE(ysrc_DH(1))
            ALLOCATE(zsrc_DH(1))
            srcStrength_DH(1) = 0.0d0
            xsrc_DH(1) = 1.0d18
            ysrc_DH(1) = 1.0d18
            zsrc_DH(1) = 1.0d18
        END IF

        IF (ttlnmbrsrc_DH > 0) THEN
            OPEN (111, file = 'Input_Source.dat', STATUS = 'OLD', IOSTAT = IOS)
            IF (IOS /= 0) THEN
                PRINT*, "'Input_Source.dat' does not exist! Please check!"
                STOP
            END IF
            DO ithprtl = 0, nmbrprtl
                IF (srcendID_DH(ithprtl) > 0) THEN
                    READ (111, *)
                    DO i = srcstaID_DH(ithprtl), srcendID_DH(ithprtl)
                        READ(111, *) j, &
                        &            xsrc_DH(i), ysrc_DH(i), zsrc_DH(i), &
                        &            srcStrength_DH(i)
                        IF (ithprtl > 0) THEN
                            ! Modification - commented out for monomer sphere problem
                            ! xsrc_DH(i) = xsrc_DH(i) + xloctn(ithprtl)
                            ! ysrc_DH(i) = ysrc_DH(i) + yloctn(ithprtl)
                            ! zsrc_DH(i) = zsrc_DH(i) + zloctn(ithprtl)
                        END IF
                    END DO
                    READ (111, *)
                END IF
            END DO
            CLOSE (111)
        END IF

        ALLOCATE(ndBCType_DH(ttlnmbrnd))
        ALLOCATE(ndBCKnown_DH(ttlnmbrnd))

        ALLOCATE(exE1x_DH(ttlnmbrnd))
        ALLOCATE(exE1y_DH(ttlnmbrnd))
        ALLOCATE(exE1z_DH(ttlnmbrnd))
        ALLOCATE(inE1x_DH(ttlnmbrnd))
        ALLOCATE(inE1y_DH(ttlnmbrnd))
        ALLOCATE(inE1z_DH(ttlnmbrnd))
        ALLOCATE(exE2x_DH(ttlnmbrnd))
        ALLOCATE(exE2y_DH(ttlnmbrnd))
        ALLOCATE(exE2z_DH(ttlnmbrnd))
        ALLOCATE(inE2x_DH(ttlnmbrnd))
        ALLOCATE(inE2y_DH(ttlnmbrnd))
        ALLOCATE(inE2z_DH(ttlnmbrnd))
        ALLOCATE(exE3x_DH(ttlnmbrnd))
        ALLOCATE(exE3y_DH(ttlnmbrnd))
        ALLOCATE(exE3z_DH(ttlnmbrnd))
        ALLOCATE(inE3x_DH(ttlnmbrnd))
        ALLOCATE(inE3y_DH(ttlnmbrnd))
        ALLOCATE(inE3z_DH(ttlnmbrnd))
        ALLOCATE(exPhi1_DH(ttlnmbrnd))
        ALLOCATE(exPhi1dn_DH(ttlnmbrnd))
        ALLOCATE(inPhi1_DH(ttlnmbrnd))
        ALLOCATE(inPhi1dn_DH(ttlnmbrnd))
        ALLOCATE(exPhi2_DH(ttlnmbrnd))
        ALLOCATE(exPhi2dn_DH(ttlnmbrnd))
        ALLOCATE(inPhi2_DH(ttlnmbrnd))
        ALLOCATE(inPhi2dn_DH(ttlnmbrnd))
        ALLOCATE(exPhi3_DH(ttlnmbrnd))
        ALLOCATE(exPhi3dn_DH(ttlnmbrnd))
        ALLOCATE(inPhi3_DH(ttlnmbrnd))
        ALLOCATE(inPhi3dn_DH(ttlnmbrnd))
        ALLOCATE(sigma_DH(ttlnmbrnd))

!$OMP PARALLEL PRIVATE (i)
!$OMP DO
        DO i = 1, ttlnmbrnd

            ndBCKnown_DH(i) = 0.0d0

            exE1x_DH(i) = 0.0d0
            exE1y_DH(i) = 0.0d0
            exE1z_DH(i) = 0.0d0
            inE1x_DH(i) = 0.0d0
            inE1y_DH(i) = 0.0d0
            inE1z_DH(i) = 0.0d0
            exE2x_DH(i) = 0.0d0
            exE2y_DH(i) = 0.0d0
            exE2z_DH(i) = 0.0d0
            inE2x_DH(i) = 0.0d0
            inE2y_DH(i) = 0.0d0
            inE2z_DH(i) = 0.0d0
            exE3x_DH(i) = 0.0d0
            exE3y_DH(i) = 0.0d0
            exE3z_DH(i) = 0.0d0
            inE3x_DH(i) = 0.0d0
            inE3y_DH(i) = 0.0d0
            inE3z_DH(i) = 0.0d0
            exPhi1_DH(i) = 0.0d0
            exPhi1dn_DH(i) = 0.0d0
            inPhi1_DH(i) = 0.0d0
            inPhi1dn_DH(i) = 0.0d0
            exPhi2_DH(i) = 0.0d0
            exPhi2dn_DH(i) = 0.0d0
            inPhi2_DH(i) = 0.0d0
            inPhi2dn_DH(i) = 0.0d0
            exPhi3_DH(i) = 0.0d0
            exPhi3dn_DH(i) = 0.0d0
            inPhi3_DH(i) = 0.0d0
            inPhi3dn_DH(i) = 0.0d0
            sigma_DH(i) = 0.0d0

        END DO
!$OMP END DO
!$OMP END PARALLEL

    END SUBROUTINE

END MODULE
