
MODULE DH_SurfCal_PhysBC

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData

    USE DH_SurfCal_GlobalData

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE GetPhysBC_DH

        INTEGER :: IOS

        INTEGER :: ithprtl, i, j, ii, jj, kk
        DOUBLE PRECISION :: tp, tp1, tp2, tp3, tp4, tp5, tp6

        DO ithprtl = 1, nmbrprtl
            xmssctr(ithprtl) = xmssctr(ithprtl) + xloctn(ithprtl)
            ymssctr(ithprtl) = ymssctr(ithprtl) + yloctn(ithprtl)
            zmssctr(ithprtl) = zmssctr(ithprtl) + zloctn(ithprtl)
        END DO

        DO ithprtl = 1, nmbrprtl
            DO i = ndstaID(ithprtl), ndendID(ithprtl)
                ndBCKnown_DH(i) = BCValue_DH(ithprtl)
            END DO
        END DO

    END SUBROUTINE


END MODULE
