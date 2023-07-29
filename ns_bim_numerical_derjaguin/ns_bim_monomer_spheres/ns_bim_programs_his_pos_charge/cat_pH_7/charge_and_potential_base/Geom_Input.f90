
MODULE Geom_Input

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE GetGeomInput_Int

        INTEGER :: IOS, ithprtl, GLQi, tpi1dGQ, tpi2dGQ

        DOUBLE PRECISION :: tp, tp1, tp2, tp3, tp4, tpcff, tpzoff
        DOUBLE PRECISION :: tp_phiincfar, tpvcmwlkw, tpwlkwdecay

        n_glqln1d = 6
        ALLOCATE (wg_glqln1d(n_glqln1d))
        ALLOCATE (xg_glqln1d(n_glqln1d))
        CALL GauLegCoeff1D(n_glqln1d,xg_glqln1d,wg_glqln1d)
        n_glqtr2d = 16
        ALLOCATE (wg_glqtr2d(n_glqtr2d))
        ALLOCATE (xg_glqtr2d(n_glqtr2d))
        ALLOCATE (yg_glqtr2d(n_glqtr2d))
        CALL GauLegUniTriCff2D(n_glqtr2d,xg_glqtr2d,yg_glqtr2d,wg_glqtr2d)
        n_glqte3d = 46
        ALLOCATE (wg_glqte3d(n_glqte3d))
        ALLOCATE (xg_glqte3d(n_glqte3d))
        ALLOCATE (yg_glqte3d(n_glqte3d))
        ALLOCATE (zg_glqte3d(n_glqte3d))
        CALL GauLegUniTetCff3D(n_glqte3d,xg_glqte3d,yg_glqte3d,zg_glqte3d,wg_glqte3d)

        OPEN (71, FILE = "Input_Geom.dat", STATUS = "OLD", IOSTAT = IOS)
        IF (IOS /= 0) THEN
            PRINT*, "'Input_Geom.dat' does not exist! Please check!"
            STOP
        END IF

        READ (71, *) !Number of the spheres, and type of mesh
        READ (71, *) nmbrprtl
        READ (71, *)

        CLOSE (71)

        ! MeshType = 'Q'
        MeshType = 'L'

    END SUBROUTINE

    SUBROUTINE GetGeomInput

        INTEGER :: i, j, k, IOS, ithprtl, id_tp

        DOUBLE PRECISION :: tp, tp1, tp2, tp3, tp4, tpcff, tpzoff
        DOUBLE PRECISION :: tp_phiincfar, tpvcmwlkw, tpwlkwdecay

        OPEN (71, FILE = "Input_Geom.dat", STATUS = "OLD", IOSTAT = IOS)
        IF (IOS /= 0) THEN
            PRINT*, "'Input_Geom.dat' does not exist! Please check!"
            STOP
        END IF

        DO i = 1, 3
            READ (71, *)
        END DO

        ALLOCATE (PrtlType(nmbrprtl))
        ALLOCATE (MeshRead(nmbrprtl))
        ALLOCATE (NrmlInOut(nmbrprtl))
        ALLOCATE (MeshGnrtn(nmbrprtl))
        ALLOCATE (Meshnlvl(nmbrprtl))
        ALLOCATE (MeshRlvlStp(nmbrprtl))
        ALLOCATE (xnrmref(nmbrprtl))
        ALLOCATE (ynrmref(nmbrprtl))
        ALLOCATE (znrmref(nmbrprtl))

        ALLOCATE (sizezoom(nmbrprtl))
        ALLOCATE (xloctn(nmbrprtl))
        ALLOCATE (yloctn(nmbrprtl))
        ALLOCATE (zloctn(nmbrprtl))

        ALLOCATE (corelnkshell(nmbrprtl))

        ALLOCATE (bvsa(nmbrprtl))
        ALLOCATE (cvsa(nmbrprtl))
        ALLOCATE (dvsa(nmbrprtl))
        ALLOCATE (bowl_a(nmbrprtl))
        ALLOCATE (bowl_b(nmbrprtl))
        ALLOCATE (dfsp_a(nmbrprtl))
        ALLOCATE (dfsp_b(nmbrprtl))
        ALLOCATE (dfsp_c(nmbrprtl))
        ALLOCATE (dfsp_l(nmbrprtl))
        ALLOCATE (dfsp_m(nmbrprtl))
        ALLOCATE (dfsp_n(nmbrprtl))
        ALLOCATE (anglecal_x(nmbrprtl))
        ALLOCATE (anglecal_y(nmbrprtl))
        ALLOCATE (anglecal_z(nmbrprtl))
        ALLOCATE (surfarea(nmbrprtl))
        ALLOCATE (volume(nmbrprtl))

        ALLOCATE (nmbrnd(nmbrprtl))
        ALLOCATE (nmbrelmnt(nmbrprtl))
        ALLOCATE (ndstaID(nmbrprtl))
        ALLOCATE (ndendID(nmbrprtl))
        ALLOCATE (elstaID(nmbrprtl))
        ALLOCATE (elendID(nmbrprtl))

        ALLOCATE (rho_den(nmbrprtl))
        ALLOCATE (xmssctr(nmbrprtl))
        ALLOCATE (ymssctr(nmbrprtl))
        ALLOCATE (zmssctr(nmbrprtl))

        PrtlType(:) = 'Sphr'
        corelnkshell(:) = 0
        MeshRead(:) = 0
        NrmlInOut(:) = 1
        xnrmref(:) = 0.0d0
        ynrmref(:) = 0.0d0
        znrmref(:) = 0.0d0
        MeshGnrtn(:) = "Icshdrl"
        MeshRlvlStp(:) = 1.0d0
        bvsa(:) = 0.0d0
        cvsa(:) = 0.0d0
        dvsa(:) = 0.0d0
        bowl_a(:) = 0.0d0
        bowl_b(:) = 0.0d0
        dfsp_a(:) = 0.0d0
        dfsp_b(:) = 0.0d0
        dfsp_c(:) = 0.0d0
        dfsp_l(:) = 0.0d0
        dfsp_m(:) = 0.0d0
        dfsp_n(:) = 0.0d0
        anglecal_x(:) = 0.0d0
        anglecal_y(:) = 0.0d0
        anglecal_z(:) = 0.0d0
        rho_den(:) = 0.0d0
        xmssctr(:) = 0.0d0
        ymssctr(:) = 0.0d0
        zmssctr(:) = 0.0d0


        DO ithprtl = 1, nmbrprtl

            READ (71, *)

            READ (71, *) !   Level of mesh on hemisphere,
            READ (71, *) Meshnlvl(ithprtl)
            READ (71, *)

            READ (71, *) !if need to zoom and move the particle
            READ (71, *) sizezoom(ithprtl), &
            &            xloctn(ithprtl),yloctn(ithprtl),zloctn(ithprtl)
            READ (71, *)

        END DO

        CLOSE (71)

    END SUBROUTINE

END MODULE
