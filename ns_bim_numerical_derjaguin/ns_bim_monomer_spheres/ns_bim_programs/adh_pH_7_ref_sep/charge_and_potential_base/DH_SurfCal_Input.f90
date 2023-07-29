
MODULE DH_SurfCal_Input

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData

    USE DH_SurfCal_GlobalData

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE GetPhysInputInt_DH

        INTEGER :: IOS, ithprtl

        DOUBLE PRECISION :: tp, tpzoff

        phi_dim =    1.602176565d0   /(8.854187817620d0) * 1.0d5              !in mV
        Efield_dim = 1.602176565d0   /(8.854187817620d0) * 1.0d5              !in mV/nm
        Energy_dim = 1.602176565d0**2/(8.854187817620d0) * 6.241509 * 1.0d1   !in eV
        Force_dim =  1.602176565d0**2/(8.854187817620d0) * 1.0d1              !in nN
        Torque_dim = 1.602176565d0**2/(8.854187817620d0) * 1.0d1              !in nN*nm
        Sigma_dim =  1.602176565d0*1.0d-1                                     !in C/m^2

        OPEN (71, FILE = "Input_Phys_DH.dat", STATUS = "OLD", IOSTAT = IOS)
        IF (IOS /= 0) THEN
            PRINT*, "'Input_Phys_DH.dat' does not exist! Please check!"
            STOP
        END IF

        READ (71, *) !Exterior medium parameters: 1/Debye length
                     !and Dielectric constant (relative permittivity)
        READ (71, *) exK_DH, exeps_DH
        READ (71, *)

        CLOSE (71)

        exnmbrsrc_DH = 0

    END SUBROUTINE


    SUBROUTINE GetPhysInput_DH

        INTEGER :: i, j, k, IOS, ithprtl

        DOUBLE PRECISION :: tp, tpzoff, tpcff

        ALLOCATE (BCType_DH(nmbrprtl))
        ALLOCATE (BCValue_DH(nmbrprtl))
        ALLOCATE (BCRead_DH(nmbrprtl))

        ALLOCATE (ink_DH(nmbrprtl))
        ALLOCATE (ineps_DH(nmbrprtl))

!        ALLOCATE (FrEngergy_DH(nmbrprtl))
        ALLOCATE (Frcx_DH(nmbrprtl))
        ALLOCATE (Frcy_DH(nmbrprtl))
        ALLOCATE (Frcz_DH(nmbrprtl))
        ALLOCATE (Trqx_DH(nmbrprtl))
        ALLOCATE (Trqy_DH(nmbrprtl))
        ALLOCATE (Trqz_DH(nmbrprtl))

        ALLOCATE (surfQ_DH(nmbrprtl))

        ALLOCATE (nmbrsrc_DH(0:nmbrprtl))
        ALLOCATE (srcstaID_DH(0:nmbrprtl))
        ALLOCATE (srcendID_DH(0:nmbrprtl))

        DO ithprtl = 0, nmbrprtl
            nmbrsrc_DH(ithprtl) = 0
            srcstaID_DH(ithprtl) = 0
            srcendID_DH(ithprtl) = 0
        END DO

        nmbrsrc_DH(0) = exnmbrsrc_DH

        IF (nmbrsrc_DH(0) > 0) THEN
            srcstaID_DH(0) = 1
            srcendID_DH(0) = nmbrsrc_DH(0)
        END IF

        ttlnmbrsrc_DH = nmbrsrc_DH(0)

        OPEN (71, FILE = "Input_Phys_DH.dat", STATUS = "OLD", IOSTAT = IOS)
        IF (IOS /= 0) THEN
            PRINT*, "'Input_Phys_DH.dat' does not exist! Please check!"
            STOP
        END IF

        DO i = 1, 3
            READ (71, *)
        END DO

        BCType_DH(:) = '2SD'
        BCValue_DH(:) = 0.0d0
        BCRead_DH(:) = 'N'

        DO ithprtl = 1, nmbrprtl

            READ (71, *)
            READ (71, *) !Interior medium parameters: 1/Debye length
                         !and Dielectric constant (relative permittivity)
            READ (71, *) ink_DH(ithprtl), ineps_DH(ithprtl)
            READ (71, *)
            READ (71, *) !number of sources
            READ (71, *) nmbrsrc_DH(ithprtl)
            READ (71, *)

            IF (nmbrsrc_DH(ithprtl) > 0) THEN
                srcstaID_DH(ithprtl) = ttlnmbrsrc_DH + 1
                srcendID_DH(ithprtl) = ttlnmbrsrc_DH + nmbrsrc_DH(ithprtl)
            END IF

            ttlnmbrsrc_DH = ttlnmbrsrc_DH + nmbrsrc_DH(ithprtl)

        END DO

        CLOSE (71)

    END SUBROUTINE

END MODULE
