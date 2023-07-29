!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!===============================================================================
!-------------------------------------------------------------------------------

!*******************************************************************************
!*******************************************************************************

!Disclaimer

!This code is written in FORTRAN 90.
!  By SUN, Qiang (qiang.sun@rmit.edu.au),
!  Australian Research Council Centre of Excellence for Nanoscale BioPhotonics
!  School of Science, RMIT University
!  Melbourne, VIC 3001, Australia

!  Version 0.1
!  Refined on Dec-2019.

!This code has been tested on
!  MAC gcc;
!  Linux gcc

!This sample code is to calculate the retation factor of a molecule next to
!  one or two ligands by using the Debye-Huckel model

!The Debye-Huckel model solver is based on
!  A robust and accurate formulation of molecular and colloidal electrostatics
!  by Qiang Sun, Evert Klaseboer, and Derek Y. C. Chan
!  The Journal of Chemical Physics 145, 054106 (2016); doi: 10.1063/1.4960033

!The retation factor is calculated based on reference -
!  Electrostatic Contributions to Protein Retention
!    in Ion-Exchange Chromatography. 1. Cytochrome c Variants
!  by Yan Yao and Abraham M. Lenhoff
!  Anal. Chem. 2004, 76, 6743-6752

!The following assumptions have been made in this code:
!  1. The quadratic triangle elements (6-node surface elements) are employed.
!  2. The quadratic interpolation functions on the element are used.
!  3. The integration over the element is mapping to the unit triangle.

!Geometric input parameters are given in file "Input_Geom.dat"
!Physical input parameters are given in file "Input_Phys_DH.dat"
!Point charge input parameters are given in file "Input_Source.dat"

!Outputs are the retation factor
!  and free energy distribution in the volume shell where the ligands locate.

!Numerical libraries, BLAS and LAPACK, are needed for solving linear matrices.

!This code supposes to be mainly used for the calculation of retation factor
!  within Lenhoff's group only.

!*******************************************************************************
!*******************************************************************************

!-------------------------------------------------------------------------------
!===============================================================================
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

INCLUDE './Pre_csvformat.f90'
INCLUDE './Pre_GQ_FDM.f90'
INCLUDE './Pre_LinearSolvers.f90'
INCLUDE './Pre_Constants.f90'



INCLUDE './Geom_GlobalData.f90'
INCLUDE './Geom_Input.f90'
INCLUDE './Geom_MeshSphereCircle.f90'
INCLUDE './Geom_Mesh.f90'
INCLUDE './Geom_NormVec.f90'



INCLUDE './BRIEFGHReal.f90'



INCLUDE './DH_SurfCal_GlobalData.f90'
INCLUDE './DH_SurfCal_Input.f90'
INCLUDE './DH_SurfCal_VariableInt.f90'
INCLUDE './DH_SurfCal_PhysBC.f90'
INCLUDE './DH_SurfCal_Solver.f90'
INCLUDE './DH_SurfCal_FreeEnergy.f90'



INCLUDE './place_2nd_charge.f90'




PROGRAM Main_DH

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData
    USE Geom_Input
    USE Geom_Mesh
    USE Geom_NormVec

    USE DH_SurfCal_GlobalData
    USE DH_SurfCal_Input
    USE DH_SurfCal_VariableInt
    USE DH_SurfCal_PhysBC
    USE DH_SurfCal_Solver
    USE DH_SurfCal_FreeEnergy

    USE place_2nd_charge

    IMPLICIT NONE


    LOGICAL :: filexist

    INTEGER :: ksmcs, Validation_Kirkwood_ion = 1
    DOUBLE PRECISION :: tempcalcs

    INTEGER :: ithprtl
    DOUBLE PRECISION :: tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8,tp9,tp10

    INTEGER, DIMENSION (8) :: LclDtTm
    CHARACTER (LEN = 12), DIMENSION (3) :: READ_REAL_CLOCK

    INTEGER :: i, j, k, id_tp, ishellnd, ishellel, GLQi, icnt
    INTEGER, ALLOCATABLE, DIMENSION (:,:) :: e_nshell
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: xshell,yshell,zshell
    DOUBLE PRECISION :: tpxi,tpyt,tpze,tpvi
    DOUBLE PRECISION :: drx_dxi,drx_dyt,drx_dze,dry_dxi,dry_dyt,dry_dze, &
    &                   drz_dxi,drz_dyt,drz_dze,JcbDtmn
    DOUBLE PRECISION :: tpendx1,tpendx2,tpendx3,tpendx4,tpendx5, &
    &                   tpendx6,tpendx7,tpendx8,tpendx9,tpendx10,&
    &                   tpendy1,tpendy2,tpendy3,tpendy4,tpendy5, &
    &                   tpendy6,tpendy7,tpendy8,tpendy9,tpendy10,&
    &                   tpendz1,tpendz2,tpendz3,tpendz4,tpendz5, &
    &                   tpendz6,tpendz7,tpendz8,tpendz9,tpendz10
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: vlmfmm_wght
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: vlmfmm_wtnd

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: ttlFrEngShell
    DOUBLE PRECISION :: ChrTetCff, sep_distance

    DOUBLE PRECISION :: tp_charge1_coor(3), tp_charge2_coor(3)

    INTEGER :: num_complete, file_loop, jk, ik
    CHARACTER (len=100) :: file_name
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: eng_res

!     CHARACTER (LEN = 50) :: x, y
!     CHARACTER (LEN = 50) :: z
!     DOUBLE PRECISION :: j_radius
!
! goto 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Record the start date and time of the simulation

    CALL DATE_AND_TIME ( READ_REAL_CLOCK (1), READ_REAL_CLOCK (2), &
                        &READ_REAL_CLOCK (3), LclDtTm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    OPEN (9001, FILE = "Rslt_PrcdSmm.dat", STATUS = "REPLACE")

    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    WRITE(9001,*) "**** This is a solver on DH problems by NSBIM ****"
    WRITE(9001,*) "**** Developed by SQ ****"
    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    WRITE(9001, "(A38,I2,A1,I2,A1,I2)") "The simulation started at (HH-MM-SS): ",&
    &                   LclDtTm(5),"-",LclDtTm(6),"-",LclDtTm(7)
    WRITE(9001, "(A17,I4,A1,I2,A1,I2,A12)") "on (YYYY-MM-DD): ",&
    &                   LclDtTm(1),"-",LclDtTm(2),"-",LclDtTm(3), " local time."
    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)

    CALL GetGeomInput_Int
    print *, 'GetGeomInput_Int OK'

    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    IF (MeshType == "L") WRITE(9001,*) "Linear triangle elements have been used."
    IF (MeshType == "Q") WRITE(9001,*) "Quadratic elements have been used."
    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    WRITE(9001,*) "The number of particles: ", nmbrprtl
    WRITE(9001,*)
    CLOSE (9001)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL GetPhysInputInt_DH
    print *, 'GetPhysInputInt_DH OK'

    OPEN (9001, FILE = "Rslt_PrcdSmm.dat", STATUS = "OLD", &
        & POSITION="APPEND", ACTION="WRITE")
    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    WRITE(9001,*) "External K (in 1/nm) and Debye Length (in nm): ", exK_DH, 1.0d0/ exK_DH
    WRITE(9001,*) "External relative eps: ", exeps_DH
    WRITE(9001,*)
    WRITE(9001,*) "**************************************************************"
    WRITE(9001,*)
    CLOSE (9001)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Get the volume mesh, GQ for integration to calculate retention coeff

    OPEN (101, FILE = 'shell_integ.msh', STATUS = "OLD")

    DO i = 1, 8
        READ (101, *)
    END DO
    READ (101, *) ishellnd

    ALLOCATE (xshell(ishellnd))
    ALLOCATE (yshell(ishellnd))
    ALLOCATE (zshell(ishellnd))
    ALLOCATE (ttlFrEngShell(ishellnd))

    CALL GetGeomInput
    tpxi = sizezoom(2) + sizezoom(1) + 0.001d0                        !smallest gap is 0.001 nm or 0.01 A
    tpyt = sizezoom(2) + sizezoom(1) + 4.0d0 * 1.0d0/exK_DH           !largest gap is 4 times of Debye length
    DEALLOCATE (PrtlType,MeshRead,NrmlInOut,MeshGnrtn,Meshnlvl,MeshRlvlStp, &
    &           xnrmref,ynrmref,znrmref)
    DEALLOCATE (sizezoom,xloctn,yloctn,zloctn)
    DEALLOCATE (corelnkshell)
    DEALLOCATE (bvsa,cvsa,dvsa,bowl_a,bowl_b, &
    &           dfsp_a,dfsp_b,dfsp_c,dfsp_l,dfsp_m,dfsp_n, &
    &           anglecal_x,anglecal_y,anglecal_z,surfarea,volume)
    DEALLOCATE (nmbrnd,nmbrelmnt,ndstaID,ndendID,elstaID,elendID)
    DEALLOCATE (rho_den,xmssctr,ymssctr,zmssctr)

    DO i = 1, ishellnd
        READ (101, *) id_tp, xshell(i), yshell(i), zshell(i)
        tp1 = xshell(i)
        tp2 = yshell(i)
        tp3 = zshell(i)
        tp4 = dsqrt(tp1**2+tp2**2+tp3**2)
        tp5 = tpxi + (tpyt - tpxi) * (tp4 - 2.0d0)/(4.0d0 - 2.0d0)      !convert mesh from 2~4 to domain needed
        xshell(i) = tp1/tp4 * tp5
        yshell(i) = tp2/tp4 * tp5
        zshell(i) = tp3/tp4 * tp5
    END DO
    READ (101, *)
    READ (101, *)
    READ (101, *) ishellel

    ALLOCATE (e_nshell(ishellel,10))

    DO k = 1, ishellel
        READ (101, *) id_tp, id_tp, id_tp, id_tp, id_tp, &
        &             e_nshell(k,1), e_nshell(k,2), e_nshell(k,3), e_nshell(k,4), &
        &             e_nshell(k,5), e_nshell(k,6), e_nshell(k,7), e_nshell(k,8), &
        &             e_nshell(k,9), e_nshell(k,10)
    END DO
    CLOSE (101)

    ALLOCATE (vlmfmm_wtnd(10,ishellel*n_glqte3d))
    ALLOCATE (vlmfmm_wght(ishellel*n_glqte3d))
    DO k = 1, ishellel

        tpendx1 = xshell(e_nshell(k, 1))
        tpendx2 = xshell(e_nshell(k, 2))
        tpendx3 = xshell(e_nshell(k, 3))
        tpendx4 = xshell(e_nshell(k, 4))
        tpendx5 = xshell(e_nshell(k, 5))
        tpendx6 = xshell(e_nshell(k, 6))
        tpendx7 = xshell(e_nshell(k, 7))
        tpendx8 = xshell(e_nshell(k, 8))
        tpendx9 = xshell(e_nshell(k, 9))
        tpendx10= xshell(e_nshell(k,10))
        tpendy1 = yshell(e_nshell(k, 1))
        tpendy2 = yshell(e_nshell(k, 2))
        tpendy3 = yshell(e_nshell(k, 3))
        tpendy4 = yshell(e_nshell(k, 4))
        tpendy5 = yshell(e_nshell(k, 5))
        tpendy6 = yshell(e_nshell(k, 6))
        tpendy7 = yshell(e_nshell(k, 7))
        tpendy8 = yshell(e_nshell(k, 8))
        tpendy9 = yshell(e_nshell(k, 9))
        tpendy10= yshell(e_nshell(k,10))
        tpendz1 = zshell(e_nshell(k, 1))
        tpendz2 = zshell(e_nshell(k, 2))
        tpendz3 = zshell(e_nshell(k, 3))
        tpendz4 = zshell(e_nshell(k, 4))
        tpendz5 = zshell(e_nshell(k, 5))
        tpendz6 = zshell(e_nshell(k, 6))
        tpendz7 = zshell(e_nshell(k, 7))
        tpendz8 = zshell(e_nshell(k, 8))
        tpendz9 = zshell(e_nshell(k, 9))
        tpendz10= zshell(e_nshell(k,10))

        DO GLQi = 1, n_glqte3d

            icnt = n_glqte3d*(k-1)+GLQi

            tpxi = xg_glqte3d(GLQi)
            tpyt = yg_glqte3d(GLQi)
            tpze = zg_glqte3d(GLQi)
            tpvi = 1.0d0 - tpxi - tpyt - tpze

            tp1 = -1.0d0*(2.0d0*tpvi - 1.0d0) + tpvi*(-2.0d0)
            tp2 =  1.0d0*(2.0d0*tpxi - 1.0d0) + tpxi*( 2.0d0)
            tp3 =  0.0d0
            tp4 =  0.0d0
            tp5 =  4.0d0*(-1.0d0)*tpxi + 4.0d0*tpvi*1.0d0
            tp6 =  4.0d0*tpyt
            tp7 =  4.0d0*(-1.0d0)*tpyt
            tp8 =  4.0d0*(-1.0d0)*tpze
            tp9 =  0.0d0
            tp10=  4.0d0*tpze

            drx_dxi  = tp1*tpendx1 + tp2*tpendx2 + tp3*tpendx3 &
            &         +tp4*tpendx4 + tp5*tpendx5 + tp6*tpendx6 &
            &         +tp7*tpendx7 + tp8*tpendx8 + tp9*tpendx9 + tp10*tpendx10
            dry_dxi  = tp1*tpendy1 + tp2*tpendy2 + tp3*tpendy3 &
            &         +tp4*tpendy4 + tp5*tpendy5 + tp6*tpendy6 &
            &         +tp7*tpendy7 + tp8*tpendy8 + tp9*tpendy9 + tp10*tpendy10
            drz_dxi  = tp1*tpendz1 + tp2*tpendz2 + tp3*tpendz3 &
            &         +tp4*tpendz4 + tp5*tpendz5 + tp6*tpendz6 &
            &         +tp7*tpendz7 + tp8*tpendz8 + tp9*tpendz9 + tp10*tpendz10

            tp1 = -1.0d0*(2.0d0*tpvi - 1.0d0) + tpvi*(-2.0d0)
            tp2 =  0.0d0
            tp3 =  1.0d0*(2.0d0*tpyt - 1.0d0) + tpyt*( 2.0d0)
            tp4 =  0.0d0
            tp5 =  4.0d0*(-1.0d0)*tpxi
            tp6 =  4.0d0*tpxi
            tp7 =  4.0d0*(-1.0d0)*tpyt + 4.0d0*tpvi*1.0d0
            tp8 =  4.0d0*(-1.0d0)*tpze
            tp9 =  4.0d0*tpze
            tp10=  0.0d0

            drx_dyt  = tp1*tpendx1 + tp2*tpendx2 + tp3*tpendx3 &
            &         +tp4*tpendx4 + tp5*tpendx5 + tp6*tpendx6 &
            &         +tp7*tpendx7 + tp8*tpendx8 + tp9*tpendx9 + tp10*tpendx10
            dry_dyt  = tp1*tpendy1 + tp2*tpendy2 + tp3*tpendy3 &
            &         +tp4*tpendy4 + tp5*tpendy5 + tp6*tpendy6 &
            &         +tp7*tpendy7 + tp8*tpendy8 + tp9*tpendy9 + tp10*tpendy10
            drz_dyt  = tp1*tpendz1 + tp2*tpendz2 + tp3*tpendz3 &
            &         +tp4*tpendz4 + tp5*tpendz5 + tp6*tpendz6 &
            &         +tp7*tpendz7 + tp8*tpendz8 + tp9*tpendz9 + tp10*tpendz10


            tp1 = -1.0d0*(2.0d0*tpvi - 1.0d0) + tpvi*(-2.0d0)
            tp2 =  0.0d0
            tp3 =  0.0d0
            tp4 =  1.0d0*(2.0d0*tpze - 1.0d0) + tpze*( 2.0d0)
            tp5 =  4.0d0*(-1.0d0)*tpxi
            tp6 =  0.0d0
            tp7 =  4.0d0*(-1.0d0)*tpyt
            tp8 =  4.0d0*(-1.0d0)*tpze + 4.0d0*tpvi*1.0d0
            tp9 =  4.0d0*tpyt
            tp10=  4.0d0*tpxi

            drx_dze  = tp1*tpendx1 + tp2*tpendx2 + tp3*tpendx3 &
            &         +tp4*tpendx4 + tp5*tpendx5 + tp6*tpendx6 &
            &         +tp7*tpendx7 + tp8*tpendx8 + tp9*tpendx9 + tp10*tpendx10
            dry_dze  = tp1*tpendy1 + tp2*tpendy2 + tp3*tpendy3 &
            &         +tp4*tpendy4 + tp5*tpendy5 + tp6*tpendy6 &
            &         +tp7*tpendy7 + tp8*tpendy8 + tp9*tpendy9 + tp10*tpendy10
            drz_dze  = tp1*tpendz1 + tp2*tpendz2 + tp3*tpendz3 &
            &         +tp4*tpendz4 + tp5*tpendz5 + tp6*tpendz6 &
            &         +tp7*tpendz7 + tp8*tpendz8 + tp9*tpendz9 + tp10*tpendz10

            JcbDtmn = drx_dxi*dry_dyt*drz_dze &
            &        +dry_dxi*drz_dyt*drx_dze &
            &        +drz_dxi*drx_dyt*dry_dze &
            &        -drx_dxi*drz_dyt*dry_dze &
            &        -drz_dxi*dry_dyt*drx_dze &
            &        -dry_dxi*drx_dyt*drz_dze


            vlmfmm_wght(icnt) = wg_glqte3d(GLQi)*JcbDtmn

            tp1 =  tpvi*(2.0d0*tpvi - 1.0d0)
            tp2 =  tpxi*(2.0d0*tpxi - 1.0d0)
            tp3 =  tpyt*(2.0d0*tpyt - 1.0d0)
            tp4 =  tpze*(2.0d0*tpze - 1.0d0)
            tp5 =  4.0d0*tpvi*tpxi
            tp6 =  4.0d0*tpxi*tpyt
            tp7 =  4.0d0*tpvi*tpyt
            tp8 =  4.0d0*tpvi*tpze
            tp9 =  4.0d0*tpyt*tpze
            tp10=  4.0d0*tpxi*tpze

            vlmfmm_wtnd( 1,icnt) = vlmfmm_wght(icnt)*tp1
            vlmfmm_wtnd( 2,icnt) = vlmfmm_wght(icnt)*tp2
            vlmfmm_wtnd( 3,icnt) = vlmfmm_wght(icnt)*tp3
            vlmfmm_wtnd( 4,icnt) = vlmfmm_wght(icnt)*tp4
            vlmfmm_wtnd( 5,icnt) = vlmfmm_wght(icnt)*tp5
            vlmfmm_wtnd( 6,icnt) = vlmfmm_wght(icnt)*tp6
            vlmfmm_wtnd( 7,icnt) = vlmfmm_wght(icnt)*tp7
            vlmfmm_wtnd( 8,icnt) = vlmfmm_wght(icnt)*tp8
            vlmfmm_wtnd( 9,icnt) = vlmfmm_wght(icnt)*tp9
            vlmfmm_wtnd(10,icnt) = vlmfmm_wght(icnt)*tp10

        END DO

    END DO

    GOTO 33

    ttlFrEngShell(:) = 1.0d0

    INQUIRE (FILE="rslt_FreeEngOrg.dat", exist=filexist)
    IF (filexist) THEN
        OPEN (111, FILE="rslt_FreeEngOrg.dat", STATUS="OLD", &
            & POSITION="APPEND", ACTION="WRITE")
    ELSE
      OPEN (111, FILE="rslt_FreeEngOrg.dat", STATUS="NEW", &
          & ACTION="WRITE")
      WRITE (111,*) 'Variables = "x in A", "y in A", "z in A", "Free Energy in eV"'
    END IF
    WRITE (111,*) 'Zone T ="Shell at Debye length ', 1.0d0/exK_DH*10.0d0, &
    &             '", n= ', ishellnd, ', e=', ishellel*8, &
    &             ', f=fepoint, et=tetrahedron'
    DO i = 1, ishellnd
        WRITE (111, *) xshell(i)*10.0d0, yshell(i)*10.0d0, zshell(i)*10.0d0, ttlFrEngShell(i)
    END DO
    DO k = 1, ishellel
        WRITE (111, *) e_nshell(k,5), e_nshell(k,2), e_nshell(k,6), e_nshell(k,10)
        WRITE (111, *) e_nshell(k,7), e_nshell(k,6), e_nshell(k,3), e_nshell(k,9)
        WRITE (111, *) e_nshell(k,8), e_nshell(k,10),e_nshell(k,9), e_nshell(k,4)
        WRITE (111, *) e_nshell(k,1), e_nshell(k,5), e_nshell(k,7), e_nshell(k,8)
        WRITE (111, *) e_nshell(k,7), e_nshell(k,10),e_nshell(k,9), e_nshell(k,8)
        WRITE (111, *) e_nshell(k,5), e_nshell(k,10),e_nshell(k,7), e_nshell(k,8)
        WRITE (111, *) e_nshell(k,7), e_nshell(k,10),e_nshell(k,9), e_nshell(k,6)
        WRITE (111, *) e_nshell(k,5), e_nshell(k,10),e_nshell(k,7), e_nshell(k,6)
    END DO
    CLOSE (111)

    tpxi = 0.0d0
    tpyt = 0.0d0
    DO k = 1, ishellel

        tp1 = dexp(-ttlFrEngShell(e_nshell(k, 1))/0.02568d0) - 1.0d0
        tp2 = dexp(-ttlFrEngShell(e_nshell(k, 2))/0.02568d0) - 1.0d0
        tp3 = dexp(-ttlFrEngShell(e_nshell(k, 3))/0.02568d0) - 1.0d0
        tp4 = dexp(-ttlFrEngShell(e_nshell(k, 4))/0.02568d0) - 1.0d0
        tp5 = dexp(-ttlFrEngShell(e_nshell(k, 5))/0.02568d0) - 1.0d0
        tp6 = dexp(-ttlFrEngShell(e_nshell(k, 6))/0.02568d0) - 1.0d0
        tp7 = dexp(-ttlFrEngShell(e_nshell(k, 7))/0.02568d0) - 1.0d0
        tp8 = dexp(-ttlFrEngShell(e_nshell(k, 8))/0.02568d0) - 1.0d0
        tp9 = dexp(-ttlFrEngShell(e_nshell(k, 9))/0.02568d0) - 1.0d0
        tp10= dexp(-ttlFrEngShell(e_nshell(k,10))/0.02568d0) - 1.0d0

        DO GLQi = 1, n_glqte3d

            icnt = n_glqte3d*(k-1)+GLQi
            tpxi = tpxi + vlmfmm_wght(icnt)

            tpyt = tpyt + vlmfmm_wtnd( 1,icnt) * tp1 &
            &           + vlmfmm_wtnd( 2,icnt) * tp2 &
            &           + vlmfmm_wtnd( 3,icnt) * tp3 &
            &           + vlmfmm_wtnd( 4,icnt) * tp4 &
            &           + vlmfmm_wtnd( 5,icnt) * tp5 &
            &           + vlmfmm_wtnd( 6,icnt) * tp6 &
            &           + vlmfmm_wtnd( 7,icnt) * tp7 &
            &           + vlmfmm_wtnd( 8,icnt) * tp8 &
            &           + vlmfmm_wtnd( 9,icnt) * tp9 &
            &           + vlmfmm_wtnd(10,icnt) * tp10

        END DO

    END DO

    tp3 = 1.9d0 + 0.24d0 + 0.001d0
    tp4 = 1.9d0 + 0.24d0 + 4.0d0 * 1.0d0/exK_DH
    tp1 = 4.0d0/3.0d0*pai*(tp4**3-tp3**3)
    tp2 = (dexp(-1.0d0/0.02568d0) - 1.0d0) * tp1
    print *, tp1, tpxi, dabs(tpxi-tp1), dabs(tpxi-tp1)/tp1
    print *, tp2, tpyt, dabs(tpyt-tp2), dabs(tpyt-tp2)/dabs(tp2)

    STOP

33  ttlFrEngShell(:) = 0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    IF (Cal_FreeEnergy_DH == 1) CALL GetPhidPhidnInf

!     DO ksmcs = 1, ishellnd, 1
! !        tempcalcs = csStartValue + DBLE(ksmcs)*csStepSize
!
!         CALL GetGeomInput
!         print *, 'GetGeomInput OK'
!         xloctn(1) = xshell(ksmcs)
!         yloctn(1) = yshell(ksmcs)
!         zloctn(1) = zshell(ksmcs)
!
!         CALL Meshediting
!         print *, 'Meshediting OK'
!
!         IF (MeshType == 'L') CALL Getndnrml
!         print *, 'Getndnrml OK'
!         IF (MeshType == 'Q') CALL GetndnrmlQdrtcLnr
!         print *, 'GetndnrmlQdrtcLnr OK'
!
! !        CALL Surf_NrmlTanCurv_Vis
! !        print *, 'Surf_NrmlTanCurv_Vis OK'
!
!         CALL GetPhysInput_DH
!         print *, 'GetPhysInput_DH OK'
!
!         CALL GetVariableInt_DH
!         print *, 'GetVariableInt_DH OK'
!
!         CALL GetPhysBC_DH
!         print *, 'GetPhysBC_DH OK'
!
!         CALL SlvPrblm_DH
!         print *, 'SlvPrblm_DH OK'
!
!         IF (Cal_FreeEnergy_DH == 1) CALL GetFreeEnergy_PhidPhidn
!         ttlFrEngShell(ksmcs) = ttlFrEngergy_DH*Energy_dim
!
!         OPEN (9001, FILE = "Rslt_PrcdSmm.dat", STATUS = "OLD", &
!             & POSITION="APPEND", ACTION="WRITE")
!         WRITE(9001,*)
!         WRITE(9001,*) "**************************************************************"
!         WRITE(9001,*)
!         WRITE(9001,*) "Case ID: ", ksmcs
!         WRITE(9001,*)
!         DO ithprtl = 1, nmbrprtl
!             WRITE(9001,*)
!             WRITE(9001,*) "------------------------------------"
!             WRITE(9001,*) "Particle ID and its core-link-shell: ", ithprtl, corelnkshell(ithprtl)
!             WRITE(9001,*) "number of nodes, node ID sta, node ID end: ", &
!             &              nmbrnd(ithprtl), ndstaID(ithprtl), ndendID(ithprtl)
!             WRITE(9001,*) "number of elements, element ID sta, element ID end: ", &
!             &              nmbrelmnt(ithprtl), elstaID(ithprtl), elendID(ithprtl)
!             WRITE(9001,*) "Particle surface area (in nm^2): ", surfarea(ithprtl)
!             WRITE(9001,*) "Particle volume (in nm^3): ", volume(ithprtl)
!             WRITE(9001,*)
!             WRITE(9001,*) "Internal K (in 1/nm) and Debye Length (in nm): ", &
!             &             inK_DH(ithprtl), 1.0d0/inK_DH(ithprtl)
!             WRITE(9001,*) "Internal relative eps: ", ineps_DH(ithprtl)
!             WRITE(9001,*) "Electrical force: ", Frcx_DH(ithprtl), Frcy_DH(ithprtl), Frcz_DH(ithprtl)
!             WRITE(9001,*) "Electrical torque: ", Trqx_DH(ithprtl), Trqy_DH(ithprtl), Trqz_DH(ithprtl)
!             WRITE(9001,*) "Total surface charge: ", surfQ_DH(ithprtl)
!             WRITE(9001,*) "------------------------------------"
!             WRITE(9001,*)
!         END DO
!         WRITE(9001,*)
!         WRITE(9001,*) "Free Energy (nondim, in eV): ", ttlFrEngergy_DH, ttlFrEngergy_DH*Energy_dim
!         WRITE(9001,*)
!         WRITE(9001,*) "**************************************************************"
!         WRITE(9001,*)
!         CLOSE (9001)
!
!         DEALLOCATE (PrtlType,MeshRead,NrmlInOut,MeshGnrtn,Meshnlvl,MeshRlvlStp, &
!         &           xnrmref,ynrmref,znrmref)
!         DEALLOCATE (sizezoom,xloctn,yloctn,zloctn)
!         DEALLOCATE (corelnkshell)
!         DEALLOCATE (bvsa,cvsa,dvsa,bowl_a,bowl_b, &
!         &           dfsp_a,dfsp_b,dfsp_c,dfsp_l,dfsp_m,dfsp_n, &
!         &           anglecal_x,anglecal_y,anglecal_z,surfarea,volume)
!         DEALLOCATE (nmbrnd,nmbrelmnt,ndstaID,ndendID,elstaID,elendID)
!         DEALLOCATE (rho_den,xmssctr,ymssctr,zmssctr)
!
!         DEALLOCATE (xnd,ynd,znd,nnx,nny,nnz,t1x,t1y,t1z,t2x,t2y,t2z)
!         DEALLOCATE (curvt1,curvt2,curvmn,curvt1th,curvt2th,curvmnth)
!         DEALLOCATE (elmntarea,nnxelmnt,nnyelmnt,nnzelmnt)
!         DEALLOCATE (srcfmm_vec,srcfmm_nrm,srcfmm_wght,srcfmm_wtnd)
!         DEALLOCATE (elmntlnknd,ndlnkelmnt, &
!         &           ndlnknd1st,ndlnknd1stslf,ndlnknd2nd,ndlnknd2ndslf)
!         IF (MeshType == "Q") DEALLOCATE (elmntlnkndlnr,ndlnkelmntlnr)
!         DEALLOCATE (d_dt1,d_dt2)
!
!         DEALLOCATE (BCType_DH,BCValue_DH,BCRead_DH)
!         DEALLOCATE (ink_DH,ineps_DH)
!         DEALLOCATE (Frcx_DH,Frcy_DH,Frcz_DH,Trqx_DH,Trqy_DH,Trqz_DH)
!         DEALLOCATE (surfQ_DH)
!         DEALLOCATE (nmbrsrc_DH,srcstaID_DH,srcendID_DH,&
!         &           srcType_DH,srcStrength_DH,xsrc_DH,ysrc_DH,zsrc_DH)
!         DEALLOCATE (ndBCType_DH,ndBCKnown_DH)
!         DEALLOCATE (exE1x_DH,exE1y_DH,exE1z_DH,inE1x_DH,inE1y_DH,inE1z_DH)
!         DEALLOCATE (exE2x_DH,exE2y_DH,exE2z_DH,inE2x_DH,inE2y_DH,inE2z_DH)
!         DEALLOCATE (exE3x_DH,exE3y_DH,exE3z_DH,inE3x_DH,inE3y_DH,inE3z_DH)
!         DEALLOCATE (exPhi1_DH,exPhi1dn_DH,inPhi1_DH,inPhi1dn_DH)
!         DEALLOCATE (exPhi2_DH,exPhi2dn_DH,inPhi2_DH,inPhi2dn_DH)
!         DEALLOCATE (exPhi3_DH,exPhi3dn_DH,inPhi3_DH,inPhi3dn_DH)
!         DEALLOCATE (sigma_DH)
!
!     END DO
!
!     ChrTetCff = 0.0d0
!     DO k = 1, ishellel
!
!         tp1 = dexp(-ttlFrEngShell(e_nshell(k, 1))/0.02568d0) - 1.0d0
!         tp2 = dexp(-ttlFrEngShell(e_nshell(k, 2))/0.02568d0) - 1.0d0
!         tp3 = dexp(-ttlFrEngShell(e_nshell(k, 3))/0.02568d0) - 1.0d0
!         tp4 = dexp(-ttlFrEngShell(e_nshell(k, 4))/0.02568d0) - 1.0d0
!         tp5 = dexp(-ttlFrEngShell(e_nshell(k, 5))/0.02568d0) - 1.0d0
!         tp6 = dexp(-ttlFrEngShell(e_nshell(k, 6))/0.02568d0) - 1.0d0
!         tp7 = dexp(-ttlFrEngShell(e_nshell(k, 7))/0.02568d0) - 1.0d0
!         tp8 = dexp(-ttlFrEngShell(e_nshell(k, 8))/0.02568d0) - 1.0d0
!         tp9 = dexp(-ttlFrEngShell(e_nshell(k, 9))/0.02568d0) - 1.0d0
!         tp10= dexp(-ttlFrEngShell(e_nshell(k,10))/0.02568d0) - 1.0d0
!
!         DO GLQi = 1, n_glqte3d
!
!             icnt = n_glqte3d*(k-1)+GLQi
!
!             ChrTetCff = ChrTetCff + vlmfmm_wtnd( 1,icnt) * tp1 &
!             &                     + vlmfmm_wtnd( 2,icnt) * tp2 &
!             &                     + vlmfmm_wtnd( 3,icnt) * tp3 &
!             &                     + vlmfmm_wtnd( 4,icnt) * tp4 &
!             &                     + vlmfmm_wtnd( 5,icnt) * tp5 &
!             &                     + vlmfmm_wtnd( 6,icnt) * tp6 &
!             &                     + vlmfmm_wtnd( 7,icnt) * tp7 &
!             &                     + vlmfmm_wtnd( 8,icnt) * tp8 &
!             &                     + vlmfmm_wtnd( 9,icnt) * tp9 &
!             &                     + vlmfmm_wtnd(10,icnt) * tp10
!
!         END DO
!
!     END DO
!
!     ChrTetCff = ChrTetCff * 5.19d-2
!
!     INQUIRE (FILE="rslt_ChrTetCff.dat", exist=filexist)
!     IF (filexist) THEN
!         OPEN (111, FILE="rslt_ChrTetCff.dat", STATUS="OLD", &
!             & POSITION="APPEND", ACTION="WRITE")
!     ELSE
!         OPEN (111, FILE="rslt_ChrTetCff.dat", STATUS="NEW", &
!             & ACTION="WRITE")
!         CALL csv_write_char(111,'Variables="Debye Length in A", "ChrTetCff"',.true.,'spc')
!     END IF
!     CALL csv_write_dble(111,1.0d0/exK_DH*10.0d0,.false.,'cmr')
!     CALL csv_write_dble(111,ChrTetCff,.true.,'spc')
!     CLOSE (111)
!
!     INQUIRE (FILE="rslt_FreeEng.dat", exist=filexist)
!     IF (filexist) THEN
!         OPEN (111, FILE="rslt_FreeEng.dat", STATUS="OLD", &
!             & POSITION="APPEND", ACTION="WRITE")
!     ELSE
!       OPEN (111, FILE="rslt_FreeEng.dat", STATUS="NEW", &
!           & ACTION="WRITE")
!       WRITE (111,*) 'Variables = "x in A", "y in A", "z in A", "Free Energy in eV"'
!     END IF
!     WRITE (111,*) 'Zone T ="Shell at Debye length ', 1.0d0/exK_DH*10.0d0, &
!     &             '", n= ', ishellnd, ', e=', ishellel*8, &
!     &             ', f=fepoint, et=tetrahedron'
!     DO i = 1, ishellnd
!         WRITE (111, *) xshell(i)*10.0d0, yshell(i)*10.0d0, zshell(i)*10.0d0, ttlFrEngShell(i)
!     END DO
!     DO k = 1, ishellel
!         WRITE (111, *) e_nshell(k,5), e_nshell(k,2), e_nshell(k,6), e_nshell(k,10)
!         WRITE (111, *) e_nshell(k,7), e_nshell(k,6), e_nshell(k,3), e_nshell(k,9)
!         WRITE (111, *) e_nshell(k,8), e_nshell(k,10),e_nshell(k,9), e_nshell(k,4)
!         WRITE (111, *) e_nshell(k,1), e_nshell(k,5), e_nshell(k,7), e_nshell(k,8)
!         WRITE (111, *) e_nshell(k,7), e_nshell(k,10),e_nshell(k,9), e_nshell(k,8)
!         WRITE (111, *) e_nshell(k,5), e_nshell(k,10),e_nshell(k,7), e_nshell(k,8)
!         WRITE (111, *) e_nshell(k,7), e_nshell(k,10),e_nshell(k,9), e_nshell(k,6)
!         WRITE (111, *) e_nshell(k,5), e_nshell(k,10),e_nshell(k,7), e_nshell(k,6)
!     END DO
!     CLOSE (111)
!
!
!     OPEN (857, FILE="rslt_FreeEng.csv", STATUS="REPLACE", ACTION="WRITE")
!         CALL csv_write_char(857,'x,y,z,eV,sep(all_in_A)',.true.,'spc')
!         DO ik = 1, ishellnd
!             CALL csv_write_dble(857,xshell(ik)*10.0d0,.false.,'cmr')
!             CALL csv_write_dble(857,yshell(ik)*10.0d0,.false.,'cmr')
!             CALL csv_write_dble(857,zshell(ik)*10.0d0,.false.,'cmr')
!             CALL csv_write_dble(857,ttlFrEngShell(ik),.false.,'cmr')
!             sep_distance = dsqrt((xshell(ik)*10.0d0)**2 + (yshell(ik)*10.0d0)**2 &
!                            & + (zshell(ik)*10.0d0)**2) - (sizezoom(2) + sizezoom(1) + 0.001d0)*10.0
!             CALL csv_write_dble(857,sep_distance,.true.,'spc')
!         END DO
!     CLOSE (857)
!
!
!     DEALLOCATE (xshell,yshell,zshell)
!     DEALLOCATE (e_nshell)
!     DEALLOCATE (vlmfmm_wght,vlmfmm_wtnd)
!     DEALLOCATE (ttlFrEngShell)
!
!     IF (Cal_FreeEnergy_DH == 1) DEALLOCATE (particle_group,phi_dphidn_inf_DH)
!     DEALLOCATE (wg_glqln1d,xg_glqln1d,wg_glqtr2d,xg_glqtr2d,yg_glqtr2d, &
!     &           wg_glqte3d,xg_glqte3d,yg_glqte3d,zg_glqte3d)
!
!     CALL DATE_AND_TIME (READ_REAL_CLOCK (1), READ_REAL_CLOCK (2), &
!     &                  READ_REAL_CLOCK (3), LclDtTm)
!
!     OPEN (9001, FILE = "Rslt_PrcdSmm.dat", STATUS = "OLD", &
!             & POSITION="APPEND", ACTION="WRITE")
!     WRITE(9001, "(A38,I2,A1,I2,A1,I2)") "The simulation ended at (HH-MM-SS): ",&
!     &                   LclDtTm(5),"-",LclDtTm(6),"-",LclDtTm(7)
!     WRITE(9001, "(A17,I4,A1,I2,A1,I2,A12)") "on (YYYY-MM-DD): ",&
!     &                   LclDtTm(1),"-",LclDtTm(2),"-",LclDtTm(3), " local time."
!     WRITE(9001,*)
!     WRITE(9001,*) "**************************************************************"
!     WRITE(9001,*)
!     WRITE(9001,*)
!     WRITE(9001,*)
!     CLOSE(9001)

! 1   x = 'Hi'
!     ! y = trim(x) // ', Chase.'
!     print*,x
!     ! print*,y
!
!     j_radius = 1.43
!     write(z, *) (INT(j_radius*100 + 0.1))
!     z = trim(x)//ADJUSTL(trim(z))
!     print*,z



END PROGRAM Main_DH
