
MODULE DH_SurfCal_GlobalData

    IMPLICIT NONE

    DOUBLE PRECISION :: distance_2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!=========================================================================================
!Domain Global
    INTEGER :: FLDCalOn
    CHARACTER (LEN = 2)  :: Postxyz_RMpln
    DOUBLE PRECISION :: PostOfstx_RMpln,PostOfsty_RMpln,PostOfstz_RMpln
    DOUBLE PRECISION :: PostEgSZx_RMpln,PostEgSZy_RMpln,PostEgSZz_RMpln
    INTEGER :: PostEgnd_RMpln
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!=========================================================================================

    DOUBLE PRECISION :: phi_dim,Efield_dim,Energy_dim,Force_dim,Torque_dim,Sigma_dim

    INTEGER :: Cal_FreeEnergy_DH = 1

    CHARACTER (LEN = 2) :: plnType_DH = 'NO'

    CHARACTER (LEN = 1) :: plnxyz_DH

    DOUBLE PRECISION :: plnPhiE0_DH

    CHARACTER (LEN = 6) :: SurfFrEngMthd_DH

    DOUBLE PRECISION :: exk_DH

    DOUBLE PRECISION :: exeps_DH

    CHARACTER (LEN = 3), ALLOCATABLE, DIMENSION (:) :: BCType_DH

    CHARACTER (LEN = 1), ALLOCATABLE, DIMENSION (:) :: BCRead_DH

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: BCValue_DH

    CHARACTER (LEN = 1), ALLOCATABLE, DIMENSION (:) :: ndBCType_DH
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: ndBCKnown_DH

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: ink_DH, ineps_DH

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: &
    &   exE1x_DH, exE1y_DH, exE1z_DH, inE1x_DH, inE1y_DH, inE1z_DH, &
    &   exE2x_DH, exE2y_DH, exE2z_DH, inE2x_DH, inE2y_DH, inE2z_DH, &
    &   exE3x_DH, exE3y_DH, exE3z_DH, inE3x_DH, inE3y_DH, inE3z_DH, &
    &   exPhi1_DH, exPhi1dn_DH, inPhi1_DH, inPhi1dn_DH, &
    &   exPhi2_DH, exPhi2dn_DH, inPhi2_DH, inPhi2dn_DH, &
    &   exPhi3_DH, exPhi3dn_DH, inPhi3_DH, inPhi3dn_DH, &
    &   plnPhi_DH, plnSigma_DH, plnEx_DH, plnEy_DH, plnEz_DH, &
    &   sigma_DH, phi_dphidn_inf_DH

    INTEGER :: ttlnmbrsrc_DH, exnmbrsrc_DH
    INTEGER, ALLOCATABLE, DIMENSION (:) :: nmbrsrc_DH, srcstaID_DH, srcendID_DH
    INTEGER, ALLOCATABLE, DIMENSION (:) :: srcType_DH
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: &
    &   srcStrength_DH, xsrc_DH, ysrc_DH, zsrc_DH

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: &
    &   FrEngergy_DH, Frcx_DH, Frcy_DH, Frcz_DH, Trqx_DH, Trqy_DH, Trqz_DH, &
    &   surfQ_DH

    DOUBLE PRECISION :: ttlFrEngergy_DH

END MODULE
