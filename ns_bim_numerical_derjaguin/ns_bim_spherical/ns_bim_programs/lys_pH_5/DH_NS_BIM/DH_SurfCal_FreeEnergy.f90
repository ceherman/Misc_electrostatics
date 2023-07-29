
MODULE DH_SurfCal_FreeEnergy

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

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE GetPhidPhidnInf

        LOGICAL :: filexist

        INTEGER :: igroup, ithprtl, jthprtl, i, id_tp, r_index, j, j_up, &
        & r_nodes, k, offset, nd_i

        DOUBLE PRECISION :: j_radius, j_step, r_1, ref_farfield = 1.0d12

        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: r_radius, r_phi, &
        & r_theta, r_x, r_y, r_z

        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: protein_phi, protein_charge
        CHARACTER (LEN=50) :: radius_string_1
        CHARACTER (LEN=3) :: radius_string_2

        j_up   = 2


          CALL GetGeomInput
          print *, 'GetGeomInput OK'

          r_1 = sizezoom(2)

          CALL Meshediting
          print *, 'Meshediting OK'

          ALLOCATE (particle_group(nmbrprtl))
          id_tp = 0
          DO  ithprtl = 1, nmbrprtl
              particle_group(ithprtl) = 0
              IF (corelnkshell(ithprtl) == 0) THEN
                  id_tp = id_tp + 1
                  particle_group(ithprtl) = id_tp
              END IF
          END DO
          ttlgroup = id_tp
          DO
              id_tp = 0
              DO ithprtl = 1, nmbrprtl
                  IF (particle_group(ithprtl) == 0) THEN
                      id_tp = 1
                      jthprtl = ithprtl
                      EXIT
                  END IF
              END DO
              IF (id_tp == 0) EXIT
              DO ithprtl = 1, nmbrprtl
                  IF (corelnkshell(jthprtl) == ithprtl) THEN
                      particle_group(jthprtl) = particle_group(ithprtl)
                      EXIT
                  END IF
              END DO
          END DO

          ALLOCATE (phi_dphidn_inf_DH(ttlnmbrnd))

          ALLOCATE (protein_phi((ndendID(2)-ndstaID(2)+2), j_up))
          ALLOCATE (protein_charge((ndendID(2)-ndstaID(2)+2), j_up))

          r_nodes = ndendID(2)-ndstaID(2)+1



          DEALLOCATE (PrtlType,MeshRead,NrmlInOut,MeshGnrtn,Meshnlvl,MeshRlvlStp, &
          &           xnrmref,ynrmref,znrmref)
          DEALLOCATE (sizezoom,xloctn,yloctn,zloctn)
          DEALLOCATE (corelnkshell)
          DEALLOCATE (bvsa,cvsa,dvsa,bowl_a,bowl_b, &
          &           dfsp_a,dfsp_b,dfsp_c,dfsp_l,dfsp_m,dfsp_n, &
          &           anglecal_x,anglecal_y,anglecal_z,surfarea,volume)
          DEALLOCATE (nmbrnd,nmbrelmnt,ndstaID,ndendID,elstaID,elendID)
          DEALLOCATE (rho_den,xmssctr,ymssctr,zmssctr)

          DEALLOCATE (xnd,ynd,znd,nnx,nny,nnz,t1x,t1y,t1z,t2x,t2y,t2z)
          DEALLOCATE (curvt1,curvt2,curvmn,curvt1th,curvt2th,curvmnth)
          DEALLOCATE (elmntarea,nnxelmnt,nnyelmnt,nnzelmnt)
          DEALLOCATE (srcfmm_vec,srcfmm_nrm,srcfmm_wght,srcfmm_wtnd)
          DEALLOCATE (elmntlnknd,ndlnkelmnt, &
          &           ndlnknd1st,ndlnknd1stslf,ndlnknd2nd,ndlnknd2ndslf)
          IF (MeshType == "Q") DEALLOCATE (elmntlnkndlnr,ndlnkelmntlnr)
          DEALLOCATE (d_dt1,d_dt2)

          DO i = 1, ttlnmbrnd
              phi_dphidn_inf_DH(i) = 0.0d0
          END DO

      !    DO igroup = 1, ttlgroup


          j_step = (5.0-r_1)/(j_up- 1)
          DO j = 1, (j_up-1)
              j_radius = r_1 + (j-1)*j_step
              protein_phi(1, j) = j_radius
              protein_charge(1, j) = j_radius
              ! print *, j_radius-0.15

              CALL GetGeomInput
              print *, 'GetGeomInput OK'

              xloctn(1) = ref_farfield
              yloctn(1) = ref_farfield
              zloctn(1) = ref_farfield

              xloctn(2) = 0.0d0
              yloctn(2) = 0.0d0
              zloctn(2) = 0.0d0

              ! Change the radius
              sizezoom(2) = j_radius

              CALL Meshediting
              print *, 'Meshediting OK'

              IF (MeshType == 'L') CALL Getndnrml
              print *, 'Getndnrml OK'
              IF (MeshType == 'Q') CALL GetndnrmlQdrtcLnr
              print *, 'GetndnrmlQdrtcLnr OK'

              CALL GetPhysInput_DH
              print *, 'GetPhysInput_DH OK'

              CALL GetVariableInt_DH
              print *, 'GetVariableInt_DH OK'

              ! Reposition charges for particle 2
              ALLOCATE(r_radius(srcendID_DH(2)-srcstaID_DH(2)+1))
              ALLOCATE(r_phi(srcendID_DH(2)-srcstaID_DH(2)+1))
              ALLOCATE(r_theta(srcendID_DH(2)-srcstaID_DH(2)+1))
              ALLOCATE(r_x(srcendID_DH(2)-srcstaID_DH(2)+1))
              ALLOCATE(r_y(srcendID_DH(2)-srcstaID_DH(2)+1))
              ALLOCATE(r_z(srcendID_DH(2)-srcstaID_DH(2)+1))

              DO i = srcstaID_DH(2), srcendID_DH(2)
                  r_index = i - srcstaID_DH(2) + 1

                  r_radius(r_index) = dsqrt(xsrc_DH(i)**2 + ysrc_DH(i)**2 + zsrc_DH(i)**2)
                  r_phi(r_index)    = atan2(ysrc_DH(i), xsrc_DH(i))
                  r_theta(r_index)  = acos(zsrc_DH(i)/r_radius(r_index))

                  xsrc_DH(i) = (j_radius-0.15)*sin(r_theta(r_index))*cos(r_phi(r_index))
                  ysrc_DH(i) = (j_radius-0.15)*sin(r_theta(r_index))*sin(r_phi(r_index))
                  zsrc_DH(i) = (j_radius-0.15)*cos(r_theta(r_index))

                  !! Check:
                  ! r_x(r_index) = r_radius(r_index)*sin(r_theta(r_index))*cos(r_phi(r_index))
                  ! r_y(r_index) = r_radius(r_index)*sin(r_theta(r_index))*sin(r_phi(r_index))
                  ! r_z(r_index) = r_radius(r_index)*cos(r_theta(r_index))
                  !
                  ! print*,''
                  ! print*, xsrc_DH(i) - r_x(r_index)
                  ! print*, ysrc_DH(i) - r_y(r_index)
                  ! print*, zsrc_DH(i) - r_z(r_index)

              END DO

              CALL GetPhysBC_DH
              print *, 'GetPhysBC_DH OK'

              CALL SlvPrblm_DH
              print *, 'SlvPrblm_DH OK'

              DO i = ndstaID(2), ndendID(2)
                  r_index = i - ndstaID(2) + 2
                  protein_phi(r_index, j) = exPhi3_DH(i)*phi_dim
                  protein_charge(r_index, j) = exPhi3dn_DH(i)*exeps_DH*Sigma_dim
              END DO

              ! Record node coordinates for the current radius; particle 2 is the protein
              write(radius_string_1, *) (INT(j_radius*100 + 0.1))
              radius_string_2 = ADJUSTL(trim(radius_string_1))
              offset = nmbrnd(1)
              OPEN (857, FILE = "rslt_node_coordinates.csv", STATUS="REPLACE", ACTION="WRITE")
                  WRITE(857,*) 'node_id,x,y,z'
                  DO i = 1, nmbrnd(2)
                      nd_i = i + offset
                      CALL csv_write_integer(857,i,.false.,'cmr')
                      CALL csv_write_dble(857,xnd(nd_i),.false.,'cmr')
                      CALL csv_write_dble(857,ynd(nd_i),.false.,'cmr')
                      CALL csv_write_dble(857,znd(nd_i),.true.,'spc')
                  END DO
              CLOSE (857)


              IF (j .eq. (j_up-1)) then
                  ! Record node linkages for quadratic surface elements
                  IF (MeshType == 'Q') THEN
                      offset = elendID(1)
                      OPEN (857, FILE="rslt_node_linkages.csv", STATUS="REPLACE", ACTION="WRITE")
                          WRITE(857,*) 'element_id,node_1_id,node_2_id,node_3_id,node_4_id,node_5_id,node_6_id'
                          DO i = elstaID(2), elendID(2)
                              nd_i = i + offset
                              CALL csv_write_integer(857,(i-offset),.false.,'cmr')
                              CALL csv_write_integer(857,elmntlnknd(i,1),.false.,'cmr')
                              CALL csv_write_integer(857,elmntlnknd(i,2),.false.,'cmr')
                              CALL csv_write_integer(857,elmntlnknd(i,3),.false.,'cmr')
                              CALL csv_write_integer(857,elmntlnknd(i,4),.false.,'cmr')
                              CALL csv_write_integer(857,elmntlnknd(i,5),.false.,'cmr')
                              CALL csv_write_integer(857,elmntlnknd(i,6),.true.,'spc')
                          END DO
                      CLOSE (857)
                  END IF

                  IF (MeshType == 'L') THEN
                      offset = elendID(1)
                      OPEN (857, FILE="rslt_node_linkages.csv", STATUS="REPLACE", ACTION="WRITE")
                          WRITE(857,*) 'element_id,node_1_id,node_2_id,node_3_id'
                          DO i = elstaID(2), elendID(2)
                              nd_i = i + offset
                              CALL csv_write_integer(857,(i-offset),.false.,'cmr')
                              CALL csv_write_integer(857,elmntlnknd(i,1),.false.,'cmr')
                              CALL csv_write_integer(857,elmntlnknd(i,2),.false.,'cmr')
                              CALL csv_write_integer(857,elmntlnknd(i,3),.true.,'spc')
                          END DO
                      CLOSE (857)
                  END IF

                  ! Record potential (mV) at each node for each radius
                  ! First column is the radius, all remaining columns are the potentials for a given node
                  protein_phi = TRANSPOSE(protein_phi)
                  OPEN (857, FILE="rslt_protein_potential_variable_radius.csv", STATUS="REPLACE", ACTION="WRITE")
                      DO i = 1, j_up
                          DO k = 1, r_nodes+1
                              IF (k .lt. r_nodes+1) THEN
                                  CALL csv_write_dble(857,protein_phi(i, k),.false.,'cmr')
                              ELSE
                                  CALL csv_write_dble(857,protein_phi(i, k),.true.,'spc')
                              END IF
                          END DO
                      END DO
                  CLOSE (857)

                  ! Record charge density (C/m^2) at each node for each radius
                  ! First column is the radius, all remaining columns are the charge densities for a given node
                  protein_charge = TRANSPOSE(protein_charge)
                  OPEN (857, FILE="rslt_protein_charge_variable_radius.csv", STATUS="REPLACE", ACTION="WRITE")
                      DO i = 1, j_up
                          DO k = 1, r_nodes+1
                              IF (k .lt. r_nodes+1) THEN
                                  CALL csv_write_dble(857,protein_charge(i, k),.false.,'cmr')
                              ELSE
                                  CALL csv_write_dble(857,protein_charge(i, k),.true.,'spc')
                              END IF
                          END DO
                      END DO
                  CLOSE (857)
              END IF

              ! DO i = 1, ttlnmbrnd
              !     phi_dphidn_inf_DH(i) = exPhi3_DH(i) * exPhi3dn_DH(i)
              ! END DO

              DEALLOCATE (PrtlType,MeshRead,NrmlInOut,MeshGnrtn,Meshnlvl,MeshRlvlStp, &
              &           xnrmref,ynrmref,znrmref)
              DEALLOCATE (sizezoom,xloctn,yloctn,zloctn)
              DEALLOCATE (corelnkshell)
              DEALLOCATE (bvsa,cvsa,dvsa,bowl_a,bowl_b, &
              &           dfsp_a,dfsp_b,dfsp_c,dfsp_l,dfsp_m,dfsp_n, &
              &           anglecal_x,anglecal_y,anglecal_z,surfarea,volume)
              DEALLOCATE (nmbrnd,nmbrelmnt,ndstaID,ndendID,elstaID,elendID)
              DEALLOCATE (rho_den,xmssctr,ymssctr,zmssctr)

              DEALLOCATE (xnd,ynd,znd,nnx,nny,nnz,t1x,t1y,t1z,t2x,t2y,t2z)
              DEALLOCATE (curvt1,curvt2,curvmn,curvt1th,curvt2th,curvmnth)
              DEALLOCATE (elmntarea,nnxelmnt,nnyelmnt,nnzelmnt)
              DEALLOCATE (srcfmm_vec,srcfmm_nrm,srcfmm_wght,srcfmm_wtnd)
              DEALLOCATE (elmntlnknd,ndlnkelmnt, &
              &           ndlnknd1st,ndlnknd1stslf,ndlnknd2nd,ndlnknd2ndslf)
              IF (MeshType == "Q") DEALLOCATE (elmntlnkndlnr,ndlnkelmntlnr)
              DEALLOCATE (d_dt1,d_dt2)

              DEALLOCATE (BCType_DH,BCValue_DH,BCRead_DH)
              DEALLOCATE (ink_DH,ineps_DH)
              DEALLOCATE (Frcx_DH,Frcy_DH,Frcz_DH,Trqx_DH,Trqy_DH,Trqz_DH)
              DEALLOCATE (surfQ_DH)
              DEALLOCATE (nmbrsrc_DH,srcstaID_DH,srcendID_DH,&
              &           srcType_DH,srcStrength_DH,xsrc_DH,ysrc_DH,zsrc_DH)
              DEALLOCATE (ndBCType_DH,ndBCKnown_DH)
              DEALLOCATE (exE1x_DH,exE1y_DH,exE1z_DH,inE1x_DH,inE1y_DH,inE1z_DH)
              DEALLOCATE (exE2x_DH,exE2y_DH,exE2z_DH,inE2x_DH,inE2y_DH,inE2z_DH)
              DEALLOCATE (exE3x_DH,exE3y_DH,exE3z_DH,inE3x_DH,inE3y_DH,inE3z_DH)
              DEALLOCATE (exPhi1_DH,exPhi1dn_DH,inPhi1_DH,inPhi1dn_DH)
              DEALLOCATE (exPhi2_DH,exPhi2dn_DH,inPhi2_DH,inPhi2dn_DH)
              DEALLOCATE (exPhi3_DH,exPhi3dn_DH,inPhi3_DH,inPhi3dn_DH)
              DEALLOCATE (sigma_DH)

              DEALLOCATE (r_radius, r_phi, r_theta, r_x, r_y, r_z)

         END DO

         DEALLOCATE (protein_phi, protein_charge)

      END SUBROUTINE


    SUBROUTINE GetFreeEnergy_PhidPhidn      !This is working that is confirmed by doing force path idea

        LOGICAL :: filexist

        INTEGER :: i,j,k,ithprtl,GLQi,icnt
        DOUBLE PRECISION :: tp1,tp2,tp3,tp4,tp5,tp6
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: tpphiphinrml

        ALLOCATE(tpphiphinrml(ttlnmbrnd))
        ttlFrEngergy_DH = 0.0d0

        DO ithprtl = 1, nmbrprtl

            IF (corelnkshell(ithprtl) == 0) THEN

                DO i = ndstaID(ithprtl), ndendID(ithprtl)

                    tpphiphinrml(i) = exPhi3_DH(i)*exPhi3dn_DH(i) &
                    &                -phi_dphidn_inf_DH(i)

                END DO

                IF (MeshType == 'L') THEN

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        tp1 = tpphiphinrml(elmntlnknd(k,1))
                        tp2 = tpphiphinrml(elmntlnknd(k,2))
                        tp3 = tpphiphinrml(elmntlnknd(k,3))

                        DO GLQi = 1, n_glqtr2d

                            icnt = n_glqtr2d*(k-1)+GLQi

                            ttlFrEngergy_DH = ttlFrEngergy_DH &
                            &                +srcfmm_wtnd(1,icnt) * tp1 &
                            &                +srcfmm_wtnd(2,icnt) * tp2 &
                            &                +srcfmm_wtnd(3,icnt) * tp3

                        END DO

                    END DO

                END IF

                IF (MeshType == 'Q') THEN

                    DO k = elstaID(ithprtl), elendID(ithprtl)

                        tp1 = tpphiphinrml(elmntlnknd(k,1))
                        tp2 = tpphiphinrml(elmntlnknd(k,2))
                        tp3 = tpphiphinrml(elmntlnknd(k,3))
                        tp4 = tpphiphinrml(elmntlnknd(k,4))
                        tp5 = tpphiphinrml(elmntlnknd(k,5))
                        tp6 = tpphiphinrml(elmntlnknd(k,6))

                        DO GLQi = 1, n_glqtr2d

                            icnt = n_glqtr2d*(k-1)+GLQi

                            ttlFrEngergy_DH = ttlFrEngergy_DH &
                            &                +srcfmm_wtnd(1,icnt) * tp1 &
                            &                +srcfmm_wtnd(2,icnt) * tp2 &
                            &                +srcfmm_wtnd(3,icnt) * tp3 &
                            &                +srcfmm_wtnd(4,icnt) * tp4 &
                            &                +srcfmm_wtnd(5,icnt) * tp5 &
                            &                +srcfmm_wtnd(6,icnt) * tp6

                        END DO

                    END DO

                END IF

            END IF

        END DO

        ttlFrEngergy_DH = 0.5d0 * ttlFrEngergy_DH * exeps_DH

        DEALLOCATE(tpphiphinrml)

    END SUBROUTINE

END MODULE
