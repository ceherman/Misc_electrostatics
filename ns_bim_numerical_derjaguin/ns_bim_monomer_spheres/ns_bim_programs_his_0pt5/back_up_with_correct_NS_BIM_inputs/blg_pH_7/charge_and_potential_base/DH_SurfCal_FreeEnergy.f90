
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

          CALL GetGeomInput
          print *, 'GetGeomInput OK'

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

              CALL GetGeomInput
              print *, 'GetGeomInput OK'

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

              CALL GetPhysBC_DH
              print *, 'GetPhysBC_DH OK'

              CALL SlvPrblm_DH
              print *, 'SlvPrblm_DH OK'

            ! Record potential (mV) at each node for each radius
              OPEN (857, FILE = "rslt_nodes_coords_and_pots.csv", STATUS="REPLACE", ACTION="WRITE")
                  WRITE(857,*) 'node_id,x_nm,y_nm,z_nm,phi_mV'
                  DO i = 1, ttlnmbrnd
                      CALL csv_write_integer(857,i,.false.,'cmr')
                      CALL csv_write_dble(857,xnd(nd_i),.false.,'cmr')
                      CALL csv_write_dble(857,ynd(nd_i),.false.,'cmr')
                      CALL csv_write_dble(857,znd(nd_i),.false.,'cmr')
                      CALL csv_write_dble(857,exPhi3_DH(i)*phi_dim,.true.,'spc')
                  END DO
              CLOSE (857)

              IF (MeshType == 'L') THEN
                  OPEN (857, FILE="rslt_node_linkages.csv", STATUS="REPLACE", ACTION="WRITE")
                      WRITE(857,*) 'element_id,node_1_id,node_2_id,node_3_id'
                      DO i = 1, ttlnmbrelmnt
                          CALL csv_write_integer(857,(i-offset),.false.,'cmr')
                          CALL csv_write_integer(857,elmntlnknd(i,1),.false.,'cmr')
                          CALL csv_write_integer(857,elmntlnknd(i,2),.false.,'cmr')
                          CALL csv_write_integer(857,elmntlnknd(i,3),.true.,'spc')
                      END DO
                  CLOSE (857)
              END IF


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
