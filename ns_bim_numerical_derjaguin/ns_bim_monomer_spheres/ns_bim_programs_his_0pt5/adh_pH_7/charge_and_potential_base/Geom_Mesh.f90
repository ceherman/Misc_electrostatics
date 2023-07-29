
MODULE Geom_Mesh

    USE omp_lib

    USE Pre_Constants
    USE Pre_csvformat

    USE Geom_GlobalData
    USE Geom_MeshSphereCircle

!****************************************************************************************!
!This Module is a long one to generate surface mesh on the sphere.
!All needed connectivity relationships between nodes and elements are obtained.
!****************************************************************************************!

    IMPLICIT NONE

    DOUBLE PRECISION :: dmMeshRlvlStp

    CONTAINS

    SUBROUTINE Meshediting

        CHARACTER (LEN=50) :: filename

        LOGICAL :: filexist

        DOUBLE PRECISION :: xrefnn, yrefnn, zrefnn

        INTEGER :: i, j, k, ithprtl, ndindex, elmntindex
        INTEGER :: ndA, ndB, ndC, ndD, ndE, ndF, IOS

        INTEGER :: nd_i, nd_j, nd_k, elmnt_i, elmnt_j, elmnt_k

        INTEGER :: ndoffset1, elmntoffset1

        INTEGER :: thti, tpswpndid

        INTEGER :: nmbrndsznd, nmbrelmntsznd

        DOUBLE PRECISION :: tpvctACx, tpvctACy, tpvctACz, &
        &                   tpvctBCx, tpvctBCy, tpvctBCz, &
        &                   tpvctrefCx, tpvctrefCy, tpvctrefCz, &
        &                   tpnrmlchck, rdbscl, tpxyscl, tpyzscl, &
        &                   tpx, tpy, tpz, tp1, tp2, tp3

        DOUBLE PRECISION :: tpswp, tpthtx, tpthty, tpthtz

        ttlnmbrnd = 0
        ttlnmbrelmnt = 0

        DO ithprtl = 1, nmbrprtl

            CALL GetMeshNdElmntIcshdrlPaper &
                        &(Meshnlvl(ithprtl),MeshRlvlStp(ithprtl),MeshType)

            OPEN (101, FILE = "Prtl_Orgnl.inp", STATUS = "OLD", IOSTAT = IOS)
            IF (IOS /= 0) THEN
                PRINT*, '"Prtl_Orgnl.inp" does not exist! Please check!'
                STOP
            END IF

            READ (101, *)
            READ (101, *)
            READ (101, *) nmbrndsznd
            READ (101, *) nmbrelmntsznd
            READ (101, *) xrefnn, yrefnn, zrefnn

            CLOSE (101)

            ttlnmbrnd = ttlnmbrnd + nmbrndsznd
            ttlnmbrelmnt = ttlnmbrelmnt + nmbrelmntsznd

        END DO

        ALLOCATE(xnd(ttlnmbrnd))
        ALLOCATE(ynd(ttlnmbrnd))
        ALLOCATE(znd(ttlnmbrnd))
        ALLOCATE(nnx(ttlnmbrnd))
        ALLOCATE(nny(ttlnmbrnd))
        ALLOCATE(nnz(ttlnmbrnd))
        ALLOCATE(t1x(ttlnmbrnd))
        ALLOCATE(t1y(ttlnmbrnd))
        ALLOCATE(t1z(ttlnmbrnd))
        ALLOCATE(t2x(ttlnmbrnd))
        ALLOCATE(t2y(ttlnmbrnd))
        ALLOCATE(t2z(ttlnmbrnd))
        ALLOCATE(curvt1(ttlnmbrnd))
        ALLOCATE(curvt2(ttlnmbrnd))
        ALLOCATE(curvmn(ttlnmbrnd))
        ALLOCATE(curvt1th(ttlnmbrnd))
        ALLOCATE(curvt2th(ttlnmbrnd))
        ALLOCATE(curvmnth(ttlnmbrnd))

!$OMP PARALLEL PRIVATE (i)
!$OMP DO
        DO i = 1, ttlnmbrnd
            xnd(i) = 0.0d0
            ynd(i) = 0.0d0
            znd(i) = 0.0d0
            nnx(i) = 0.0d0
            nny(i) = 0.0d0
            nnz(i) = 0.0d0
            t1x(i) = 0.0d0
            t1y(i) = 0.0d0
            t1z(i) = 0.0d0
            t2x(i) = 0.0d0
            t2y(i) = 0.0d0
            t2z(i) = 0.0d0
            curvt1(i) = 0.0d0
            curvt2(i) = 0.0d0
            curvmn(i) = 0.0d0
            curvt1th(i) = 0.0d0
            curvt2th(i) = 0.0d0
            curvmnth(i) = 0.0d0
        END DO
!$OMP END DO
!$OMP END PARALLEL

        ALLOCATE(elmntarea(ttlnmbrelmnt))
        ALLOCATE(nnxelmnt (ttlnmbrelmnt))
        ALLOCATE(nnyelmnt (ttlnmbrelmnt))
        ALLOCATE(nnzelmnt (ttlnmbrelmnt))

!$OMP PARALLEL PRIVATE (i)
!$OMP DO
        DO i = 1, ttlnmbrelmnt
            elmntarea(i) = 0.0d0
            nnxelmnt(i)  = 0.0d0
            nnyelmnt(i)  = 0.0d0
            nnzelmnt(i)  = 0.0d0
        END DO
!$OMP END DO
!$OMP END PARALLEL

        ttlsrcfmm = ttlnmbrelmnt*n_glqtr2d
        ALLOCATE(srcfmm_vec(3,ttlsrcfmm))
        ALLOCATE(srcfmm_nrm(3,ttlsrcfmm))
        ALLOCATE(srcfmm_wght(ttlsrcfmm))
        IF (MeshType == "L") THEN
            ALLOCATE(elmntlnknd(ttlnmbrelmnt,3))
            ALLOCATE(srcfmm_wtnd(3,ttlsrcfmm))
        END IF
        IF (MeshType == "Q") THEN
            ALLOCATE(elmntlnknd(ttlnmbrelmnt,6))
            ALLOCATE(elmntlnkndlnr(4*ttlnmbrelmnt,3))
            elmntlnkndlnr(:,:) = 0
            ALLOCATE(srcfmm_wtnd(6,ttlsrcfmm))
        END IF
        elmntlnknd(:,:) = 0



        ndoffset1 = 0
        elmntoffset1 = 0
        DO ithprtl = 1, nmbrprtl
                
            CALL GetMeshNdElmntIcshdrlPaper &
                        &(Meshnlvl(ithprtl),MeshRlvlStp(ithprtl),MeshType)

            OPEN (101, FILE = "Prtl_Orgnl.inp", STATUS = "OLD", IOSTAT = IOS)
            IF (IOS /= 0) THEN
                PRINT*, '"Prtl_Orgnl.inp" does not exist! Please check!'
                STOP
            END IF

            READ (101, *)
            READ (101, *)
            READ (101, *) nmbrnd(ithprtl)
            READ (101, *) nmbrelmnt(ithprtl)
            READ (101, *) xrefnn, yrefnn, zrefnn

            CLOSE (101, STATUS = "DELETE")

            OPEN (111, FILE = "Prtl_Orgnl.vrt", STATUS = "OLD", IOSTAT = IOS)
            IF (IOS /= 0) THEN
                PRINT*, "'Prtl_Orgnl.vrt' does not exist! Please check!"
                STOP
            END IF

            DO i = 1, nmbrnd(ithprtl)
                nd_i = i + ndoffset1
                READ (111, *) ndindex, xnd(nd_i), ynd(nd_i), znd(nd_i)
            END DO

            CLOSE (111, STATUS = "DELETE")

            OPEN (121, FILE = "Prtl_Orgnl.cel", STATUS = "OLD", IOSTAT = IOS)
            IF (IOS /= 0) THEN
                PRINT*, '"Prtl_Orgnl.cel" does not exist! Please check!'
                STOP
            END IF

            DO k = 1, nmbrelmnt(ithprtl)

                elmnt_k = k + elmntoffset1

                IF (MeshType == "L") THEN

                    READ (121, *) elmntindex, ndA, ndB, ndC

                    ndA = ndA + ndoffset1
                    ndB = ndB + ndoffset1
                    ndC = ndC + ndoffset1

                    tpvctACx = xnd(ndA) - xnd(ndC)
                    tpvctACy = ynd(ndA) - ynd(ndC)
                    tpvctACz = znd(ndA) - znd(ndC)

                    tpvctBCx = xnd(ndB) - xnd(ndC)
                    tpvctBCy = ynd(ndB) - ynd(ndC)
                    tpvctBCz = znd(ndB) - znd(ndC)

                    tpvctRefCx = xrefnn - xnd(ndC)
                    tpvctRefCy = yrefnn - ynd(ndC)
                    tpvctRefCz = zrefnn - znd(ndC)

                    tpnrmlchck =tpvctRefCx * (tpvctACy*tpvctBCz - tpvctACz*tpvctBCy)&
                            & + tpvctRefCy * (tpvctACz*tpvctBCx - tpvctACx*tpvctBCz)&
                            & + tpvctRefCz * (tpvctACx*tpvctBCy - tpvctACy*tpvctBCx)

                    IF (NrmlInOut(ithprtl) == 1) THEN
                        IF (tpnrmlchck > 0.0d0) THEN
                            elmntlnknd(elmnt_k,1) = ndA
                            elmntlnknd(elmnt_k,2) = ndB
                            elmntlnknd(elmnt_k,3) = ndC
                        END IF

                        IF (tpnrmlchck < 0.0d0) THEN
                            elmntlnknd(elmnt_k,1) = ndB
                            elmntlnknd(elmnt_k,2) = ndA
                            elmntlnknd(elmnt_k,3) = ndC
                        END IF
                    END IF

                    IF (NrmlInOut(ithprtl) ==-1) THEN
                        IF (tpnrmlchck > 0.0d0) THEN
                            elmntlnknd(elmnt_k,1) = ndB
                            elmntlnknd(elmnt_k,2) = ndA
                            elmntlnknd(elmnt_k,3) = ndC
                        END IF

                        IF (tpnrmlchck < 0.0d0) THEN
                            elmntlnknd(elmnt_k,1) = ndA
                            elmntlnknd(elmnt_k,2) = ndB
                            elmntlnknd(elmnt_k,3) = ndC
                        END IF
                    END IF

                END IF

                IF (MeshType == "Q") THEN

                    READ (121, *) elmntindex, ndA, ndB, ndC, ndD, ndE, ndF

                    ndA = ndA + ndoffset1
                    ndB = ndB + ndoffset1
                    ndC = ndC + ndoffset1
                    ndD = ndD + ndoffset1
                    ndE = ndE + ndoffset1
                    ndF = ndF + ndoffset1

                    tpvctACx = xnd(ndA) - xnd(ndC)
                    tpvctACy = ynd(ndA) - ynd(ndC)
                    tpvctACz = znd(ndA) - znd(ndC)

                    tpvctBCx = xnd(ndB) - xnd(ndC)
                    tpvctBCy = ynd(ndB) - ynd(ndC)
                    tpvctBCz = znd(ndB) - znd(ndC)

                    tpvctRefCx = xrefnn - xnd(ndC)
                    tpvctRefCy = yrefnn - ynd(ndC)
                    tpvctRefCz = zrefnn - znd(ndC)

                    tpnrmlchck =tpvctRefCx * (tpvctACy*tpvctBCz - tpvctACz*tpvctBCy)&
                            & + tpvctRefCy * (tpvctACz*tpvctBCx - tpvctACx*tpvctBCz)&
                            & + tpvctRefCz * (tpvctACx*tpvctBCy - tpvctACy*tpvctBCx)

                    IF (NrmlInOut(ithprtl) == 1) THEN
                        IF (tpnrmlchck > 0.0d0) THEN
                            elmntlnknd(elmnt_k,1) = ndA
                            elmntlnknd(elmnt_k,2) = ndB
                            elmntlnknd(elmnt_k,3) = ndC
                            elmntlnknd(elmnt_k,4) = ndD
                            elmntlnknd(elmnt_k,5) = ndE
                            elmntlnknd(elmnt_k,6) = ndF
                        END IF

                        IF (tpnrmlchck < 0.0d0) THEN
                            elmntlnknd(elmnt_k,1) = ndB
                            elmntlnknd(elmnt_k,2) = ndA
                            elmntlnknd(elmnt_k,3) = ndC
                            elmntlnknd(elmnt_k,4) = ndD
                            elmntlnknd(elmnt_k,5) = ndF
                            elmntlnknd(elmnt_k,6) = ndE
                        END IF
                    END IF

                    IF (NrmlInOut(ithprtl) ==-1) THEN
                        IF (tpnrmlchck < 0.0d0) THEN
                            elmntlnknd(elmnt_k,1) = ndA
                            elmntlnknd(elmnt_k,2) = ndB
                            elmntlnknd(elmnt_k,3) = ndC
                            elmntlnknd(elmnt_k,4) = ndD
                            elmntlnknd(elmnt_k,5) = ndE
                            elmntlnknd(elmnt_k,6) = ndF
                        END IF

                        IF (tpnrmlchck > 0.0d0) THEN
                            elmntlnknd(elmnt_k,1) = ndB
                            elmntlnknd(elmnt_k,2) = ndA
                            elmntlnknd(elmnt_k,3) = ndC
                            elmntlnknd(elmnt_k,4) = ndD
                            elmntlnknd(elmnt_k,5) = ndF
                            elmntlnknd(elmnt_k,6) = ndE
                        END IF
                    END IF

                    elmntlnkndlnr(4*(elmnt_k-1)+1,1)=elmntlnknd(elmnt_k,1)
                    elmntlnkndlnr(4*(elmnt_k-1)+1,2)=elmntlnknd(elmnt_k,4)
                    elmntlnkndlnr(4*(elmnt_k-1)+1,3)=elmntlnknd(elmnt_k,6)

                    elmntlnkndlnr(4*(elmnt_k-1)+2,1)=elmntlnknd(elmnt_k,4)
                    elmntlnkndlnr(4*(elmnt_k-1)+2,2)=elmntlnknd(elmnt_k,2)
                    elmntlnkndlnr(4*(elmnt_k-1)+2,3)=elmntlnknd(elmnt_k,5)

                    elmntlnkndlnr(4*(elmnt_k-1)+3,1)=elmntlnknd(elmnt_k,5)
                    elmntlnkndlnr(4*(elmnt_k-1)+3,2)=elmntlnknd(elmnt_k,3)
                    elmntlnkndlnr(4*(elmnt_k-1)+3,3)=elmntlnknd(elmnt_k,6)

                    elmntlnkndlnr(4*(elmnt_k-1)+4,1)=elmntlnknd(elmnt_k,4)
                    elmntlnkndlnr(4*(elmnt_k-1)+4,2)=elmntlnknd(elmnt_k,5)
                    elmntlnkndlnr(4*(elmnt_k-1)+4,3)=elmntlnknd(elmnt_k,6)

                END IF

            END DO

            CLOSE (121, STATUS = "DELETE")

            ndoffset1 = ndoffset1 + nmbrnd(ithprtl)
            elmntoffset1 = elmntoffset1 + nmbrelmnt(ithprtl)

        END DO

        ndoffset1 = 0
        elmntoffset1 = 0
        DO ithprtl = 1, nmbrprtl

            ndstaID(ithprtl) = ndoffset1 + 1
            elstaID(ithprtl) = elmntoffset1 + 1

            ndoffset1 = ndoffset1 + nmbrnd(ithprtl)
            elmntoffset1 = elmntoffset1 + nmbrelmnt(ithprtl)

            ndendID(ithprtl) = ndoffset1
            elendID(ithprtl) = elmntoffset1

        END DO


!scale
        ndoffset1 = 0
        DO ithprtl = 1, nmbrprtl

!$OMP PARALLEL PRIVATE (i,nd_i)
!$OMP DO
            DO i = 1, nmbrnd(ithprtl)

                nd_i = i + ndoffset1

                xnd(nd_i) = xnd(nd_i)*sizezoom(ithprtl)
                ynd(nd_i) = ynd(nd_i)*sizezoom(ithprtl)
                znd(nd_i) = znd(nd_i)*sizezoom(ithprtl)

            END DO
!$OMP END DO
!$OMP END PARALLEL

!rotation and shift
!$OMP PARALLEL PRIVATE (i,nd_i,tpx,tpy,tpz,tp1,tp2,tp3)
!$OMP DO
            DO i = 1, nmbrnd(ithprtl)

                nd_i = i + ndoffset1

                tpx = xnd(nd_i)
                tpy = ynd(nd_i)
                tpz = znd(nd_i)

                tp1 = DCOS(anglecal_y(ithprtl))*DCOS(anglecal_z(ithprtl))
                tp2 =-DCOS(anglecal_y(ithprtl))*DSIN(anglecal_z(ithprtl))
                tp3 = DSIN(anglecal_y(ithprtl))
                xnd(nd_i) = tpx*tp1+tpy*tp2+tpz*tp3
                xnd(nd_i) = xnd(nd_i) + xloctn(ithprtl)

                tp1 = DSIN(anglecal_x(ithprtl))*DSIN(anglecal_y(ithprtl))&
                    &*DCOS(anglecal_z(ithprtl)) &
                    &+DCOS(anglecal_x(ithprtl))*DSIN(anglecal_z(ithprtl))
                tp2 =-DSIN(anglecal_x(ithprtl))*DSIN(anglecal_y(ithprtl))&
                    &*DSIN(anglecal_z(ithprtl)) &
                    &+DCOS(anglecal_x(ithprtl))*DCOS(anglecal_z(ithprtl))
                tp3 =-DSIN(anglecal_x(ithprtl))*DCOS(anglecal_y(ithprtl))
                ynd(nd_i) = tpx*tp1+tpy*tp2+tpz*tp3
                ynd(nd_i) = ynd(nd_i) + yloctn(ithprtl)

                tp1 =-DCOS(anglecal_x(ithprtl))*DSIN(anglecal_y(ithprtl))&
                    &*DCOS(anglecal_z(ithprtl)) &
                    &+DSIN(anglecal_x(ithprtl))*DSIN(anglecal_z(ithprtl))
                tp2 = DCOS(anglecal_x(ithprtl))*DSIN(anglecal_y(ithprtl))&
                    &*DSIN(anglecal_z(ithprtl)) &
                    &+DSIN(anglecal_x(ithprtl))*DCOS(anglecal_z(ithprtl))
                tp3 = DCOS(anglecal_x(ithprtl))*DCOS(anglecal_y(ithprtl))
                znd(nd_i) = tpx*tp1+tpy*tp2+tpz*tp3
                znd(nd_i) = znd(nd_i) + zloctn(ithprtl)

            END DO
!$OMP END DO
!$OMP END PARALLEL

            ndoffset1 = ndoffset1 + nmbrnd(ithprtl)

        END DO

        CALL Getndlnkelmnt
        CALL Getndlnknd

        IF (MeshType == 'L') THEN
!           ttlddtnd = mxnmbrndlnknd1stslf
            ttlddtnd = mxnmbrndlnknd2ndslf
        END IF

        IF (MeshType == 'Q') THEN
            ttlddtnd = mxnmbrndlnknd2ndslf
        END IF

        ALLOCATE(d_dt1(ttlnmbrnd,ttlddtnd))
        ALLOCATE(d_dt2(ttlnmbrnd,ttlddtnd))

!$OMP PARALLEL PRIVATE (i,j)
!$OMP DO
        DO i = 1, ttlnmbrnd
            DO j = 1, ttlddtnd
                d_dt1(i,j) = 0.0d0
                d_dt2(i,j) = 0.0d0
            END DO
        END DO
!$OMP END DO
!$OMP END PARALLEL

    END SUBROUTINE MeshEditing

    SUBROUTINE Getndlnkelmnt

        INTEGER :: ithprtl, rec_kmx, rec_k, k, i, j

        mxnmbrndlnkelmnt = 0

        rec_kmx = 0

        !work out the maximun number of the elements connecting with the node
        DO k = 1, ttlnmbrnd

            rec_k = 0

            DO i = 1, ttlnmbrelmnt

                IF (MeshType == "L") THEN
                    DO j = 1, 3
                        IF (k == elmntlnknd(i,j)) THEN
                            rec_k = rec_k+1
                        END IF
                    END DO
                END IF
                IF (MeshType == "Q") THEN
                    DO j = 1, 6
                        IF (k == elmntlnknd(i,j)) THEN
                            rec_k = rec_k+1
                        END IF
                    END DO
                END IF
            END DO

            IF (rec_k>rec_kmx) rec_kmx = rec_k

        END DO

        !locate the elements connecting with host node
        IF (mxnmbrndlnkelmnt < rec_kmx) mxnmbrndlnkelmnt = rec_kmx

        ALLOCATE (ndlnkelmnt(ttlnmbrnd, mxnmbrndlnkelmnt))
        ndlnkelmnt(:,:) = 0

        DO k = 1, ttlnmbrnd

            rec_k = 0

            DO i = 1, ttlnmbrelmnt

                IF (MeshType == "L") THEN
                    DO j = 1, 3
                        IF (k == elmntlnknd(i,j)) THEN
                            rec_k = rec_k+1
                            ndlnkelmnt(k, rec_k) = i
                        END IF
                    END DO
                END IF
                IF (MeshType == "Q") THEN
                    DO j = 1, 6
                        IF (k == elmntlnknd(i,j)) THEN
                            rec_k = rec_k+1
                            ndlnkelmnt(k, rec_k) = i
                        END IF
                    END DO
                END IF
            END DO

        END DO

        IF (MeshType == "Q") THEN

            mxnmbrndlnkelmntlnr = 0

            rec_kmx = 0

            !work out the maximun number of the elements connecting with the node
            DO k = 1, ttlnmbrnd

                rec_k = 0

                DO i = 1, 4*ttlnmbrelmnt

                    DO j = 1, 3
                        IF (k == elmntlnkndlnr(i,j)) THEN
                            rec_k = rec_k+1
                        END IF
                    END DO

                END DO

                IF (rec_k>rec_kmx) rec_kmx = rec_k

            END DO

            !locate the elements connecting with host node
            IF (mxnmbrndlnkelmntlnr < rec_kmx) mxnmbrndlnkelmntlnr = rec_kmx


            ALLOCATE (ndlnkelmntlnr(ttlnmbrnd,mxnmbrndlnkelmntlnr))
            ndlnkelmntlnr(:,:) = 0

            DO k = 1, ttlnmbrnd

                rec_k = 0

                DO i = 1, 4*ttlnmbrelmnt

                    DO j = 1, 3
                        IF (k == elmntlnkndlnr(i,j)) THEN
                            rec_k = rec_k+1
                            ndlnkelmntlnr(k, rec_k) = i
                        END IF
                    END DO

                END DO

            END DO

        END IF

!       OPEN (302, FILE = "particle_NdLnkElmnt.dat", STATUS = "REPLACE")
!        DO i = 1, ttlnmbrnd
!            IF (MeshType == "L") THEN
!               WRITE (302, *) i, ndlnkelmnt(i,:)
!           END IF
!            IF (MeshType == "Q") THEN
!               WRITE (302, *) i, ndlnkelmntlnr(i,:)
!           END IF
!        END DO
!        CLOSE (302)

    END SUBROUTINE Getndlnkelmnt
    
    SUBROUTINE Getndlnknd

        INTEGER :: rec_i, rec_imx, rec_imxthd, k, i, j, rdj, rdmxnmbrndlnknd, ithprtl
        INTEGER, ALLOCATABLE,  DIMENSION (:,:) :: rdndlnknd

        rdmxnmbrndlnknd = 100
        ALLOCATE (rdndlnknd(ttlnmbrnd, rdmxnmbrndlnknd))
        rdndlnknd(:,:) = 0

        mxnmbrndlnknd1st = 0

!Get the connected nodes of the first ring of the host node

        rec_imx = 0
        rec_imxthd = 0
        DO i = 1, ttlnmbrnd

!record all the nodes on the elements connecting with the host node
            rec_i = 0
            IF (MeshType == "L") THEN
                DO k = 1, mxnmbrndlnkelmnt
                    IF (ndlnkelmnt(i,k) /= 0) THEN
                        IF (elmntlnknd(ndlnkelmnt(i,k),1) == i) THEN
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnknd(ndlnkelmnt(i,k),2)
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnknd(ndlnkelmnt(i,k),3)
                        END IF
                        IF (elmntlnknd(ndlnkelmnt(i,k),2) == i) THEN
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnknd(ndlnkelmnt(i,k),3)
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnknd(ndlnkelmnt(i,k),1)
                        END IF
                        IF (elmntlnknd(ndlnkelmnt(i,k),3) == i) THEN
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnknd(ndlnkelmnt(i,k),1)
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnknd(ndlnkelmnt(i,k),2)
                        END IF
                    END IF
                END DO
            END IF

            IF (MeshType == "Q") THEN
                DO k = 1, mxnmbrndlnkelmntlnr
                    IF (ndlnkelmntlnr(i,k) /= 0) THEN
                        IF (elmntlnkndlnr(ndlnkelmntlnr(i,k),1) == i) THEN
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnkndlnr(ndlnkelmntlnr(i,k),2)
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnkndlnr(ndlnkelmntlnr(i,k),3)
                        END IF
                        IF (elmntlnkndlnr(ndlnkelmntlnr(i,k),2) == i) THEN
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnkndlnr(ndlnkelmntlnr(i,k),3)
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnkndlnr(ndlnkelmntlnr(i,k),1)
                        END IF
                        IF (elmntlnkndlnr(ndlnkelmntlnr(i,k),3) == i) THEN
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnkndlnr(ndlnkelmntlnr(i,k),1)
                            rec_i = rec_i + 1
                            rdndlnknd(i,rec_i) &
                        & = elmntlnkndlnr(ndlnkelmntlnr(i,k),2)
                        END IF
                    END IF
                END DO
            END IF

!eliminate the repeated nodes
            DO k = 1, rdmxnmbrndlnknd-1
                DO j = k+1, rdmxnmbrndlnknd
                    IF (rdndlnknd(i,k) == rdndlnknd(i,j)) &
                    & rdndlnknd(i,j) = 0
                END DO
            END DO

!get the maximum number of nodes connecting to one node
            rec_i = 0
            DO k = 1, rdmxnmbrndlnknd
                IF (rdndlnknd(i,k) /= 0) THEN
                    rec_i = rec_i + 1
                END IF
            END DO
            rec_imxthd=MAX(rec_imxthd,rec_i)

        END DO

        rec_imx=MAX(rec_imx,rec_imxthd)

        mxnmbrndlnknd1st = MAX(mxnmbrndlnknd1st,rec_imx)

!copy the nodes connecting with the host node to array NdprtlLnkNd
        ALLOCATE (ndlnknd1st(ttlnmbrnd, mxnmbrndlnknd1st))
        ndlnknd1st(:,:) = 0

        DO i = 1, ttlnmbrnd

            rec_i = 0

            DO k = 1, rdmxnmbrndlnknd
                IF (rdndlnknd(i,k) /= 0) THEN
                    rec_i = rec_i + 1
                    ndlnknd1st(i,rec_i) = rdndlnknd(i,k)
                END IF
            END DO
        END DO

!        OPEN (302, FILE = "particle_NdLnkNd1st.dat", STATUS = "REPLACE")
!        DO i = 1, ttlnmbrnd
!            WRITE (302, *) i, ndlnknd1st(i,:)
!        END DO
!        CLOSE (302)

        mxnmbrndlnknd1stslf = mxnmbrndlnknd1st+1
        ALLOCATE (ndlnknd1stslf(ttlnmbrnd, mxnmbrndlnknd1stslf))
        ndlnknd1stslf(:,:) = 0

        DO i = 1, ttlnmbrnd

            ndlnknd1stslf(i,1) = i
            DO k = 1, mxnmbrndlnknd1st
                ndlnknd1stslf(i,k+1) = ndlnknd1st(i,k)
            END DO
        END DO

!==================
!second ring

        rdndlnknd(:,:) = 0

        mxnmbrndlnknd2nd = 0

        rec_imx = 0
        rec_imxthd = 0
        DO i = 1, ttlnmbrnd

!record the first ring
            rdj = 0
            DO j = 1, mxnmbrndlnknd1st
                IF (ndlnknd1st(i,j)/=0) THEN
                    rdj = rdj + 1
                    rdndlnknd(i,rdj) = ndlnknd1st(i,j)
                END IF
            END DO

!record the connected nodes of the first ring nodes of the host node
            DO j = 1, mxnmbrndlnknd1st
                IF (ndlnknd1st(i,j)/=0) THEN
                    DO k = 1, mxnmbrndlnknd1st
                        IF (ndlnknd1st(ndlnknd1st(i,j),k)/=0) THEN
                            rdj = rdj + 1
                            rdndlnknd(i,rdj) &
                        & = ndlnknd1st(ndlnknd1st(i,j),k)
                        END IF
                    END DO
                END IF
            END DO

!eliminate the host node
            DO k = 1, rdmxnmbrndlnknd
                IF (rdndlnknd(i,k)==i) rdndlnknd(i,k) = 0
            END DO

!eliminate the repeated nodes
            DO k = 1, rdmxnmbrndlnknd-1
                DO j = k+1, rdmxnmbrndlnknd
                    IF (rdndlnknd(i,j)==rdndlnknd(i,k)) &
                    & rdndlnknd(i,j) = 0
                END DO
            END DO

!get the max number of the connected nodes of the first and second rings of the host node
            rec_i = 0
            DO k = 1, rdmxnmbrndlnknd
                IF (rdndlnknd(i,k)/=0) rec_i = rec_i+1
            END DO
            rec_imxthd=MAX(rec_imxthd,rec_i)

        END DO

        rec_imx=MAX(rec_imx,rec_imxthd)

        mxnmbrndlnknd2nd = MAX(mxnmbrndlnknd2nd, rec_imx)


!copy the connected nodes of the first and second rings of the host node
        ALLOCATE (ndlnknd2nd(ttlnmbrnd, mxnmbrndlnknd2nd))
        ndlnknd2nd(:,:) = 0

        DO i = 1, ttlnmbrnd
            rec_i = 0
            DO k = 1, rdmxnmbrndlnknd
                IF (rdndlnknd(i,k)/=0) THEN
                    rec_i = rec_i+1
                    ndlnknd2nd(i,rec_i) = rdndlnknd(i,k)
                END IF
            END DO
        END DO

!       OPEN (302, FILE = "particle_NdLnkNd2cd.dat", STATUS = "REPLACE")
!        DO i = 1, ttlnmbrnd
!            WRITE (302, *) i, ndlnknd2nd(i,:)
!        END DO
!        CLOSE (302)

        mxnmbrndlnknd2ndslf = mxnmbrndlnknd2nd+1
        ALLOCATE (ndlnknd2ndslf(ttlnmbrnd, mxnmbrndlnknd2ndslf))
        ndlnknd2ndslf(:,:) = 0

        DO i = 1, ttlnmbrnd

            ndlnknd2ndslf(i,1) = i
            DO k = 1, mxnmbrndlnknd2nd
                ndlnknd2ndslf(i,k+1) = ndlnknd2nd(i,k)
            END DO
        END DO

        DEALLOCATE (rdndlnknd)

    END SUBROUTINE Getndlnknd

END MODULE

