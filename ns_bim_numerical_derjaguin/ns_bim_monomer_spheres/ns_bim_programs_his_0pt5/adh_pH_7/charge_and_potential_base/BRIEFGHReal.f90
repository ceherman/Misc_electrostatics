
MODULE BRIEFGHReal

    USE Pre_Constants

    USE Geom_GlobalData

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE CalGHLnrBRIEFLnrREAL(dmKwn, dmelmnid, &
    &                               dmp0x, dmp0y, dmp0z, dmp0nnx, dmp0nny, dmp0nnz, &
    &                               dmga, dmgb, dmgc, &
    &                               dmha, dmhb, dmhc, &
    &                               dmgnsbim, dmhnsbim )

        DOUBLE PRECISION, INTENT(IN) ::  dmKwn
        INTEGER, INTENT(IN) ::  dmelmnid
        DOUBLE PRECISION, INTENT(IN) ::  dmp0x, dmp0y, dmp0z, dmp0nnx, dmp0nny, dmp0nnz

        DOUBLE PRECISION, INTENT(OUT) :: dmga, dmgb, dmgc

        DOUBLE PRECISION, INTENT(OUT) :: dmha, dmhb, dmhc

        DOUBLE PRECISION, INTENT(OUT) :: dmgnsbim, dmhnsbim

        DOUBLE PRECISION :: tpga, tpgb, tpgc, tpha, tphb, tphc, tpgnsbim, tphnsbim
        DOUBLE PRECISION :: tpGrnK, tpdGKdn, tpGrn0, tpdG0dn
        DOUBLE PRECISION :: tpkwni, tpExpkwnr
        DOUBLE PRECISION :: tpn0xx0, tpn0En

        DOUBLE PRECISION :: tpEnx, tpEny, tpEnz

        DOUBLE PRECISION :: tprx, tpry, tprz
        DOUBLE PRECISION :: tprr0x, tprr0y, tprr0z, mdl_r0r, &
        &                   over_r0r1, over_r0r2, over_r0r3, tpdp

        INTEGER :: GLQi, icnt

        tpkwni = (-1.0d0)*dmKwn

        tpga = 0.0d0
        tpgb = 0.0d0
        tpgc = 0.0d0

        tpha = 0.0d0
        tphb = 0.0d0
        tphc = 0.0d0

        tpgnsbim = 0.0d0
        tphnsbim = 0.0d0

        DO GLQi = 1, n_glqtr2d

            icnt = n_glqtr2d*(dmelmnid-1)+GLQi

            tprx = srcfmm_vec(1, icnt)
            tpry = srcfmm_vec(2, icnt)
            tprz = srcfmm_vec(3, icnt)

            tpEnx = srcfmm_nrm(1, icnt)
            tpEny = srcfmm_nrm(2, icnt)
            tpEnz = srcfmm_nrm(3, icnt)

            tprr0x = tprx - dmp0x
            tprr0y = tpry - dmp0y
            tprr0z = tprz - dmp0z

            mdl_r0r = DSQRT(tprr0x**2 + tprr0y**2 + tprr0z**2)
            over_r0r1 = 1.0d0/mdl_r0r
            over_r0r2 = over_r0r1*over_r0r1
            over_r0r3 = over_r0r2*over_r0r1
            tpdp = -(tpEnx*tprr0x + tpEny*tprr0y + tpEnz*tprr0z)

            tpExpkwnr = DEXP(tpkwni*mdl_r0r)

            tpGrnK = over_r0r1*tpExpkwnr
            tpdGKdn = tpdp*(over_r0r1-tpkwni)*tpExpkwnr*over_r0r2
            tpGrn0 = over_r0r1
            tpdG0dn = tpdp*over_r0r3

            tpga = tpga + srcfmm_wtnd(1,icnt) * tpGrnK
            tpgb = tpgb + srcfmm_wtnd(2,icnt) * tpGrnK
            tpgc = tpgc + srcfmm_wtnd(3,icnt) * tpGrnK

            tpn0xx0 = dmp0nnx*tprr0x+dmp0nny*tprr0y+dmp0nnz*tprr0z
            tpn0En  = dmp0nnx*tpEnx+dmp0nny*tpEny+dmp0nnz*tpEnz

            tpgnsbim = tpgnsbim + srcfmm_wght(icnt) * ( tpn0xx0*tpdG0dn -tpn0En*tpGrn0 )

            tpha = tpha + srcfmm_wtnd(1,icnt) * tpdGKdn
            tphb = tphb + srcfmm_wtnd(2,icnt) * tpdGKdn
            tphc = tphc + srcfmm_wtnd(3,icnt) * tpdGKdn

            tphnsbim = tphnsbim + srcfmm_wght(icnt) * (-tpdG0dn)


        END DO

        dmga = tpga
        dmgb = tpgb
        dmgc = tpgc

        dmha = tpha
        dmhb = tphb
        dmhc = tphc

        dmgnsbim = tpgnsbim
        dmhnsbim = tphnsbim

    END SUBROUTINE

    SUBROUTINE CalGHBMLnrBRIEFLnrREAL(dmKwn, dmelmnid, &
    &                                 dmp0x, dmp0y, dmp0z, dmp0nnx, dmp0nny, dmp0nnz, &
    &                                 dmga, dmgb, dmgc, &
    &                                 dmha, dmhb, dmhc, &
    &                                 dmgnsbim, dmhnsbim, &
    &                                 dmaxga, dmaxgb, dmaxgc, dmaxgnsbim, &
    &                                 dmaxg0a, dmaxg0b, dmaxg0c, &
    &                                 dmaxh0a, dmaxh0b, dmaxh0c, dmaxg0nsbim, dmaxh0nsbim )

        DOUBLE PRECISION, INTENT(IN) ::  dmKwn
        INTEGER, INTENT(IN) ::  dmelmnid
        DOUBLE PRECISION, INTENT(IN) ::  dmp0x, dmp0y, dmp0z, dmp0nnx, dmp0nny, dmp0nnz

        DOUBLE PRECISION, INTENT(OUT) :: dmga, dmgb, dmgc, dmha, dmhb, dmhc, &
        &   dmgnsbim, dmhnsbim, &
        &   dmaxga, dmaxgb, dmaxgc, dmaxgnsbim, &
        &   dmaxg0a, dmaxg0b, dmaxg0c, dmaxh0a, dmaxh0b, dmaxh0c, dmaxg0nsbim, dmaxh0nsbim

        DOUBLE PRECISION :: tpga, tpgb, tpgc, tpha, tphb, tphc, &
        &   tpgnsbim, tphnsbim, &
        &   tpaxga, tpaxgb, tpaxgc, tpaxgnsbim, &
        &   tpaxg0a, tpaxg0b, tpaxg0c, tpaxh0a, tpaxh0b, tpaxh0c, tpaxg0nsbim, tpaxh0nsbim
        DOUBLE PRECISION :: tpGrnK, tpdGKdn, tpGrn0, tpdG0dn, &
        &   tpddG0dndn0, tpddGKdndn0, tpdG0dn0, tpdGKdn0
        DOUBLE PRECISION :: tpkwni, tpExpkwnr
        DOUBLE PRECISION :: ztp1, ztp2, ztp3, ztp4, ztp5, &
        &   ztp_n0x, ztp_n0y, ztp_n0z
        DOUBLE PRECISION :: tp, tp1, tp2, tp3, tp4, tp5, tp_n0x, tp_n0y, tp_n0z
        DOUBLE PRECISION :: tpn0xx0, tpn0En

        DOUBLE PRECISION :: tpEnx, tpEny, tpEnz

        DOUBLE PRECISION :: tprx, tpry, tprz
        DOUBLE PRECISION :: tprr0x, tprr0y, tprr0z, &
        &   mdl_r0r, over_r0r1, over_r0r2, over_r0r3, over_r0r4, over_r0r5, tpdp

        INTEGER :: GLQi, icnt

        tpkwni = (-1.0d0)*dmKwn

        tpga = 0.0d0
        tpgb = 0.0d0
        tpgc = 0.0d0

        tpha = 0.0d0
        tphb = 0.0d0
        tphc = 0.0d0

        tpgnsbim = 0.0d0
        tphnsbim = 0.0d0

        tpaxga = 0.0d0
        tpaxgb = 0.0d0
        tpaxgc = 0.0d0

        tpaxgnsbim = 0.0d0

        tpaxg0a = 0.0d0
        tpaxg0b = 0.0d0
        tpaxg0c = 0.0d0

        tpaxh0a = 0.0d0
        tpaxh0b = 0.0d0
        tpaxh0c = 0.0d0

        tpaxg0nsbim = 0.0d0
        tpaxh0nsbim = 0.0d0

        DO GLQi = 1, n_glqtr2d

            icnt = n_glqtr2d*(dmelmnid-1)+GLQi

            tprx = srcfmm_vec(1, icnt)
            tpry = srcfmm_vec(2, icnt)
            tprz = srcfmm_vec(3, icnt)

            tpEnx = srcfmm_nrm(1, icnt)
            tpEny = srcfmm_nrm(2, icnt)
            tpEnz = srcfmm_nrm(3, icnt)

            tprr0x = tprx - dmp0x
            tprr0y = tpry - dmp0y
            tprr0z = tprz - dmp0z

            mdl_r0r = DSQRT(tprr0x**2 + tprr0y**2 + tprr0z**2)
            over_r0r1 = 1.0d0/mdl_r0r
            over_r0r2 = over_r0r1*over_r0r1
            over_r0r3 = over_r0r2*over_r0r1
            over_r0r4 = over_r0r3*over_r0r1
            over_r0r5 = over_r0r4*over_r0r1

            tpExpkwnr = DEXP(tpkwni*mdl_r0r)

            ztp5 = 3.0d0*over_r0r2
            ztp4 = 3.0d0*(-1.0d0)*over_r0r1
            ztp3 = 1.0d0
            ztp2 = (-1.0d0)*mdl_r0r
            ztp1 = -ztp5+ztp4*dmKwn+ztp3*dmKwn*dmKwn

            ztp_n0x= &
            & tpEnx*(tprr0x*tprr0x*ztp1+ztp3-ztp2*dmKwn) &
            &+tpEny*(tprr0x*tprr0y*ztp1) &
            &+tpEnz*(tprr0x*tprr0z*ztp1)

            ztp_n0y= &
            & tpEnx*(tprr0y*tprr0x*ztp1) &
            &+tpEny*(tprr0y*tprr0y*ztp1+ztp3-ztp2*dmKwn) &
            &+tpEnz*(tprr0y*tprr0z*ztp1)

            ztp_n0z= &
            & tpEnx*(tprr0z*tprr0x*ztp1) &
            &+tpEny*(tprr0z*tprr0y*ztp1) &
            &+tpEnz*(tprr0z*tprr0z*ztp1+ztp3-ztp2*dmKwn)

            tp_n0x= &
            & tpEnx*(-over_r0r2*3.0d0*tprr0x*tprr0x+1.0d0) &
            &+tpEny*(-over_r0r2*3.0d0*tprr0x*tprr0y) &
            &+tpEnz*(-over_r0r2*3.0d0*tprr0x*tprr0z)

            tp_n0y= &
            & tpEnx*(-over_r0r2*3.0d0*tprr0y*tprr0x) &
            &+tpEny*(-over_r0r2*3.0d0*tprr0y*tprr0y+1.0d0) &
            &+tpEnz*(-over_r0r2*3.0d0*tprr0y*tprr0z)

            tp_n0z= &
            & tpEnx*(-over_r0r2*3.0d0*tprr0z*tprr0x) &
            &+tpEny*(-over_r0r2*3.0d0*tprr0z*tprr0y) &
            &+tpEnz*(-over_r0r2*3.0d0*tprr0z*tprr0z+1.0d0)

            tpdp = -(tpEnx*tprr0x + tpEny*tprr0y + tpEnz*tprr0z)
            tpGrnK = over_r0r1*tpExpkwnr
            tpdGKdn = tpdp*(over_r0r1-tpkwni)*tpExpkwnr*over_r0r2
            tpGrn0 = over_r0r1
            tpdG0dn = tpdp*over_r0r3
            tpddGKdndn0 = (dmp0nnx*ztp_n0x+dmp0nny*ztp_n0y+dmp0nnz*ztp_n0z)
            tpddGKdndn0 = tpddGKdndn0 * over_r0r3 * tpExpkwnr
            tp = dmp0nnx*tp_n0x+dmp0nny*tp_n0y+dmp0nnz*tp_n0z
            tpddG0dndn0 = tp
            tpddG0dndn0 = tpddG0dndn0 * over_r0r3
            tpdp = (dmp0nnx*tprr0x + dmp0nny*tprr0y + dmp0nnz*tprr0z)
            tpdGKdn0 = tpdp*(over_r0r1-tpkwni)*tpExpkwnr*over_r0r2
            tpdG0dn0 = tpdp*over_r0r3

            tpn0xx0 = dmp0nnx*tprr0x+dmp0nny*tprr0y+dmp0nnz*tprr0z
            tpn0En  = dmp0nnx*tpEnx+dmp0nny*tpEny+dmp0nnz*tpEnz

            tpga = tpga + srcfmm_wtnd(1,icnt) * tpdGKdn0
            tpgb = tpgb + srcfmm_wtnd(2,icnt) * tpdGKdn0
            tpgc = tpgc + srcfmm_wtnd(3,icnt) * tpdGKdn0

            tpgnsbim = tpgnsbim + srcfmm_wght(icnt) * ( tpdG0dn )

            tpha = tpha + srcfmm_wtnd(1,icnt) * (tpddGKdndn0-tpddG0dndn0)
            tphb = tphb + srcfmm_wtnd(2,icnt) * (tpddGKdndn0-tpddG0dndn0)
            tphc = tphc + srcfmm_wtnd(3,icnt) * (tpddGKdndn0-tpddG0dndn0)

            tpgnsbim = tpgnsbim + srcfmm_wght(icnt) * ( tpdG0dn )
            tphnsbim = tphnsbim + srcfmm_wght(icnt) * &
            &                    (tpn0xx0*tpdG0dn-tpn0En*tpGrn0)*0.5d0*dmKwn*dmKwn

            tpaxga = tpaxga + srcfmm_wtnd(1,icnt) * (-tpdG0dn0)
            tpaxgb = tpaxgb + srcfmm_wtnd(2,icnt) * (-tpdG0dn0)
            tpaxgc = tpaxgc + srcfmm_wtnd(3,icnt) * (-tpdG0dn0)

            tpaxgnsbim = tpaxgnsbim + srcfmm_wght(icnt) * (-tpdG0dn)

            tpaxg0a = tpaxg0a + srcfmm_wtnd(1,icnt) * tpGrn0
            tpaxg0b = tpaxg0b + srcfmm_wtnd(2,icnt) * tpGrn0
            tpaxg0c = tpaxg0c + srcfmm_wtnd(3,icnt) * tpGrn0

            tpaxg0nsbim = tpaxg0nsbim + srcfmm_wght(icnt) &
            &                          *( tpn0xx0*tpdG0dn -tpn0En*tpGrn0 )

            tpaxh0a = tpaxh0a + srcfmm_wtnd(1,icnt) * tpdG0dn
            tpaxh0b = tpaxh0b + srcfmm_wtnd(2,icnt) * tpdG0dn
            tpaxh0c = tpaxh0c + srcfmm_wtnd(3,icnt) * tpdG0dn

            tpaxh0nsbim = tpaxh0nsbim + srcfmm_wght(icnt) * (-tpdG0dn)


        END DO

        dmga = tpga
        dmgb = tpgb
        dmgc = tpgc

        dmha = tpha
        dmhb = tphb
        dmhc = tphc

        dmgnsbim = tpgnsbim
        dmhnsbim = tphnsbim

        dmaxga = tpaxga
        dmaxgb = tpaxgb
        dmaxgc = tpaxgc

        dmaxgnsbim = tpaxgnsbim

        dmaxg0a = tpaxg0a
        dmaxg0b = tpaxg0b
        dmaxg0c = tpaxg0c

        dmaxh0a = tpaxh0a
        dmaxh0b = tpaxh0b
        dmaxh0c = tpaxh0c

        dmaxg0nsbim = tpaxg0nsbim
        dmaxh0nsbim = tpaxh0nsbim


    END SUBROUTINE

    SUBROUTINE CalGHQdrBRIEFLnrREAL(dmKwn, dmelmnid, &
    &                               dmp0x, dmp0y, dmp0z, dmp0nnx, dmp0nny, dmp0nnz,  &
    &                               dmga, dmgb, dmgc, dmgd, dmge, dmgf, &
    &                               dmha, dmhb, dmhc, dmhd, dmhe, dmhf, &
    &                               dmgnsbim, dmhnsbim )

        DOUBLE PRECISION, INTENT(IN) ::  dmKwn
        INTEGER, INTENT(IN) ::  dmelmnid
        DOUBLE PRECISION, INTENT(IN) ::  dmp0x, dmp0y, dmp0z, dmp0nnx, dmp0nny, dmp0nnz

        DOUBLE PRECISION, INTENT(OUT) :: dmga, dmgb, dmgc, dmgd, dmge, dmgf

        DOUBLE PRECISION, INTENT(OUT) :: dmha, dmhb, dmhc, dmhd, dmhe, dmhf

        DOUBLE PRECISION, INTENT(OUT) :: dmgnsbim, dmhnsbim

        DOUBLE PRECISION :: tpga, tpgb, tpgc, tpgd, tpge, tpgf
        DOUBLE PRECISION :: tpha, tphb, tphc, tphd, tphe, tphf
        DOUBLE PRECISION :: tpgnsbim, tphnsbim
        DOUBLE PRECISION :: tpGrnK, tpdGKdn, tpGrn0, tpdG0dn
        DOUBLE PRECISION :: tpkwni, tpExpkwnr
        DOUBLE PRECISION :: tpn0xx0, tpn0En

        DOUBLE PRECISION :: tprx, tpry, tprz
        DOUBLE PRECISION :: tprr0x, tprr0y, tprr0z, mdl_r0r, &
        &                   over_r0r1, over_r0r2, over_r0r3, tpdp

        DOUBLE PRECISION :: drx_deps,drx_dyet,dry_deps,dry_dyet,drz_deps,drz_dyet
        DOUBLE PRECISION :: tpEnx, tpEny, tpEnz, JcbDtmn
        DOUBLE PRECISION :: tpxieta

        INTEGER :: GLQi, icnt

        tpkwni = (-1.0d0)*dmKwn

        tpga = 0.0d0
        tpgb = 0.0d0
        tpgc = 0.0d0
        tpgd = 0.0d0
        tpge = 0.0d0
        tpgf = 0.0d0

        tpha = 0.0d0
        tphb = 0.0d0
        tphc = 0.0d0
        tphd = 0.0d0
        tphe = 0.0d0
        tphf = 0.0d0

        tpgnsbim = 0.0d0
        tphnsbim = 0.0d0

        DO GLQi = 1, n_glqtr2d

            icnt = n_glqtr2d*(dmelmnid-1)+GLQi

            tprx = srcfmm_vec(1, icnt)
            tpry = srcfmm_vec(2, icnt)
            tprz = srcfmm_vec(3, icnt)

            tpEnx = srcfmm_nrm(1, icnt)
            tpEny = srcfmm_nrm(2, icnt)
            tpEnz = srcfmm_nrm(3, icnt)

            tprr0x = tprx - dmp0x
            tprr0y = tpry - dmp0y
            tprr0z = tprz - dmp0z

            mdl_r0r = DSQRT(tprr0x**2 + tprr0y**2 + tprr0z**2)
            over_r0r1 = 1.0d0/mdl_r0r
            over_r0r2 = over_r0r1*over_r0r1
            over_r0r3 = over_r0r1*over_r0r2
            tpdp = -(tpEnx*tprr0x + tpEny*tprr0y + tpEnz*tprr0z)

            tpExpkwnr = DEXP(tpkwni*mdl_r0r)

            tpGrnK = over_r0r1*tpExpkwnr
            tpdGKdn = tpdp*(over_r0r1-tpkwni)*tpExpkwnr*over_r0r2
            tpGrn0 = over_r0r1
            tpdG0dn = tpdp*over_r0r3

            tpn0xx0 = dmp0nnx*tprr0x+dmp0nny*tprr0y+dmp0nnz*tprr0z
            tpn0En  = dmp0nnx*tpEnx+dmp0nny*tpEny+dmp0nnz*tpEnz


            tpga = tpga + srcfmm_wtnd(1,icnt) * tpGrnK
            tpgb = tpgb + srcfmm_wtnd(2,icnt) * tpGrnK
            tpgc = tpgc + srcfmm_wtnd(3,icnt) * tpGrnK
            tpgd = tpgd + srcfmm_wtnd(4,icnt) * tpGrnK
            tpge = tpge + srcfmm_wtnd(5,icnt) * tpGrnK
            tpgf = tpgf + srcfmm_wtnd(6,icnt) * tpGrnK

            tpha = tpha + srcfmm_wtnd(1,icnt) * tpdGKdn
            tphb = tphb + srcfmm_wtnd(2,icnt) * tpdGKdn
            tphc = tphc + srcfmm_wtnd(3,icnt) * tpdGKdn
            tphd = tphd + srcfmm_wtnd(4,icnt) * tpdGKdn
            tphe = tphe + srcfmm_wtnd(5,icnt) * tpdGKdn
            tphf = tphf + srcfmm_wtnd(6,icnt) * tpdGKdn

            tpgnsbim = tpgnsbim + srcfmm_wght(icnt) * ( tpn0xx0*tpdG0dn -tpn0En*tpGrn0 )

            tphnsbim = tphnsbim + srcfmm_wght(icnt) * (-tpdG0dn)


        END DO

        dmga = tpga
        dmgb = tpgb
        dmgc = tpgc
        dmgd = tpgd
        dmge = tpge
        dmgf = tpgf

        dmha = tpha
        dmhb = tphb
        dmhc = tphc
        dmhd = tphd
        dmhe = tphe
        dmhf = tphf

        dmgnsbim = tpgnsbim
        dmhnsbim = tphnsbim

    END SUBROUTINE

    SUBROUTINE CalGHBMQdrBRIEFLnrREAL(dmKwn, dmelmnid, &
    &                                 dmp0x, dmp0y, dmp0z, dmp0nnx, dmp0nny, dmp0nnz, &
    &                                 dmga, dmgb, dmgc, dmgd, dmge, dmgf, &
    &                                 dmha, dmhb, dmhc, dmhd, dmhe, dmhf, &
    &                                 dmgnsbim, dmhnsbim, &
    &                                 dmaxga, dmaxgb, dmaxgc, dmaxgd, dmaxge, dmaxgf, &
    &                                 dmaxgnsbim, &
    &                                 dmaxg0a, dmaxg0b, dmaxg0c, &
    &                                 dmaxg0d, dmaxg0e, dmaxg0f, &
    &                                 dmaxh0a, dmaxh0b, dmaxh0c, &
    &                                 dmaxh0d, dmaxh0e, dmaxh0f, &
    &                                 dmaxg0nsbim, dmaxh0nsbim )


        DOUBLE PRECISION, INTENT(IN) ::  dmKwn
        INTEGER, INTENT(IN) ::  dmelmnid
        DOUBLE PRECISION, INTENT(IN) ::  dmp0x, dmp0y, dmp0z, dmp0nnx, dmp0nny, dmp0nnz

        DOUBLE PRECISION, INTENT(OUT) :: dmga, dmgb, dmgc, dmgd, dmge, dmgf, &
        &   dmha, dmhb, dmhc, dmhd, dmhe, dmhf, &
        &   dmgnsbim, dmhnsbim, &
        &   dmaxga, dmaxgb, dmaxgc, dmaxgd, dmaxge, dmaxgf, &
        &   dmaxgnsbim, &
        &   dmaxg0a, dmaxg0b, dmaxg0c, dmaxg0d, dmaxg0e, dmaxg0f, &
        &   dmaxh0a, dmaxh0b, dmaxh0c, dmaxh0d, dmaxh0e, dmaxh0f, &
        &   dmaxg0nsbim, dmaxh0nsbim

        DOUBLE PRECISION :: tpga, tpgb, tpgc, tpgd, tpge, tpgf, &
        &   tpha, tphb, tphc, tphd, tphe, tphf, tpgnsbim, tphnsbim
        DOUBLE PRECISION :: tpaxga, tpaxgb, tpaxgc, tpaxgd, tpaxge, tpaxgf, &
        &   tpaxgnsbim
        DOUBLE PRECISION :: tpaxg0a, tpaxg0b, tpaxg0c, tpaxg0d, tpaxg0e, tpaxg0f
        DOUBLE PRECISION :: tpaxh0a, tpaxh0b, tpaxh0c, tpaxh0d, tpaxh0e, tpaxh0f,&
        &   tpaxg0nsbim, tpaxh0nsbim
        DOUBLE PRECISION :: tpGrnK, tpdGKdn, tpGrn0, tpdG0dn, &
        &   tpddG0dndn0, tpddGKdndn0, tpdG0dn0, tpdGKdn0
        DOUBLE PRECISION :: tpkwni, tpExpkwnr
        DOUBLE PRECISION :: ztp1, ztp2, ztp3, ztp4, ztp5, &
        &   ztp_n0x, ztp_n0y, ztp_n0z
        DOUBLE PRECISION :: tp, tp1, tp2, tp3, tp4, tp5, tp_n0x, tp_n0y, tp_n0z
        DOUBLE PRECISION :: tpn0xx0, tpn0En

        DOUBLE PRECISION :: dmax, dmay, dmaz, dmbx, dmby, dmbz, dmcx, dmcy, dmcz, &
                            &dmEarea, dmEnx, dmEny, dmEnz

        DOUBLE PRECISION :: tprx, tpry, tprz
        DOUBLE PRECISION :: tprr0x, tprr0y, tprr0z, &
        &   mdl_r0r, over_r0r1, over_r0r2, over_r0r3, over_r0r4, over_r0r5, tpdp

        DOUBLE PRECISION :: drx_deps,drx_dyet,dry_deps,dry_dyet,drz_deps,drz_dyet
        DOUBLE PRECISION :: tpEnx, tpEny, tpEnz, JcbDtmn
        DOUBLE PRECISION :: tpxieta

        INTEGER :: GLQi, icnt

        tpkwni = (-1.0d0)*dmKwn

        tpga = 0.0d0
        tpgb = 0.0d0
        tpgc = 0.0d0
        tpgd = 0.0d0
        tpge = 0.0d0
        tpgf = 0.0d0

        tpha = 0.0d0
        tphb = 0.0d0
        tphc = 0.0d0
        tphd = 0.0d0
        tphe = 0.0d0
        tphf = 0.0d0

        tpgnsbim = 0.0d0
        tphnsbim = 0.0d0

        tpaxga = 0.0d0
        tpaxgb = 0.0d0
        tpaxgc = 0.0d0
        tpaxgd = 0.0d0
        tpaxge = 0.0d0
        tpaxgf = 0.0d0

        tpaxgnsbim = 0.0d0

        tpaxg0a = 0.0d0
        tpaxg0b = 0.0d0
        tpaxg0c = 0.0d0
        tpaxg0d = 0.0d0
        tpaxg0e = 0.0d0
        tpaxg0f = 0.0d0

        tpaxh0a = 0.0d0
        tpaxh0b = 0.0d0
        tpaxh0c = 0.0d0
        tpaxh0d = 0.0d0
        tpaxh0e = 0.0d0
        tpaxh0f = 0.0d0

        tpaxg0nsbim = 0.0d0
        tpaxh0nsbim = 0.0d0

        DO GLQi = 1, n_glqtr2d

            icnt = n_glqtr2d*(dmelmnid-1)+GLQi

            tprx = srcfmm_vec(1, icnt)
            tpry = srcfmm_vec(2, icnt)
            tprz = srcfmm_vec(3, icnt)

            tpEnx = srcfmm_nrm(1, icnt)
            tpEny = srcfmm_nrm(2, icnt)
            tpEnz = srcfmm_nrm(3, icnt)

            tprr0x = tprx - dmp0x
            tprr0y = tpry - dmp0y
            tprr0z = tprz - dmp0z

            mdl_r0r = DSQRT(tprr0x**2 + tprr0y**2 + tprr0z**2)
            over_r0r1 = 1.0d0/mdl_r0r
            over_r0r2 = over_r0r1*over_r0r1
            over_r0r3 = over_r0r2*over_r0r1
            over_r0r4 = over_r0r3*over_r0r1
            over_r0r5 = over_r0r4*over_r0r1

            tpExpkwnr = DEXP(tpkwni*mdl_r0r)

            ztp5 = 3.0d0*over_r0r2
            ztp4 = 3.0d0*(-1.0d0)*over_r0r1
            ztp3 = 1.0d0
            ztp2 = (-1.0d0)*mdl_r0r
            ztp1 = -ztp5+ztp4*dmKwn+ztp3*dmKwn*dmKwn

            ztp_n0x= &
            & tpEnx*(tprr0x*tprr0x*ztp1+ztp3-ztp2*dmKwn) &
            &+tpEny*(tprr0x*tprr0y*ztp1) &
            &+tpEnz*(tprr0x*tprr0z*ztp1)

            ztp_n0y= &
            & tpEnx*(tprr0y*tprr0x*ztp1) &
            &+tpEny*(tprr0y*tprr0y*ztp1+ztp3-ztp2*dmKwn) &
            &+tpEnz*(tprr0y*tprr0z*ztp1)

            ztp_n0z= &
            & tpEnx*(tprr0z*tprr0x*ztp1) &
            &+tpEny*(tprr0z*tprr0y*ztp1) &
            &+tpEnz*(tprr0z*tprr0z*ztp1+ztp3-ztp2*dmKwn)

            tp_n0x= &
            & tpEnx*(-over_r0r2*3.0d0*tprr0x*tprr0x+1.0d0) &
            &+tpEny*(-over_r0r2*3.0d0*tprr0x*tprr0y) &
            &+tpEnz*(-over_r0r2*3.0d0*tprr0x*tprr0z)

            tp_n0y= &
            & tpEnx*(-over_r0r2*3.0d0*tprr0y*tprr0x) &
            &+tpEny*(-over_r0r2*3.0d0*tprr0y*tprr0y+1.0d0) &
            &+tpEnz*(-over_r0r2*3.0d0*tprr0y*tprr0z)

            tp_n0z= &
            & tpEnx*(-over_r0r2*3.0d0*tprr0z*tprr0x) &
            &+tpEny*(-over_r0r2*3.0d0*tprr0z*tprr0y) &
            &+tpEnz*(-over_r0r2*3.0d0*tprr0z*tprr0z+1.0d0)

            tpdp = -(tpEnx*tprr0x + tpEny*tprr0y + tpEnz*tprr0z)
            tpGrnK = over_r0r1*tpExpkwnr
            tpdGKdn = tpdp*(over_r0r1-tpkwni)*tpExpkwnr*over_r0r2
            tpGrn0 = over_r0r1
            tpdG0dn = tpdp*over_r0r3
            tpddGKdndn0 = (dmp0nnx*ztp_n0x+dmp0nny*ztp_n0y+dmp0nnz*ztp_n0z)
            tpddGKdndn0 = tpddGKdndn0 * over_r0r3 * tpExpkwnr
            tp = dmp0nnx*tp_n0x+dmp0nny*tp_n0y+dmp0nnz*tp_n0z
            tpddG0dndn0 = tp
            tpddG0dndn0 = tpddG0dndn0 * over_r0r3
            tpdp = (dmp0nnx*tprr0x + dmp0nny*tprr0y + dmp0nnz*tprr0z)
            tpdGKdn0 = tpdp*(over_r0r1-tpkwni)*tpExpkwnr*over_r0r2
            tpdG0dn0 = tpdp*over_r0r3

            tpn0xx0 = dmp0nnx*tprr0x+dmp0nny*tprr0y+dmp0nnz*tprr0z
            tpn0En  = dmp0nnx*tpEnx+dmp0nny*tpEny+dmp0nnz*tpEnz

            tpga = tpga + srcfmm_wtnd(1,icnt) * tpdGKdn0
            tpgb = tpgb + srcfmm_wtnd(2,icnt) * tpdGKdn0
            tpgc = tpgc + srcfmm_wtnd(3,icnt) * tpdGKdn0
            tpgd = tpgd + srcfmm_wtnd(4,icnt) * tpdGKdn0
            tpge = tpge + srcfmm_wtnd(5,icnt) * tpdGKdn0
            tpgf = tpgf + srcfmm_wtnd(6,icnt) * tpdGKdn0

            tpha = tpha + srcfmm_wtnd(1,icnt) * (tpddGKdndn0-tpddG0dndn0)
            tphb = tphb + srcfmm_wtnd(2,icnt) * (tpddGKdndn0-tpddG0dndn0)
            tphc = tphc + srcfmm_wtnd(3,icnt) * (tpddGKdndn0-tpddG0dndn0)
            tphd = tphd + srcfmm_wtnd(4,icnt) * (tpddGKdndn0-tpddG0dndn0)
            tphe = tphe + srcfmm_wtnd(5,icnt) * (tpddGKdndn0-tpddG0dndn0)
            tphf = tphf + srcfmm_wtnd(6,icnt) * (tpddGKdndn0-tpddG0dndn0)

            tpgnsbim = tpgnsbim + srcfmm_wght(icnt) * ( tpdG0dn )
            tphnsbim = tphnsbim + srcfmm_wght(icnt) * &
            &                    (tpn0xx0*tpdG0dn-tpn0En*tpGrn0)*0.5d0*dmKwn*dmKwn

            tpaxga = tpaxga + srcfmm_wtnd(1,icnt) * (-tpdG0dn0)
            tpaxgb = tpaxgb + srcfmm_wtnd(2,icnt) * (-tpdG0dn0)
            tpaxgc = tpaxgc + srcfmm_wtnd(3,icnt) * (-tpdG0dn0)
            tpaxgd = tpaxgd + srcfmm_wtnd(4,icnt) * (-tpdG0dn0)
            tpaxge = tpaxge + srcfmm_wtnd(5,icnt) * (-tpdG0dn0)
            tpaxgf = tpaxgf + srcfmm_wtnd(6,icnt) * (-tpdG0dn0)

            tpaxgnsbim = tpaxgnsbim + srcfmm_wght(icnt) * (-tpdG0dn)

            tpaxg0a = tpaxg0a + srcfmm_wtnd(1,icnt) * tpGrn0
            tpaxg0b = tpaxg0b + srcfmm_wtnd(2,icnt) * tpGrn0
            tpaxg0c = tpaxg0c + srcfmm_wtnd(3,icnt) * tpGrn0
            tpaxg0d = tpaxg0d + srcfmm_wtnd(4,icnt) * tpGrn0
            tpaxg0e = tpaxg0e + srcfmm_wtnd(5,icnt) * tpGrn0
            tpaxg0f = tpaxg0f + srcfmm_wtnd(6,icnt) * tpGrn0

            tpaxg0nsbim = tpaxg0nsbim + srcfmm_wght(icnt) * &
            &                           ( tpn0xx0*tpdG0dn - tpn0En*tpGrn0 )

            tpaxh0a = tpaxh0a + srcfmm_wtnd(1,icnt) * tpdG0dn
            tpaxh0b = tpaxh0b + srcfmm_wtnd(2,icnt) * tpdG0dn
            tpaxh0c = tpaxh0c + srcfmm_wtnd(3,icnt) * tpdG0dn
            tpaxh0d = tpaxh0d + srcfmm_wtnd(4,icnt) * tpdG0dn
            tpaxh0e = tpaxh0e + srcfmm_wtnd(5,icnt) * tpdG0dn
            tpaxh0f = tpaxh0f + srcfmm_wtnd(6,icnt) * tpdG0dn

            tpaxh0nsbim = tpaxh0nsbim + srcfmm_wght(icnt) * (-tpdG0dn)


        END DO

        dmga = tpga
        dmgb = tpgb
        dmgc = tpgc
        dmgd = tpgd
        dmge = tpge
        dmgf = tpgf

        dmha = tpha
        dmhb = tphb
        dmhc = tphc
        dmhd = tphd
        dmhe = tphe
        dmhf = tphf

        dmgnsbim = tpgnsbim
        dmhnsbim = tphnsbim

        dmaxga = tpaxga
        dmaxgb = tpaxgb
        dmaxgc = tpaxgc
        dmaxgd = tpaxgd
        dmaxge = tpaxge
        dmaxgf = tpaxgf

        dmaxgnsbim = tpaxgnsbim

        dmaxg0a = tpaxg0a
        dmaxg0b = tpaxg0b
        dmaxg0c = tpaxg0c
        dmaxg0d = tpaxg0d
        dmaxg0e = tpaxg0e
        dmaxg0f = tpaxg0f

        dmaxh0a = tpaxh0a
        dmaxh0b = tpaxh0b
        dmaxh0c = tpaxh0c
        dmaxh0d = tpaxh0d
        dmaxh0e = tpaxh0e
        dmaxh0f = tpaxh0f

        dmaxg0nsbim = tpaxg0nsbim
        dmaxh0nsbim = tpaxh0nsbim

    END SUBROUTINE

END MODULE
