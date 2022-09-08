      subroutine mosq2(vds,vbs,vgs,vdsat,vth,vbin,gamasd,cox,phi,
     1   qg,qc,qb,cggb,cgdb,cgsb,cbgb,cbdb,cbsb)
      implicit double precision (a-h,o-z)
c
c     initialize charges;
c     change reference voltages for charge computation
c
      qg=0.0d0
      qb=0.0d0
      vbd=vbs-vds
      vgb=vgs-vbs
      vd=dmax1(phi-vbd,1.0d-8)
      vs=dmax1(phi-vbs,1.0d-8)
      vg=vgb-vbin+phi
      vsp5=dsqrt(vs)
c
c     determine operating region
c
      if (vgs.le.vth) go to 1100
c
c     compute charges for "on" region
c
 1020 vsat=vdsat+vs
      vs2=vs*vs
      vs3=vs2*vs
      vs5=vs3*vs2
      vs1p5=vs*vsp5
      vs2p5=vs1p5*vs
 1025 if (vd.ge.vsat) go to 1035
      ve=vd
 1030 dvedvd=1.0d0
      dvedvg=0.0d0
      go to 1040
 1035 ve=vsat
      dvedvd=0.0d0
      dvedvg=0.0d0
 1040 ve2=ve*ve
      ve3=ve2*ve
      ve5=ve2*ve3
      vep5=dsqrt(ve)
      ve1p5=ve*vep5
      ve2p5=ve1p5*ve
      term0=ve+vs
      term1=vep5+vsp5
      term2=vep5*vsp5
      term3=ve2+vs2
      term4=ve*vs
      term5=term0*term1
      term6=(term3+term4)+term2*term0
      term7=(term3+term4)*term1
      term10=vep5+0.5d0*vsp5
      term11=1.5d0*ve+vsp5*term10
      term12=2.0d0*ve1p5+vsp5*term11
      term20=0.5d0*vep5+vsp5
      term21=1.5d0*vs+vep5*term20
      term22=2.0d0*vs1p5+vep5*term21
      argn=0.5d0*vg*term5-0.4d0*gamasd*term6-term7/3.0d0
      argd=vg*term1-gamasd*(term0+term2)/1.5d0-0.5d0*term1*term0
      argd2=argd*argd
      qg=cox*(vg-argn/argd)
      dgndve=0.5d0*vg*term11-0.4d0*gamasd*term12-
     1   (2.5d0*ve2+vsp5*term12)/3.0d0
      dddve=0.5d0*vg-gamasd*term10/1.5d0-0.5d0*term11
      dqgdve=-cox/argd*(dgndve-(vg-qg/cox)*dddve)
      dgndvs=0.5d0*vg*term21-0.4d0*gamasd*term22-
     1   (2.5d0*vs2+vep5*term22)/3.0d0
      dddvs=0.5d0*vg-gamasd*term20/1.5d0-0.5d0*term21
      cgdb=-cox/(argd*vep5)*(dgndve-(vg-qg/cox)*dddve)*dvedvd
      cgsb=-cox/(argd*vsp5)*(dgndvs-(vg-qg/cox)*dddvs)
      cggb=cox*(1.0d0-term1/argd*(0.5d0*term0-vg+qg/cox))
      argn=vg*(term0+term2)/1.5d0-0.5d0*gamasd*term5-0.4d0*term6
      dgndve=vg*term10/1.5d0-0.5d0*gamasd*term11-0.4d0*term12
      dgndvs=vg*term20/1.5d0-0.5d0*gamasd*term21-0.4d0*term22
      qb=-gamasd*cox*argn/argd
      cbdb=-cox/(vep5*argd)*(qb/cox*dddve+gamasd*dgndve)*dvedvd
      cbsb=-cox/(vsp5*argd)*(qb/cox*dddvs+gamasd*dgndvs)
      cbgb=-cox/argd*(gamasd*(term0+term2)/1.5d0+qb/cox*term1)
      go to 2000
c
c  finish special cases
c
 1100 if (vg.gt.0.0d0) go to 1110
      qg=cox*vg
      cggb=cox
      go to 1120
 1110 gamma2=gamasd*0.5d0
      sqarg=dsqrt(gamma2*gamma2+vg)
      qg=gamasd*cox*(sqarg-gamma2)
      cggb=0.5d0*cox*gamasd/sqarg
 1120 qb=-qg
      cbgb=-cggb
      cgdb=0.0d0
      cgsb=0.0d0
      cbdb=0.0d0
      cbsb=0.0d0
c
c  finished
c
 2000 qc=-(qg+qb)
 2050 return
      end
