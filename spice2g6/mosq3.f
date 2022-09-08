      subroutine mosq3(vds,vbs,vgs,vdsat,vth,vbin,gamasd,cox,phi,
     1   qg,qc,qb,cggb,cgdb,cgsb,cbgb,cbdb,cbsb)
      implicit double precision (a-h,o-z)
c
      equivalence (xlamda,alpha),(vbp,theta),(uexp,eta),(utra,xkappa)
c
c     charge equations are referenced to bulk
c
      vgb=vgs-vbs
      vfb=vbi-phi
      onxl=1.0d0/xl
      phibs=sqphbs*sqphbs
c
c     body effect
c
      gammas=gamma*fshort
      fbodys=gammas/(sqphbs+sqphbs)*0.5d0
      fbody=fbodys+fnarrw
      onfbdy=1.0d0/(1.0d0+fbody)
      dfbdvb=-fbodys*dsqdvb/sqphbs+fbodys*dfsdvb/fshort
      qbonco=gammas*sqphbs+fnarrw*phibs
      dqbdvb=gammas*dsqdvb+gamma*dfsdvb*sqphbs-fnarrw
c
c.....static feedback effect
c
      vbix=vbi-eta*vds
c
c.....threshold voltage
c
      vth=vbix+qbonco
      dvtdvd=-eta
      dvtdvb=dqbdvb
c
c     branch according to region of operation
c
      if (vgs.le.vth) go to 800
      vgsx=dmax1(vgs,von)
c
c     branch on vds=0.0d0
c
      vdsx=dmin1(vds,vdsat)
      if ( vdsx.eq.0.0d0 ) go to 900
      cdo=vgsx-vth-0.5d0*(1.0d0+fbody)*vdsx
      dcodvg=1.0d0
      if (vds.lt.vdsat) dcodvd=-dvtdvd-0.5d0*(1.0d0+fbody)
      dcodvb=-dvtdvb-0.5d0*dfbdvb*vdsx
c
c.....charge terms
c
420   continue
      arga=(1.0d0+fbody)*vdsx*vdsx/(12.0d0*cdo)
      dadco=-arga/cdo
      if (vds.lt.vdsat) dadvd=arga/vdsx
      dadfb=arga*onfbdy
c
c.....gate charge
c
      qg=cox*(vgs-vbix-0.5d0*vdsx+arga)
      cggb=cox*(1.0d0+dadco*dcodvg)
      if (vds.lt.vdsat) cgdb=cox*(-dvtdvd-0.5d0+dadvd+dadco*dcodvd)
      cgsb=-cggb-cgdb-cox*(dadco*dcodvb+dadfb*dfbdvb)
c
c.....bulk charge
c
      arga=arga*fbody
      dadco=dadco*fbody
      if (vds.lt.vdsat) dadvd=dadvd*fbody
      dadfb=dadfb*(1.0d0+fbody+fbody)
c
      qb=-cox*(qbonco+0.5d0*fbody*vdsx-arga)
      cbgb=cox*dadco*dcodvg
      if (vds.lt.vdsat) cbdb=-cox*(0.5d0*fbody-dadvd-dadco*dcodvd)
      cbsb=-cbgb-cbdb
     1          +cox*(dqbdvb+(0.5d0*vdsx-dadfb)*dfbdvb-dadco*dcodvb)
      go to 1000
c
c.....charge terms of vgs<vth
c
800   continue
      if ( vgb.gt.vfb ) go to 810
      qg=cox*(vgb-vfb)
      cggb=cox
      go to  820
810   continue
      gamma2=gammas*0.5d0
      arga=dsqrt(gamma2*gamma2+(vgb-vfb))
      qg=gammas*cox*(arga-gamma2)
      cggb=0.5d0*cox*gammas/arga
820   continue
      qb=-qg
      cbgb=-cggb
      cgdb=0.0d0
      cgsb=0.0d0
      cbdb=0.0d0
      cbsb=0.0d0
      go to 1000
c
c     special case vds=0.0d0
c
  900 qg=cox*(vgs-vbi)
      qb=-cox*qbonco
      cggb=cox
      cgdb=-cox*(0.5d0+dvtdvd)
      cgsb=-cox*(0.5d0-dvtdvb)
      cbgb=0.0d0
      cbdb=-0.5d0*cox*fbody
      cbsb=cox*(dqbdvb+0.5d0*fbody)
c
c     done
c
 1000 qc=-(qg+qb)
      return
      end
