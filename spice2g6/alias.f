      double precision function alias(anam)
      implicit double precision (a-h,o-z)
      dimension anam1(14),anam2(14)
      data anam1 /3hva ,3hvb ,3hccs,3hns ,3hc2 ,3hpt ,3hc4 ,
     1            3hpe ,3hme ,3hpc ,3hmc ,3hps ,3hms ,3hik /
      data anam2 /3hvaf,3hvar,3hcjs,3hnss,3hise,3hxti,3hisc,
     1            3hvje,3hmje,3hvjc,3hmjc,3hvjs,3hmjs,3hikf/
c
c  this function returns the mgp equivalent of the gp parameters
c  (those which apply)
c
      iknt=0
      alias=anam
   10 iknt=iknt+1
      if(iknt.gt.14) return
      if(anam1(iknt).ne.anam) go to 10
      alias=anam2(iknt)
      return
      end
