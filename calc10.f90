!------------------------------------------------------------------------------------------
! ��������ե�����(2016.6.27)
! ��StationaryMac.dat
! x1,x2,̩��,ή®1,ή®2,����,����
! ��StationaryPN.dat
! x1,x2,�ƥ����γ�ҿ�
! ��PorousData(pt,dpdx,porosity).dat
! ���ϸ���,����Ψ,������α�ư�̤���,��¦�����μ���ή��,��¦�����μ���ή��,¿������ɽ����
! ��porous.dat
! ¿����������,Knudsen1,Knudsen2,¿�������γ�ҿ��γ��,�����Ѥ�����¿����������Ψ
! ��ave_dis_fun1.dat
! ��ave_dis_fun2.dat
! x1,x2������®��ʬ�۴ؿ��������֤Ǥ�ʿ���͡�
!------------------------------------------------------------------------------------------
! �ɤ߹���ե����� (2016.6.27)
! ��Macroscopic(Rho1,Tau1,pt).dat
! ��Parameter(CN1,CN2).dat
! ��Flux(Rho1,Tau1,pt).dat
! ��DistributionFunc1(Rho1,Tau1,pt).dat
! ��DistributionFunc2(Rho1,Tau1,pt).dat
! ��KnudsenData(CN1,CN2,pt).dat
! ��DataRaw(Rho1,Tau1,pt).dat
!------------------------------------------------------------------------------------------
! ��������
! 2016.1.22 ������ɽ���Ѥ�׻�������ʬ���ɲ�
! 2016.5.31 ��ư�̤���뼰��Ƴ�� 
! 2016.6.26 �ɤ߼��ե�����̾���ѹ�
! 2016.6.27 ��������
! 2016.6.27 module���̥ե�����˰�ư(module.f90)
!   �¹Ԥ� gfortran -o calc8 module.f90 calc8.f90 -O3 -fopenmp 
! 2016.6.28 PorousData(pt,dpdx,porosity).dat�ˤۤȤ�ɤΥǡ�������Ϥ���褦���ѹ�
!   (grad~.dat��porous.dat�����)
! 2016.6.30 PorousData(pt,dpdx,porosity).dat��ʪ�Τ�Ư���Ϥξ�������
! 2016.7.3 ή��ή�ж����ˤ����������֤Ǥΰ���ʿ���ͤ����(���Ϥ�PorousData().dat)
! 2016.7.15 test3�Ѥ��ɤ߹���Knudsen��񤭴���
! 2016.7.15 ®�٤Υǡ�������(test3��)
! 2016.7.15 ʪ�β��̤�ή��ή��γ�ҿ��η׻�
! 2016.7.27 ���ϥե������iPattern��4����ѹ�
!------------------------------------------------------------------------------------------

program calc
  use commn
  use crf
  implicit none
  integer incx,incy,iincx,iincy
  integer, parameter :: samp_step = (iTotalStep - iFinalStep) / iSampleStep
  real(8), parameter :: P0 = 1.0d0 * Rho0
  real(8), parameter :: P1 = Tau1 * Rho1
  
  real(8) xx1,xx2

  real(8) sta_rho(iCellNumber1,iCellNumber2),sta_v(3,iCellNumber1,iCellNumber2)
  real(8) sta_t(iCellNumber1,iCellNumber2),sta_p(iCellNumber1,iCellNumber2)
  real(8) sta_pp(3,3,iCellNumber1,iCellNumber2)
  real(8) sta_m0(iCellNumber1,iCellNumber2),sta_m1(iCellNumber1,iCellNumber2)
  real(8) sta_m2(iCellNumber1,iCellNumber2),sta_m3(iCellNumber1,iCellNumber2)
  real(8) sta_m4(iCellNumber1,iCellNumber2),sta_m5(iCellNumber1,iCellNumber2)
  real(8) sta_m6(iCellNumber1,iCellNumber2),sta_m7(iCellNumber1,iCellNumber2)
  real(8) sta_m8(iCellNumber1,iCellNumber2),sta_m9(iCellNumber1,iCellNumber2)
  real(8) sta_m10(iCellNumber1,iCellNumber2),sta_m11(iCellNumber1,iCellNumber2)
  real(8) sta_m12(iCellNumber1,iCellNumber2),sta_m13(iCellNumber1,iCellNumber2)

  real(8) sta_Fm,sta_Fp !��¦��������¦�����Ǥμ���ή��
  real(8) sta_F !������α�ư�̤���
  real(8) total_n,total_kp1,total_kp2

  character filename*128

  real(8) por_area !¿������ɽ����
  real(8), parameter :: grad_p = (P1 - P0)/2.0d0/L1 !���ϸ���

  !--- ��ư�̷׻� ---
  real(8) mom_l,mom_r,mom_t,mom_b !�����岼�α�ư��
  real(8) samp_data !���§�׻��ѤΥǡ�����Ǽ�ѿ�
  real(8) p_l,p_r,p_t,p_b
  real(8) rho_r,t_r

  !--- ®��ʬ�۴ؿ��η׻� ---
  integer izeta

  real(8) ave_f1(iCellNumber1,iCellNumber2,-nzeta1:nzeta1)
  real(8) ave_f2(iCellNumber1,iCellNumber2,-nzeta2:nzeta2)
  integer zeta1(-nzeta1:nzeta1),zeta2(-nzeta2:nzeta2)
  real(8) dbzeta1(-nzeta1:nzeta1),dbzeta2(-nzeta2:nzeta2)

  !--- ����ʿ�ѷ׻��� ---
  real(8) average_p_left,average_p_right

  !--- �����Τؤ�ή��ή��γ�ҷ׻��� ---
  integer iNBIn_ts0,iNBOut_ts0,iNBIn_ts1,iNBOut_ts1

  write(filepara,'("Parameter(CN1_",i2,",CN2_",i2,").dat")') &
       iCellNumber1,iCellNumber2

  !--- ep.dat�ե����뤫�����Ψ���ɤ߼�� ---
  open(20,file=filepara)
  read(20,'(i3)') KnudsenPt
  read(20,'(i3)') iParallel
  read(20,'(i4)') iPattern
  read(20,'(f8.4)') porosity
  close(20)

  write(fileinidata,'("KnudsenData(CN1_",i2,",CN2_",i2,",pt",i4,").dat")') &
       iCellNumber1,iCellNumber2,iPattern
  write(fileflu,'("Flux(Rho1_",f4.2,",Tau1_",f4.2,",pt",i4,").dat")') &
       Rho1,Tau1,iPattern
  write(filemac,'("Macroscopic(Rho1_",f4.2,",Tau1_",f4.2,",pt",i4,").dat")') &
       Rho1,Tau1,iPattern
  write(filedis1,'("DistributionFunc1(Rho1_",f4.2,",Tau1_",f4.2,",pt",i4,").dat")') &
       Rho1,Tau1,iPattern
  write(filedis2,'("DistributionFunc2(Rho1_",f4.2,",Tau1_",f4.2,",pt",i4,").dat")') &
       Rho1,Tau1,iPattern
  write(filedataraw,'("DataRaw(Rho1_",f4.2,",Tau1_",f4.2,",pt",i4,").dat")') &
       Rho1,Tau1,iPattern
  ! write(filebound,'("BoundsNumber(Rho1_",f4.2,",Tau1_",f4.2,",pt",i4,").dat")') &
  !      Rho1,Tau1,iPattern

  open(50,file=fileflu)
  open(70,file=filemac,access="stream",form="unformatted")
  open(80,file=filedis1)
  open(85,file=filedis2)

  !--- �ʻ���������Ȥ��Ƴ����Ѥ�׻� ---
  dv(:,:) = dx1 * dx2 * L3

  open(30,file='StationaryMac.dat')
  open(60,file='StationaryPN.dat')

 !--- ʪ�ΤΥѥ����󡤰��ϸ��ۡ�����Ψ�˴ؤ���ǡ���
  write(filename,'("PorousData(pt",i4,",dpdx_",f4.2,",poro_",f4.2").dat")') &
       iPattern,grad_p,porosity
  open(40,file=filename,status='replace')

  por_area = 0.0d0

  open(20,file=fileinidata)
  do incy=1,iCellNumber2 
     do incx=1,iCellNumber1
        read(20,*) iincx,iincy,Knudsen(incx,incy)
     end do
     read(20,*)
  end do
  close(20)
  ! Knudsen(:,:) = Knudsen1

  ! do incy=1,iCellNumber2 
  !    do incx=1,iCellNumber1
  !       if(incy<=iCellNumber2/4 .or. incy>iCellNumber2/4*3) then
  !          Knudsen(incx,incy) = Knudsen2
  !       else
  !          Knudsen(incx,incy) = Knudsen1
  !       end if
  !    end do
  ! end do

  open(20,file=filedataraw)
  read(20,*) delta
  read(20,*) dzeta1
  read(20,*) dzeta2 
  close(20)

  !--- ������ʪ�Τؤ�ή��ή��γ�ҿ� ---
  ! iNBIn_ts0 = 0
  ! iNBOut_ts0 = 0 
  ! iNBIn_ts1 = 0 
  ! iNBOut_ts1 = 0
  ! open(90,file=filebound)
  ! do incx=1,iCellNumber1
  !    read(90,'((f8.4),2(i12))') x1(incx),iNumberBoundsIn1(incx),iNumberBoundsOut1(incx)
  !    if(incx<9 .or. incx>24) then
  !       iNBIn_ts0  = iNBIn_ts0  + iNumberBoundsIn1(incx)
  !       iNBOut_ts0 = iNBOut_ts0 + iNumberBoundsOut1(incx)
  !    else
  !       iNBIn_ts1  = iNBIn_ts1  + iNumberBoundsIn1(incx)
  !       iNBOut_ts1 = iNBOut_ts1 + iNumberBoundsOut1(incx)
  !    end if
  ! end do
  ! write(40,*) "iNBIn_ts0=" ,iNBIn_ts0*delta/(iTotalStep - iFinalStep)/dt/2.0d0
  ! write(40,*) "iNBOut_ts0=",iNBOut_ts0*delta/(iTotalStep - iFinalStep)/dt/2.0d0
  ! write(40,*) "iNBIn_ts0-iNBOut_ts0=",(iNBIn_ts0*delta-iNBOut_ts0*delta)/(iTotalStep - iFinalStep)/dt/2.0d0
  ! write(40,*) "iNBIn_ts1=" ,iNBIn_ts1*delta/(iTotalStep - iFinalStep)/dt/2.0d0
  ! write(40,*) "iNBOut_ts1=",iNBOut_ts1*delta/(iTotalStep - iFinalStep)/dt/2.0d0
  ! write(40,*) "iNBIn_ts1-iNBOut_ts1=",(iNBIn_ts1*delta-iNBOut_ts1*delta)/(iTotalStep - iFinalStep)/dt/2.0d0
  ! write(40,*)
  ! close(90)

  !--- ¿������ɽ���Ѥ�׻� ---
  do incy=1,iCellNumber2
     do incx=1,iCellNumber1
        if(incx<iCellNumber1 .and. Knudsen(incx,incy)/=Knudsen(incx+1,incy)) por_area = por_area + dx2
        if(incy<iCellNumber2 .and. Knudsen(incx,incy)/=Knudsen(incx,incy+1)) por_area = por_area + dx1
        if(incx==1              .and. Knudsen(incx,incy)==Knudsen1) por_area = por_area + dx2
        if(incx==iCellNumber1   .and. Knudsen(incx,incy)==Knudsen1) por_area = por_area + dx2
        if(incy==1              .and. Knudsen(incx,incy)/=Knudsen(incx,iCellNumber2)) por_area = por_area + dx1
     end do
  end do

  close(20)

  sta_rho(:,:) = 0.0d0
  sta_v(:,:,:) = 0.0d0
  sta_t(:,:) = 0.0d0
  sta_p(:,:) = 0.0d0
  sta_pp(:,:,:,:) = 0.0d0

  sta_m0(:,:) = 0.0d0
  sta_m1(:,:) = 0.0d0
  sta_m2(:,:) = 0.0d0
  sta_m3(:,:) = 0.0d0
  sta_m4(:,:) = 0.0d0
  sta_m5(:,:) = 0.0d0
  sta_m6(:,:) = 0.0d0
  sta_m7(:,:) = 0.0d0
  sta_m8(:,:) = 0.0d0
  sta_m9(:,:) = 0.0d0
  sta_m10(:,:) = 0.0d0
  sta_m11(:,:) = 0.0d0
  sta_m12(:,:) = 0.0d0
  sta_m13(:,:) = 0.0d0
  total_n = 0.0d0

  sta_Fp = 0.0d0
  sta_Fm = 0.0d0
  sta_F = 0.0d0

  total_n   = 0.0d0
  total_kp1 = 0.0d0
  total_kp2 = 0.0d0

  mom_l = 0.0d0
  mom_r = 0.0d0
  mom_t = 0.0d0
  mom_b = 0.0d0

  p_l = 0.0d0
  p_r = 0.0d0
  p_t = 0.0d0
  p_b = 0.0d0
  
  rho_r = 0.0d0
  t_r = 0.0d0

  average_p_left = 0.0d0
  average_p_right = 0.0d0

  do izeta = -nzeta1,nzeta1
     ave_f1(:,:,izeta) = 0.0d0
  end do
  do izeta = -nzeta2,nzeta2
     ave_f2(:,:,izeta) = 0.0d0
  end do

!---x��ɸ,y��ɸ������---
  do incy = 1,iCellNumber2
     do incx = 1,iCellNumber1
           x1(incx) = -L1 + dx1*dble(incx) - dx1/2.0d0
           x2(incy) = -L2 + dx2*dble(incy) - dx2/2.0d0
     end do
  end do

!---mac.dat���ɤ߼�ä�ʿ���ͤ����---
  do irep=1,samp_step
     do incy=1,iCellNumber2
        do incx=1,iCellNumber1
           read(70) &
                time,xx1,xx2,m0(incx,incy),m1(incx,incy),m2(incx,incy),m3(incx,incy),m4(incx,incy),&
                m5(incx,incy),m6(incx,incy),m7(incx,incy),m8(incx,incy),m9(incx,incy),&
                m10(incx,incy),m11(incx,incy),m12(incx,incy),m13(incx,incy)

           sta_m0(incx,incy) = sta_m0(incx,incy) + m0(incx,incy) 
           sta_m1(incx,incy) = sta_m1(incx,incy) + m1(incx,incy)
           sta_m2(incx,incy) = sta_m2(incx,incy) + m2(incx,incy)
           sta_m3(incx,incy) = sta_m3(incx,incy) + m3(incx,incy)
           sta_m4(incx,incy) = sta_m4(incx,incy) + m4(incx,incy)
           sta_m5(incx,incy) = sta_m5(incx,incy) + m5(incx,incy)
           sta_m6(incx,incy) = sta_m6(incx,incy) + m6(incx,incy)
           sta_m7(incx,incy) = sta_m7(incx,incy) + m7(incx,incy)
           sta_m8(incx,incy) = sta_m8(incx,incy) + m8(incx,incy)
           sta_m9(incx,incy) = sta_m9(incx,incy) + m9(incx,incy)
           sta_m10(incx,incy) = sta_m10(incx,incy) + m10(incx,incy)
           sta_m11(incx,incy) = sta_m11(incx,incy) + m11(incx,incy)
           sta_m12(incx,incy) = sta_m12(incx,incy) + m12(incx,incy)
           sta_m13(incx,incy) = sta_m13(incx,incy) + m13(incx,incy)

          
           ! mac_q(1,incx,incy) = delta/dv(incx,incy) &
           !      *( m11(incx,incy) - m4(incx,incy)*m1(incx,incy)/m0(incx,incy) &
           !      - 2.0d0*(m5(incx,incy)*m1(incx,incy)&
           !      +m8(incx,incy)*m2(incx,incy)&
           !      +m9(incx,incy)*m3(incx,incy))/m0(incx,incy) &
           !      + m1(incx,incy)*(m1(incx,incy)*m1(incx,incy)&
           !      +m2(incx,incy)*m2(incx,incy)&
           !      +m3(incx,incy)*m3(incx,incy))/m0(incx,incy)/m0(incx,incy) &
           !      + 2.0d0*m1(incx,incy)*(m1(incx,incy)*m1(incx,incy)&
           !      +m2(incx,incy)*m2(incx,incy)&
           !      +m3(incx,incy)*m3(incx,incy))/m0(incx,incy)/m0(incx,incy) &
           !      - m1(incx,incy)*(m1(incx,incy)*m1(incx,incy)&
           !      +m2(incx,incy)*m2(incx,incy)&
           !      +m3(incx,incy)*m3(incx,incy))/m0(incx,incy)/m0(incx,incy)/m0(incx,incy) )
        
           !--- ��γ�ҿ��η׻� ---
           total_n = total_n + m0(incx,incy)

           !--- ¿������γ�ҿ��η׻� ---
           if(Knudsen(incx,incy)<Knudsen1) then
              total_kp1 = total_kp1 + m0(incx,incy)
              !--- ¿��������γ�ҿ� ---
           else 
              total_kp2 = total_kp2 + m0(incx,incy)
           end if

           sta_F = sta_F + m1(incx,incy) * delta / dv(incx,incy)

           !------------------------------------------------------------------------------------
           ! ®��ʬ�۴ؿ��η׻� 
           !------------------------------------------------------------------------------------
           
           if(mod(irep,100)==0) then
              do izeta = -nzeta1,nzeta1

                 read(80,*) time,xx1,xx2,zeta1(izeta),f1(izeta,incx,incy)

                 ave_f1(incx,incy,izeta) = ave_f1(incx,incy,izeta) &
                      + f1(izeta,incx,incy)*delta/dv(incx,incy)/dzeta1
              end do
              read(80,*)
              do izeta = -nzeta2,nzeta2

                 read(85,*) time,xx1,xx2,zeta2(izeta),f2(izeta,incx,incy)

                 ave_f2(incx,incy,izeta) = ave_f2(incx,incy,izeta) &
                      + f2(izeta,incx,incy)*delta/dv(incx,incy)/dzeta2
              end do
              read(85,*)
           end if

        end do !incx
        !read(10)
        !read(70,*)
        read(70)
        
     end do !incy

     if(mod(irep,100)==0) write(*,*) irep,"/",samp_step

     !--- flux(rhi1_,iPattern).dat���ɤ߼�ä�sam_step�ĤΥǡ�����­����碌�� ---
     read(50,'(i8,i8,i8,i8)') flux_inp,flux_inm,flux_outp,flux_outm
     sta_Fp  = sta_Fp + dble(flux_outp) - dble(flux_inp) 
     sta_Fm  = sta_Fm + dble(flux_inm)  - dble(flux_outm) 
  end do !irep

  close(80)
  close(85)

  sta_m0(:,:) = sta_m0(:,:) / dble(samp_step)
  sta_m1(:,:) = sta_m1(:,:) / dble(samp_step)
  sta_m2(:,:) = sta_m2(:,:) / dble(samp_step)
  sta_m3(:,:) = sta_m3(:,:) / dble(samp_step)
  sta_m4(:,:) = sta_m4(:,:) / dble(samp_step)
  sta_m5(:,:) = sta_m5(:,:) / dble(samp_step)
  sta_m6(:,:) = sta_m6(:,:) / dble(samp_step)
  sta_m7(:,:) = sta_m7(:,:) / dble(samp_step)
  sta_m8(:,:) = sta_m8(:,:) / dble(samp_step)
  sta_m9(:,:) = sta_m9(:,:) / dble(samp_step)
  sta_m10(:,:) = sta_m10(:,:) / dble(samp_step)
  sta_m11(:,:) = sta_m11(:,:) / dble(samp_step)
  sta_m12(:,:) = sta_m12(:,:) / dble(samp_step)
  sta_m13(:,:) = sta_m13(:,:) / dble(samp_step)

  sta_rho(:,:) = sta_m0(:,:) * delta / dv(:,:) 

  sta_v(1,:,:) = sta_m1(:,:) / sta_m0(:,:)

  sta_v(2,:,:) = sta_m2(:,:) / sta_m0(:,:) 

  sta_v(3,:,:) = sta_m3(:,:) / sta_m0(:,:) 

  sta_t(:,:) = 2.0d0/3.0d0 &
                * ( sta_m4(:,:)/sta_m0(:,:) &
                - (sta_m1(:,:)*sta_m1(:,:)&
                +sta_m2(:,:)*sta_m2(:,:)&
                +sta_m3(:,:)*sta_m3(:,:))&
                /(sta_m0(:,:)*sta_m0(:,:))) 

  sta_p(:,:) = sta_rho(:,:) * sta_t(:,:)

  sta_pp(1,1,:,:) = 2.0d0 &
       * delta / dv(:,:) &
       * ( sta_m5(:,:) - sta_m1(:,:)*sta_m1(:,:)/sta_m0(:,:) ) 

  sta_pp(1,2,:,:) = 2.0d0 &
       * delta / dv(:,:) &
       * ( sta_m8(:,:) - sta_m1(:,:)*sta_m2(:,:)/sta_m0(:,:) ) 

  sta_pp(1,3,:,:) = 2.0d0 &
       * delta / dv(:,:) &
       * ( sta_m9(:,:) - sta_m1(:,:)*sta_m2(:,:)/sta_m0(:,:) ) 
 
  sta_pp(2,1,:,:) = sta_pp(1,2,:,:)

  sta_pp(2,2,:,:) = 2.0d0 &
       * delta / dv(:,:) &
       * ( sta_m6(:,:) - sta_m2(:,:)*sta_m2(:,:)/sta_m0(:,:) )

  sta_pp(2,3,:,:) = 2.0d0 &
       * delta / dv(:,:) &
       * ( sta_m10(:,:) - sta_m2(:,:)*sta_m3(:,:)/sta_m0(:,:) )

  sta_pp(3,1,:,:) = sta_pp(1,3,:,:)

  sta_pp(3,2,:,:) = sta_pp(2,3,:,:)

  sta_pp(3,3,:,:) = 2.0d0 &
       * delta / dv(:,:) &
       * ( sta_m7(:,:) - sta_m3(:,:)*sta_m3(:,:)/sta_m0(:,:) ) 

  !--- γ�ҿ��γ������
  open(20,file="number.dat")
  write(20,*) delta*total_kp1/dble(samp_step)/100.0d0,&
       delta*total_kp2/dble(samp_step)/100.0d0,&
       delta*total_n/dble(samp_step)/100.0d0,&
       delta*sta_Fp/dble(samp_step)/100.0d0/dt/2.0d0,&
       delta*sta_Fm/dble(samp_step)/100.0d0/dt/2.0d0 
  close(20)
  write(*,*) "porous_num",total_kp1/total_n,"porous_area",1-porosity

  !--- ʪ���ΰ�����������Ǥ���γ�Ҥγ������ ---
  write(40,*) "Knudsen1=",Knudsen1
  write(40,*) "Knudsen2=",Knudsen2
  write(40,*) "PN in the object=",total_kp1/total_n
  write(40,*) "1-Porosity=",1-porosity
  
  sta_F = sta_F / dble(samp_step)
  sta_Fp = delta * sta_Fp / dble(samp_step) / dt 
  sta_Fm = delta * sta_Fm / dble(samp_step) / dt 

  write(40,*) "dp/dx=",grad_p
  write(40,*) "Porosity=",porosity
  write(40,*) "Momentum=",sta_F
  write(40,*) "OutflowFlux(x1=L1)",sta_Fp
  write(40,*) "InflowFlux(x1=-L1)",sta_Fm
  write(40,*) "Surface area of porous",por_area
  
  do incy=1,iCellNumber2
     do incx=1,iCellNumber1
        write(30,'(13(f14.8))') &
             x1(incx),x2(incy),sta_rho(incx,incy),&
             sta_v(1,incx,incy),sta_v(2,incx,incy),&
             sta_t(incx,incy),sta_p(incx,incy),&
             sta_pp(1,1,incx,incy),sta_pp(1,2,incx,incy),&
             sta_pp(1,3,incx,incy),sta_pp(2,2,incx,incy),&
             sta_pp(2,3,incx,incy),sta_pp(3,3,incx,incy)
        write(60,*) x1(incx),x2(incy),sta_m0(incx,incy),sta_m1(incx,incy)
     
        !--- ���§���Ѥ��ƺǤ⳰¦�Υ���ˤ����뱿ư��ή��ή�Ф�׻� ---

        !-- ������α�ư�̤ȱ��Ϸ׻� --
        if(incx==1) then
           !-- ��ư�̤η׻�
           samp_data = -delta*sta_m1(incx,incy)/dv(incx,incy) &
                * sta_v(1,incx,incy)
           if(incy==1 .or. incy==iCellNumber2) then
              mom_l = mom_l + samp_data
           else
              mom_l = mom_l + 2.0d0*samp_data
           end if

           !-- ���Ϥη׻�
           if(incy==1 .or. incy==iCellNumber2) then
              p_l = p_l + sta_pp(1,1,incx,incy)
           else
              p_l = p_l + 2.0d0*sta_pp(1,1,incx,incy)
           end if

           average_p_left = average_p_left + sta_pp(1,1,incx,incy)

        end if

        !-- ������α�ư�̤ȱ��Ϸ׻� --
        if(incx==iCellNumber1) then
           samp_data = delta*sta_m1(incx,incy)/dv(incx,incy) &
                * sta_v(1,incx,incy)
           if(incy==1 .or. incy==iCellNumber2) then
              mom_r = mom_r + samp_data
           else
              mom_r = mom_r + 2.0d0*samp_data
           end if

           !-- ���Ϥη׻�
           if(incy==1 .or. incy==iCellNumber2) then
              p_r = p_r - sta_pp(1,1,incx,incy)
              rho_r = rho_r + sta_rho(incx,incy)
              t_r = t_r + sta_t(incx,incy)
           else
              p_r = p_r - 2.0d0*sta_pp(1,1,incx,incy)
              rho_r = rho_r + 2.0d0*sta_rho(incx,incy)
              t_r = t_r + 2.0d0*sta_t(incx,incy)
           end if

           average_p_right = average_p_right + sta_pp(1,1,incx,incy)

        end if

         !-- �夫��α�ư�̤ȱ��Ϸ׻� --
        if(incy==iCellNumber2) then
           samp_data = delta*sta_m1(incx,incy)/dv(incx,incy) &
                * sta_v(2,incx,incy)
           if(incx==1 .or. incx==iCellNumber1) then
              mom_t = mom_t + samp_data
           else
              mom_t = mom_t + 2.0d0*samp_data
           end if

           !-- ���Ϥη׻�
           if(incx==1 .or. incx==iCellNumber1) then
              p_t = p_t + sta_pp(1,2,incx,incy)
           else
              p_t = p_t + 2.0d0*sta_pp(1,2,incx,incy)
           end if
        end if

        !-- ������α�ư�̤ȱ��Ϸ׻� --
        if(incy==1) then
           samp_data = -delta*sta_m1(incx,incy)/dv(incx,incy) &
                * sta_v(2,incx,incy)
           if(incx==1 .or. incx==iCellNumber1) then
              mom_b = mom_b + samp_data
           else
              mom_b = mom_b + 2.0d0*samp_data
           end if

           !-- ���Ϥη׻�
           if(incx==1 .or. incx==iCellNumber1) then
              p_b = p_b + sta_pp(1,2,incx,incy)
           else
              p_b = p_b + 2.0d0*sta_pp(1,2,incx,incy)
           end if
        end if     

     end do
     write(30,*)
     write(60,*)
 
  end do

  !--- ���§�η׻���³�� --- 
  mom_l = mom_l * dx2 * L3 / 2.0d0
  mom_r = mom_r * dx2 * L3 / 2.0d0
  mom_t = mom_t * dx1 * L3 / 2.0d0
  mom_b = mom_b * dx1 * L3 / 2.0d0

  p_l = p_l * dx2 * L3 / 2.0d0
  p_r = p_r * dx2 * L3 / 2.0d0
  p_t = p_t * dx1 * L3 / 2.0d0
  p_b = p_b * dx1 * L3 / 2.0d0

  rho_r = rho_r * dx2 * L3 / 2.0d0
  t_r = t_r * dx2 * L3 / 2.0d0

  !--- ʿ�Ѱ��Ϥη׻�1 ---
  average_p_left  = average_p_left  / dble(iCellNumber2)
  average_p_right = average_p_right / dble(iCellNumber2)

!------------------------------------------------------------------------------------
! ®��ʬ�۴ؿ������
!------------------------------------------------------------------------------------
  open(80,file="ave_dis_fun1.dat")
  open(85,file="ave_dis_fun2.dat")

  do incy = 1,iCellNumber2
     do incx = 1,iCellNumber1
        do izeta = -nzeta1,nzeta1
           dbzeta1(izeta) = (dble(zeta1(izeta)) + 0.5d0) * dzeta1
           write(80,'(4(f15.4))') &
                x1(incx),x2(incy),dbzeta1(izeta),ave_f1(incx,incy,izeta)/samp_step*100.0d0
        end do
        write(80,*)
        do izeta = -nzeta2,nzeta2
           dbzeta2(izeta) = (dble(zeta2(izeta)) + 0.5d0) * dzeta2
           write(85,'(4(f15.4))') &
                x1(incx),x2(incy),dbzeta2(izeta),ave_f2(incx,incy,izeta)/samp_step*100.0d0
        end do
        write(85,*)
     end do
  end do
  close(80)
  close(85)
!------------------------------------------------------------------------------------
! ʪ�Τ�Ư���Ϥη׻��ǡ�������
!------------------------------------------------------------------------------------
  write(40,*)
  write(40,*) "momentum_l=",mom_l
  write(40,*) "momentum_r=",mom_r
  write(40,*) "momentum_t=",mom_t
  write(40,*) "momentum_b=",mom_b
  write(40,*) "x_momentum=",mom_l + mom_r + mom_t + mom_b
  write(40,*) "p_l=",p_l
  write(40,*) "p_r=",p_r
  write(40,*) "p_l-p_r=",p_l-p_r
  write(40,*) "p_t=",p_t
  write(40,*) "p_b=",p_b
  write(40,*) "x_p=",p_l + p_r + p_t + p_b 
  write(40,*) "x_totalforce=",mom_l + mom_r + mom_t + mom_b - p_l - p_r - p_t - p_b

  write(40,*) "rho_r=",rho_r
  write(40,*) "t_r=",t_r

  write(40,*) "average_p_left =",average_p_left
  write(40,*) "average_p_right=",average_p_right
  write(40,*) "dp=",average_p_right - average_p_left

!------------------------------------------------------------------------------------
! ή®�ǡ����ν���
!------------------------------------------------------------------------------------
  write(40,*)
  write(40,'(4(f12.6))') &
       x1(iCellNumber1/4),x2(iCellNumber2/4),&
       sta_v(1,iCellNumber1/4,iCellNumber2/4),sta_v(2,iCellNumber1/4,iCellNumber2/4)
  write(40,'(4(f12.6))') &
       x1(iCellNumber1/4+1),x2(iCellNumber2/4),&
       sta_v(1,iCellNumber1/4+1,iCellNumber2/4),sta_v(2,iCellNumber1/4+1,iCellNumber2/4)
   write(40,'(4(f12.6))') &
       x1(iCellNumber1/4*3),x2(iCellNumber2/4),&
       sta_v(1,iCellNumber1/4*3,iCellNumber2/4),sta_v(2,iCellNumber1/4*3,iCellNumber2/4)
  write(40,'(4(f12.6))') &
       x1(iCellNumber1/4*3+1),x2(iCellNumber2/4),&
       sta_v(1,iCellNumber1/4*3+1,iCellNumber2/4),sta_v(2,iCellNumber1/4*3+1,iCellNumber2/4)
  write(40,'(4(f12.6))') &
       x1(iCellNumber1/4),x2(iCellNumber2/4*3+1),&
       sta_v(1,iCellNumber1/4,iCellNumber2/4*3+1),sta_v(2,iCellNumber1/4,iCellNumber2/4*3+1)
  write(40,'(4(f12.6))') &
       x1(iCellNumber1/4+1),x2(iCellNumber2/4*3+1),&
       sta_v(1,iCellNumber1/4+1,iCellNumber2/4*3+1),sta_v(2,iCellNumber1/4+1,iCellNumber2/4*3+1)
  write(40,'(4(f12.6))') &
       x1(iCellNumber1/4*3),x2(iCellNumber2/4*3+1),&
       sta_v(1,iCellNumber1/4*3,iCellNumber2/4*3+1),sta_v(2,iCellNumber1/4*3,iCellNumber2/4*3+1)
  write(40,'(4(f12.6))') &
       x1(iCellNumber1/4*3+1),x2(iCellNumber2/4*3+1),&
       sta_v(1,iCellNumber1/4*3+1,iCellNumber2/4*3+1),sta_v(2,iCellNumber1/4*3+1,iCellNumber2/4*3+1)

end program calc
