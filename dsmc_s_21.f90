!--- �����ɤγ��� -----------------------------------------------------------------------
!
! �����ϴ����Ž̶������岼�ϼ����������Σ������������ΰ�ˤĤ��ƹͤ��롥
! ʬ��Ʊ�Τξ��ͤϤʤ���ʬ�Ҥϻ����ΤȾ��ͤ��ʤ�����־���ư���롥
! ���֤򥻥��ʬ�䤷���ƥ��뤴�Ȥ˻����ΤȤξ������٤��б����륯�̡��������Ϳ���뤳�Ȥ��Ǥ��롥
! 
!
! ���Ƚ�ꤢ��
!
! �ɼ�ǡ���
! ��KnudsenData(CN1,CN2,iPattern).dat���ƥ���Υ��̡���������ɤ߼�롥
! ��Parameter(CN1,CN2).dat���������ʪ�Υѥ����󡤶���Ψ�Υǡ������ɤ߼��
!
! ���ϥǡ���
! ������Ƚ�ꤵ��Ƥ�����ϡ�
! ��Macroscopic(Rho1,Tau1,iPattern).dat������̤���뤿����ǥǡ���
!       Rho1:��¦�δ����Ž̶������ˤ�����̩��
!       Tau1:��¦�δ����Ž̶������ˤ����벹��
!       iPattern:�ΰ�����֤���ʪ�ΤΥѥ�����ʤɤ��¿�����ΤʤΤ����ɤ�ʲ���ʬ�ۤ�Ϳ����ʪ�ΤʤΤ��ʤɡ�
! ��DistributionFunc1(Rho1,Tau1,iPattern).dat��1������®��ʬ�۴ؿ�
! ��DistributionFunc2(Rho1,Tau1,iPattern).dat��2������®��ʬ�۴ؿ�
! ��Flux(Rho1,Tau1,iPattern).dat�������ζ�������ή��ή�Ф���γ�ҤθĿ��򥫥����
! ��KnudsenConf(Rho1,Tau1,iPattern).dat����ɸ��Knudsen���Υǡ���
! ��Data(Rho1,Tau1,iPattern).dat����������ꤷ���ѥ�᡼���䡤�׻��ˤ�ä���������ǡ�����
! ��DataRaw(Rho1,Tau1,iPattern).dat���׻��ˤ�ä���������ǡ���
! ��BoundsNumber(Rho1,Tau1,iPattern).dat ������ʪ�Τ�ή��ή�Ф���γ�Ҥο��������
!
!--- �ѹ����� --------------------------------------------------------------------------
!
! 2016.6.2 ����̤ε�����ѹ� 
! 2016.6.6 ®��ʬ�۴ؿ������(1�����Τ�) 
! 2016.6.13 ®��ʬ�۴ؿ������(2������) 
! 2016.6.14 ���ͳ�Ψ100%��0%�Ѥ˽񤭴���
! 2016.6.27 ���Ϥ���ե�����̾���ѹ�
! 2016.6.27 module���̥ե�����ˤ��롥
!     �¹Ԥ� gfortran -o dsmc_s_18 module.f90 dsmc_s_18.f90 -O3 -fopenmp 
! 2016.6.27 DataRaw(Rho1,Tau1,iPattern).dat�ե�����κ���
! 2016.6.27 �����Τα�¦��¦�ǲ��٤��ѹ����Ƥ�����ʬ��̵������
! 2016.6.28 Data(Rho1,Tau1,iPattern).dat�˶�������̩�١����٤������ͤ����
! 2016.7.11 ���Ƚ�����ʬ��ifʸ����
! 2016.7.11 γ�Ҥν�����֤��濴�����֤��줿������ʪ�Τμ�������֤���ʤ��褦�ˤʤäƤ����Τ�
!           ʪ�Τΰ��ִط��ʤ�γ�Ҥ���������������֤�����褦�ˤ��롥
!           �����ȼ�ä�delta�η׻����ѹ�
! 2016.7.11 sum()��minval()���Ѥ��ƽ����γ�ҿ���Knusen�κǾ��ͤη׻���Ԥ�
! 2016.7.11 ������֤γƥ���ˤ�����γ�ҤθĿ��η׻�����ʬ���ѹ�,���֤ΰܤ��ؤ�����ʬ���ѹ�
! 2016.7.11 ��껻��ݤ������ѹ�
! 2016.7.11 zmax�ν��Ϥ���Ȥ������
! 2016.7.12 dis_fun��dzeta�������ʬ���ѹ���zmax���礭���ʤꤹ���Ƥ��������н��
! 2016.7.15 false sharing���θ����count_rf��ä���
! 2016.7.15 �̾��Knudsen�Υѥ�����(0)�ȡ����ͳ�Ψ0%������ѥ�����(1)�ξ��ʬ��
! 2016.7.15 ���ꤷ���ΰ�ˤ�����ή��ή��γ�Ҥο�������륳���ɤ��ɲ�
!
!---------------------------------------------------------------------------------------

!------------------------------------------------------------------
!  main
!------------------------------------------------------------------
program dsmc
  use commn
  use crf
  implicit none
  integer :: incx,incy,iincx,iincy,j,ip,jran,jp,jr,jl
  real(8) dBoundaryPNLeft,dBoundaryPNRight
  real(8) dBPNSubLeft,dBPNSubRight 
  integer iBoundaryPNLeft,iBoundaryPNRight !������������ή������γ�ҿ�
  integer np_max_re
  real(8) djrandom

  !-- ������ʪ�Τ�ή��ή�Ф���γ�Ҥη׻� --
  real(8) boundx

  real(8), parameter :: Redx1 = 1.0d0 / dx1
  real(8), parameter :: Redx2 = 1.0d0 / dx2

  integer, allocatable :: icellPN_sub(:,:,:)
  real(8), allocatable :: xx(:,:,:,:,:)
  real(8), allocatable :: zzeta(:,:,:,:,:)

  !--- Parameter.dat�ե����뤫���������ʪ�Υѥ����󡤶���Ψ���ɤ߼�� ---
  write(filepara,'("Parameter(CN1_",i2,",CN2_",i2,").dat")') &
       iCellNumber1,iCellNumber2

  open(10,file=filepara)
  read(10,*) KnudsenPt
  read(10,*) iParallel
  read(10,*) iPattern
  read(10,*) porosity
  close(10)

  !write(*,*) iParallel

  allocate (icellPN_sub(iCellNumber1,iCellNumber2,0:iParallel))
  allocate (xx(2,iAssumptionPN,iCellNumber1,iCellNumber2,0:iParallel))
  allocate (zzeta(3,iAssumptionPN,iCellNumber1,iCellNumber2,0:iParallel))

  !--- Knudsen����ƥ���˿���ʬ���� ---
  ! write(fileinidata,'("KnudsenData(CN1_",i2,",CN2_",i2,",pt",i3,").dat")') &
  !      iCellNumber1,iCellNumber2,iPattern

  ! open(10,file=fileinidata) 

  ! do incy=1,iCellNumber2 
  !    do incx=1,iCellNumber1
  !       read(10,*) iincx,iincy,Knudsen(incx,incy)
  !    end do
  !    read(10,*)
  ! end do

  ! close(10)

  !Knudsen(:,:) = Knudsen1

  do incy=1,iCellNumber2 
     do incx=1,iCellNumber1
        if(incy<=iCellNumber2/4 .or. incy>iCellNumber2/4*3) then
           Knudsen(incx,incy) = Knudsen2
        else
           Knudsen(incx,incy) = Knudsen1
        end if
     end do
  end do

  write(fileknu,'("KnudsenConf(Rho1_",f4.2,",Tau1_",f4.2,",pt",i3,").dat")') &
       Rho1,Tau1,iPattern
  write(fileflu,'("Flux(Rho1_",f4.2,",Tau1_",f4.2,",pt",i3,").dat")') &
       Rho1,Tau1,iPattern
  write(filemac,'("Macroscopic(Rho1_",f4.2,",Tau1_",f4.2,",pt",i3,").dat")') &
       Rho1,Tau1,iPattern
  write(filedis1,'("DistributionFunc1(Rho1_",f4.2,",Tau1_",f4.2,",pt",i3,").dat")') &
       Rho1,Tau1,iPattern
  write(filedis2,'("DistributionFunc2(Rho1_",f4.2,",Tau1_",f4.2,",pt",i3,").dat")') &
       Rho1,Tau1,iPattern

  open(18,file=fileknu)
  open(20,file=fileflu)
  open(40,file=filemac,access="stream",form="unformatted",status="replace")
  open(50,file=filedis1)
  open(55,file=filedis2)

!----- ������ -----
  time = 0.0d0
  irep = 0
  ip = 0
  j_ste = iTotalStep
  x(:,:,:,:) = 0.0d0
  xx(:,:,:,:,:) = 0.0d0
  zeta(:,:,:,:) = 0.0d0
  zzeta(:,:,:,:,:) = 0.0d0
  ste_mac_m0(:,:) = 0.0d0
  bste_mac_m0(:,:) = 0.0d0
  count_rf(:) = 0
  count_rff(:) = 0
  total_rf = 0
  KnudsenMin = 10000
  zmax1 = 0
  zmax2 = 0
  nfun = 0

  flux_inp = 0
  flux_inm = 0
  flux_outp = 0
  flux_outm = 0
  flux_outp1 = 0
  flux_outm1 = 0
  flux_outp2 = 0
  flux_outm2 = 0
  flux_outp3 = 0
  flux_outm3 = 0
 
  ste(:,:)  = 0
  ste2(:,:) = iTotalStep

  icellPN (:,:) = 0
  icellPN_sub(:,:,:) = 0

  iNumberBoundsIn1(:) = 0
  iNumberBoundsOut1(:) = 0
  
  !--- �ƥ�������� ---
  dv(:,:) = dx1*dx2*L3

  !--- �ƥ����γ�ҿ������ ---
  icellPN(:,:) = int(dble(iParticleNumber)/dble(iCellNumber1)/dble(iCellNumber2))

  !--- �������γ�ҿ� ---
  np_max_re = sum(icellPN)

  !--- ���ߥ�졼�����γ�Ҥμ��� ---
  delta = L1*2.0d0*L2*2.0d0*L3 * RhoInitial / np_max_re 

  !--- Knudsen�κǾ��ͤ�׻� ---
  KnudsenMin = minval(Knudsen)

  !--- �����Ž̶������Ǳ������������äƤ���γ�ҿ��η׻�      ---
  !--- �¿��������ˤ���Ȥ��˳�ΨdBPNSubRight��dBPNSubLeft���Ѥ�������Ū�������ͤ˹礦�褦�ˤ���  ---
  dBoundaryPNLeft = Rho0 / 2.0d0 * sqrt(1.0d0/pi) * L2*2.0d0*L3 * dt / delta  
  dBPNSubLeft = dBoundaryPNLeft - dble(int(dBoundaryPNLeft))
  dBoundaryPNRight = Rho1 / 2.0d0 * sqrt(Tau1/pi) * L2*2.0d0*L3 * dt / delta  
  dBPNSubRight = dBoundaryPNRight - dble(int(dBoundaryPNRight))

  allocate (x_boundary_right(1:2,1:int(dBoundaryPNRight)+1))
  allocate (zeta_boundary_right(1:3,1:int(dBoundaryPNRight)+1))
  allocate (x_boundary_left(1:2,1:int(dBoundaryPNLeft)+1))
  allocate (zeta_boundary_left(1:3,1:int(dBoundaryPNLeft)+1))

  !--- ������ ---
  do incy = 1,iCellNumber2
     do incx = 1,iCellNumber1

        !-- ���Υ����γ�Ҥ�0�Ǥ���з׻���Ԥ�ʤ� --
        if(icellPN(incx,incy)>0) then
           do j=1,icellPN(incx,incy)
        
              !--- ���֡�®�٤����������� ---
              x(1,j,incx,incy) = -L1 + 2.0d0 * L1 * rf(iSeed,0)
              x(2,j,incx,incy) = -L2 + 2.0d0 * L2 * rf(iSeed,0)
              
              !count_rf(0) = count_rf(0) + 2

              !--- ®��ʬ�ۤ���� ---
              call makevv(1.0d0,dtemporary,zeta(1,j,incx,incy),0)
              call makevv(1.0d0,zeta(2,j,incx,incy),zeta(3,j,incx,incy),0)

              !count_rf(0) = count_rf(0) + 4

              !-- ����ʬ�� --
                 !-- �ΰ��Ф�γ�Ҥ�̵�� --
                 if(x(1,j,incx,incy)<=L1 .and. x(1,j,incx,incy)>=-L1) then

                    if(x(1,j,incx,incy)==L1) then
                       iincx = iCellNumber1
                    else
                       iincx = int((x(1,j,incx,incy)+L1)*Redx1) + 1   
                    end if

                    if(x(2,j,incx,incy)==L2) then
                       iincy = iCellNumber2
                    else 
                       iincy = int((x(2,j,incx,incy)+L2)*Redx2) + 1     
                    end if

                    !-- ����������ˤ�����γ�ҿ� --
                    icellPN_sub(iincx,iincy,ip) = icellPN_sub(iincx,iincy,ip) + 1

                    !-- ������̾�����ѹ�
                    xx(1,icellPN_sub(iincx,iincy,ip),iincx,iincy,ip) = x(1,j,incx,incy)
                    xx(2,icellPN_sub(iincx,iincy,ip),iincx,iincy,ip) = x(2,j,incx,incy)
                    zzeta(:,icellPN_sub(iincx,iincy,ip),iincx,iincy,ip) = zeta(:,j,incx,incy)
                 end if

           end do
        end if

        !--- x,y��ɸ�η׻� ---
        x1(incx) = -L1 + dx1*dble(incx) - dx1/2.0d0
        x2(incy) = -L2 + dx2*dble(incy) - dx2/2.0d0

        write(18,'(f15.8,f15.8,f8.2)') x1(incx),x2(incy),Knudsen(incx,incy)

     end do
     write(18,*)
  end do

  icellPN(:,:) = icellPN_sub(:,:,ip)

  !--- ��γ�Ҥξ��������x�˰ܤ��ؤ���(iAssumptionPN����礭��j���ͤ�0�ˤ�����) ---
  x(:,:,:,:) = 0.0d0
  zeta(:,:,:,:) = 0.0d0

  do incy=1,iCellNumber2
     do incx=1,iCellNumber1        
        if(icellPN(incx,incy)>0) then
           do j=1,icellPN(incx,incy)                   
              x(1,j,incx,incy) = xx(1,j,incx,incy,ip) 
              x(2,j,incx,incy) = xx(2,j,incx,incy,ip) 
              zeta(:,j,incx,incy) = zzeta(:,j,incx,incy,ip)
           end do
        end if !icellPN
     end do !incx
  end do !incy

     !----- �����֤��׻� -----
  do irep = 1,iTotalStep
     time = time + dt

     icellPN_sub(:,:,:) = 0

     min_ste = 10000
     max_ste = 1

     zeta_boundary_left(:,:) = 0.0d0
     zeta_boundary_right(:,:) = 0.0d0
     x_boundary_left(:,:) = 0.0d0
     x_boundary_right(:,:) = 0.0d0

     !--- �����Ž̶������Ǻ����ζ����������äƤ���γ�ҿ�������Ū�˸��������ͤȹ礦�褦�˷׻� ---
     if(dBPNSubRight>rf(iSeed,0)) then
        iBoundaryPNRight = int( dBoundaryPNRight ) + 1
     else
        iBoundaryPNRight = int( dBoundaryPNRight ) 
     end if
     if(dBPNSubLeft>rf(iSeed,0)) then
        iBoundaryPNLeft = int( dBoundaryPNLeft ) + 1
     else
        iBoundaryPNLeft = int( dBoundaryPNLeft ) 
     end if

     !count_rf(0) = count_rf(0) + 2

     !-- �ΰ�����äƤ���γ�Ҥο�������� --
     flux_inp = flux_inp + iBoundaryPNRight
     flux_inm = flux_inm + iBoundaryPNLeft

     !$omp parallel private(j,ip,iincx,iincy,icollisionCPN,jran,dtsub,dcollisionCPN,djrandom) &
     !$omp private(itemporary,dtemporary,dtemporary2,boundx) &
     !$omp num_threads(iParallel+1)
     !$ ip = omp_get_thread_num()
     !$omp do reduction(+:flux_outp1,flux_outm1) &
     !$omp reduction(+:iNumberBoundsIn1,iNumberBoundsOut1,iNumberBoundsIn2,iNumberBoundsOut2) collapse(2)
     do incx = 1,iCellNumber1
        do incy = 1,iCellNumber2
           if(icellPN(incx,incy)>0) then
              do j=1,icellPN(incx,incy)

                 !-- ������ʪ�Τؤ�ή��ή�Ф���γ�Ҥη׻��� --
                 buf_x(1,j,incx,incy) = x(1,j,incx,incy)
                 buf_x(2,j,incx,incy) = x(2,j,incx,incy)

                 !-- γ�Ҥΰ�ή --
                 x(1,j,incx,incy) = x(1,j,incx,incy) + dt * zeta(1,j,incx,incy)
                 x(2,j,incx,incy) = x(2,j,incx,incy) + dt * zeta(2,j,incx,incy)

                 !--- �����Τؤ�ή��ή��γ�Ҥη׻� ------------------------------------------
                 ! incy=iCellNumber2/4��incy=iCellNumber2/4+1�δ֤ˤ���
                 ! �����ˤ�����γ�Ҥ�ή��ή�Ф�Ĵ�٤� 
                 ! dtemporary��γ�Ҥ����ζ������̲ᤷ������Ƚ�Ǥ��뤿��Τ�Ρ���ʤ��̲��
                 ! boundx��objx2_b���̲ᤷ������x��ɸ
                 !-------------------------------------------------------------------------

                 if( irep>j_ste .or. irep>iFinalStep ) then

                    if(incy==iCellNumber2/4 .or. incy==iCellNumber2/4+1) then

                       dtemporary = (objx2_b-x(2,j,incx,incy)) * (objx2_b-buf_x(2,j,incx,incy))

                       if(dtemporary < 0.0d0) then
                          boundx = x(1,j,incx,incy) &
                               + (x(1,j,incx,incy)-buf_x(1,j,incx,incy)) / (x(2,j,incx,incy)-buf_x(2,j,incx,incy)) &
                               * (objx2_b-buf_x(2,j,incx,incy))

                          if(boundx<=L1 .and. boundx>=-L1) then
                             if(boundx==L1) then
                                iincx = iCellNumber1
                             else
                                iincx = int((boundx+L1)/dx1) + 1   
                             end if

                             if(x(2,j,incx,incy)>objx2_b) then
                                iNumberBoundsIn1(iincx) = iNumberBoundsIn1(iincx) + 1
                             else
                                iNumberBoundsOut1(iincx) = iNumberBoundsOut1(iincx) + 1
                             end if
                          end if

                       end if

                    end if

                 end if

                 !-- �岼�μ���������� --
400              continue
                
                 if(x(2,j,incx,incy)>L2) then
                    x(2,j,incx,incy) = x(2,j,incx,incy) - 2.0d0 * L2
                 else if(x(2,j,incx,incy)<-L2) then
                    x(2,j,incx,incy) = x(2,j,incx,incy) + 2.0d0 * L2
                 end if

                 if(x(2,j,incx,incy)>L2 .or. x(2,j,incx,incy)<-L2) go to 400

                 !-- �ΰ��ФƹԤ�γ�Ҥο�������� --
                 if(x(1,j,incx,incy)> L1) flux_outp1 = flux_outp1 + 1
                 if(x(1,j,incx,incy)<-L1) flux_outm1 = flux_outm1 + 1

                 !-- ����ʬ�� --
                 !-- �ΰ��Ф�γ�Ҥ�̵�� --
                 if(x(1,j,incx,incy)<=L1 .and. x(1,j,incx,incy)>=-L1) then

                    if(x(1,j,incx,incy)==L1) then
                       iincx = iCellNumber1
                    else
                       iincx = int((x(1,j,incx,incy)+L1)*Redx1) + 1   
                    end if

                    if(x(2,j,incx,incy)==L2) then
                       iincy = iCellNumber2
                    else 
                       iincy = int((x(2,j,incx,incy)+L2)*Redx2) + 1     
                    end if

                    !-- ����������ˤ�����γ�ҿ� --
                    icellPN_sub(iincx,iincy,ip) = icellPN_sub(iincx,iincy,ip) + 1
                    itemporary =  icellPN_sub(iincx,iincy,ip) 

                    !-- ������̾�����ѹ�
                    xx(1,itemporary,iincx,iincy,ip) = x(1,j,incx,incy)
                    xx(2,itemporary,iincx,iincy,ip) = x(2,j,incx,incy)
                    zzeta(:,itemporary,iincx,iincy,ip) = zeta(:,j,incx,incy)

                 end if

              end do !j
           end if
        end do !incy
     end do !incx
     !$omp end do

     !-- �����Ž̶������ --    
     if(iBoundaryPNLeft>0) then

     !$omp do reduction(+:flux_outp3,flux_outm3)
     
        do jl=1,iBoundaryPNLeft

           !-- ��������®��ʬ�ۤ˽��ä�®�٤���� ---
           call makev(Rho0,zeta_boundary_left(1,jl),ip)
           call makevv(Rho0,zeta_boundary_left(2,jl),zeta_boundary_left(3,jl),ip)

           !count_rf(ip) = count_rf(ip) + 3

           !--- �����������äƤ���γ�Ҥΰ��֡�®�٤���� ---
           dtsub = dt * rf(iSeed,ip)
           dtemporary2 = -L2 + rf(iSeed,ip) * L2 * 2.0d0
           x_boundary_left(1,jl) = -L1 + zeta_boundary_left(1,jl) * dtsub
           x_boundary_left(2,jl) = dtemporary2 + zeta_boundary_left(2,jl) * dtsub

           ! !--- �����Τؤ�ή��ή��γ�Ҥη׻� ------------------------------------------
           ! ! incy=iCellNumber2/4��incy=iCellNumber2/4+1�δ֤ˤ���
           ! ! �����ˤ�����γ�Ҥ�ή��ή�Ф�Ĵ�٤� 
           ! ! dtemporary��γ�Ҥ����ζ������̲ᤷ������Ƚ�Ǥ��뤿��Τ�Ρ���ʤ��̲��
           ! ! boundx��objx2_b���̲ᤷ������x��ɸ
           ! !-------------------------------------------------------------------------
           
           ! if( irep>j_ste .or. irep>iFinalStep ) then

           !    if(incy==iCellNumber2/4 .or. incy==iCellNumber2/4+1) then

           !       dtemporary = (objx2_b-x(2,j,incx,incy)) * (objx2_b-dtemporary2)

           !       if(dtemporary < 0.0d0) then
           !          boundx = x(1,j,incx,incy) &
           !               + (x(1,j,incx,incy)-buf_x(1,j,incx,incy)) / (x(2,j,incx,incy)-buf_x(2,j,incx,incy)) &
           !               * (objx2_b-buf_x(2,j,incx,incy))

           !          if(boundx<=L1 .and. boundx>=-L1) then
           !             if(boundx==L1) then
           !                iincx = iCellNumber1
           !             else
           !                iincx = int((boundx+L1)/dx1) + 1   
           !             end if

           !             if(x(2,j,incx,incy)>objx2_b) then
           !                iNumberBoundsIn1(iincx) = iNumberBoundsIn1(iincx) + 1
           !             else
           !                iNumberBoundsOut1(iincx) = iNumberBoundsOut1(iincx) + 1
           !             end if
           !          end if

           !       end if

           !    end if

           ! end if


           !count_rf(ip) = count_rf(ip) + 2

           !--- �岼�μ���������� ---
600        continue
          
           if(x_boundary_left(2,jl)>L2) then
              x_boundary_left(2,jl) = x_boundary_left(2,jl) - 2.0d0 * L2
           else if(x_boundary_left(2,jl)<-L2) then
              x_boundary_left(2,jl) = x_boundary_left(2,jl) + 2.0d0 * L2
           end if

           if(x_boundary_left(2,jl)>L2 .or. x_boundary_left(2,jl)<-L2) go to 600

           !-- �ΰ��ФƹԤ�γ�Ҥο�������� --
           if(x_boundary_left(1,jl)> L1) flux_outp3 = flux_outp3 + 1
           if(x_boundary_left(1,jl)<-L1) flux_outm3 = flux_outm3 + 1

           !--- γ�Ҥο����������ֹ����� ---
           if(x_boundary_left(1,jl)<=L1 .and. x_boundary_left(1,jl)>=-L1) then
              if(x_boundary_left(1,jl)==L1) then
                 iincx = iCellNumber1
              else 
                 iincx = int((x_boundary_left(1,jl)+L1)*Redx1) + 1
              end if

              if(x_boundary_left(2,jl)==L2) then
                 iincy = iCellNumber2              
              else
                 iincy = int((x_boundary_left(2,jl)+L2)*Redx2) + 1                    
              end if

              !--- ����������Ǥ�γ�ҿ��򥫥���� ---
              icellPN_sub(iincx,iincy,ip) = icellPN_sub(iincx,iincy,ip) + 1
              itemporary = icellPN_sub(iincx,iincy,ip) 
              !--- ������γ�Ҥ�̾������� ---
              xx(1,itemporary,iincx,iincy,ip) = x_boundary_left(1,jl)
              xx(2,itemporary,iincx,iincy,ip) = x_boundary_left(2,jl)
              zzeta(:,itemporary,iincx,iincy,ip) = zeta_boundary_left(:,jl)

           end if

        end do !jl

     !$omp end do

     end if


     if(iBoundaryPNRight>0) then

        !$omp do reduction(+:flux_outp2,flux_outm2)

        do jr=1,iBoundaryPNRight

           !-- ��������®��ʬ�ۤ˽��ä�®�٤���� ---
           call makev(1.0d0,zeta_boundary_right(1,jr),ip)
           call makevv(1.0d0,zeta_boundary_right(2,jr),zeta_boundary_right(3,jr),ip)

           !count_rf(ip) = count_rf(ip) + 3

           !--- �����������äƤ���γ�Ҥΰ��֡�®�٤���� ---
           zeta_boundary_right(1,jr) = -zeta_boundary_right(1,jr)
           dtsub = dt * rf(iSeed,ip)
           x_boundary_right(1,jr) = L1 + zeta_boundary_right(1,jr) * dtsub
           x_boundary_right(2,jr) = -L2 + rf(iSeed,ip) * L2 * 2.0d0 + zeta_boundary_right(2,jr) * dtsub

           !count_rf(ip) = count_rf(ip) + 2

           !--- �岼�μ���������� ---
500        continue
          
           if(x_boundary_right(2,jr)>L2) then
              x_boundary_right(2,jr) = x_boundary_right(2,jr) - 2.0d0 * L2
           else if(x_boundary_right(2,jr)<-L2) then
              x_boundary_right(2,jr) = x_boundary_right(2,jr) + 2.0d0 * L2
           end if

           if(x_boundary_right(2,jr)>L2 .or. x_boundary_right(2,jr)<-L2) go to 500

           !-- �ΰ��ФƹԤ�γ�Ҥο�������� --
           if(x_boundary_right(1,jr)> L1) flux_outp2 = flux_outp2 + 1
           if(x_boundary_right(1,jr)<-L1) flux_outm2 = flux_outm2 + 1

           !--- γ�Ҥο����������ֹ����� ---
           if(x_boundary_right(1,jr)<=L1 .and. x_boundary_right(1,jr)>=-L1) then
              if(x_boundary_right(1,jr)==L1) then
                 iincx = iCellNumber1              
              else 
                 iincx = int((x_boundary_right(1,jr)+L1)*Redx1) + 1                    
              end if

              if(x_boundary_right(2,jr)==L2) then
                 iincy = iCellNumber2              
              else 
                 iincy = int((x_boundary_right(2,jr)+L2)*Redx2) + 1                      
              end if

              !--- ����������Ǥ�γ�ҿ��򥫥���� ---
              icellPN_sub(iincx,iincy,ip) = icellPN_sub(iincx,iincy,ip) + 1
              itemporary = icellPN_sub(iincx,iincy,ip)
              !--- ������γ�Ҥ�̾������� ---
              xx(1,itemporary,iincx,iincy,ip) = x_boundary_right(1,jr)
              xx(2,itemporary,iincx,iincy,ip) = x_boundary_right(2,jr)  
              zzeta(:,itemporary,iincx,iincy,ip) = zeta_boundary_right(:,jr)

           end if

        end do !jr

        !$omp end do
     end if

     !$omp single

     flux_outp = flux_outp + flux_outp1 + flux_outp2 + flux_outp3
     flux_outm = flux_outm + flux_outm1 + flux_outm2 + flux_outm3

     !--- �ƥ����γ�ҿ�������icellPN�˰ܤ��ؤ��� ---
     !���� 1
     do incy=1,iCellNumber2
        do incx=1,iCellNumber1
           icellPN(incx,incy) = 0
           do ip=0,iParallel
              icellPN(incx,incy) = icellPN(incx,incy) + icellPN_sub(incx,incy,ip)
           end do
        end do
     end do

     !--- ��γ�Ҥξ��������x�˰ܤ��ؤ��� ---
     x(:,:,:,:) = 0.0d0
     zeta(:,:,:,:) = 0.0d0

     !���� 2
     do incy=1,iCellNumber2
        do incx=1,iCellNumber1
           jp = 0
           do ip=0,iParallel
              if(icellPN_sub(incx,incy,ip)>0) then
                 do j=1,icellPN_sub(incx,incy,ip)
                    jp = jp + 1
                    x(1,jp,incx,incy) = xx(1,j,incx,incy,ip) 
                    x(2,jp,incx,incy) = xx(2,j,incx,incy,ip) 
                    zeta(:,jp,incx,incy) = zzeta(:,j,incx,incy,ip)
                 end do
              end if !icellPN_sub
           end do !ip
        end do !incx
     end do !incy

     !$omp end single

     !--- ���ͤη׻� ---
     !$ ip = omp_get_thread_num() 
     !$omp do collapse(2)
     
     do incx=1,iCellNumber1
        do incy=1,iCellNumber2

           !--- �̾��Knudsen�Υѥ�����(0)�ȡ����ͳ�Ψ0%������ѥ�����(1)�ξ��ʬ�� ---
           if(KnudsenPt == 0) then

              !--- dt�δ֤˾��ͤ���γ�ҿ��η׻� ---
              dcollisionCPN = dble(icellPN(incx,incy))*(1.0d0-exp(-dt/Knudsen(incx,incy)))

              !--- Knudsen=Knudsen2�ʤ���ͤ��ʤ� ---
              if(knudsen(incx,incy)==Knudsen2) dcollisionCPN = 0.0d0

              dtemporary = dcollisionCPN - dble(int(dcollisionCPN)) 
              !-- ���ͤ���γ�ҿ� --
              if(dtemporary>rf(iSeed,ip)) then
                 icollisionCPN = int(dcollisionCPN) + 1
              else
                 icollisionCPN = int(dcollisionCPN)
              end if

           else if(KnudsenPt == 1) then

              if(knudsen(incx,incy)==Knudsen2) then
                 icollisionCPN = 0
              else
                 icollisionCPN = icellPN(incx,incy)
              end if

           end if

           !--- 1���뤢����icollisionCPN�Ĥ�γ�Ҥ����ͤ��� ---
           if(icollisionCPN>0) then
              do j=1,icollisionCPN
                 !--------------------------------------------------------------------------- 
                 !  j�Ͽ������������,����ή����γ�ҡ������������äƤ���γ�ҡ��ˤʤäƤ��뤿�ᡤ 
                 !  ���Τޤ�j��1����icollisionCPN�ޤǷ׻�����ȡ����ͤ���γ�Ҥξ��֤ˤ�����꤬�����롥
                 !  ��äơ�j��������˥Х餱��褦��������Ѥ��ƾ��ͤ���γ�Ҥ����
                 !--------------------------------------------------------------------------- 
               
                 !--- �̾��Knudsen�Υѥ�����(0)�ȡ����ͳ�Ψ0%������ѥ�����(1)�ξ��ʬ�� ---
                 if(KnudsenPt == 0) then

300                 continue

                    djrandom = (dble(icellPN(incx,incy))+1.0d0) * rf(iSeed,ip)

                    jran = int(djrandom) + 1

                    !count_rf(ip) = count_rf(ip) + 2

                    if(jran==0 .or. jran==icellPN(incx,incy)+1) go to 300

                    call makevv(TauObject0,zeta(1,jran,incx,incy),dtemporary,ip)
                    call makevv(TauObject0,zeta(2,jran,incx,incy),zeta(3,jran,incx,incy),ip)

                 else if(KnudsenPt == 1) then

                    call makevv(TauObject0,zeta(1,j,incx,incy),dtemporary,ip)
                    call makevv(TauObject0,zeta(2,j,incx,incy),zeta(3,j,incx,incy),ip)

                 end if
                
                 !count_rf(ip) = count_rf(ip) + 4

              end do
           end if !icollisionCPN

        end do
     end do
     !$omp end do

!----- ®�٤κ����ͤ���� -----
     !$omp do reduction(max:zmax1,zmax2)
     do incx=1,iCellNumber1
        do incy=1,iCellNumber2
           if(icellPN(incx,incy)>0) then
              do j=1,icellPN(incx,incy)
                 if(dabs(zeta(1,j,incx,incy))>zmax1) zmax1 = dabs(zeta(1,j,incx,incy))
                 if(dabs(zeta(2,j,incx,incy))>zmax2) zmax2 = dabs(zeta(2,j,incx,incy))
              end do
           end if
        end do
     end do
     !$omp end do

     !$omp end parallel 

 !----- irep>iInitialStep �ΤȤ�������̤�׻����� -----
     if(irep>iInitialStep) then
        do incy=1,iCellNumber2
           do incx=1,iCellNumber1

              !--- ����̤�׻� ---
              call macro

              if( irep>j_ste .or. irep>iFinalStep ) then
                 !-- ®��ʬ�۴ؿ��η׻��Ƚ��� --
                 call dis_fun
              end if

              !--- ���ȤʤäƤ���iSampleStep���Ȥ˹Ԥ����� ---
              if(irep>j_ste .or. irep>iFinalStep) then
                 if( mod(irep,iSampleStep)==0) then
                    !-- ����̤ν���
                    call out_put
                 end if
              end if

              !--- ����̤�®��ʬ�۴ؿ��ν���� ---
              if(mod(irep,iSampleStep)==0) then
                 m0(incx,incy) = 0.0d0
                 m1(incx,incy) = 0.0d0
                 m2(incx,incy) = 0.0d0
                 m3(incx,incy) = 0.0d0
                 m4(incx,incy) = 0.0d0
                 m5(incx,incy) = 0.0d0
                 m6(incx,incy) = 0.0d0
                 m7(incx,incy) = 0.0d0
                 m8(incx,incy) = 0.0d0
                 m9(incx,incy) = 0.0d0
                 m10(incx,incy) = 0.0d0
                 m11(incx,incy) = 0.0d0
                 m12(incx,incy) = 0.0d0
                 m13(incx,incy) = 0.0d0
              end if
     
           end do   !incx
           if(irep>j_ste .or. irep>iFinalStep) then
              if( mod(irep,iSampleStep)==0) then
                
                 write(40)
                 
              end if
           end if
        end do   !incy

        if(min_ste==3 .and. max_ste/=iTotalStep) j_ste = irep

        !--- flux�ν��� ---
        if(irep>j_ste .or. irep>iFinalStep) then
           if( mod(irep,iSampleStep)==0) then
              write(20,'(i8,i8,i8,i8)') flux_inp,flux_inm,flux_outp,flux_outm    
           end if
        end if

     end if

     !--- �����ߥʥ�˽��� ---
     if(irep<=iInitialStep .and. mod(irep,1000)==0) then
        write(*,'("step:",i8,"/",i8,", time:,",f8.2,"s",", zmax1:", f8.3," m/s^2")') &
             irep,iTotalStep,time,zmax1
     else if (irep>iInitialStep .and. mod(irep,1000)==0) then
         write(*,'("step:",i8,"/",i8,", time:,",f8.2,"s",", m0(10,10):", f8.3,", m4(15,15):", f8.3,", steady:",i8)') &
             irep,iTotalStep,time,ste_jud_m0(1,1),ste_jud_m4(2,1),min_ste
     end if

     !---�����ȯ������η׻�---
     ! !���� 3
     ! do ip=0,iParallel
     !    if(count_rf(ip)>1000000000) then
     !       count_rf(ip) = 0
     !       count_rff(ip) = count_rff(ip) + 1
     !    end if
     ! end do

     !--- ���Ƚ�ꤵ��Ƥ���100000���ƥå׷в��˷׻���λ ---
     if(irep==j_ste+200000) exit

     !--- γ�Ҥ�x����®�ٺ�����zmax1�ν����(iSampleStep���ƥåפ��Ȥκ�����) ---
     if(mod(irep,iSampleStep)==0) then
        zmax1 = 0
        zmax2 = 0
     end if

     !--- flux������ ---
     if(mod(irep,iSampleStep)==0) then
        flux_inp   = 0
        flux_inm   = 0
        flux_outp  = 0
        flux_outm  = 0
     end if
     flux_outp1 = 0 
     flux_outm1 = 0
     flux_outp2 = 0
     flux_outm2 = 0
     flux_outp3 = 0
     flux_outm3 = 0

  end do

  !--- �ǡ����ν��� ---
  write(filedata,'("Data(Rho1_",f4.2,",Tau1_",f4.2,",pt",i3,").dat")') &
       Rho1,Tau1,iPattern
  open(15,file=filedata)
  write(15,*) "Rho1=",Rho1
  write(15,*) "Tau1=",Tau1
  write(15,*) "iParticleNumber=",iParticleNumber
  write(15,*) "iCellNumber1=",iCellNumber1,"iCellNumber2=",iCellNumber2
  write(15,*) "iSeed=",iSeed
  write(15,*) "iTotalStep=",iTotalStep
  write(15,*) "iFinalStep=",iFinalStep
  write(15,*) "dt=",dt
  write(15,*) "delta=",delta
  write(15,*) "final_step=",irep
  write(15,*) "dBoundaryPNLeft=",dBoundaryPNLeft,"dBoundaryPNRight=",dBoundaryPNRight
  write(15,*) "dzeta1=",dzeta1,"dzeta2=",dzeta2 !®�٤ι��������
  write(15,*) "zmax1=",dzeta1*dble(nzeta1),"zmax2=",dzeta2*dble(nzeta2) !®�٤κ����ͽ���
    !���� 4 �Ǹ�
  ! do ip=0,iParallel
  !   write(15,*) "count_rf",ip,count_rff(ip),count_rf(ip)
  ! end do
  ! close(15)
  
  !--- �ǡ����ν��� ---
  write(filedataraw,'("DataRaw(Rho1_",f4.2,",Tau1_",f4.2,",pt",i3,").dat")') &
       Rho1,Tau1,iPattern
  open(15,file=filedataraw)
  write(15,*) delta
  write(15,*) dzeta1
  write(15,*) dzeta2 !®�٤ι��������
  close(15)

  !--- ������ʪ�Τˤ�����ή��ή��γ�ҿ��η׻� ---
  write(filebound,'("BoundsNumber(Rho1_",f4.2,",Tau1_",f4.2,",pt",i3,").dat")') &
       Rho1,Tau1,iPattern
  open(15,file=filebound)
  do incx=1,iCellNumber1
     write(15,'((f8.4),2(i12))') x1(incx),iNumberBoundsIn1(incx),iNumberBoundsOut1(incx)
  end do
  close(15)
  
  write(*,*) "program finish"

contains
  !------------------------------------------------------------------
  !   function rf
  !------------------------------------------------------------------

  real(8) function rf(iidum,ip)
    ! use iPara
    implicit none

    integer, intent(in) :: iidum,ip
    integer idum
    integer, save       :: ma(55,0:31),inext(0:31),inextp(0:31),iff(0:31)=0
    integer, parameter  :: mbig=1000000000,mseed=161803398,mz=0
    real(8), parameter  :: fac=1.d-9
    integer             :: mj,mk,i,ii,k

    idum = iidum + ip

    if((idum<0) .or. (iff(ip)==0)) then
       iff(ip) = 1
       mj = mseed-iabs(idum)
       mj = mod(mj,mbig)
       ma(55,ip) = mj
       mk = 1

       do i=1,54
          ii = mod(21*i,55)
          ma(ii,ip) = mk
          mk = mj - mk
          if(mk<mz) mk=mk+mbig
          mj = ma(ii,ip)
       end do

       do k=1,4
          do i=1,55
             ma(i,ip) = ma(i,ip) - ma(1+mod(i+30,55),ip)
             if(ma(i,ip)<mz) ma(i,ip) = ma(i,ip) + mbig
          end do
       end do

       inext(ip) = 0 
       inextp(ip) = 31
    end if

200 inext(ip)=inext(ip)+1
    if(inext(ip)==56) inext(ip)=1
    inextp(ip) = inextp(ip) + 1
    if(inextp(ip)==56) inextp(ip)=1
    mj = ma(inext(ip),ip) - ma(inextp(ip),ip)
    if(mj<mz) mj=mj+mbig
    ma(inext(ip),ip) = mj
    rf = mj * fac
    if((1.0d-8<rf) .and. (rf<0.99999999)) return
    go to 200

  end function rf

!------------------------------------------------------------------
!   subroutine makev f(x)=(2/a)*x*exp(-(x*x/a)) x>0
!------------------------------------------------------------------

  subroutine makev(a,v1,ip)
    implicit none
    real(8), intent(in)  :: a
    integer, intent(in)  :: ip 
    real(8), intent(out) :: v1
    real(4)              :: r1

    r1 = rf(iSeed,ip)
 
    v1 = sqrt(-a*log(r1)) 
    
    return
  end subroutine makev

!------------------------------------------------------------------
!   subroutine makevv f(x)={ 1/[(pi*a)^(1/2)] }*exp(-x*x/a)
!------------------------------------------------------------------
 
  subroutine makevv(a,v2,v3,ip)
    implicit none
    real(8), intent(in)  :: a 
    integer, intent(in)  :: ip
    real(8), intent(out) :: v2, v3
    real(4)              :: r2, r3

    r2 = rf(iSeed,ip)
    r3 = rf(iSeed,ip)

    v2 = sqrt(-a*log(r2)) * cos(2.0d0*pi*r3)
    v3 = sqrt(-a*log(r2)) * sin(2.0d0*pi*r3)
    
    return
  end subroutine makevv

!------------------------------------------------------------------
!   calc macroscopic quantities
!------------------------------------------------------------------
subroutine macro
  use commn
  implicit none
  integer i,j
  real(8) sqm4

  if(irep==iInitialStep+1) then
     m0(:,:) = 0.0d0
     m1(:,:) = 0.0d0
     m2(:,:) = 0.0d0
     m3(:,:) = 0.0d0
     m4(:,:) = 0.0d0
     m5(:,:) = 0.0d0
     m6(:,:) = 0.0d0
     m7(:,:) = 0.0d0
     m8(:,:) = 0.0d0
     m9(:,:) = 0.0d0
     m10(:,:) = 0.0d0
     m11(:,:) = 0.0d0
     m12(:,:) = 0.0d0
     m13(:,:) = 0.0d0
  end if

  m0(incx,incy) = m0(incx,incy) + icellPN(incx,incy)
  
  if(icellPN(incx,incy)>0) then
     do j=1,icellPN(incx,incy)
        sqm4 = zeta(1,j,incx,incy)*zeta(1,j,incx,incy) + zeta(2,j,incx,incy)*zeta(2,j,incx,incy) &
             + zeta(3,j,incx,incy)*zeta(3,j,incx,incy)

        m1(incx,incy) = m1(incx,incy) + zeta(1,j,incx,incy)
        m2(incx,incy) = m2(incx,incy) + zeta(2,j,incx,incy)
        m3(incx,incy) = m3(incx,incy) + zeta(3,j,incx,incy)

        m4(incx,incy) = m4(incx,incy) + sqm4

        m5(incx,incy) = m5(incx,incy) + zeta(1,j,incx,incy)*zeta(1,j,incx,incy)
        m6(incx,incy) = m6(incx,incy) + zeta(2,j,incx,incy)*zeta(2,j,incx,incy)
        m7(incx,incy) = m7(incx,incy) + zeta(3,j,incx,incy)*zeta(3,j,incx,incy)
        m8(incx,incy) = m8(incx,incy) + zeta(1,j,incx,incy)*zeta(2,j,incx,incy)
        m9(incx,incy) = m9(incx,incy) + zeta(1,j,incx,incy)*zeta(3,j,incx,incy)
        m10(incx,incy) = m10(incx,incy) + zeta(2,j,incx,incy)*zeta(3,j,incx,incy)

        m11(incx,incy) = m11(incx,incy) + zeta(1,j,incx,incy)*sqm4
        m12(incx,incy) = m12(incx,incy) + zeta(2,j,incx,incy)*sqm4
        m13(incx,incy) = m13(incx,incy) + zeta(3,j,incx,incy)*sqm4
     end do
  end if!icellPN

  !--- iSampleStep��˵���̤���� ---
  if(mod(irep,iSampleStep)==0) then

     m0(incx,incy) = m0(incx,incy) / dble(iSampleStep)
     m1(incx,incy) = m1(incx,incy) / dble(iSampleStep)
     m2(incx,incy) = m2(incx,incy) / dble(iSampleStep)
     m3(incx,incy) = m3(incx,incy) / dble(iSampleStep)
     m4(incx,incy) = m4(incx,incy) / dble(iSampleStep)

     m5(incx,incy) = m5(incx,incy) / dble(iSampleStep)
     m6(incx,incy) = m6(incx,incy) / dble(iSampleStep)
     m7(incx,incy) = m7(incx,incy) / dble(iSampleStep)
     m8(incx,incy) = m8(incx,incy) / dble(iSampleStep)
     m9(incx,incy) = m9(incx,incy) / dble(iSampleStep)
     m10(incx,incy) = m10(incx,incy) / dble(iSampleStep)

     m11(incx,incy) = m11(incx,incy) / dble(iSampleStep)
     m12(incx,incy) = m12(incx,incy) / dble(iSampleStep)
     m13(incx,incy) = m13(incx,incy) / dble(iSampleStep)

     !--- ���Ƚ�� ---
     !--- ���Ƚ��ϡ��ƥ���ˤ����벹�٤ˤĤ���500��ʬʿ�Ѥ�Ȥꡤ ---
     
     ste_mac_m0(incx,incy) = ste_mac_m0(incx,incy) + m0(incx,incy)
     ste_mac_m4(incx,incy) = ste_mac_m4(incx,incy) + m4(incx,incy)

     !----------------------------------------------------------------------
     ! iSampleStep��˵�᤿����̤򡤤����500��ʬʿ�Ѥ���
     ! ����500��ʬ��ʿ���ͤ���Ӥ����Ƚ�ꤹ��
     !----------------------------------------------------------------------
     
     if(mod(irep,iSampleStep*500)==0) then
        !--- ���ʤ��礭���ͤˤʤ�ʤ��褦�˰��500�ǳ�äȤ� ---
        ste_mac_m0(incx,incy) = ste_mac_m0(incx,incy) / 500.0d0
        ste_mac_m4(incx,incy) = ste_mac_m4(incx,incy) / 500.0d0
        !ste_jud = dabs((ste_mac_t-bste_mac_t)/ste_mac_t)
        ste_jud_m0(incx,incy) = (ste_mac_m0(incx,incy)-bste_mac_m0(incx,incy))/ste_mac_m0(incx,incy)
        ste_jud_m4(incx,incy) = (ste_mac_m4(incx,incy)-bste_mac_m4(incx,incy))/ste_mac_m4(incx,incy)

        !--- ste_jud��1.4d-2����礭����Фޤ����Ǥʤ��Ȥ���ste_mac_m0,m4��ꥻ�å� ---
        !--- ��������3��Ϣ³1.4d-2̤���ˤʤä����Ϥ���ʹ�1.4d-2����礭���ʤä��Ȥ��Ƥ����Ȥ��� ---
        !---------------------------------------------------------------------------------- 
        ! ����Ƚ�ꤷ�Ƥ���10�����ƥåפˤ����Ƶ���̤���Ƥ��ꡤ���ε���̤�0.01�ʲ���
        ! ���ˤʤäƤۤ����Τǡ�5�����ƥåפˤ����ƤϤ��΢�2�ܤ�0.014�ʲ��ˤʤ�Ф褤�Ȥ��롥
        !----------------------------------------------------------------------------------

        if(dabs(ste_jud_m0(incx,incy))>5.0d-2 .or. dabs(ste_jud_m4(incx,incy))>1.4d-2 .and. ste2(incx,incy)==iTotalStep) then
           bste_mac_m0(incx,incy) = ste_mac_m0(incx,incy)
           ste_mac_m0(incx,incy)  = 0.0d0
           bste_mac_m4(incx,incy) = ste_mac_m4(incx,incy)
           ste_mac_m4(incx,incy)  = 0.0d0
           ste(incx,incy) = 0
        else
           ste(incx,incy) = ste(incx,incy) + 1

           if(ste(incx,incy)==3) then
              ste2(incx,incy)=0
           end if
        end if

        !-- 3��Ϣ³�����������Ф��Υ����ste2��1�������䤹 --
        if(ste2(incx,incy)<iTotalStep) ste2(incx,incy) = ste2(incx,incy) + 1

        !-- ������ˤ����롤ste2�κǾ��ͤ���� --
        if(ste2(incx,incy)<min_ste) min_ste = ste2(incx,incy)

        !-- ������ˤ����롤ste2�κ����ͤ���� --
        if(ste2(incx,incy)>max_ste) max_ste = ste2(incx,incy)

     end if

  end if

end subroutine macro

!------------------------------------------------------------------
!   out put data
!------------------------------------------------------------------
  subroutine out_put
    use commn
    implicit none
    integer i,j
    integer, parameter :: irep_sta =10000

    !--- iSampleStep���ƥåפ��ȤΤ��줾����̤ι���ͤ���� ---
    write(40) time,x1(incx),x2(incy),m0(incx,incy),m1(incx,incy),m2(incx,incy),&
         m3(incx,incy),m4(incx,incy),m5(incx,incy),m6(incx,incy),m7(incx,incy),&
         m8(incx,incy),m9(incx,incy),m10(incx,incy),m11(incx,incy),m12(incx,incy),m13(incx,incy)

  end subroutine out_put

!------------------------------------------------------------------
!   subroutine distribution function
!------------------------------------------------------------------
  subroutine dis_fun
    use commn
    implicit none
    
    integer j,izeta1,izeta2

    nfun = nfun + 1 !���֥롼�����ƤӽФ�����򥫥����
    
    !---------------------------------------------------------------------------------- 
    ! ľ��iSampleStep*500���ƥå״֤ˤ�����zeta1�κ����ͤ���®��ʬ�۴ؿ���®�ٹ�ߤ����
    !----------------------------------------------------------------------------------
    if(nfun==1) then
       dzeta1 = zmax1 / dble(nzeta1)
       dzeta2 = zmax2 / dble(nzeta2)
       
       !--- zmax���礭��������� ---
       if(dzeta1>0.5d0) dzeta1 = 0.2d0
       if(dzeta2>0.5d0) dzeta2 = 0.2d0

       f1(:,:,:) = 0.0d0
       f2(:,:,:) = 0.0d0
    end if

    do j=1,icellPN(incx,incy)
      
       if(zeta(1,j,incx,incy)>=0)then
          izeta1 = int(zeta(1,j,incx,incy)/dzeta1) 
       else
          izeta1 = int(zeta(1,j,incx,incy)/dzeta1) - 1
       end if

       if(zeta(2,j,incx,incy)>=0)then
          izeta2 = int(zeta(2,j,incx,incy)/dzeta2) 
       else
          izeta2 = int(zeta(2,j,incx,incy)/dzeta2) - 1
       end if

       f1(izeta1,incx,incy) = f1(izeta1,incx,incy) + 1.0d0
       f2(izeta2,incx,incy) = f2(izeta2,incx,incy) + 1.0d0

    end do


    !--- ��̤ν��� ---
    if(mod(irep,iSampleStep*100)==0 ) then
       do izeta1 = -nzeta1,nzeta1
          f1(izeta1,incx,incy) = f1(izeta1,incx,incy) / dble(iSampleStep) / 100
          write(50,'((f15.4),3(i8),(f15.2))') time,incx,incy,izeta1,f1(izeta1,incx,incy)
       end do

       do izeta2 = -nzeta2,nzeta2
          f2(izeta2,incx,incy) = f2(izeta2,incx,incy) / dble(iSampleStep) / 100
          write(55,'((f15.4),3(i8),(f15.2))') time,incx,incy,izeta2,f2(izeta2,incx,incy)
       end do

       write(50,*)
       write(55,*)

       !--- f1,f2������ ---
       do izeta1 = -nzeta1,nzeta1
          f1(izeta1,incx,incy) = 0.0d0
       end do

       do izeta2 = -nzeta2,nzeta2
          f2(izeta2,incx,incy) = 0.0d0
       end do
    end if

       
  end subroutine dis_fun
 
end program dsmc
