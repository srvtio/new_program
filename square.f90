program kptest4
  implicit none
  integer, parameter :: ncx = 32
  integer, parameter :: ncy = 32
  real(8), parameter :: kp1 = 0.0
  real(8), parameter :: kp2 = 100

  integer, parameter :: seed = 41
  integer, parameter :: npor = 800  !¿�������ΤθĿ�
  integer, parameter :: porsize = 1
  real(8) ep !����Ψ

  integer i,j
  integer irep
  integer incx,incy

  integer xcell1,xcell2
  integer ycell1,ycell2
  integer xcell1p,xcell2m
  integer ycell1p,ycell2m

  integer count
  integer sumreg
  integer j0or1(ncx,ncy)

  integer fin_n !�������줿¿�������ΤθĿ�

  integer KnudsenPt
  integer iPattern  !¿�������֥ѥ�����
  integer iParallel !�����

  real(8) kp(ncx,ncy)
  real(8) lbotx(1:npor),lboty(1:npor) !¿�������Τ����������ΰ�κ����Υ���ΰ���

  character filedata*128
  character filepara*128

  ! write(*,*) "KnudsenPt="
  ! read(*,*) KnudsenPt
  KnudsenPt = 1
  
  write(*,*) "iPattern="
  read(*,*) iPattern

  ! write(*,*) "iParallel="
  ! read(*,*) iParallel
  iParallel = 1

  write(filedata,'("KnudsenData(CN1_",i2,",CN2_",i2,",pt",i3,").dat")') &
       ncx,ncy,iPattern
  write(filepara,'("Parameter(CN1_",i2,",CN2_",i2,").dat")') &
       ncx,ncy

  open(10,file=filedata) !kp������
  open(20,file=filepara)

  write(20,'(i3)') KnudsenPt
  write(20,'(i3)') iParallel-1
  write(20,'(i3)') iPattern

  count = 0
  fin_n = 0
  j0or1(:,:) = 0
  kp(:,:) = kp2

  do irep=1,npor

     sumreg = 0
     count = count + 1

     !--------------------------------------------------------------
     ! ¿�������Τκ����ΰ��֤������˷���
     ! x�����˴ؤ��Ƥ�ξü2����ʬ������¿���������֤���
     !--------------------------------------------------------------
     !lbotx(irep) = 2 + ceiling(dble(ncx-4-porsize+1)*rf(seed,0))
     lbotx(irep) = ceiling(dble(ncx-porsize+1)*rf(seed,0))
     lboty(irep) = ceiling(dble(ncy)*rf(seed,0))
     !write(*,*) lbotx(irep),lboty(irep)
     !--------------------------------------------------------------
     !¿���������֤����������ΰ�λͶ��ΰ��֤����
     !--------------------------------------------------------------
     xcell1 = lbotx(irep)
     ycell1 = lboty(irep)
     xcell2 = xcell1 + porsize - 1
     ycell2 = ycell1 + porsize - 1

     !------------------------------------------------------------------
     ! �������ΰ���󤹤٤�kp1�ˤ���
     !------------------------------------------------------------------
     
     !---¿�����ΰ褬x����y����������ξ����ۤ������---
     if(xcell2>ncx .and. ycell2>ncy) then 

        xcell2 = xcell2 - ncx
        ycell2 = ycell2 - ncy

        !--------------------------------------------------------------------------
        ! j0or1�ϳƥ���Ǥ��Ǥ�¿�������Τ�����Ƥ����1,�ޤ�����Ƥ��ʤ����0
        ! ��ɽ�������ߤΥ��ƥåפǺ��������������ΰ�ˡ�����ޤǺ�������¿������
        ! ���֤�ʤ�����sumreg��Ƚ�ꤹ��
        !--------------------------------------------------------------------------
        do i=xcell1,ncx
           do j=1,ycell2
              sumreg = sumreg + j0or1(i,j)
           end do
           do j=ycell1,ncy
              sumreg = sumreg + j0or1(i,j)
           end do
        end do
        do i=1,xcell2
           do j=1,ycell2
              sumreg = sumreg + j0or1(i,j)
           end do
           do j=ycell1,ncy
              sumreg = sumreg + j0or1(i,j)
           end do
        end do

        !-- 100�����������֤��Ƥ�¿���������֤äƤ��ޤ����Ϸ׻���λ --
        if(count>500) then
           write(*,*) "error"
           exit
        end if
        !-- ¿���������֤����������ΰ褬���֤�ʤ��褦�ˤ��� --
        if(sumreg>0) cycle

        kp(xcell1:ncx,ycell1:ncy) = kp1
        kp(1:xcell2,ycell1:ncy) = kp1
        kp(1:xcell2,1:ycell2) = kp1
        kp(xcell1:ncx,1:ycell2) = kp1

        !-------------------------------------------------------------------
        ! j0or1��¿���������������ΰ褬���֤�ͤ��褦��Ƚ�ꤹ�뤿��Τ��
        ! ¿�������Τ��������ΰ��1���ɤꡤ���Υ��ƥåפ�j0or1�ι����sumreg
        ! ��׻����Ƥ��줬0�����礭���ʤ�Ф���ޤǤΥ��ƥåפǤ����ΰ��
        ! ¿�������Τ����줿�Ȥ������Ȥ��狼��
        !-------------------------------------------------------------------
        j0or1(xcell1:ncx,ycell1:ncy) = 1
        j0or1(1:xcell2,ycell1:ncy) = 1
        j0or1(1:xcell2,1:ycell2) = 1
        j0or1(xcell1:ncx,1:ycell2) = 1

        !---¿�����ΰ褬x����������ۤ������---
     else if(xcell2>ncx .and. ycell2<=ncy) then 

        xcell2 = xcell2 - ncx

        do j=ycell1,ycell2
           do i=1,xcell2
              sumreg = sumreg + j0or1(i,j)
           end do
           do i=xcell1,ncx
              sumreg = sumreg + j0or1(i,j)
           end do
        end do

        if(count>500) then
           write(*,*) "error"
           exit
        end if
        if(sumreg>0) cycle

        kp(xcell1:ncx,ycell1:ycell2) = kp1
        kp(1:xcell2,ycell1:ycell2) = kp1

        !---¿�������Τ��������ΰ��1���ɤ�---
        j0or1(xcell1:ncx,ycell1:ycell2) = 1
        j0or1(1:xcell2,ycell1:ycell2) = 1

        !---¿�����ΰ褬y����������ۤ������---
     else if(xcell2<=ncx .and. ycell2>ncy) then 

        ycell2 = ycell2 - ncy

        do i=xcell1,xcell2
           do j=1,ycell2
              sumreg = sumreg + j0or1(i,j)
           end do
           do j=ycell1,ncy
              sumreg = sumreg + j0or1(i,j)
           end do
        end do

        if(count>500) then
           write(*,*) "error"
           exit
        end if
        if(sumreg>0) cycle

        kp(xcell1:xcell2,ycell1:ncy) = kp1
        kp(xcell1:xcell2,1:ycell2) = kp1

        !---¿�������Τ��������ΰ��1���ɤ�---
        j0or1(xcell1:xcell2,ycell1:ncy) = 1
        j0or1(xcell1:xcell2,1:ycell2) = 1

        !---¿�����ΰ褬������ۤ��ʤ��ä����---
     else if(xcell2<=ncx .and. ycell2<=ncy) then 

        do j=ycell1,ycell2
           do i=xcell1,xcell2
              sumreg = sumreg + j0or1(i,j)
           end do
        end do

        if(count>500) then
           write(*,*) "error"
           exit
        end if
        if(sumreg>0) cycle

        kp(xcell1:xcell2,ycell1:ycell2) = kp1

        !---¿�������Τ��������ΰ��1���ɤ�---
        j0or1(xcell1:xcell2,ycell1:ycell2) = 1


     end if


     !---�������ΰ褫��Ѥ�ݤ��---
     xcell1p = xcell1 + 1
     xcell2m = xcell2 - 1
     ycell1p = ycell1 + 1
     ycell2m = ycell2 - 1 

     !---�����������---
     if(xcell1p>ncx) xcell1p = 1
     if(xcell2m<1   ) xcell2m = ncx
     if(ycell1p>ncy) ycell1p = 1
     if(ycell2m<1   ) ycell2m = ncy

     fin_n = fin_n + 1
  end do

  !------------------------------------------------------------------------
  ! ����Ψ�η׻�
  ! ¿�������������ΰ�porsize*porsize����Ͷ��ζ���0�������¿������Ĥ�
  ! �ΰ�Ȥ��Ʒ׻�
  !------------------------------------------------------------------------
  ep = 1 - dble((porsize*porsize-0)*fin_n)/dble(ncx*ncy)

  write(*,*) fin_n,ep
  write(20,'(f8.4)') ep

  !------------------------------------------------------------------------

  do incy=1,ncy
     do incx=1,ncx
        !---kp().dat�˽���---
        write(10,*) incx,incy,kp(incx,incy)
        
        ! !--- 1/19�ɲ� ̤���� ---
        ! if(kp(incx,incy)==kp1 .and. kp(incx+1,incy)==kp2) then
        !    por_area = por_area + dx2
        ! else if(kp(incx,incy)==kp2 .and. kp(incx+1,incy)==kp1) then
        !    por_area = por_area + dx2
        ! end if
     end do
     write(10,*)
  end do

contains

  !------------------------------------------------------------------
  !   function rf
  !------------------------------------------------------------------

  real(8) function rf(iidum,ip)

    implicit none
    integer, intent(in) :: iidum,ip
    integer idum
    integer, save       :: ma(55,0:7),inext(0:7),inextp(0:7),iff(0:7)=0
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

end program kptest4
