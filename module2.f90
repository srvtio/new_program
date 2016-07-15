module commn
  !$ use omp_lib
  implicit none
  integer, parameter :: iParticleNumber = 400000
  integer, parameter :: iCellNumber1 = 32
  integer, parameter :: iCellNumber2 = 32
  integer, parameter :: iAssumptionPN = 3*iParticleNumber/iCellNumber1/iCellNumber2
  integer, parameter :: iSeed = 1
  integer, parameter :: iSampleStep = 100
  integer, parameter :: iTotalStep = 1000000
  integer, parameter :: iFinalStep = 800000
  integer, parameter :: iInitialStep = 500000
 
  real(8), parameter :: L1 = 1.0d0
  real(8), parameter :: L2 = 1.0d0
  real(8), parameter :: L3 = 2.0d0
  real(8), parameter :: Knudsen1 = 1.0d0
  real(8), parameter :: Knudsen2 = 100.0d0
  real(8), parameter :: Tau1 = 1.0d0                !temperature of wall1
  real(8), parameter :: TauObject0 = 1.0d0          !temperature of porous left
  real(8), parameter :: TauObject1 = 1.0d0          !temperature of porous right
  real(8), parameter :: RhoInitial = 1.0d0
  real(8), parameter :: Rho0 = 0.0d0
  real(8), parameter :: Rho1 = 1.0d0

  real(8), parameter :: dx1 = 2*L1/dble(iCellNumber1)
  real(8), parameter :: dx2 = 2*L2/dble(iCellNumber2)
  real(8), parameter :: dt = 0.0005d0
  real(8), parameter :: pi = 4.0d0 * atan(1.0d0)
  
  integer itemporary !���Ū���Ѥ����ѿ���������
  real(8) dtemporary !���Ū���Ѥ����ѿ��ʼ¿���

  real(8) delta
  real(8) porosity
  integer iPattern  !¿�������֥ѥ�����
  integer iParallel !����� 

  !--- ��������ФƤ���γ�Ҥΰ��֤�®�٤Υǡ��� ---
  real(8), allocatable :: zeta_boundary_right(:,:)
  real(8), allocatable :: x_boundary_right(:,:)
  real(8), allocatable :: zeta_boundary_left(:,:)
  real(8), allocatable :: x_boundary_left(:,:)

  integer icellPN(iCellNumber1,iCellNumber2) !�ƥ������γ�ҿ�

  integer icollisionCPN !�ƥ���ˤ����ơ������ΤȾ��ͤ���γ�ҿ���������
  real(8) dcollisionCPN

  integer irep

  real(8) dv(iCellNumber1,iCellNumber2)

  real(8) x(2,iAssumptionPN,iCellNumber1,iCellNumber2)
  real(8) x1(iCellNumber1)
  real(8) x2(iCellNumber2)
  real(8) zeta(3,iAssumptionPN,iCellNumber1,iCellNumber2)

  real(8) m0(iCellNumber1,iCellNumber2)
  real(8) m1(iCellNumber1,iCellNumber2)
  real(8) m2(iCellNumber1,iCellNumber2)
  real(8) m3(iCellNumber1,iCellNumber2)
  real(8) m4(iCellNumber1,iCellNumber2)
  real(8) m5(iCellNumber1,iCellNumber2)
  real(8) m6(iCellNumber1,iCellNumber2)
  real(8) m7(iCellNumber1,iCellNumber2)
  real(8) m8(iCellNumber1,iCellNumber2)
  real(8) m9(iCellNumber1,iCellNumber2)
  real(8) m10(iCellNumber1,iCellNumber2)
  real(8) m11(iCellNumber1,iCellNumber2)
  real(8) m12(iCellNumber1,iCellNumber2)
  real(8) m13(iCellNumber1,iCellNumber2)

  real(8) Knudsen(iCellNumber1,iCellNumber2)  !�����ޤ��ѹ��Ѥߡ�2016.6.16��
  real(8) KnudsenMin !Knudsen���κǾ���

  real(8) dtsub !������狼�����äƤ���γ�ҿ���dt�֤˶����˥Х餱��褦�ˤ���ݤ��Ѥ���
  real(8) time

  !--- ���Ƚ�� ---
  integer ste(iCellNumber1,iCellNumber2)
  integer ste2(iCellNumber1,iCellNumber2)
  integer min_ste,max_ste
  real(8) ste_mac_m0(iCellNumber1,iCellNumber2)
  real(8) bste_mac_m0(iCellNumber1,iCellNumber2)
  real(8) ste_mac_m4(iCellNumber1,iCellNumber2)
  real(8) bste_mac_m4(iCellNumber1,iCellNumber2)
  real(8) ste_jud_m0(iCellNumber1,iCellNumber2)
  real(8) ste_jud_m4(iCellNumber1,iCellNumber2)

  real(8) j_ste 

  !--- γ�ҿ����Ѥ������Ƚ�� ---
  ! integer s_fluxm !�����鼡��iSampleStep*10���ƥå״֤�ή��������̣��γ�ҿ� 
  ! integer s_fluxp !�����鼡��iSampleStep*10���ƥå״֤�ή��������̣��γ�ҿ�
  ! integer b_fluxm !����������iSampleStep*10���ƥå״֤�ή��������̣��γ�ҿ� 
  ! integer b_fluxp !����������iSampleStep*10���ƥå״֤�ή��������̣��γ�ҿ�
  ! integer s_flux  !����iSampleStep*10���ƥå�ή��γ�Ҥ�ή��γ�Ҥκ�
  ! integer b_flux  !����iSampleStep*10���ƥå�ή��γ�Ҥ�ή��γ�Ҥκ�

  !--- ή«�׻� ---
  integer flux_inp,flux_inm     !�����������äƤ���γ�ҿ�
  integer flux_outp ,flux_outm  !��������ФƹԤ�γ�ҿ�
  integer flux_outp1,flux_outm1
  integer flux_outp2,flux_outm2
  integer flux_outp3,flux_outm3

  ! !--- ʪ�Τ��ΰ������ ---
  ! real(8), parameter :: objx1_l = -0.5d0*L1
  ! real(8), parameter :: objx1_r = 0.5d0*L1
  ! real(8), parameter :: objx2_b = -0.5d0*L2
  ! real(8), parameter :: objx2_u = 0.5d0*L2

  !--- ®��ʬ�۴ؿ� ---
  integer nfun !���֥롼����򲿲�ƤӽФ����򥫥���� ---
  real(8) dzeta1,dzeta2 !®�٤ι����
  integer, parameter :: nzeta1 = 20
  integer, parameter :: nzeta2 = 20
  real(8) f1(-nzeta1:nzeta1,iCellNumber1,iCellNumber2)
  real(8) f2(-nzeta2:nzeta2,iCellNumber1,iCellNumber2)
  real(8) zmax1,zmax2

  !---���ե�����̾��� ---
  character filepara*128
  character fileinidata*128
  character fileknu*128
  character fileflu*128
  character filemac*128
  character filedis1*128
  character filedis2*128
  character filedata*128
  character filedataraw*128

end module commn

module crf
  implicit none
  integer,save :: count_rf(0:31)
  integer,save :: count_rff(0:31)
  integer total_rf

end module crf
