!--- コードの概要 -----------------------------------------------------------------------
!
! 左右は完全凝縮境界条件，上下は周期境界条件の２次元正方形領域について考える．
! 分子同士の衝突はなく，分子は散乱体と衝突しながら空間上を移動する．
! 空間をセルに分割し，各セルごとに散乱体との衝突頻度に対応するクヌーセン数を与えることができる．
! 
!
! 定常判定あり
!
! 読取データ
! ・KnudsenData(CN1,CN2,iPattern).dat：各セルのクヌーセン数を読み取る．
! ・Parameter(CN1,CN2).dat：並列数，物体パターン，空隙率のデータを読み取る
!
! 出力データ
! （定常と判定されてから出力）
! ・Macroscopic(Rho1,Tau1,iPattern).dat：巨視量を求めるための素データ
!       Rho1:右側の完全凝縮境界条件における密度
!       Tau1:右側の完全凝縮境界条件における温度
!       iPattern:領域内に置いた物体のパターン（どんな多孔質体なのか，どんな温度分布を与えた物体なのかなど）
! ・DistributionFunc1(Rho1,Tau1,iPattern).dat：1方向の速度分布関数
! ・DistributionFunc2(Rho1,Tau1,iPattern).dat：2方向の速度分布関数
! ・Flux(Rho1,Tau1,iPattern).dat：左右の境界から流入流出する粒子の個数をカウント
! ・KnudsenConf(Rho1,Tau1,iPattern).dat：座標とKnudsen数のデータ
! ・Data(Rho1,Tau1,iPattern).dat：初期で設定したパラメータや，計算によって得られる条件データ等
! ・DataRaw(Rho1,Tau1,iPattern).dat：計算によって得られる条件データ
! ・BoundsNumber(Rho1,Tau1,iPattern).dat 散乱体物体に流入流出する粒子の数を数える
!
!--- 変更履歴 --------------------------------------------------------------------------
!
! 2016.6.2 巨視量の求め方変更 
! 2016.6.6 速度分布関数を求める(1方向のみ) 
! 2016.6.13 速度分布関数を求める(2方向も) 
! 2016.6.14 衝突確率100%と0%用に書き換え
! 2016.6.27 出力するファイル名の変更
! 2016.6.27 moduleを別ファイルにする．
!     実行は gfortran -o dsmc_s_18 module.f90 dsmc_s_18.f90 -O3 -fopenmp 
! 2016.6.27 DataRaw(Rho1,Tau1,iPattern).datファイルの作成
! 2016.6.27 散乱体の右側左側で温度を変更していた部分を無くした
! 2016.6.28 Data(Rho1,Tau1,iPattern).datに境界条件の密度，温度の設定値を出力
! 2016.7.11 定常判定の部分のif文修正
! 2016.7.11 粒子の初期配置が中心に配置された正方形物体の周りに配置されないようになっていたので
!           物体の位置関係なく粒子をランダムに全域に配置させるようにする．
!           それに伴ってdeltaの計算も変更
! 2016.7.11 sum()とminval()を用いて初期の粒子数やKnusenの最小値の計算を行う
! 2016.7.11 初期配置の各セルにおける粒子の個数の計算の部分を変更,位置の移し替えの部分も変更
! 2016.7.11 割り算を掛け算に変更
! 2016.7.11 zmaxの出力するところを修正
! 2016.7.12 dis_funのdzetaを求める部分を変更（zmaxが大きくなりすぎていた場合の対処）
! 2016.7.15 false sharingを考慮してcount_rfを消した
! 2016.7.15 通常のKnudsenのパターン(0)と，衝突確率0%があるパターン(1)の場合分け
! 2016.7.15 指定した領域における流入流出粒子の数を数えるコードの追加
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
  integer iBoundaryPNLeft,iBoundaryPNRight !左右境界から流入する粒子数
  integer np_max_re
  real(8) djrandom

  !-- 散乱体物体に流入流出する粒子の計算 --
  real(8) boundx

  real(8), parameter :: Redx1 = 1.0d0 / dx1
  real(8), parameter :: Redx2 = 1.0d0 / dx2

  integer, allocatable :: icellPN_sub(:,:,:)
  real(8), allocatable :: xx(:,:,:,:,:)
  real(8), allocatable :: zzeta(:,:,:,:,:)

  !--- Parameter.datファイルから並列数，物体パターン，空隙率を読み取る ---
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

  !--- Knudsen数を各セルに振り分ける ---
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

!----- 初期条件 -----
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
  
  !--- 各セルの体積 ---
  dv(:,:) = dx1*dx2*L3

  !--- 各セルの粒子数を決める ---
  icellPN(:,:) = int(dble(iParticleNumber)/dble(iCellNumber1)/dble(iCellNumber2))

  !--- 初期の全粒子数 ---
  np_max_re = sum(icellPN)

  !--- シミュレーション粒子の質量 ---
  delta = L1*2.0d0*L2*2.0d0*L3 * RhoInitial / np_max_re 

  !--- Knudsenの最小値を計算 ---
  KnudsenMin = minval(Knudsen)

  !--- 完全凝縮境界条件で右，左から入ってくる粒子数の計算      ---
  !--- 実数を整数にするときに確率dBPNSubRightとdBPNSubLeftを用いて統計的に理論値に合うようにする  ---
  dBoundaryPNLeft = Rho0 / 2.0d0 * sqrt(1.0d0/pi) * L2*2.0d0*L3 * dt / delta  
  dBPNSubLeft = dBoundaryPNLeft - dble(int(dBoundaryPNLeft))
  dBoundaryPNRight = Rho1 / 2.0d0 * sqrt(Tau1/pi) * L2*2.0d0*L3 * dt / delta  
  dBPNSubRight = dBoundaryPNRight - dble(int(dBoundaryPNRight))

  allocate (x_boundary_right(1:2,1:int(dBoundaryPNRight)+1))
  allocate (zeta_boundary_right(1:3,1:int(dBoundaryPNRight)+1))
  allocate (x_boundary_left(1:2,1:int(dBoundaryPNLeft)+1))
  allocate (zeta_boundary_left(1:3,1:int(dBoundaryPNLeft)+1))

  !--- 初期条件 ---
  do incy = 1,iCellNumber2
     do incx = 1,iCellNumber1

        !-- そのセルに粒子が0であれば計算を行わない --
        if(icellPN(incx,incy)>0) then
           do j=1,icellPN(incx,incy)
        
              !--- 位置，速度をランダムに配置 ---
              x(1,j,incx,incy) = -L1 + 2.0d0 * L1 * rf(iSeed,0)
              x(2,j,incx,incy) = -L2 + 2.0d0 * L2 * rf(iSeed,0)
              
              !count_rf(0) = count_rf(0) + 2

              !--- 速度分布を求める ---
              call makevv(1.0d0,dtemporary,zeta(1,j,incx,incy),0)
              call makevv(1.0d0,zeta(2,j,incx,incy),zeta(3,j,incx,incy),0)

              !count_rf(0) = count_rf(0) + 4

              !-- セル分割 --
                 !-- 領域を出た粒子は無視 --
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

                    !-- 新しいセルにおける粒子数 --
                    icellPN_sub(iincx,iincy,ip) = icellPN_sub(iincx,iincy,ip) + 1

                    !-- 新しい名前に変更
                    xx(1,icellPN_sub(iincx,iincy,ip),iincx,iincy,ip) = x(1,j,incx,incy)
                    xx(2,icellPN_sub(iincx,iincy,ip),iincx,iincy,ip) = x(2,j,incx,incy)
                    zzeta(:,icellPN_sub(iincx,iincy,ip),iincx,iincy,ip) = zeta(:,j,incx,incy)
                 end if

           end do
        end if

        !--- x,y座標の計算 ---
        x1(incx) = -L1 + dx1*dble(incx) - dx1/2.0d0
        x2(incy) = -L2 + dx2*dble(incy) - dx2/2.0d0

        write(18,'(f15.8,f15.8,f8.2)') x1(incx),x2(incy),Knudsen(incx,incy)

     end do
     write(18,*)
  end do

  icellPN(:,:) = icellPN_sub(:,:,ip)

  !--- 各粒子の情報を配列xに移し替える(iAssumptionPNより大きいjの値は0にしたい) ---
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

     !----- 繰り返し計算 -----
  do irep = 1,iTotalStep
     time = time + dt

     icellPN_sub(:,:,:) = 0

     min_ste = 10000
     max_ste = 1

     zeta_boundary_left(:,:) = 0.0d0
     zeta_boundary_right(:,:) = 0.0d0
     x_boundary_left(:,:) = 0.0d0
     x_boundary_right(:,:) = 0.0d0

     !--- 完全凝縮境界条件で左右の境界から入ってくる粒子数を，統計的に見て理論値と合うように計算 ---
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

     !-- 領域に入ってくる粒子の数を数える --
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

                 !-- 散乱体物体への流入流出する粒子の計算用 --
                 buf_x(1,j,incx,incy) = x(1,j,incx,incy)
                 buf_x(2,j,incx,incy) = x(2,j,incx,incy)

                 !-- 粒子の移流 --
                 x(1,j,incx,incy) = x(1,j,incx,incy) + dt * zeta(1,j,incx,incy)
                 x(2,j,incx,incy) = x(2,j,incx,incy) + dt * zeta(2,j,incx,incy)

                 !--- 散乱体への流入流出粒子の計算 ------------------------------------------
                 ! incy=iCellNumber2/4とincy=iCellNumber2/4+1の間にある
                 ! 境界における粒子の流入流出を調べる 
                 ! dtemporaryは粒子がその境界を通過したかを判断するためのもの（負なら通過）
                 ! boundxはobjx2_bを通過した点のx座標
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

                 !-- 上下の周期境界条件 --
400              continue
                
                 if(x(2,j,incx,incy)>L2) then
                    x(2,j,incx,incy) = x(2,j,incx,incy) - 2.0d0 * L2
                 else if(x(2,j,incx,incy)<-L2) then
                    x(2,j,incx,incy) = x(2,j,incx,incy) + 2.0d0 * L2
                 end if

                 if(x(2,j,incx,incy)>L2 .or. x(2,j,incx,incy)<-L2) go to 400

                 !-- 領域を出て行く粒子の数を数える --
                 if(x(1,j,incx,incy)> L1) flux_outp1 = flux_outp1 + 1
                 if(x(1,j,incx,incy)<-L1) flux_outm1 = flux_outm1 + 1

                 !-- セル分割 --
                 !-- 領域を出た粒子は無視 --
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

                    !-- 新しいセルにおける粒子数 --
                    icellPN_sub(iincx,iincy,ip) = icellPN_sub(iincx,iincy,ip) + 1
                    itemporary =  icellPN_sub(iincx,iincy,ip) 

                    !-- 新しい名前に変更
                    xx(1,itemporary,iincx,iincy,ip) = x(1,j,incx,incy)
                    xx(2,itemporary,iincx,iincy,ip) = x(2,j,incx,incy)
                    zzeta(:,itemporary,iincx,iincy,ip) = zeta(:,j,incx,incy)

                 end if

              end do !j
           end if
        end do !incy
     end do !incx
     !$omp end do

     !-- 完全凝縮境界条件 --    
     if(iBoundaryPNLeft>0) then

     !$omp do reduction(+:flux_outp3,flux_outm3)
     
        do jl=1,iBoundaryPNLeft

           !-- 境界条件の速度分布に従った速度を決める ---
           call makev(Rho0,zeta_boundary_left(1,jl),ip)
           call makevv(Rho0,zeta_boundary_left(2,jl),zeta_boundary_left(3,jl),ip)

           !count_rf(ip) = count_rf(ip) + 3

           !--- 境界から入ってきた粒子の位置，速度を求める ---
           dtsub = dt * rf(iSeed,ip)
           dtemporary2 = -L2 + rf(iSeed,ip) * L2 * 2.0d0
           x_boundary_left(1,jl) = -L1 + zeta_boundary_left(1,jl) * dtsub
           x_boundary_left(2,jl) = dtemporary2 + zeta_boundary_left(2,jl) * dtsub

           ! !--- 散乱体への流入流出粒子の計算 ------------------------------------------
           ! ! incy=iCellNumber2/4とincy=iCellNumber2/4+1の間にある
           ! ! 境界における粒子の流入流出を調べる 
           ! ! dtemporaryは粒子がその境界を通過したかを判断するためのもの（負なら通過）
           ! ! boundxはobjx2_bを通過した点のx座標
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

           !--- 上下の周期境界条件 ---
600        continue
          
           if(x_boundary_left(2,jl)>L2) then
              x_boundary_left(2,jl) = x_boundary_left(2,jl) - 2.0d0 * L2
           else if(x_boundary_left(2,jl)<-L2) then
              x_boundary_left(2,jl) = x_boundary_left(2,jl) + 2.0d0 * L2
           end if

           if(x_boundary_left(2,jl)>L2 .or. x_boundary_left(2,jl)<-L2) go to 600

           !-- 領域を出て行く粒子の数を数える --
           if(x_boundary_left(1,jl)> L1) flux_outp3 = flux_outp3 + 1
           if(x_boundary_left(1,jl)<-L1) flux_outm3 = flux_outm3 + 1

           !--- 粒子の新しいセル番号を求める ---
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

              !--- 新しいセルでの粒子数をカウント ---
              icellPN_sub(iincx,iincy,ip) = icellPN_sub(iincx,iincy,ip) + 1
              itemporary = icellPN_sub(iincx,iincy,ip) 
              !--- 新しい粒子の名前を決める ---
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

           !-- 境界条件の速度分布に従った速度を決める ---
           call makev(1.0d0,zeta_boundary_right(1,jr),ip)
           call makevv(1.0d0,zeta_boundary_right(2,jr),zeta_boundary_right(3,jr),ip)

           !count_rf(ip) = count_rf(ip) + 3

           !--- 境界から入ってきた粒子の位置，速度を求める ---
           zeta_boundary_right(1,jr) = -zeta_boundary_right(1,jr)
           dtsub = dt * rf(iSeed,ip)
           x_boundary_right(1,jr) = L1 + zeta_boundary_right(1,jr) * dtsub
           x_boundary_right(2,jr) = -L2 + rf(iSeed,ip) * L2 * 2.0d0 + zeta_boundary_right(2,jr) * dtsub

           !count_rf(ip) = count_rf(ip) + 2

           !--- 上下の周期境界条件 ---
500        continue
          
           if(x_boundary_right(2,jr)>L2) then
              x_boundary_right(2,jr) = x_boundary_right(2,jr) - 2.0d0 * L2
           else if(x_boundary_right(2,jr)<-L2) then
              x_boundary_right(2,jr) = x_boundary_right(2,jr) + 2.0d0 * L2
           end if

           if(x_boundary_right(2,jr)>L2 .or. x_boundary_right(2,jr)<-L2) go to 500

           !-- 領域を出て行く粒子の数を数える --
           if(x_boundary_right(1,jr)> L1) flux_outp2 = flux_outp2 + 1
           if(x_boundary_right(1,jr)<-L1) flux_outm2 = flux_outm2 + 1

           !--- 粒子の新しいセル番号を求める ---
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

              !--- 新しいセルでの粒子数をカウント ---
              icellPN_sub(iincx,iincy,ip) = icellPN_sub(iincx,iincy,ip) + 1
              itemporary = icellPN_sub(iincx,iincy,ip)
              !--- 新しい粒子の名前を決める ---
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

     !--- 各セルの粒子数を配列icellPNに移し替える ---
     !並列 1
     do incy=1,iCellNumber2
        do incx=1,iCellNumber1
           icellPN(incx,incy) = 0
           do ip=0,iParallel
              icellPN(incx,incy) = icellPN(incx,incy) + icellPN_sub(incx,incy,ip)
           end do
        end do
     end do

     !--- 各粒子の情報を配列xに移し替える ---
     x(:,:,:,:) = 0.0d0
     zeta(:,:,:,:) = 0.0d0

     !並列 2
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

     !--- 衝突の計算 ---
     !$ ip = omp_get_thread_num() 
     !$omp do collapse(2)
     
     do incx=1,iCellNumber1
        do incy=1,iCellNumber2

           !--- 通常のKnudsenのパターン(0)と，衝突確率0%があるパターン(1)の場合分け ---
           if(KnudsenPt == 0) then

              !--- dtの間に衝突する粒子数の計算 ---
              dcollisionCPN = dble(icellPN(incx,incy))*(1.0d0-exp(-dt/Knudsen(incx,incy)))

              !--- Knudsen=Knudsen2なら衝突しない ---
              if(knudsen(incx,incy)==Knudsen2) dcollisionCPN = 0.0d0

              dtemporary = dcollisionCPN - dble(int(dcollisionCPN)) 
              !-- 衝突する粒子数 --
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

           !--- 1セルあたりicollisionCPN個の粒子が衝突する ---
           if(icollisionCPN>0) then
              do j=1,icollisionCPN
                 !--------------------------------------------------------------------------- 
                 !  jは数が小さい順に,　移流した粒子，境界から入ってきた粒子　になっているため， 
                 !  そのままjが1からicollisionCPNまで計算すると，衝突する粒子の状態にかたよりが生じる．
                 !  よって，jがランダムにバラけるように乱数を用いて衝突する粒子を決める
                 !--------------------------------------------------------------------------- 
               
                 !--- 通常のKnudsenのパターン(0)と，衝突確率0%があるパターン(1)の場合分け ---
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

!----- 速度の最大値を求める -----
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

 !----- irep>iInitialStep のとき，巨視量を計算する -----
     if(irep>iInitialStep) then
        do incy=1,iCellNumber2
           do incx=1,iCellNumber1

              !--- 巨視量を計算 ---
              call macro

              if( irep>j_ste .or. irep>iFinalStep ) then
                 !-- 速度分布関数の計算と出力 --
                 call dis_fun
              end if

              !--- 定常となってからiSampleStepごとに行う処理 ---
              if(irep>j_ste .or. irep>iFinalStep) then
                 if( mod(irep,iSampleStep)==0) then
                    !-- 巨視量の出力
                    call out_put
                 end if
              end if

              !--- 巨視量と速度分布関数の初期化 ---
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

        !--- fluxの出力 ---
        if(irep>j_ste .or. irep>iFinalStep) then
           if( mod(irep,iSampleStep)==0) then
              write(20,'(i8,i8,i8,i8)') flux_inp,flux_inm,flux_outp,flux_outm    
           end if
        end if

     end if

     !--- ターミナルに出力 ---
     if(irep<=iInitialStep .and. mod(irep,1000)==0) then
        write(*,'("step:",i8,"/",i8,", time:,",f8.2,"s",", zmax1:", f8.3," m/s^2")') &
             irep,iTotalStep,time,zmax1
     else if (irep>iInitialStep .and. mod(irep,1000)==0) then
         write(*,'("step:",i8,"/",i8,", time:,",f8.2,"s",", m0(10,10):", f8.3,", m4(15,15):", f8.3,", steady:",i8)') &
             irep,iTotalStep,time,ste_jud_m0(1,1),ste_jud_m4(2,1),min_ste
     end if

     !---乱数の発生回数の計算---
     ! !並列 3
     ! do ip=0,iParallel
     !    if(count_rf(ip)>1000000000) then
     !       count_rf(ip) = 0
     !       count_rff(ip) = count_rff(ip) + 1
     !    end if
     ! end do

     !--- 定常判定されてから100000ステップ経過後に計算終了 ---
     if(irep==j_ste+200000) exit

     !--- 粒子のx方向速度最大値zmax1の初期化(iSampleStepステップごとの最大値) ---
     if(mod(irep,iSampleStep)==0) then
        zmax1 = 0
        zmax2 = 0
     end if

     !--- fluxを初期化 ---
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

  !--- データの出力 ---
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
  write(15,*) "dzeta1=",dzeta1,"dzeta2=",dzeta2 !速度の刻み幅出力
  write(15,*) "zmax1=",dzeta1*dble(nzeta1),"zmax2=",dzeta2*dble(nzeta2) !速度の最大値出力
    !並列 4 最後
  ! do ip=0,iParallel
  !   write(15,*) "count_rf",ip,count_rff(ip),count_rf(ip)
  ! end do
  ! close(15)
  
  !--- データの出力 ---
  write(filedataraw,'("DataRaw(Rho1_",f4.2,",Tau1_",f4.2,",pt",i3,").dat")') &
       Rho1,Tau1,iPattern
  open(15,file=filedataraw)
  write(15,*) delta
  write(15,*) dzeta1
  write(15,*) dzeta2 !速度の刻み幅出力
  close(15)

  !--- 散乱体物体における流入流出粒子数の計算 ---
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

  !--- iSampleStep毎に巨視量を求める ---
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

     !--- 定常判定 ---
     !--- 定常判定は，各セルにおける温度について500回分平均をとり， ---
     
     ste_mac_m0(incx,incy) = ste_mac_m0(incx,incy) + m0(incx,incy)
     ste_mac_m4(incx,incy) = ste_mac_m4(incx,incy) + m4(incx,incy)

     !----------------------------------------------------------------------
     ! iSampleStep毎に求めた巨視量を，さらに500回分平均し，
     ! 次の500回分の平均値と比較し定常判定する
     !----------------------------------------------------------------------
     
     if(mod(irep,iSampleStep*500)==0) then
        !--- かなり大きい値にならないように一応500で割っとく ---
        ste_mac_m0(incx,incy) = ste_mac_m0(incx,incy) / 500.0d0
        ste_mac_m4(incx,incy) = ste_mac_m4(incx,incy) / 500.0d0
        !ste_jud = dabs((ste_mac_t-bste_mac_t)/ste_mac_t)
        ste_jud_m0(incx,incy) = (ste_mac_m0(incx,incy)-bste_mac_m0(incx,incy))/ste_mac_m0(incx,incy)
        ste_jud_m4(incx,incy) = (ste_mac_m4(incx,incy)-bste_mac_m4(incx,incy))/ste_mac_m4(incx,incy)

        !--- ste_judが1.4d-2より大きければまだ定常でないとしてste_mac_m0,m4をリセット ---
        !--- ただし，3回連続1.4d-2未満になった場合はそれ以降1.4d-2より大きくなったとしても定常とする ---
        !---------------------------------------------------------------------------------- 
        ! 定常と判定してから10万ステップにおいて巨視量を求めており，その巨視量が0.01以下の
        ! 誤差になってほしいので，5万ステップにおいてはその√2倍の0.014以下になればよいとする．
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

        !-- 3回連続条件を満たせばそのセルのste2を1ずつ増やす --
        if(ste2(incx,incy)<iTotalStep) ste2(incx,incy) = ste2(incx,incy) + 1

        !-- 全セルにおける，ste2の最小値を求める --
        if(ste2(incx,incy)<min_ste) min_ste = ste2(incx,incy)

        !-- 全セルにおける，ste2の最大値を求める --
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

    !--- iSampleStepステップごとのそれぞれの量の合計値を出力 ---
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

    nfun = nfun + 1 !サブルーチンを呼び出す回数をカウント
    
    !---------------------------------------------------------------------------------- 
    ! 直前iSampleStep*500ステップ間におけるzeta1の最大値から速度分布関数の速度刻みを決める
    !----------------------------------------------------------------------------------
    if(nfun==1) then
       dzeta1 = zmax1 / dble(nzeta1)
       dzeta2 = zmax2 / dble(nzeta2)
       
       !--- zmaxが大きすぎた場合 ---
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


    !--- 結果の出力 ---
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

       !--- f1,f2を初期化 ---
       do izeta1 = -nzeta1,nzeta1
          f1(izeta1,incx,incy) = 0.0d0
       end do

       do izeta2 = -nzeta2,nzeta2
          f2(izeta2,incx,incy) = 0.0d0
       end do
    end if

       
  end subroutine dis_fun
 
end program dsmc
