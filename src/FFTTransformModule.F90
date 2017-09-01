 module FFTTransformModule
   use KindParamModule, only : IntKind
   use KindParamModule, only : RealKind
   use KindParamModule, only : CmplxKind
   use MathParamModule, only : ZERO, ONE, CZERO, CONE, TWO, HALF
   use MathParamModule, only : PI, PI4, sqrtm1, TEN2m8
   use ErrorHandlerModule, only : ErrorHandler, StopHandler, WarningHandler
!
#ifdef FFTW
   include 'fftw3.f'
#endif
!
public ::              &
   initFFTTransform,   &
   endFFTTransform,    &
   setFFTStruct,       &
   delFFTStruct,       &
   calFFTTransform,    &
   calFFTRadialProj,   &
   getTransformDirect, &
   calPseudoPotAtSite, &
   FFT_struct
!
   interface initFFTTransform
      module procedure initFFTTransform0, initFFTTransform1
   end interface initFFTTransform
!
   interface calFFTRadialProj
      module procedure calRadialProj,         &
                       calRadialProj_LeastSqChebFit
!                       calRadialProj_LeastSqFit
   end interface calFFTRadialProj
!
private
   integer (kind=IntKind) :: na, nb, nc, nn, out_na
#ifdef FFTW
   integer (kind=8) :: plan
#endif
!
   integer (kind=IntKind) :: nK_cut = 1
!
   real (kind=RealKind) :: rbox(3,3)
   real (kind=RealKind) :: kbox(3,3)
   real (kind=RealKind), allocatable :: kvec(:,:)
!
   real (kind=RealKind), allocatable, target :: fftw_in_r(:)
   real (kind=RealKind), allocatable, target :: fftw_out(:,:,:)
!
   complex (kind=CmplxKind), allocatable, target :: fftw_in_c(:)
!   complex (kind=CmplxKind), allocatable, target :: v_jl(:)
!
   complex (kind=CmplxKind), allocatable, target :: fftw_tmp(:)
!
   logical :: Initialized = .false.
!
   character (len=50) :: stop_routine
!
   integer (kind=IntKind) :: print_level
!
   real (kind=RealKind), allocatable :: wks_kvec(:)
!
   type FFT_struct
      integer (kind=IntKind) :: nga
      integer (kind=IntKind) :: ngb
      integer (kind=IntKind) :: ngc
      integer (kind=IntKind) :: ng
      real (kind=RealKind), pointer :: kvec(:,:)
!
      real (kind=RealKind), pointer :: fft_r(:)
#ifndef FFTW
      complex (kind=CmplxKind), pointer :: ft(:)
#endif
   end type FFT_struct
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initFFTTransform0(nga,ngb,ngc,a0)
!  ===================================================================
!
   implicit none
!
   integer(kind=IntKind), intent(in) :: nga, ngb, ngc
   real(kind=RealKind), intent(in) :: a0
!
   rbox(1:3,1:3) = ZERO
   rbox(1,1) = a0
   rbox(2,2) = a0
   rbox(3,3) = a0
   kbox(1:3,1:3) = ZERO
   kbox(1,1) = ONE/a0
   kbox(2,2) = ONE/a0
   kbox(3,3) = ONE/a0
   na = nga
   nb = ngb
   nc = ngc
   nn = na*nb*nc
!
   end subroutine initFFTTransform0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initFFTTransform1(nga,ngb,ngc,rcell)
!  ===================================================================
!
   implicit none
!
   integer(kind=IntKind), intent(in) :: nga, ngb, ngc
   real(kind=RealKind), intent(in) :: rcell(3,3)
   real(kind=RealKind) :: vol
!
   rbox(1:3,1:3) = rcell(1:3,1:3)  ! Unit cell box
   vol = (rcell(2,1)*rcell(3,2)-rcell(3,1)*rcell(2,2))*rcell(1,3)+    &
         (rcell(3,1)*rcell(1,2)-rcell(1,1)*rcell(3,2))*rcell(2,3)+    &
         (rcell(1,1)*rcell(2,2)-rcell(2,1)*rcell(1,2))*rcell(3,3)
   kbox(1,1) = (rcell(2,2)*rcell(3,3)-rcell(3,2)*rcell(2,3))/vol
   kbox(2,1) = (rcell(3,2)*rcell(1,3)-rcell(1,2)*rcell(3,3))/vol
   kbox(3,1) = (rcell(1,2)*rcell(2,3)-rcell(2,2)*rcell(1,3))/vol
   kbox(1,2) = (rcell(2,3)*rcell(3,1)-rcell(3,3)*rcell(2,1))/vol
   kbox(2,2) = (rcell(3,3)*rcell(1,1)-rcell(1,3)*rcell(3,1))/vol
   kbox(3,2) = (rcell(1,3)*rcell(2,1)-rcell(2,3)*rcell(1,1))/vol
   kbox(1,3) = (rcell(2,1)*rcell(3,2)-rcell(3,1)*rcell(2,2))/vol
   kbox(2,3) = (rcell(3,1)*rcell(1,2)-rcell(1,1)*rcell(3,2))/vol
   kbox(3,3) = (rcell(1,1)*rcell(2,2)-rcell(2,1)*rcell(1,2))/vol
   na = nga
   nb = ngb
   nc = ngc
   nn = na*nb*nc
!
   end subroutine initFFTTransform1
!  ===================================================================
!
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endFFTTransform()
!  ===================================================================
   implicit none
!
#ifdef FFTW
   if ( allocated(fftw_in_r) ) then
      deallocate(fftw_in_r)
   endif
   if ( allocated(fftw_in_c) ) then
      deallocate(fftw_in_c)
   endif
#endif
   if ( allocated(fftw_out) ) then
      deallocate(fftw_out)
   endif
!   if ( allocated(v_jl) ) then
!      deallocate(v_jl)
!   endif
!
   end subroutine endFFTTransform
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setFFTStruct( nga, ngb, ngc, FFT_struct_in )
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nga, ngb, ngc
   type (FFT_struct), intent(inout), target :: FFT_struct_in
!
   FFT_struct_in%nga = nga
   FFT_struct_in%ngb = ngb
   FFT_struct_in%ngc = ngc
   FFT_struct_in%ng  = nga*ngb*ngc
!
!   allocate(FFT_struct_in%kvec(3,nga*ngb*ngc) )
   if ( .not. allocated(wks_kvec) ) then
      allocate(wks_kvec(3*nga*ngb*ngc))
   endif
   FFT_struct_in%kvec => aliasArray2_r(wks_kvec,3,FFT_struct_in%ng)
!
#ifdef FFTW
   allocate( FFT_struct_in%fft_r( (nga+2)*ngb*ngc ) )
#else
   nullify( FFT_struct_in%fft_r )
   allocate( FFT_struct_in%ft( ngb*ngc ) )
#endif
!
   end subroutine setFFTStruct
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine delFFTStruct(FFT_struct_in)
!  ===================================================================
   implicit none
!
   type ( FFT_struct ), intent(inout), target :: FFT_struct_in
!
   FFT_struct_in%nga = -1
   FFT_struct_in%ngb = -1
   FFT_struct_in%ngc = -1
   FFT_struct_in%ng  = -1
!
   nullify( FFT_struct_in%kvec )
   if ( allocated(wks_kvec) ) then
      deallocate( wks_kvec )
   endif
!
#ifdef FFTW
   deallocate( FFT_struct_in%fft_r )
#else
   nullify( FFT_struct_in%fft_r )
   deallocate(FFT_struct_in%ft)
#endif
!
   end subroutine delFFTStruct
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTransformDirect(nga,ngb,ngc,den_in)    result(den_out)
!  ===================================================================
   use DataServiceCenterModule, only : getDataStorage, &
                                       RealMark
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: nga, ngb, ngc
!
   real (kind=RealKind), target  :: den_in(nga,ngb,ngc)
!
   integer (kind=IntKind) :: i, j, k, ng, ind_n
!
   real (kind=RealKind) :: dummy, ksqr
   real (kind=RealKind) :: kpr
   real (kind=RealKind) :: k_nk, k_nj, k_ni
   real (kind=RealKind), pointer :: den_out(:,:,:)
!
   integer (kind=IntKind) :: ir, jr, kr, ind_nr
   real (kind=RealKind)   :: rvec(3), rna, rnb, rnc
!
!  Assign the values to the array
!
   ng = nga*ngb*ngc
   rna = real(nga,kind=RealKind)
   rnb = real(ngb,kind=RealKind)
   rnc = real(ngc,kind=RealKind)
!
   allocate(fftw_tmp(1:ng))
   allocate(fftw_in_c(1:ng))
   if ( .not.allocated(fftw_out) ) then
      allocate( fftw_out(nga,ngb,ngc) )
   else if ( ng > size(fftw_out) ) then
      deallocate(fftw_out)
      allocate( fftw_out(nga,ngb,ngc) )
   endif
!
   ind_n = 0
   dummy = TWO*PI ! Mar 27, 2017 - Yang
   do k = 1, ngc
      k_nk = getKIndex(k,ngc)
      do j = 1, ngb
         k_nj = getKIndex(j,ngb)
         do i = 1, nga
            ind_n = ind_n+1
            k_ni = getKIndex(i,nga)
!
            fftw_in_c(ind_n) = cmplx( den_in(i, j, k),ZERO,kind=CmplxKind)
            kvec(1,ind_n) = dummy*(k_ni*kbox(1,1)+k_nj*kbox(1,2)+k_nk*kbox(1,3))
            kvec(2,ind_n) = dummy*(k_ni*kbox(2,1)+k_nj*kbox(2,2)+k_nk*kbox(2,3))
            kvec(3,ind_n) = dummy*(k_ni*kbox(3,1)+k_nj*kbox(3,2)+k_nk*kbox(3,3))
         enddo
      enddo
   enddo
!
!**    Test for direct sum
!
   fftw_tmp = CZERO

   ind_nr = 0
   do kr = 1, ngc
      do jr = 1, ngb
         do ir = 1, nga
            ind_nr = ind_nr + 1
!
            rvec(1) = (ir-1)*rbox(1,1)/rna+(jr-1)*rbox(1,2)/rnb+(kr-1)*rbox(1,3)/rnc
            rvec(2) = (ir-1)*rbox(2,1)/rna+(jr-1)*rbox(2,2)/rnb+(kr-1)*rbox(2,3)/rnc
            rvec(3) = (ir-1)*rbox(3,1)/rna+(jr-1)*rbox(3,2)/rnb+(kr-1)*rbox(3,3)/rnc
!
            do ind_n = 1, ng
               kpr = kvec(1,ind_n)*rvec(1)+kvec(2,ind_n)*rvec(2)+kvec(3,ind_n)*rvec(3)
               fftw_tmp(ind_n) = fftw_tmp(ind_n)+fftw_in_c(ind_nr)*exp(-sqrtm1*kpr)
            enddo
!
         enddo
      enddo
   enddo

!
   do ind_n = 2, ng
      ksqr = kvec(1,ind_n)**2 + kvec(2,ind_n)**2 + kvec(3,ind_n)**2
      fftw_in_c(ind_n) = fftw_tmp(ind_n)/ksqr
   enddo
!
   fftw_tmp = CZERO
   dummy = TWO*PI4/real(ng,kind=RealKind)
!
   ind_nr = 0
   do kr = 1, ngc
      do jr = 1, ngb
         do ir = 1, nga
            ind_nr = ind_nr + 1
            rvec(1) = (ir-1)*rbox(1,1)/rna+(jr-1)*rbox(1,2)/rnb+(kr-1)*rbox(1,3)/rnc
            rvec(2) = (ir-1)*rbox(2,1)/rna+(jr-1)*rbox(2,2)/rnb+(kr-1)*rbox(2,3)/rnc
            rvec(3) = (ir-1)*rbox(3,1)/rna+(jr-1)*rbox(3,2)/rnb+(kr-1)*rbox(3,3)/rnc

            do ind_n = 2, ng
               kpr = kvec(1,ind_n)*rvec(1)+kvec(2,ind_n)*rvec(2)+kvec(3,ind_n)*rvec(3)
               fftw_tmp(ind_nr) = fftw_tmp(ind_nr)+fftw_in_c(ind_n)*exp(sqrtm1*kpr)
            enddo

            fftw_out(ir,jr,kr) = dummy*real(fftw_tmp(ind_nr),kind=RealKind)
       enddo
      enddo
   enddo
!
   den_out => fftw_out(1:nga,1:ngb,1:ngc)
!
   deallocate(fftw_tmp, fftw_in_c)
!
   end function getTransformDirect
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calFFTTransform( den_in, fft_out, fwd )
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), optional :: fwd
   real (kind=RealKind), target  :: den_in(:,:,:)
   type (FFT_struct), target  :: fft_out
!
   integer (kind=IntKind) :: ng_a, ng_b, ng_c, ng
   integer (kind=IntKind) :: i, j, k, ind_n, k_ni, k_nj, k_nk, ind_m
!
#ifdef TEST_FFT
   integer (kind=IntKind) :: ia, ib, ic, jconj, kconj, ind_k, ind_bc
   real (kind=RealKind) :: rvec(3), rna, rnb, rnc, den_comp
   real (kind=RealKind), allocatable :: den(:,:,:)
   complex (kind=CmplxKind) :: theta_k, expikR
#endif
!
   real (kind=RealKind) :: dummy
   real (kind=RealKind), pointer :: p_kvec(:,:)
   real (kind=RealKind), pointer :: p_fft_r(:)
#ifndef FFTW
   complex (kind=CmplxKind), pointer :: p_ft(:)
#endif
!
   ng_a = fft_out%nga
   ng_b = fft_out%ngb
   ng_c = fft_out%ngc
   ng   = fft_out%ng
!
#ifdef TEST_FFT
!  --------------------------------------------------------------------
!  For the test purpose
   allocate(den(ng_a,ng_b,ng_c))
   den = den_in
!  --------------------------------------------------------------------
#endif
!
   p_kvec => fft_out%kvec(1:3,1:ng)
#ifndef FFTW
   fft_out%fft_r => aliasArray1_r(den_in,ng)
   p_ft => fft_out%ft(1:ng_b*ng_c)
   p_fft_r => fft_out%fft_r(1:ng)
#else
   p_fft_r => fft_out%fft_r(1:(ng_a+2)*ng_b*ng_c)
#endif

   dummy = TWO*PI  ! Mar 27, 2017 -Yang.
   ind_n = 0
   do k = 1, ng_c
      k_nk = getKIndex(k,ng_c)
      do j = 1, ng_b
         k_nj = getKIndex(j,ng_b)
         do i = 1, ng_a
            ind_n = ind_n+1
            ind_m = i + (j-1)*(ng_a+2) + (k-1)*ng_b*(ng_a+2)
            k_ni = getKIndex(i,ng_a)
!
#ifdef FFTW
! *** Complex ***
!            fftw_in_c(ind_n) = cmplx( den_in(i, j, k),ZERO,kind=CmplxKind)
! *** Real in-place ***
            p_fft_r(ind_m) = den_in(i, j, k)
#endif
            p_kvec(1,ind_n) = dummy*(k_ni*kbox(1,1)+k_nj*kbox(1,2)+k_nk*kbox(1,3))
            p_kvec(2,ind_n) = dummy*(k_ni*kbox(2,1)+k_nj*kbox(2,2)+k_nk*kbox(2,3))
            p_kvec(3,ind_n) = dummy*(k_ni*kbox(3,1)+k_nj*kbox(3,2)+k_nk*kbox(3,3))
!
         enddo
      enddo
   enddo
!
#ifdef FFTW
!
! *** Real in-place ***
! Note: We follow the NR convention, where fwd==1 signifies a forward transform,
!       and fwd==-1 signifies an inverse transform.
   if (present(fwd)) then
      if ( fwd==1 ) then
         call dfftw_plan_dft_r2c_3d( plan, ng_a, ng_b, ng_c, p_fft_r, p_fft_r, &
                                     FFTW_ESTIMATE )
      else if ( fwd==-1) then
         call dfftw_plan_dft_c2r_3d( plan, ng_a, ng_b, ng_c, p_fft_r, p_fft_r, &
                                     FFTW_ESTIMATE )
      else
         call ErrorHandler("calFFTTransform","Bad parameter: fwd", fwd)
      endif
   else
      call dfftw_plan_dft_r2c_3d( plan, ng_a, ng_b, ng_c, p_fft_r, p_fft_r, &
                                  FFTW_ESTIMATE )
   endif
! *** Complex ***
!   call dfftw_plan_dft_3d( plan, ng_a, ng_b, ng_c, fftw_in_c, fftw_in_c,        &
!                           FFTW_FORWARD, FFTW_ESTIMATE )
!
   call dfftw_execute( plan )
   call dfftw_destroy_plan( plan )
#else
   if (present(fwd)) then
      if ( fwd==1 .or. fwd==-1 ) then
         call rlft3(p_fft_r,p_ft,ng_a,ng_b,ng_c,fwd)
      else
         call ErrorHandler("calFFTTransform","Bad parameter: fwd", fwd)
      endif
   else
      call rlft3(p_fft_r,p_ft,ng_a,ng_b,ng_c,1)
   endif
#endif
!
#ifdef TEST_FFT
   ia = 2; ib = 3; ic = 4
   rna = real(ng_a,kind=RealKind)
   rnb = real(ng_b,kind=RealKind)
   rnc = real(ng_c,kind=RealKind)
   rvec(1) = (ia-1)*rbox(1,1)/rna+(ib-1)*rbox(1,2)/rnb+(ic-1)*rbox(1,3)/rnc
   rvec(2) = (ia-1)*rbox(2,1)/rna+(ib-1)*rbox(2,2)/rnb+(ic-1)*rbox(2,3)/rnc
   rvec(3) = (ia-1)*rbox(3,1)/rna+(ib-1)*rbox(3,2)/rnb+(ic-1)*rbox(3,3)/rnc
   den_comp = CZERO
   ind_bc = 0
   ind_n = 0
   do k = 1, ng_c
      if (k == 1) then
         kconj = 1
      else
         kconj = ng_c + 2 - k
      endif
!
      do j = 1, ng_b
         if (j == 1) then
            jconj = 1
         else
            jconj = ng_b + 2 - j
         endif
!
         ind_bc = ind_bc + 1
         do i = 1, ng_a
            ind_n = ind_n + 1
            if (i == ng_a/2+1) then
               theta_k = p_ft(ind_bc)
            else if (i < ng_a/2+1) then
               ind_k = 2*i + (j-1)*ng_a + (k-1)*ng_b*ng_a
               theta_k = cmplx(p_fft_r(ind_k-1),p_fft_r(ind_k),kind=CmplxKind)
            else
               ind_k = 2*(ng_a+2-i) + (jconj-1)*ng_a + (kconj-1)*ng_b*ng_a
               theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)
            endif
            expikR  = exp( SQRTm1*(p_kvec(1,ind_n)*rvec(1) + &
                                   p_kvec(2,ind_n)*rvec(2) + &
                                   p_kvec(3,ind_n)*rvec(3)) )
            den_comp = den_comp + real(theta_k*expikR,kind=RealKind)
         enddo
      enddo
   enddo
   den_comp = den_comp/real(ng,kind=RealKind)
   write(6,'(a,d20.10,a,d20.10)')'den, FFT of Theta = ',den(ia,ib,ic),', ',den_comp
   deallocate(den)
#endif
!
   end subroutine calFFTTransform
!  ===================================================================
!
!  *******************************************************************
!
#ifndef DEBUG_BesselProj
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calRadialProj( lmax, n_rmesh, rmesh, posi_in, fft_in, kpow, pv_jl)
!  ===================================================================
   use BesselModule, only : SphericalBessel
   use SphericalHarmonicsModule, only : calYlmConjg
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax, n_rmesh, kpow
!
   real (kind=RealKind), target :: rmesh(n_rmesh)
   real (kind=RealKind), intent(in) :: posi_in(3)
   complex(kind=CmplxKind), target :: pv_jl(:,:)
!
   type (FFT_struct), intent(in) :: fft_in
!
   integer (kind=IntKind) :: i, j, k, ind_n, ng, nr
   integer (kind=IntKind) :: k0, k1, j0, j1, i0, i1, i1n, nstep
   integer (kind=IntKind) :: jmax, kmax, kl, jl, l, m
   integer (kind=IntKind) :: ng_a, ng_b, ng_c
!
   integer (kind=IntKind) :: jconj, kconj, ind_bc, ind_k, n0
!
   real (kind=RealKind) :: dummy, k2, k2inv, posi(3), kfact_r, kfact_i
   real (kind=RealKind), pointer :: p_kvec(:,:)
   real (kind=RealKind), allocatable ::  kri(:), Bj_l(:,:)
!
   real (kind=RealKind), pointer :: p_fft_r(:)
!
#ifndef FFTW
   complex (kind=CmplxKind), pointer :: p_ft(:)
!
#endif
   complex (kind=CmplxKind) :: expikr, kfact, i2l, theta_k
   complex (kind=CmplxKind), allocatable :: Ylm(:)
!
!  Assign local variables
!
!
   ng_a = fft_in%nga
   ng_b = fft_in%ngb
   ng_c = fft_in%ngc
   ng   = fft_in%ng
!
!  write(6,'(a,4i8)') 'calRadialProj:: ', ng_a,ng_b,ng_c,ng
!
   jmax = (lmax+1)*(lmax+2)/2
   kmax = (lmax+1)*(lmax+1)
   allocate( Ylm(kmax) )
   allocate( Bj_l(1:n_rmesh,0:lmax) )
   allocate( kri(n_rmesh) )
!
   posi = posi_in
!
   p_kvec => fft_in%kvec(1:3,1:ng)
#ifndef FFTW
   p_fft_r => fft_in%fft_r(1:ng)
   p_ft => fft_in%ft(1:ng_b*ng_c)
#else
   p_fft_r => fft_in%fft_r(1:(ng_a+2)*ng_b*ng_c)
#endif
!
!  Back FFT to compute the radial jl components of potential
!
   pv_jl(1:n_rmesh,1:jmax) = CZERO
!
   ind_bc = 0
   if ( kpow == -2 ) then
      ind_bc = ng_c*ng_b+1
      ind_n = ng+1
      n0 = 2 ! for Poisson Solver skip the k=0 point
      k0 = ng_c
      k1 = 1
      j0 = ng_b
      j1 = 1
      i0 = ng_a
      i1 = 1
      nstep = -1
   else
      ind_bc = 0
      ind_n = 0
      n0 = 1
      k0 = 1
      k1 = ng_c
      j0 = 1
      j1 = ng_b
      i0 = 1
      i1 = ng_a
      nstep = 1
   endif
   do k = k0,k1,nstep
      if (k == 1) then
         kconj = 1
      else
         kconj = ng_c + 2 - k
      endif
!
      do j = j0,j1,nstep
         i1n = i1
         if (j == 1) then
            jconj = 1
            if ( k==1 .and. kpow==-2 ) i1n=n0
         else
            jconj = ng_b + 2 - j
         endif
!
         ind_bc = ind_bc + nstep
         do i = i0,i1n,nstep
            ind_n = ind_n + nstep
#ifdef FFTW
! *** Complex ***
!           theta_k = fftw_in_c(ind_n)
!
! *** Real in-place ***
            if (i <= ng_a/2+1) then
               ind_k = 2*i + (j-1)*(ng_a+2) + (k-1)*ng_b*(ng_a+2)
    !          theta_k = cmplx(p_fft_r(ind_k-1),p_fft_r(ind_k),kind=CmplxKind)
               theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)
            else
               ind_k = 2*(ng_a+2-i) + (jconj-1)*(ng_a+2) + (kconj-1)*ng_b*(ng_a+2)
    !          theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)
               theta_k = cmplx(p_fft_r(ind_k-1),p_fft_r(ind_k),kind=CmplxKind)
            endif
!
#else
            if (i == ng_a/2+1) then
               theta_k = p_ft(ind_bc)
            else if (i < ng_a/2+1) then
               ind_k = 2*i + (j-1)*ng_a + (k-1)*ng_b*ng_a
    !          theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)
               theta_k = cmplx(p_fft_r(ind_k-1),p_fft_r(ind_k),kind=CmplxKind)
            else
               ind_k = 2*(ng_a+2-i) + (jconj-1)*ng_a + (kconj-1)*ng_b*ng_a
    !          theta_k = cmplx(p_fft_r(ind_k-1),p_fft_r(ind_k),kind=CmplxKind)
               theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)
            endif
#endif
!
            k2 = p_kvec(1,ind_n)**2 + p_kvec(2,ind_n)**2 + p_kvec(3,ind_n)**2
            expikR  = exp( -( p_kvec(1,ind_n)*posi(1) + &
                              p_kvec(2,ind_n)*posi(2) + &
                              p_kvec(3,ind_n)*posi(3) )*sqrtm1 )
            do nr = 1,n_rmesh
               kri(nr) = sqrt(k2)*rmesh(nr)
            enddo
!
            if ( kpow==-2 ) then
               kfact = (expikR*theta_k)/k2
            else
               kfact = (expikR*theta_k)
            endif
            kfact_r = real(kfact,kind=RealKind)
            kfact_i = -real(sqrtm1*kfact,kind=RealKind)
!
            call calYlmConjg( p_kvec(1:3,ind_n), lmax, Ylm )
            call SphericalBessel( lmax, n_rmesh, kri(1:n_rmesh), &
                                     Bj_L(1:n_rmesh,0:lmax) )
!
            do l = 0,lmax
               if (mod(l,2) == 0) then
                  do m = 0,l
                     kl = (l+1)*(l+1) -l + m
                     jl = ((l+1)*(l+2))/2 -l + m
                     pv_jl(1:n_rmesh,jl) = pv_jl(1:n_rmesh,jl)+ &
                                           kfact_r*Ylm(kl)*Bj_L(1:n_rmesh,l)
                  enddo
               else
                  do m = 0,l
                     kl = (l+1)*(l+1) -l + m
                     jl = ((l+1)*(l+2))/2 -l + m
                     pv_jl(1:n_rmesh,jl) = pv_jl(1:n_rmesh,jl)+ &
                                           sqrtm1*kfact_i*Ylm(kl)*Bj_L(1:n_rmesh,l)
                  enddo
               endif
            enddo
!
         enddo
      enddo
   enddo
!
   if ( kpow == -2) then
      dummy =TWO*PI4*PI4/real(ng,kind=RealKind)
   else
      dummy = PI4/real(ng,kind=RealKind)
   endif
   i2l = CONE
   do l = 0,lmax
      do m = 0,l
         jl = ((l+1)*(l+2))/2 - l + m
         pv_jl(1:n_rmesh,jl) = dummy*i2l*pv_jl(1:n_rmesh,jl)
      enddo
      i2l = i2l*(-sqrtm1)
   enddo
!
   deallocate( Ylm, Bj_L, kri )
!
   end subroutine calRadialProj
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calRadialProj_LeastSqFit( lmax, n_rmesh, rmesh, posi_in,    &
                                     fft_in, ninterp, pv_jl )
!  ===================================================================
   use IntegerFactorsModule, only : lofj
   use BesselModule, only : SphericalBessel
   use SphericalHarmonicsModule, only : calYlm
   use InterpolationModule, only : LeastSqFitInterp
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax, n_rmesh, ninterp
!
   real (kind=RealKind), intent(in) :: posi_in(3)
   real (kind=RealKind), target :: rmesh(n_rmesh)
!
   complex(kind=CmplxKind), target :: pv_jl(:,:)
!
   type (FFT_struct), intent(in) :: fft_in
!
   integer (kind=IntKind) :: i, j, k, ind_n, ng, nr
   integer (kind=IntKind) :: jmax, kmax, kl, jl, l, m
   integer (kind=IntKind) :: ng_a, ng_b, ng_c
!
   integer (kind=IntKind) :: jconj, kconj, ind_bc, ind_k, n0
!
   real (kind=RealKind) :: dummy, k2, k2inv, posi(3)
   real (kind=RealKind), pointer ::  p_kvec(:,:)
   real (kind=RealKind), allocatable ::  kri(:), Bj_l(:,:)
!
   real (kind=RealKind), pointer :: p_fft_r(:)
!
#ifndef FFTW
   complex (kind=CmplxKind), pointer :: p_ft(:)
!
#endif
   complex (kind=CmplxKind) :: expikr, kfact, i2l, theta_k
   complex(kind=CmplxKind), pointer :: pv_interp(:,:)
   complex (kind=CmplxKind), allocatable :: Ylm(:)
!
!   complex (kind=CmplxKind), pointer :: den_out(:,:)
!
!  Interpolation Data
!
   integer (kind=IntKind) :: n_interp
   real (kind=RealKind) ::  d_pv,d_r
   real (kind=RealKind), allocatable :: r_interp(:)
   complex (kind=CmplxKind), allocatable, target :: v_interp(:)
!
!  Assign local variables
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
   ng_a = fft_in%nga
   ng_b = fft_in%ngb
   ng_c = fft_in%ngc
   ng   = fft_in%ng
!
!  write(6,'(a,4i8)') 'calRadialProj:: ', ng_a,ng_b,ng_c,ng
!
   jmax = (lmax+1)*(lmax+2)/2
   kmax = (lmax+1)*(lmax+1)
!
   n_interp = ninterp
!
   allocate( Ylm(kmax) )
   allocate( Bj_l(1:n_interp,0:lmax) )
   allocate( kri(n_interp) )
!
   allocate( r_interp(n_interp) )
   allocate( v_interp(1:n_interp*jmax) )
!
   pv_interp => aliasArray2_c(v_interp,n_interp,jmax)
!
   posi = posi_in
!
   p_kvec => fft_in%kvec(1:3,1:ng)
#ifndef FFTW
   p_fft_r => fft_in%fft_r(1:ng)
   p_ft => fft_in%ft(1:ng_b*ng_c)
#else
   p_fft_r => fft_in%fft_r(1:(ng_a+2)*ng_b*ng_c)
#endif
!
!  Back FFT to compute the radial jl components of potential
!
   pv_jl(1:n_rmesh,1:jmax) = CZERO
   pv_interp(1:n_interp,1:jmax) = CZERO
   j = 1
   d_r = (rmesh(n_rmesh) - 0.80d0)/real(n_interp,kind=RealKind)
   d_pv = 0.80d0
   Loop_I: do i = 1, n_interp
      r_interp(i) = d_pv
      d_pv = d_pv + d_r
   enddo Loop_I
   r_interp(n_interp) = rmesh(n_rmesh)
!
   ind_n = 1
   ind_bc = 0
   n0    = 2
   do k = 1,ng_c
      if (k == 1) then
         kconj = 1
      else
         kconj = ng_c + 2 - k
      endif
!
      do j = 1,ng_b
         if (j == 1) then
            jconj = 1
         else
            jconj = ng_b + 2 - j
         endif
!
         ind_bc = ind_bc + 1
         do i = n0,ng_a
            ind_n = ind_n + 1
#ifdef FFTW
! *** Complex ***
!           theta_k = fftw_in_c(ind_n)
!
! *** Real in-place ***
            if (i <= ng_a/2+1) then
               ind_k = 2*i + (j-1)*(ng_a+2) + (k-1)*ng_b*(ng_a+2)
    !          theta_k = cmplx(p_fft_r(ind_k-1),p_fft_r(ind_k),kind=CmplxKind)
               theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)
            else
               ind_k = 2*(ng_a+2-i) + (jconj-1)*(ng_a+2) + (kconj-1)*ng_b*(ng_a+2)
    !          theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)
               theta_k = cmplx(p_fft_r(ind_k-1),p_fft_r(ind_k),kind=CmplxKind)
            endif
!
#else
            if (i == ng_a/2+1) then
               theta_k = p_ft(ind_bc)
            else if (i < ng_a/2+1) then
               ind_k = 2*i + (j-1)*ng_a + (k-1)*ng_b*ng_a
    !          theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)
               theta_k = cmplx(p_fft_r(ind_k-1),p_fft_r(ind_k),kind=CmplxKind)
            else
               ind_k = 2*(ng_a+2-i) + (jconj-1)*ng_a + (kconj-1)*ng_b*ng_a
    !          theta_k = cmplx(p_fft_r(ind_k-1),p_fft_r(ind_k),kind=CmplxKind)
               theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)
            endif
#endif
!
      k2    = ( p_kvec(1,ind_n)**2 + p_kvec(2,ind_n)**2 + p_kvec(3,ind_n)**2)
      k2inv = ONE/k2
      expikR  = exp( -( p_kvec(1,ind_n)*posi(1) + &
                        p_kvec(2,ind_n)*posi(2) + &
                        p_kvec(3,ind_n)*posi(3) )*sqrtm1 )
      do nr = 1,n_interp
         kri(nr) = sqrt(k2)*r_interp(nr)
      enddo
!
      kfact = expikR*k2inv*theta_k
!
      call calYlm( p_kvec(1:3,ind_n), lmax, Ylm )
      call SphericalBessel( lmax, n_interp, kri(1:n_interp), Bj_L(1:n_interp,0:lmax) )
!
      do l = 0,lmax
         do m = 0,l
            kl = (l+1)*(l+1) -l + m
            jl = (l+1)*(l+2)/2 -l + m
            pv_interp(1:n_interp,jl) = pv_interp(1:n_interp,jl)+ &
                                  kfact*Ylm(kl)*Bj_L(1:n_interp,l)
         enddo
      enddo
!
         enddo
         n0 = 1
      enddo
   enddo
!
   dummy =TWO*PI4*PI4/real(ng,kind=RealKind)
   i2l = CONE
   do l = 0,lmax
      do m = 0,l
         jl = (l+1)*(l+2)/2 - l + m
         do i =1,n_interp
            pv_interp(i,jl) = dummy*i2l*pv_interp(i,jl)/(r_interp(i)**l)
         enddo
      enddo
      i2l = i2l*(-sqrtm1)
   enddo
!
   LoopI: do i = 1, n_rmesh
      do jl = 1,jmax
         l = lofj(jl)
         call LeastSqFitInterp( n_interp, r_interp, &
              pv_interp(1:n_interp,jl), rmesh(i), pv_jl(i,jl))
         pv_jl(i,jl) = pv_jl(i,jl)*(rmesh(i)**l)
      enddo
   enddo LOOPI
!
   nullify( pv_interp )
   deallocate( Ylm, Bj_L, kri )
   deallocate( r_interp, v_interp )
!
   end subroutine calRadialProj_LeastSqFit
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  This is the subroutine actually used.... Yang @ Mar 25, 2017.......
   subroutine calRadialProj_LeastSqChebFit( lmax, n_rmesh, rmesh, rmt_in, &
                            posi_in, fft_in, kpow, ninterp, pv_jl )
!  ===================================================================
   use IntegerFactorsModule, only : lofj
   use BesselModule, only : SphericalBessel
   use SphericalHarmonicsModule, only : calYlmConjg
   use ChebyshevModule, only : ChebyshevStruct,     &
                               initChebyshevSeries, &
                               ChebyshevEval,       &
                               calChebGrid,         &
                               endChebyshevSeries
   use InterpolationModule, only : LeastSqFitInterp, &
                                   FitInterp
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax, n_rmesh, kpow, ninterp
!
   real (kind=RealKind), intent(in) :: posi_in(3)
   real (kind=RealKind), target :: rmesh(n_rmesh)
   real (kind=RealKind), optional :: rmt_in
   complex(kind=CmplxKind), target :: pv_jl(:,:)
!
   type (FFT_struct), intent(in) :: fft_in
!
   integer (kind=IntKind) :: i, j, k, ind_n, ng, nr
   integer (kind=IntKind) :: jmax, kmax, kl, jl, l, m
   integer (kind=IntKind) :: ng_a, ng_b, ng_c
!
   integer (kind=IntKind) :: jconj, kconj, ind_bc, ind_k, n0
!
   real (kind=RealKind) :: dummy, k2, k2pow, posi(3), rmt, p_vec(3), kfact_r, kfact_i
   real (kind=RealKind), pointer ::  p_kvec(:,:)
   real (kind=RealKind), allocatable ::  kri(:), Bj_l(:,:)
!
   real (kind=RealKind), pointer :: p_fft_r(:)
!
#ifndef FFTW
   complex (kind=CmplxKind), pointer :: p_ft(:)
!
#endif
   complex(kind=CmplxKind) :: expikr, kfact, i2l, theta_k
   complex(kind=CmplxKind), pointer :: pv_interp(:,:)
   complex(kind=CmplxKind), pointer :: pv_interp0(:,:)
   complex(kind=CmplxKind), pointer :: pv_interp1(:,:)
   complex(kind=CmplxKind), allocatable :: Ylm(:)
!
!   complex (kind=CmplxKind), pointer :: den_out(:,:)
!
!  Interpolation Data
!
   integer (kind=IntKind) ::  k0, k1, i0, i1, i1n, j0, j1, nstep
   integer (kind=IntKind) :: n_sz
   integer (kind=IntKind) :: n_interp, ir_lsq, ir_lsqm1, n_interp_in, nr_int
   integer (kind=IntKind), allocatable :: ind_interp(:), d_ir(:)
   real (kind=RealKind) :: r0
   real (kind=RealKind), allocatable :: r_interp(:), dr_int(:)
   complex(kind=CmplxKind), pointer :: pv(:), pv0(:)
   complex(kind=CmplxKind), allocatable, target :: v_interp(:)
   complex(kind=CmplxKind), allocatable, target :: v_interp0(:)
   complex(kind=CmplxKind), allocatable, target :: v_interp1(:)
   type (ChebyshevStruct), allocatable, target :: chebv_struct(:)
!
!  Assign local variables
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
!  print *, " calRadialProj_LeastSqChebFit:: Enter "
!
   ng_a = fft_in%nga
   ng_b = fft_in%ngb
   ng_c = fft_in%ngc
   ng   = fft_in%ng
!
   jmax = (lmax+1)*(lmax+2)/2
   kmax = (lmax+1)*(lmax+1)
!
!  With increase of l_max, Gibs fenomena occurs at larger r.
!  For different values of l_max define ranges to apply the Chebyshev fit,
!  and use a linear square fit( order 3 polynomial )
!  to extrapolate the Chebyshev fit to small r
!
   n_interp  = ninterp
   if ( present(rmt_in) ) then
      rmt = rmt_in
   else
      rmt = TWO
   endif
   if ( lmax>8) then
      nr_int = 3
      allocate( dr_int(nr_int), d_ir(nr_int) )
      dr_int(1) = rmesh(1)
      dr_int(2) = HALF*rmt
      dr_int(3) = rmt
   else if ( lmax>4 ) then
      nr_int = 2
      allocate( dr_int(nr_int), d_ir(nr_int) )
      dr_int(1) = rmesh(1)
      dr_int(2) = HALF*rmt
   else
      nr_int = 1
      allocate( dr_int(nr_int), d_ir(nr_int) )
      dr_int(1) = rmesh(1)
   endif
!
   allocate( Ylm(kmax) )
   allocate( Bj_l(1:n_interp*nr_int,0:lmax) )
   allocate( kri(n_interp*nr_int) )
   allocate( ind_interp(1:n_interp*nr_int) )
   allocate( r_interp(n_interp*nr_int) )
   allocate( v_interp(1:n_interp*nr_int*jmax) )
   allocate( v_interp0(1:n_interp*nr_int*jmax) )
   allocate( v_interp1(1:n_interp*nr_int*jmax) )
   allocate( chebv_struct(jmax) )
!
   do i = 1,nr_int
      r0 = dr_int(i)
      call hunt(n_rmesh,rmesh,r0,ir_lsq)
      d_ir(i) = ir_lsq
      dr_int(i) = rmesh(ir_lsq)
      call calChebGrid( rmesh(ir_lsq), rmesh(n_rmesh), n_interp, &
                        r_interp( (i-1)*n_interp+1:i*n_interp ) )
      if ( n_rmesh-ir_lsq+1 <= 10 ) then
         n_interp_in = n_rmesh-ir_lsq+1
      else
         n_interp_in = 10
      endif
   enddo
!
   n_sz = n_interp*nr_int*jmax
   pv_interp  => aliasArray2_c(v_interp,n_interp*nr_int,jmax)
   pv_interp0 => aliasArray2_c(v_interp0,n_interp*nr_int,jmax)
   pv_interp1 => aliasArray2_c(v_interp1,n_interp*nr_int,jmax)
!
   posi = posi_in
!
   p_kvec => fft_in%kvec(1:3,1:ng)
#ifndef FFTW
   p_fft_r => fft_in%fft_r(1:ng)
   p_ft => fft_in%ft(1:ng_b*ng_c)
#else
   p_fft_r => fft_in%fft_r(1:(ng_a+2)*ng_b*ng_c)
#endif
!
!  Back FFT to compute the radial jl components of potential
!
   pv_jl(1:n_rmesh,1:jmax) = CZERO
   pv_interp(1:n_interp*nr_int,1:jmax) = CZERO
!
!  print *, " calRadialProj_LeastSqChebFit:: Init done "
!
   if ( kpow == -2 ) then
      ind_bc = ng_c*ng_b+1
      ind_n = ng+1
      n0 = 2 ! for Poisson Solver skip the k=0 point
      k0 = ng_c
      k1 = 1
      j0 = ng_b
      j1 = 1
      i0 = ng_a
      i1 = 1
      nstep = -1
   else
      ind_bc = 0
      ind_n = 0
      n0 = 1
      k0 = 1
      k1 = ng_c
      j0 = 1
      j1 = ng_b
      i0 = 1
      i1 = ng_a
      nstep = 1
   endif
!
   do k = k0,k1,nstep
      if (k == 1) then
         kconj = 1
      else
         kconj = ng_c + 2 - k
      endif
      pv_interp1(1:n_interp*nr_int,1:jmax) = CZERO
      do j = j0,j1,nstep
         i1n = i1
         if (j == 1) then
            jconj = 1
            if ( k==1 .and. kpow==-2 ) i1n=n0
         else
            jconj = ng_b + 2 - j
         endif
!
         pv_interp0(1:n_interp*nr_int,1:jmax) = CZERO
         ind_bc = ind_bc + nstep
         do i = i0,i1n,nstep
            ind_n = ind_n + nstep
#ifdef FFTW
! *** Complex ***
!           theta_k = fftw_in_c(ind_n)
!
! *** Real in-place ***
            if (i <= ng_a/2+1) then
               ind_k = 2*i + (j-1)*(ng_a+2) + (k-1)*ng_b*(ng_a+2)
               theta_k = cmplx(p_fft_r(ind_k-1),p_fft_r(ind_k),kind=CmplxKind)
    !          theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)  ! Mar 27, 2017 - Yang
            else
               ind_k = 2*(ng_a+2-i) + (jconj-1)*(ng_a+2) + (kconj-1)*ng_b*(ng_a+2)
               theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)
    !          theta_k = cmplx(p_fft_r(ind_k-1),p_fft_r(ind_k),kind=CmplxKind)   ! Mar 27, 2017 - Yang
            endif
!
#else
            if (i == ng_a/2+1) then
               theta_k = p_ft(ind_bc)
            else if (i < ng_a/2+1) then
               ind_k = 2*i + (j-1)*ng_a + (k-1)*ng_b*ng_a
               theta_k = cmplx(p_fft_r(ind_k-1),p_fft_r(ind_k),kind=CmplxKind)
    !          theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)
            else
               ind_k = 2*(ng_a+2-i) + (jconj-1)*ng_a + (kconj-1)*ng_b*ng_a
               theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)
    !          theta_k = cmplx(p_fft_r(ind_k-1),p_fft_r(ind_k),kind=CmplxKind)
            endif
#endif
!
            k2    = p_kvec(1,ind_n)**2 + p_kvec(2,ind_n)**2 + p_kvec(3,ind_n)**2
   !        expikR  = exp( -( p_kvec(1,ind_n)*posi(1) + & ! Mar 27, 2017 - Yang
            expikR  = exp(  ( p_kvec(1,ind_n)*posi(1) + &
                              p_kvec(2,ind_n)*posi(2) + &
                              p_kvec(3,ind_n)*posi(3) )*sqrtm1 )
            do nr = 1,n_interp*nr_int
               kri(nr) = sqrt(k2)*r_interp(nr)
            enddo
!
            if (kpow==-2) then
               kfact = (expikR*theta_k)/k2
            else
               kfact = expikR*theta_k
            endif
            kfact_r = real(kfact,kind=RealKind)
            kfact_i = -real(sqrtm1*kfact,kind=RealKind)
!
            p_vec(1:3) = p_kvec(1:3,ind_n) !/sqrt(k2)
            call calYlmConjg( p_vec, lmax, Ylm )
            call SphericalBessel( lmax, n_interp*nr_int, kri(1:n_interp*nr_int), Bj_L(1:n_interp*nr_int,0:lmax) )
!
            do l = 0,lmax
               if (mod(l,2) == 0) then
                  do m = 0,l
                     kl = (l+1)*(l+1) -l + m
                     jl = ((l+1)*(l+2))/2 -l + m
                     pv_interp0(1:n_interp*nr_int,jl) = pv_interp0(1:n_interp*nr_int,jl) + &
                                                        kfact_r*Ylm(kl)*Bj_L(1:n_interp*nr_int,l)
                  enddo
               else
                  do m = 0,l
                     kl = (l+1)*(l+1) -l + m
                     jl = ((l+1)*(l+2))/2 -l + m
                     pv_interp0(1:n_interp*nr_int,jl) = pv_interp0(1:n_interp*nr_int,jl) + &
                                                        sqrtm1*kfact_i*Ylm(kl)*Bj_L(1:n_interp*nr_int,l)
                  enddo
               endif
            enddo
         enddo
         pv_interp1(1:n_interp*nr_int,1:jmax) = pv_interp1(1:n_interp*nr_int,1:jmax) + &
                                                pv_interp0(1:n_interp*nr_int,1:jmax)
      enddo
      pv_interp(1:n_interp*nr_int,1:jmax) = pv_interp(1:n_interp*nr_int,1:jmax) + &
                                            pv_interp1(1:n_interp*nr_int,1:jmax)
   enddo
!
!  print *, " calRadialProj_LeastSqChebFit:: Loop K "
!
   if ( kpow == -2) then
      dummy = TWO*PI4*PI4/real(ng,kind=RealKind)
   else
      dummy = PI4/real(ng,kind=RealKind)
   endif
   i2l = dummy
   do l = 0,lmax
      do m = 0,l
         jl = ((l+1)*(l+2))/2 - l + m
         do i =1,n_interp*nr_int
            pv_interp(i,jl) = i2l*pv_interp(i,jl)/(r_interp(i)**l)
         enddo
      enddo
!     i2l = i2l*(-sqrtm1)  ! Mar 27, 2017 - Yang
      i2l = i2l*sqrtm1
   enddo
!
   do jl = 1,jmax
      l = lofj(jl)
      if ( l<=4 ) then
         pv => pv_interp(1:n_interp,jl)
         r0 = rmesh(d_ir(1))
         ir_lsq = d_ir(1)
      else if ( l>4 .and. l<=8 ) then
         pv => pv_interp(n_interp+1:2*n_interp,jl)
         r0 = rmesh(d_ir(2))
         ir_lsq = d_ir(2)
      else
         pv => pv_interp(2*n_interp+1:3*n_interp,jl)
         r0 = rmesh(d_ir(3))
         ir_lsq = d_ir(3)
      endif
      call initChebyshevSeries( n_interp, r0, rmesh(n_rmesh), &
                               pv, chebv_struct(jl) )
      pv => pv_jl(ir_lsq:n_rmesh,jl)
      call ChebyshevEval( chebv_struct(jl), n_interp, n_rmesh-ir_lsq+1, &
                          rmesh(ir_lsq:n_rmesh), pv )
      call endChebyshevSeries(chebv_struct(jl))
   enddo
!
!  print *, " calRadialProj_LeastSqChebFit:: Chebyshev Factors "
!
   if ( nr_int>1 ) then
      do jl = 1,jmax
         l = lofj(jl)
         if ( l<=4 ) then
            ir_lsq = d_ir(1)
            ir_lsqm1 = ir_lsq - 1
         else if ( l>4 .and. l<=8 ) then
            ir_lsq = d_ir(2)
            ir_lsqm1 = ir_lsq - 1
         else
            ir_lsq = d_ir(3)
            ir_lsqm1 = ir_lsq - 1
         endif
         pv => pv_jl(1:ir_lsqm1,jl)
         pv0 => pv_jl(ir_lsq:ir_lsq+n_interp_in-1,jl)
         call LeastSqFitInterp( n_interp_in, rmesh(ir_lsq:n_rmesh), &
                                pv0, ir_lsqm1, rmesh(1:ir_lsqm1), pv )
         do i = 1, n_rmesh
            pv_jl(i,jl) = pv_jl(i,jl)*(rmesh(i)**l)
         enddo
      enddo
   else
      do jl = 1,jmax
         l = lofj(jl)
         do i = 1, n_rmesh
            pv_jl(i,jl) = pv_jl(i,jl)*(rmesh(i)**l)
         enddo
      enddo
   endif
!
!  print *, " calRadialProj_LeastSqChebFit:: Deallocate "
!
   nullify( pv_interp, pv_interp0, pv_interp1 )
   deallocate( Ylm, Bj_L, kri )
   deallocate( ind_interp, r_interp, v_interp, v_interp0, v_interp1 )
   deallocate( chebv_struct )
   deallocate( dr_int, d_ir )
!
!  print *, " calRadialProj_LeastSqChebFit:: Exit "
!
   end subroutine calRadialProj_LeastSqChebFit
#else
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calRadialProj( lmax, n_rmesh, rmesh, posi_in, fft_in, pv_jl)
!  ===================================================================
   use IntegerFactorsModule, only : lofj, mofj
   use BesselModule, only : SphericalBessel
   use SphericalHarmonicsModule, only : calYlm
   use MPPModule, only : MyPE, syncAllPEs
#ifdef TIMING
   use TimerModule, only : getTime
#endif
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax, n_rmesh
!
   real (kind=RealKind), intent(in) :: posi_in(3)
   real (kind=RealKind), target :: rmesh(n_rmesh)
!
   type (FFT_struct), intent(in) :: fft_in
!
   integer (kind=IntKind) :: i, j, k, ind_n, ng, nr
   integer (kind=IntKind) :: ir
   integer (kind=IntKind) :: jmax, kmax, kl, jl, l, m
   integer (kind=IntKind) :: ng_a, ng_b, ng_c
!
   integer (kind=IntKind) :: jconj, kconj, ind_bc, ind_k, n0
   integer (kind=IntKind), allocatable :: ind_k_ord(:), ind_n_ord(:)
   integer (kind=IntKind), allocatable :: sign_Im_k(:)
!
   real (kind=RealKind) :: dummy, k2, k2inv, k2_prev, k2_max, posi(3)
   real (kind=RealKind) :: k_nk, k_nj, k_ni
   real (kind=RealKind), pointer ::  p_kvec(:,:)
   real (kind=RealKind), allocatable ::  kri(:), Bj_l(:,:), k2_ord(:)
!
   real (kind=RealKind), pointer :: p_fft_r(:)
!
#ifndef FFTW
   complex (kind=CmplxKind), pointer :: p_ft(:)
!
#endif
   complex (kind=CmplxKind) :: expikr, kfact, i2l, theta_k
   complex(kind=CmplxKind), pointer :: pv_jl(:,:)
   complex (kind=CmplxKind), allocatable :: Ylm(:)
!
!   complex (kind=CmplxKind), pointer :: den_out(:,:)
!
#ifdef TIMING
   real (kind=RealKind) :: t0, t1
#endif
!
!  Bessel debug
!
   character(len=9), allocatable :: fname(:)
   character(len=9) :: fname_tmp
   character(len=7), allocatable :: l_name(:)
!
   integer (kind=IntKind) :: nr_print, njl_print, offset_l = 100
   integer (kind=IntKind), allocatable :: rind_print(:), jl_print(:)
   integer (kind=IntKind), allocatable :: funit(:)
!
   real (kind=RealKind), allocatable :: r_print(:)
!
   complex (kind=CmplxKind), allocatable :: kfact_y(:)
!
!  Assign local variables
!
!
   ng_a = fft_in%nga
   ng_b = fft_in%ngb
   ng_c = fft_in%ngc
   ng   = fft_in%ng
!
!  write(6,'(a,4i8)') 'calRadialProj:: ', ng_a,ng_b,ng_c,ng
!
   jmax = (lmax+1)*(lmax+2)/2
   kmax = (lmax+1)*(lmax+1)
   allocate( Ylm(kmax) )
   allocate( Bj_l(1:n_rmesh,0:lmax) )
   allocate( kri(n_rmesh) )
!
!   if ( .not.allocated(v_jl) ) then
!      allocate(v_jl(1:n_rmesh*jmax))
!   else if ( size(v_jl)/=n_rmesh*jmax )  then
!      deallocate(v_jl)
!      allocate(v_jl(1:n_rmesh*jmax))
!   endif
!
   allocate( ind_k_ord(ng), ind_n_ord(ng), sign_Im_k(ng) )
   allocate( k2_ord(ng) )
!
!   pv_jl => aliasArray2_c(v_jl,n_rmesh,jmax)
!
   posi = posi_in
!
   kvec => fft_in%kvec(1:3,1:ng)
#ifndef FFTW
   p_fft_r => fft_in%fft_r(1:ng)
   p_ft => fft_in%ft(1:ng_b*ng_c)
#else
   p_fft_r => fft_in%fft_r(1:(ng_a+2)*ng_b*ng_c)
#endif
!
!
!  Bessel debug
!
   njl_print = 4
   nr_print = 5
!
   allocate( rind_print(nr_print), jl_print(njl_print) )
   allocate( funit(nr_print), fname(nr_print) )
   allocate( r_print(nr_print), kfact_y(njl_print) )
   allocate( l_name(njl_print) )
!
   jl_print(1) = 1
   jl_print(2) = (4+1)*(4+2)/2-4
   jl_print(3) = (6+1)*(6+2)/2-6
   jl_print(4) = (8+1)*(8+2)/2-8
   do i = 1,njl_print
      write(fname_tmp,'(a1,i3,i3)') "L",100+lofj(jl_print(i)),100+mofj(jl_print(i))
      fname_tmp(2:2) = '_'
      fname_tmp(5:5) = '_'
      l_name(i) = fname_tmp(1:7)
   enddo
!
   r_print(1) = rmesh(1)
   r_print(2) = 0.5d0
   r_print(3) = 1.0d0
   r_print(4) = 2.0d0
   r_print(5) = 2.50d0
!
   rind_print(1) = 1
   do i = 2,nr_print
!     -------------------------------------------------------------
      call hunt(n_rmesh,rmesh(1:n_rmesh),r_print(i), ir )
!     -------------------------------------------------------------
      rind_print(i) = ir
   enddo
!
   if (MyPE==1) then
   do i = 1,nr_print
      funit(i) = 20+i
      write(fname(i),'(a6,i3)') "Bessel",100+i
      fname_tmp = fname(i)
      fname_tmp(7:7) = '_'
      fname(i) = fname_tmp
      open(unit=funit(i),file=trim(fname(i)),status='unknown')
      write(funit(i),'(a27,d16.8)') '# Bessel functions for r = ', rmesh(rind_print(i))
      write(funit(i),'(a22,a23,$)') "#         |r|         ","           |k|         "
      do j = 1, njl_print
         write(funit(i),'(5x,a6,a4,10x,a9,a7,6x,a9,a7,2x,a14,a7,a1,1x,a14,a7,a1,$)') &
              "j(kr)_",l_name(j), &
              "Re[F(k)]_",l_name(j),"Im[F(k)]_",l_name(j), &
              "Sum{Re[W(kr)]_",l_name(j),"}","Sum{Im[W(kr)]_",l_name(j),"}"
      enddo
      write(funit(i),'(a1)') " "
   enddo
   endif
!
   p_kvec => fft_in%kvec(1:3,1:ng)
!  dummy = TWO*PI
   dummy = -TWO*PI  ! since we are using "forward" FFT that has a "-" sign in the exponent
   k2_max = ZERO
   ind_n = 0
   do k = 1, ng_c
      k_nk = getKIndex(k,ng_c)
      do j = 1, ng_b
         k_nj = getKIndex(j,ng_b)
         do i = 1, ng_a
            ind_n = ind_n+1
            k_ni = getKIndex(i,ng_a)
!
#ifdef FFTW
! *** Complex ***
!            fftw_in_c(ind_n) = cmplx( den_in(i, j, k),ZERO,kind=CmplxKind)
! *** Real in-place ***
            fftw_in_r(ind_n) = den_in(i, j, k)
#endif
            p_kvec(1,ind_n) = dummy*(k_ni*kbox(1,1)+k_nj*kbox(1,2)+k_nk*kbox(1,3))
            p_kvec(2,ind_n) = dummy*(k_ni*kbox(2,1)+k_nj*kbox(2,2)+k_nk*kbox(2,3))
            p_kvec(3,ind_n) = dummy*(k_ni*kbox(3,1)+k_nj*kbox(3,2)+k_nk*kbox(3,3))
            k2 = (p_kvec(1,ind_n)**2 + p_kvec(2,ind_n)**2 + p_kvec(3,ind_n)**2)
            k2_max = max(k2_max,k2)
!
         enddo
      enddo
   enddo
!
#ifdef TIMING
   t0 = getTime()
#endif
!
   ind_n_ord = -1
   ind_n_ord(1) = 1
   k2_ord(1) = ZERO
   ind_n = 1
   ind_bc = 0
   n0    = 2
   do k = 1,ng_c
      if (k /= 1) then
         kconj = ng_c + 2 - k
      else
         kconj = 1
      endif
!
      do j = 1,ng_b
         if (j /= 1) then
            jconj = ng_b + 2 - j
         else
            jconj = 1
         endif
!
         ind_bc = ind_bc + 1
!
         do i = n0,ng_a/2
            ind_n = ind_n + 1
            ind_n_ord(ind_n) = ind_n
            k2_ord(ind_n) = (p_kvec(1,ind_n)**2 + p_kvec(2,ind_n)**2 + &
                             p_kvec(3,ind_n)**2)
#ifdef FFTW
            ind_k = 2*i + (j-1)*(ng_a+2) + (k-1)*ng_b*(ng_a+2)
#else
            ind_k = 2*i + (j-1)*ng_a + (k-1)*ng_b*ng_a
#endif
            sign_Im_k(ind_n) = 1
            ind_k_ord(ind_n) = ind_k
         enddo
         i = ng_a/2+1
         ind_n = ind_n + 1
         ind_n_ord(ind_n) = ind_n
         k2_ord(ind_n) = (p_kvec(1,ind_n)**2 + p_kvec(2,ind_n)**2 + &
                          p_kvec(3,ind_n)**2)
#ifdef FFTW
         ind_k = 2*i + (j-1)*(ng_a+2) + (k-1)*ng_b*(ng_a+2)
#else
         ind_k = ind_bc
#endif
         sign_Im_k(ind_n) = 0
         ind_k_ord(ind_n) = ind_k
         do i = ng_a/2+2,ng_a
            ind_n = ind_n + 1
            ind_n_ord(ind_n) = ind_n
            k2_ord(ind_n) = (p_kvec(1,ind_n)**2 + p_kvec(2,ind_n)**2 + &
                             p_kvec(3,ind_n)**2)
#ifdef FFTW
            ind_k = 2*(ng_a+2-i) + (jconj-1)*(ng_a+2) + (kconj-1)*ng_b*(ng_a+2)
#else
            ind_k = 2*(ng_a+2-i) + (jconj-1)*ng_a + (kconj-1)*ng_b*ng_a
#endif
            sign_Im_k(ind_n) = -1
            ind_k_ord(ind_n) = ind_k
         enddo
         n0 = 1
      enddo
   enddo
   call sort_kmin(ng, ind_n_ord, ind_k_ord, sign_Im_k, k2_ord)
   deallocate( k2_ord )
!
!  Back FFT to compute the radial jl components of potential
!
   pv_jl(1:n_rmesh,1:jmax) = CZERO
!
   kfact_y = CZERO
   k2 = ZERO
   k2_prev = ZERO
   do ind_n = 2,ng
      ind_k = ind_k_ord(ind_n)
      if ( sign_Im_k(ind_n)/=0 ) then
#ifdef FFTW
! *** Complex ***
!           theta_k = fftw_in_c(ind_n)
!
! *** Real in-place ***
         theta_k = cmplx(p_fft_r(ind_k-1), &
                         -sign_Im_k(ind_n)*p_fft_r(ind_k),kind=CmplxKind)
!
#else
         theta_k = cmplx(p_fft_r(ind_k-1), &
                        sign_Im_k(ind_n)*p_fft_r(ind_k),kind=CmplxKind)
#endif
      else
#ifdef FFTW
!        theta_k = cmplx(p_fft_r(ind_k-1),p_fft_r(ind_k),kind=CmplxKind)
         theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)
!
#else
         theta_k = p_ft(ind_k)
#endif
      endif
!
      k2 = ( p_kvec(1,ind_n_ord(ind_n))**2 + &
             p_kvec(2,ind_n_ord(ind_n))**2 + &
             p_kvec(3,ind_n_ord(ind_n))**2)
      k2inv = ONE/k2
!     expikR = exp( ( p_kvec(1,ind_n_ord(ind_n))*posi(1) + &
      expikR = exp(-( p_kvec(1,ind_n_ord(ind_n))*posi(1) + &
                      p_kvec(2,ind_n_ord(ind_n))*posi(2) + &
                      p_kvec(3,ind_n_ord(ind_n))*posi(3) )*sqrtm1 )
!
      kfact = expikR*k2inv*theta_k
!
      call calYlm( p_kvec(1:3,ind_n_ord(ind_n)), lmax, Ylm )
!
      if ( k2/=k2_prev ) then
         if ( MyPE==1 .and. k2_prev/=ZERO ) then
            do i = 1,nr_print
               write(funit(i),'(2(1x,d21.13),$)') rmesh(rind_print(i)), sqrt(k2_prev)
               do j = 1,njl_print
                  jl = jl_print(j)
                  l = lofj(jl)
                  write(funit(i),'(5(1x,d21.13),$)') &
                   Bj_l(rind_print(i),lofj(jl)), kfact_y(j), &
                   pv_jl(rind_print(i),jl)/(rmesh(rind_print(i))**l)
               enddo
               write(funit(i),'(a)') " "
            enddo
         endif
         kfact_y = CZERO
         kri(1:n_rmesh) = sqrt(k2)*rmesh(1:n_rmesh)
         call SphericalBessel( lmax, n_rmesh, kri(1:n_rmesh), Bj_L(1:n_rmesh,0:lmax) )
      endif
!
      i2l = CONE
      do l = 0,lmax
         do m = 0,l
            kl = (l+1)*(l+1) -l + m
            jl = (l+1)*(l+2)/2 -l + m
            if ( MyPE==1 ) then
               do j = 1,njl_print
                  if (jl == jl_print(j)) then
                     kfact_y(j) = kfact_y(j) + kfact*Ylm(kl)
                  endif
               enddo
            endif
            pv_jl(1:n_rmesh,jl) = pv_jl(1:n_rmesh,jl)+ &
                                  TWO*PI4*PI4/real(ng,kind=RealKind)* &
                                  kfact*Ylm(kl)*Bj_L(1:n_rmesh,l)
         enddo
         i2l = i2l*(-sqrtm1)
      enddo
      k2_prev = k2
!
   enddo
!
   if ( MyPE==1 ) then
      do i = 1,nr_print
         write(funit(i),'(2(1x,d21.13),$)') rmesh(rind_print(i)),sqrt(k2_prev)
         do j = 1,njl_print
            jl = jl_print(j)
            l = lofj(jl)
            write(funit(i),'(5(1x,d21.13),$)') &
            Bj_l(rind_print(i),lofj(jl)), kfact_y(j), &
            pv_jl(rind_print(i),jl)/(rmesh(rind_print(i))**l)
         enddo
         write(funit(i),'(a)') ' '
         close( funit(i) )
      enddo
   endif
   call syncAllPEs()
!
#ifdef TIMING
   t0 = getTime()
   t1 = t0-t1
   write(6,*) "getFFTTransformRadial:: Time in back l-projection: ", t1
#endif
   deallocate( Ylm, Bj_L, kri )
!
   deallocate( ind_k_ord, ind_n_ord, sign_Im_k )
!
   deallocate( rind_print, jl_print )
   deallocate( funit, fname )
   deallocate( r_print, kfact_y )
!
   end subroutine calRadialProj
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine sort_kmin(ng, ind_n, ind_k, sign_k, k2)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ng
   integer (kind=IntKind), intent(inout) :: ind_n(ng), ind_k(ng), sign_k(ng)
   real (kind=RealKind), intent(inout) :: k2(ng)
!
   integer (kind=IntKind) :: ind_tmp, i, j
   real (kind=RealKind) :: k2_tmp
!
   do i = 2,ng-1
      do j = i+1,ng
         if ( k2(i) > k2(j) ) then
            k2_tmp = k2(j)
            k2(j) = k2(i)
            k2(i) = k2_tmp
            ind_tmp = ind_n(j)
            ind_n(j) = ind_n(i)
            ind_n(i) = ind_tmp
            ind_tmp = ind_k(j)
            ind_k(j) = ind_k(i)
            ind_k(i) = ind_tmp
            ind_tmp = sign_k(j)
            sign_k(j) = sign_k(i)
            sign_k(i) = ind_tmp
         endif
      enddo
   enddo
!
   end subroutine sort_kmin
!
#endif
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKIndex(n,nfft)                             result(g_ind)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n, nfft
   integer (kind=IntKind) :: g_ind
!
   if (n<nfft/2+1 .or. nfft==1) then
      g_ind = n-1
   else
      g_ind = n-1-nfft
   endif
!
   end function getKIndex
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calPseudoPotAtSite(posi,fft_in,v,dv)
!  ===================================================================
   use MathParamModule, only : THIRD
   use SphericalHarmonicsModule, only : calYlmConjg
   implicit none
!
   real (kind=RealKind), intent(in) :: posi(3)
!
   type (FFT_struct), intent(in) :: fft_in
!
   complex (kind=CmplxKind), intent(out) :: v, dv(3)
!
   integer (kind=IntKind) :: ng_a, ng_b, ng_c, ng
   integer (kind=IntKind) :: i, i1, j, k, ind_n
   integer (kind=IntKind) :: jconj, kconj, ind_bc, ind_k
!
   real (kind=RealKind) :: kfact
   real (kind=RealKind), pointer ::  kvec(:,:), p_fft_r(:)
!
   complex (kind=CmplxKind) :: theta_k, expiKR, cfac, Ylm(4)
!
#ifndef FFTW
   complex (kind=CmplxKind), pointer :: p_ft(:)
#endif
!
   ng_a = fft_in%nga
   ng_b = fft_in%ngb
   ng_c = fft_in%ngc
   ng   = fft_in%ng
!
   kvec => fft_in%kvec(1:3,1:ng)
#ifndef FFTW
   p_fft_r => fft_in%fft_r(1:ng)
   p_ft => fft_in%ft(1:ng_b*ng_c)
#else
   p_fft_r => fft_in%fft_r(1:(ng_a+2)*ng_b*ng_c)
#endif
!
   ind_bc = ng_c*ng_b+1
   ind_n = ng+1
!
   v = CZERO; dv = CZERO
   do k = ng_c, 1, -1
      if (k == 1) then
         kconj = 1
      else
         kconj = ng_c + 2 - k
      endif
!
      do j = ng_b, 1, -1
         if (j == 1) then
            jconj = 1
         else
            jconj = ng_b + 2 - j
         endif
!
         if ( j == 1 .and. k == 1 ) then
            i1 = 2 ! for Poisson Solver skip the k=0 point
         else
            i1 = 1
         endif
!
         ind_bc = ind_bc - 1
         do i = ng_a, i1, -1
            ind_n = ind_n - 1
#ifdef FFTW
! *** Complex ***
!           theta_k = fftw_in_c(ind_n)
!
! *** Real in-place ***
            if (i <= ng_a/2+1) then
               ind_k = 2*i + (j-1)*(ng_a+2) + (k-1)*ng_b*(ng_a+2)
!              theta_k = cmplx(p_fft_r(ind_k-1),p_fft_r(ind_k),kind=CmplxKind)
               theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)
            else
               ind_k = 2*(ng_a+2-i) + (jconj-1)*(ng_a+2) + (kconj-1)*ng_b*(ng_a+2)
!              theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)
               theta_k = cmplx(p_fft_r(ind_k-1), p_fft_r(ind_k),kind=CmplxKind)
            endif
!
#else
            if (i == ng_a/2+1) then
               theta_k = p_ft(ind_bc)
            else if (i < ng_a/2+1) then
               ind_k = 2*i + (j-1)*ng_a + (k-1)*ng_b*ng_a
               theta_k = cmplx(p_fft_r(ind_k-1),p_fft_r(ind_k),kind=CmplxKind)
            else
               ind_k = 2*(ng_a+2-i) + (jconj-1)*ng_a + (kconj-1)*ng_b*ng_a
               theta_k = cmplx(p_fft_r(ind_k-1),-p_fft_r(ind_k),kind=CmplxKind)
            endif
#endif
!
!           expikR  = exp( ( kvec(1,ind_n)*posi(1) + kvec(2,ind_n)*posi(2) + &
            expikR  = exp(-( kvec(1,ind_n)*posi(1) + kvec(2,ind_n)*posi(2) + &
                             kvec(3,ind_n)*posi(3) )*sqrtm1 )
            call calYlmConjg( kvec(1:3,ind_n), 1, Ylm )
            kfact = kvec(1,ind_n)**2 + kvec(2,ind_n)**2 + kvec(3,ind_n)**2
!
            v = v + expikR*theta_k/kfact
!
            cfac = expikR*theta_k/sqrt(kfact)
            dv(1) = dv(1) + cfac*Ylm(2)
            dv(2) = dv(2) + cfac*Ylm(3)
            dv(3) = dv(3) + cfac*Ylm(4)
!
         enddo
      enddo
   enddo
   cfac = TWO*PI4*sqrt(PI4)/real(ng,kind=RealKind)
   v = cfac*v
   cfac = TWO*PI4**2*THIRD*sqrtm1/real(ng,kind=RealKind)
   dv = cfac*dv
!
   end subroutine calPseudoPotAtSite
!  ===================================================================
end module FFTTransformModule
