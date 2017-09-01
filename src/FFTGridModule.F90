module FFTGridModule
   use KindParamModule, only : IntKind
   use KindParamModule, only : RealKind
   use PublicTypeDefinitionsModule, only : FFTGridStruct, FFTGAtomInfoStruct
   use Atom2ProcModule, only : getMaxLocalNumAtoms, getLocalNumAtoms
   use PolyhedraModule, only: isExternalPoint, isSurfacePoint
   use ErrorHandlerModule, only : ErrorHandler, StopHandler, WarningHandler
!
public ::            &
   initFFTGrid,      &
   endFFTGrid,       &
   setFFTGrid,       &
   printFFTGrid,     &
   getFFTGrid,       &
   getNumGridA,      &
   getNumGridB,      &
   getNumGridC,      &
   getNumGridPoints, &
   getGridStepA,     &
   getGridStepB,     &
   getGridStepC,     &
   getInfoPoint,     &
   getGPAtomsI,      &
   getGPAtomsN,      &
   getGPAtomsR,      &
   isVP_Point,       &
   getGPAtomsR_VP
!
private
   type (FFTGridStruct), pointer :: FFTGrid
   type (FFTGAtomInfoStruct), allocatable, target :: FFTGAInfo(:)
!
   logical :: Initialized = .false.
   logical :: DefRad = .false.
!
   character (len=50) :: stop_routine
!
   integer (kind=IntKind) :: print_level
   integer (kind=IntKind) :: FFTGAinfo_size
   integer (kind=IntKind) :: nAtomsGP_size
!
   integer (kind=IntKind), allocatable :: nAtomsGP(:)
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initFFTGrid(istop,iprint)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: iprint
!
   if (Initialized) then
      call WarningHandler('initFFTGrid','module has already been initialized')
   endif
!
   stop_routine = istop
   print_level = iprint
!
   allocate( FFTGrid )
!
   Initialized = .true.
!
   end subroutine initFFTGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endFFTGrid()
!  ===================================================================
   implicit none
!
   if (.not.Initialized) then
      call WarningHandler('endFFTGrid','module is not initialized')
      return
   endif
!
   call delete_FFTGAInfoList()
   deallocate(FFTGAInfo)
   deallocate(nAtomsGP)
   if ( DefRad ) then
      DefRad = .false.
   endif
!
   deallocate( FFTGrid )
!
   Initialized = .false.
!
   end subroutine endFFTGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine delete_FFTGAInfoList()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: i, ng
   type (FFTGAtomInfoStruct), pointer :: p_FFTGAInfo, p_FFTGAInfo_tmp
!
   if (.not.Initialized) then
      call WarningHandler('delete_FFTGAInfoList', &
                          'module FFTGrid is not initialized')
      return
   endif
!
   do i = 1,FFTGAInfo_size
      p_FFTGAInfo => FFTGAInfo(i)
      if ( nAtomsGP(i) > 1 ) then 
         p_FFTGAInfo => p_FFTGAInfo%nextAtom
         LoopNG: do ng = 2,nAtomsGP(i)
             p_FFTGAInfo_tmp => p_FFTGAInfo%nextAtom
             deallocate(p_FFTGAInfo)
             p_FFTGAInfo => p_FFTGAInfo_tmp
         enddo LoopNG
      else
         nullify(p_FFTGAInfo%nextAtom)
      endif
   enddo
   nullify( p_FFTGAInfo,  p_FFTGAInfo_tmp )
!
   end subroutine delete_FFTGAInfoList
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setFFTGrid(na,nb,nc,unit_cell, radius)
!  ===================================================================
   use VectorModule, only : getVecLength
   use MathParamModule, only : ZERO, TEN2m12
   use SystemModule, only : getBravaisLattice
   use AtomModule, only : getLocalAtomPosition
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: na, nb, nc
   real (kind=RealKind), intent(in) :: unit_cell(3,3)
   real (kind=RealKind), optional :: radius
!
   integer (kind=IntKind) :: gCounter, ic, jc, kc
   integer (kind=IntKind) :: ia, nLocal, lp, mp, np
!
   real (kind=RealKind) :: length, sqrtr2
   real (kind=RealKind) :: r(3), rg(3), bravais(3,3)
!
   type (FFTGAtomInfoStruct), pointer :: currentFFTGAInfo
   type (FFTGAtomInfoStruct), pointer :: previousFFTGAInfo
!
   if (.not.Initialized) then
      call ErrorHandler('setFFTGrid','module not initialized')
   else if (na < 1) then
      call ErrorHandler('setFFTGrid','na < 1',na)
   else if (nb < 1) then
      call ErrorHandler('setFFTGrid','nb < 1',nb)
   else if (nc < 1) then
      call ErrorHandler('setFFTGrid','nc < 1',nc)
   endif
!
   FFTGrid%nga = na
   FFTGrid%ngb = nb
   FFTGrid%ngc = nc
   FFTGrid%ng = na*nb*nc
!
   length = getVecLength(3,unit_cell(1:3,1))
   FFTGrid%vec_a(1:3) = unit_cell(1:3,1)/length
   FFTGrid%step_a = length/real(na,kind=RealKind)

!
   length = getVecLength(3,unit_cell(1:3,2))
   FFTGrid%vec_b(1:3) = unit_cell(1:3,2)/length
   FFTGrid%step_b = length/real(nb,kind=RealKind)
!
   length = getVecLength(3,unit_cell(1:3,3))
   FFTGrid%vec_c(1:3) = unit_cell(1:3,3)/length
   FFTGrid%step_c = length/real(nc,kind=RealKind)
!
   FFTGAInfo_size = FFTGrid%ng
   if ( present(radius) ) then
      DefRad = .true.
   else
      DefRad = .false.
   endif
!
   allocate( FFTGAInfo(1:FFTGrid%ng) )
   allocate( nAtomsGP(1:FFTGrid%ng) )
   do ic = 1, FFTGrid%ng
      nullify(FFTGAInfo(ic)%nextAtom)
      FFTGAInfo(ic)%r = ZERO
      FFTGAInfo(ic)%indexAtom = 0
      FFTGAInfo(ic)%VP_point = -1
   enddo
!
   nAtomsGP = 0
   nLocal = getLocalNumAtoms()
   bravais(1:3,1:3) = getBravaisLattice()
!  ===================================================================
!  Allocate the arrays holding atomic indices
!  Procedure to modify the search for atoms contributing to the grid point
!  Loop over grid points
!  ===================================================================
   if ( DefRad ) then
      gCounter = 0
      do kc = 0, FFTGrid%ngc-1
         do jc = 0, FFTGrid%ngb-1
            do ic = 0, FFTGrid%nga-1
               gCounter = gCounter + 1
               rg(1:3) = FFTGrid%step_a * FFTGrid%vec_a * ic + &
                         FFTGrid%step_b * FFTGrid%vec_b * jc + &
                         FFTGrid%step_c * FFTGrid%vec_c * kc
!           ===================================
!           Loop over local atoms
!           ===================================
            currentFFTGAInfo => FFTGAInfo(gCounter)
            previousFFTGAInfo => currentFFTGAInfo
            do lp = -5, 5
               do mp = -5, 5
                  do np = -5, 5
!
                     do ia = 1, nLocal
                        r(1:3) = getLocalAtomPosition(ia)+lp*bravais(1:3,1)+ &
                                 mp*bravais(1:3,2)+np*bravais(1:3,3)
                        sqrtr2 = sqrt( (rg(1)-r(1))**2 + (rg(2)-r(2))**2 + &
                                 (rg(3)-r(3))**2 )
                        if ( radius >= sqrtr2 ) then
                           if (nAtomsGP(gCounter) == 0) then
                              FFTGAInfo(gCounter)%indexAtom = ia
                              FFTGAInfo(gCounter)%r = rg(1:3) - r(1:3)
                              allocate(FFTGAInfo(gCounter)%nextAtom)
                              previousFFTGAInfo => currentFFTGAInfo
                              currentFFTGAInfo => FFTGAInfo(gCounter)%nextAtom
                              currentFFTGAInfo%VP_point = -1
                              nullify(currentFFTGAInfo%nextAtom)
                           else
                              currentFFTGAInfo%indexAtom = ia
                              currentFFTGAInfo%r = rg(1:3) - r(1:3)
                              allocate(currentFFTGAInfo%nextAtom)
                              previousFFTGAInfo => currentFFTGAInfo
                              currentFFTGAInfo => currentFFTGAInfo%nextAtom
                              currentFFTGAInfo%VP_point = -1
                              nullify(currentFFTGAInfo%nextAtom)
                           endif
                           nAtomsGP(gCounter) = nAtomsGP(gCounter)+1
                           if ( isExternalPoint(ia, rg(1)-r(1), &
                                rg(2)-r(2), rg(3)-r(3)) ) then
                              previousFFTGAInfo%VP_point = -1
                           else if ( isSurfacePoint(ia, rg(1)-r(1), &
                                rg(2)-r(2), rg(3)-r(3)) ) then
                              previousFFTGAInfo%VP_point = 0
                           else
                              previousFFTGAInfo%VP_point = 1
                           endif
                        endif
                     enddo
!
                  enddo
               enddo
            enddo
!
            enddo
         enddo
      enddo
   else
      gCounter = 0
      do kc = 0, FFTGrid%ngc-1
         do jc = 0, FFTGrid%ngb-1
            do ic = 0, FFTGrid%nga-1
               gCounter = gCounter + 1
               rg(1:3) = FFTGrid%step_a * FFTGrid%vec_a * ic + &
                         FFTGrid%step_b * FFTGrid%vec_b * jc + &
                         FFTGrid%step_c * FFTGrid%vec_c * kc
!           ===================================
!           Loop over local atoms
!           ===================================
            currentFFTGAInfo => FFTGAInfo(gCounter)
            do lp = -1, 1
               do mp = -1, 1
                  do np = -1, 1
!
                     do ia = 1, nLocal
                        r(1:3) = getLocalAtomPosition(ia)+lp*bravais(1:3,1)+ &
                                 mp*bravais(1:3,2)+np*bravais(1:3,3)
                        if (.not.(isExternalPoint(ia, rg(1)-r(1), &
                                       rg(2)-r(2), rg(3)-r(3)))) then
                           if (nAtomsGP(gCounter) == 0) then
                              FFTGAInfo(gCounter)%indexAtom = ia
                              FFTGAInfo(gCounter)%r = rg(1:3) - r(1:3)
                              allocate(FFTGAInfo(gCounter)%nextAtom)
                              previousFFTGAInfo => currentFFTGAInfo
                              currentFFTGAInfo => FFTGAInfo(gCounter)%nextAtom
                              currentFFTGAInfo%VP_point = -1
                              nullify(currentFFTGAInfo%nextAtom)
                           else
                              currentFFTGAInfo%indexAtom = ia
                              currentFFTGAInfo%r = rg(1:3) - r(1:3)
                              allocate(currentFFTGAInfo%nextAtom)
                              previousFFTGAInfo => currentFFTGAInfo
                              currentFFTGAInfo => currentFFTGAInfo%nextAtom
                              currentFFTGAInfo%VP_point = -1
                              nullify(currentFFTGAInfo%nextAtom)
                           endif
                           nAtomsGP(gCounter) = nAtomsGP(gCounter)+1
                           if ( isSurfacePoint(ia, rg(1)-r(1), &
                                rg(2)-r(2), rg(3)-r(3)) ) then
                              previousFFTGAInfo%VP_point = 0         ! on the surface of VP
                           else
                              previousFFTGAInfo%VP_point = 1         ! inside VP
                           endif
                        endif
                     enddo
!
                  enddo
               enddo
            enddo
!
            enddo
         enddo
      enddo
   endif
!
   end subroutine setFFTGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printFFTGrid()
!  ===================================================================
   implicit none
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('printFFTGrid','need to call initFFTGrid first')
!     ----------------------------------------------------------------
   endif
!
   write(6,'(/,80(''-''))')
   write(6,'(/,24x,a)')' ******************************'
   write(6,'(24x,a)')  ' *  Output from printFFTGrid  *'
   write(6,'(24x,a,/)')' ******************************'
   write(6,'(a)')'============================================================'
   write(6,'(  ''cell side a '',t40,''='',3d13.6)')FFTGrid%vec_a(1:3)
   write(6,'(  ''cell side b '',t40,''='',3d13.6)')FFTGrid%vec_b(1:3)
   write(6,'(  ''cell side c '',t40,''='',3d13.6)')FFTGrid%vec_c(1:3)
   write(6,'(  ''nga        '',t40,''='',i5)')     FFTGrid%nga
   write(6,'(  ''ngb        '',t40,''='',i5)')     FFTGrid%ngb 
   write(6,'(  ''ngc        '',t40,''='',i5)')     FFTGrid%ngc
   write(6,'(  ''number of grid points'',t40,''='',i10)') FFTGrid%ng
   write(6,'(  ''step_a     '',t40,''='',d13.6)')  FFTGrid%step_a
   write(6,'(  ''step_b     '',t40,''='',d13.6)')  FFTGrid%step_b
   write(6,'(  ''step_c     '',t40,''='',d13.6)')  FFTGrid%step_c
   write(6,'(a)')'============================================================'
!
   end subroutine printFFTGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getFFTGrid() result(gp)
!  ===================================================================
   implicit none
   type (FFTGridStruct), pointer :: gp
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getFFTGrid','need to call initFFTGrid first')
!     ----------------------------------------------------------------
   endif
!
   gp=>FFTGrid
   end function getFFTGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumGridA() result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getNumGridA','need to call initFFTGrid first')
!     ----------------------------------------------------------------
   endif
!
   n = FFTGrid%nga
   end function getNumGridA
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumGridB() result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getNumGridB','need to call initFFTGrid first')
!     ----------------------------------------------------------------
   endif
!
   n = FFTGrid%ngb
   end function getNumGridB
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumGridC() result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getNumGridC','need to call initFFTGrid first')
!     ----------------------------------------------------------------
   endif
!
   n = FFTGrid%ngc
   end function getNumGridC
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumGridPoints() result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getNumGridPoints','need to call initFFTGrid first')
!     ----------------------------------------------------------------
   endif
!
   n = FFTGrid%ng
!
   end function getNumGridPoints
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGridStepA() result(h)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: h
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGridStepA','need to call initFFTGrid first')
!     ----------------------------------------------------------------
   endif
!
   h = FFTGrid%step_a
!
   end function getGridStepA
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGridStepB() result(h)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: h
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGridStepB','need to call initFFTGrid first')
!     ----------------------------------------------------------------
   endif
!
   h = FFTGrid%step_b
!
   end function getGridStepB
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGridStepC() result(h)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: h
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGridStepC','need to call initFFTGrid first')
!     ----------------------------------------------------------------
   endif
!
   h = FFTGrid%step_c
!
   end function getGridStepC
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getInfoPoint( gCounter, ig )               result(fftInfo)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: gCounter
   integer (kind=IntKind), intent(out) :: ig
!
   type (FFTGAtomInfoStruct), pointer :: fftInfo
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGPAtomsI','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   ig = nAtomsGP(gCounter)
   fftInfo => FFTGAInfo(gCounter)
!
   end function getInfoPoint
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGPAtomsI(gCounter, ig) result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: gCounter, ig
   integer (kind=IntKind) :: n
!
   integer (kind=IntKind) :: ia
   type (FFTGAtomInfoStruct), pointer :: currentFFTGAInfo
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGPAtomsI','need to call initFFTGrid first')
!     ----------------------------------------------------------------
   endif
!
   currentFFTGAInfo => FFTGAInfo(gCounter)
   la: do ia = 1, ig
      n = currentFFTGAInfo%indexAtom
      if (associated(currentFFTGAInfo%nextAtom)) then
         currentFFTGAInfo => currentFFTGAInfo%nextAtom
      else
         exit la
      endif
   enddo la
!
   end function getGPAtomsI
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGPAtomsN(ig) result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: ig
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGPAtomsN','need to call initFFTGrid first')
!     ----------------------------------------------------------------
   endif
!
   n = nAtomsGP(ig)
!
   end function getGPAtomsN
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGPAtomsR(gCounter, ig) result(r)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: gCounter, ig
   real (kind=RealKind) :: r(3)
!
   integer (kind=IntKind) :: ia
   type (FFTGAtomInfoStruct), pointer :: currentFFTGAInfo
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGPAtomsR','need to call initFFTGrid first')
!     ----------------------------------------------------------------
   endif
!
   currentFFTGAInfo => FFTGAInfo(gCounter)
   la:do ia = 1, ig
      r = currentFFTGAInfo%r
      if (associated(currentFFTGAInfo%nextAtom)) then
         currentFFTGAInfo => currentFFTGAInfo%nextAtom
      else
         exit la
      endif
   enddo la
!
   end function getGPAtomsR
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isVP_Point(gCounter, ig) result(isVP)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: gCounter, ig
!
   integer (kind=IntKind), pointer :: isVP
!
   integer (kind=IntKind) :: ia
   type (FFTGAtomInfoStruct), pointer :: currentFFTGAInfo
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGPAtomsR','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   currentFFTGAInfo => FFTGAInfo(gCounter)
   la:do ia = 1, ig
      isVP => currentFFTGAInfo%VP_point
      if (associated(currentFFTGAInfo%nextAtom)) then
         currentFFTGAInfo => currentFFTGAInfo%nextAtom
      else
         exit la
      endif
   enddo la
!
   end function isVP_Point
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGPAtomsR_VP(gCounter, ig) result(r)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: gCounter, ig
   real (kind=RealKind) :: r(3)
!
   integer (kind=IntKind) :: ia, ing, ng
   type (FFTGAtomInfoStruct), pointer :: currentFFTGAInfo
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGPAtomsR','need to call initFFTGrid first')
!     ----------------------------------------------------------------
   endif
!
   currentFFTGAInfo => FFTGAInfo(gCounter)
   if ( DefRad ) then
      ng = nAtomsGP(gCounter)
   else
      ng = ig
   endif
   ing = 0
   la:do ia = 1, ng
      if ( currentFFTGAInfo%VP_point >=0 ) then
         ing = ing+1
         if (ing==ig) then
            r = currentFFTGAInfo%r
            exit la
         else
            currentFFTGAInfo => currentFFTGAInfo%nextAtom
            cycle la
         endif
      endif
      if (associated(currentFFTGAInfo%nextAtom)) then
         currentFFTGAInfo => currentFFTGAInfo%nextAtom
      else
         exit la
      endif
   enddo la
!
   end function getGPAtomsR_VP
!  ===================================================================
end module FFTGridModule
