module VisualGridModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler, StopHandler, WarningHandler
   use PublicTypeDefinitionsModule, only : VisualGridStruct, VGAtomInfoStruct
   use Atom2ProcModule, only : getLocalNumAtoms
   use AtomModule, only : getLocalAtomPosition
   use PolyhedraModule, only: isExternalPoint, isSurfacePoint
   use SystemModule, only : getBravaisLattice
!
public ::            &
   initVisualGrid,   &
   endVisualGrid,    &
   setVisualGrid,    &
   printVisualGrid,  &
   getVisualGrid,    &
   getVisualDim,     &
   getNumGridA,      &
   getNumGridB,      &
   getNumGridC,      &
   getNumGridPoints, &
   getGridStepA,     &
   getGridStepB,     &
   getGridStepC,     &
   getInfoPoint,     &
   getGPAtomsN,      &
   getGPAtomsR,      &
   getGPAtomsR_VP,   &
   getGPAtomsI,      &
   isVP_Point,       &
   printGPAtoms
!
private
   type (VisualGridStruct), target :: VisualGrid
   type (VGAtomInfoStruct), allocatable, target :: VGAInfo(:)
!
   logical :: Initialized = .false.
   logical :: DefRad = .false.
!
   character (len=50) :: stop_routine
!
   integer (kind=IntKind) :: print_level
   integer (kind=IntKind) :: VGAinfo_size
   integer (kind=IntKind) :: nAtomsGP_size
!
   integer (kind=IntKind), allocatable :: nAtomsGP(:)
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initVisualGrid(istop,iprint,radius)
!  ===================================================================
   use MathParamModule, only : ZERO
   use ScfDataModule, only : TableID
   use InputModule, only : getKeyValue, getNumKeyValues
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: iprint
   real (kind=RealKind), intent(in), optional :: radius
!
   character (len=80) :: svalue
!
   integer (kind=IntKind) :: vg_dim, vg_na, vg_nb, vg_nc
   real (kind=RealKind) :: vg_v0(3), vg_vec(3,3), a0_vgrid
!
   if (Initialized) then
      call WarningHandler('initVisualGrid','module has already been initialized')
   endif
!
   stop_routine = istop
   print_level = iprint
!
   Initialized = .true.
!
   if ( getKeyValue(TableID,'Visual Grid Type (0<D<4)', svalue) == 0 ) then
      read(svalue,*)vg_dim
   else
      call ErrorHandler('initVisualGrid','Visual Grid Type (0<D<4)','Not exist')
   endif
!
   if ( getKeyValue(TableID,'Origin Grid Vector', svalue) == 0 ) then
      read(svalue,*)vg_v0
   else
      call ErrorHandler('initVisualGrid','Origin Grid Vector','Not exist')
   endif
!
   if ( getNumKeyValues(TableID,'Grid Vector') < vg_dim ) then
      call ErrorHandler('initVisualGrid','number of grid vectors < dimension')
   else if (vg_dim > 3) then
      call ErrorHandler('initVisualGrid','dimension > 3',vg_dim)
   endif
!
   vg_vec = ZERO
   if ( getKeyValue(TableID,'Grid Vector', 3, vg_vec, vg_dim) /= 0 ) then
      call ErrorHandler('initVisualGrid','Grid Vector','Data is flawed')
   endif
!
   if ( getKeyValue(TableID,'Grid Points', svalue) == 0 ) then
      read(svalue,*)vg_na, vg_nb, vg_nc
   else
      call ErrorHandler('initVisualGrid','Visual Grid Type (0<D<4)','Not exist')
   endif
!
   if ( getKeyValue(TableID,'Grid Scale', a0_vgrid) == 0 ) then
      if ( a0_vgrid > ZERO ) then
         vg_v0  = a0_vgrid*vg_v0
         vg_vec = a0_vgrid*vg_vec
      endif
   endif
!
   if ( present(radius) ) then
      call setVisualGrid(vg_dim,vg_na,vg_nb,vg_nc,vg_v0,vg_vec(1:3,1),vg_vec(1:3,2),vg_vec(1:3,3),radius)
   else
      call setVisualGrid(vg_dim,vg_na,vg_nb,vg_nc,vg_v0,vg_vec(1:3,1),vg_vec(1:3,2),vg_vec(1:3,3))
   endif
!
   end subroutine initVisualGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endVisualGrid()
!  ===================================================================
   implicit none
!
   if (.not.Initialized) then
      call WarningHandler('endVisualGrid','module is not initialized')
      return
   endif
!
   call delete_VGAInfoList()
   deallocate(VGAInfo)
   deallocate(nAtomsGP)
   DefRad = .false.
!
   Initialized = .false.
!
   end subroutine endVisualGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine delete_VGAInfoList()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: i, ng
   type (VGAtomInfoStruct), pointer :: p_VGAInfo, p_VGAInfo_tmp
!
   if (.not.Initialized) then
      call WarningHandler('endVisualGrid','module is not initialized')
      return
   endif
!
   do i = 1,VGAInfo_size
      p_VGAInfo => VGAInfo(i)
      if ( nAtomsGP(i) > 1 ) then 
         p_VGAInfo => p_VGAInfo%nextAtom
         LoopNG: do ng = 2,nAtomsGP(i)
             p_VGAInfo_tmp => p_VGAInfo%nextAtom
             deallocate(p_VGAInfo)
             p_VGAInfo => p_VGAInfo_tmp
         enddo LoopNG
      else
         nullify(p_VGAInfo%nextAtom)
      endif
   enddo
   nullify( p_VGAInfo,  p_VGAInfo_tmp )
!
   end subroutine delete_VGAInfoList
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setVisualGrid(dim, na, nb, nc, vb0, vb1, vb2, vb3, radius)
!  ===================================================================
   use VectorModule, only : getVecLength
   use MathParamModule, only : ZERO, TEN2m12
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: na, nb, nc, dim
!
   real (kind=RealKind), intent(in) :: vb0(3), vb1(3), vb2(3), vb3(3)
   real (kind=RealKind), optional :: radius
!
   integer (kind=IntKind) :: ia, nLocal, lp, mp, np
   integer (kind=IntKind) :: gCounter, ic, jc, kc
!
   real (kind=RealKind) :: length, sqrtr2
   real (kind=RealKind) :: r(3), rg(3), bravais(3,3)
!
   type (VGAtomInfoStruct), pointer :: currentVGAInfo
   type (VGAtomInfoStruct), pointer :: previousVGAInfo
!

   if (.not.Initialized) then
      call ErrorHandler('setVisualGrid','module not initialized')
   else if (na < 1) then
      call ErrorHandler('setVisualGrid','na < 1',na)
   else if (nb < 1) then
      call ErrorHandler('setVisualGrid','nb < 1',nb)
   else if (nc < 1) then
      call ErrorHandler('setVisualGrid','nc < 1',nc)
   else if (dim < 1) then
      call ErrorHandler('setVisualGrid','dim < 1',dim)
   else if (dim > 3) then
      call ErrorHandler('setVisualGrid','dim > 3',dim)
   endif
!
   VisualGrid%nga = 1
   VisualGrid%ngb = 1
   VisualGrid%ngc = 1
   VisualGrid%vec_o = vb0
   VisualGrid%vec_a = ZERO
   VisualGrid%vec_b = ZERO
   VisualGrid%vec_c = ZERO
   VisualGrid%step_a = ZERO
   VisualGrid%step_b = ZERO
   VisualGrid%step_c = ZERO
   VisualGrid%dim = dim
!
   VisualGrid%ng = na+1
   if ( VisualGrid%ng <= 1 ) then 
      call ErrorHandler('setVisualGrid','Number of grid points',   &
                        VisualGrid%ng )
   endif
!
   length = getVecLength(3,vb1(1:3) - vb0(1:3))
   if ( length <= TEN2m12 ) then 
      call ErrorHandler('setVisualGrid','length =< 0', length )
   endif
!
   VisualGrid%nga = na+1
   VisualGrid%step_a=length/real(na,kind=RealKind)
   VisualGrid%vec_a(1:3) = (vb1(1:3) - vb0(1:3))/length
!
   if ( VisualGrid%dim > 1 ) then
!
      VisualGrid%ng = (na+1)*(nb+1)
      if ( VisualGrid%ng <= 1 ) then 
         call ErrorHandler('setVisualGrid','Number of grid points',   &
                        VisualGrid%ng )
      endif
!
      length = getVecLength(3,vb2(1:3) - vb0(1:3))
      if ( length <= TEN2m12 ) then 
         call ErrorHandler('setVisualGrid','length =< 0', length )
      endif
!
      VisualGrid%ngb = nb+1
      VisualGrid%step_b=length/real(nb,kind=RealKind)
      VisualGrid%vec_b(1:3) = (vb2(1:3) - vb0(1:3))/length
!
   endif
   if ( VisualGrid%dim > 2 ) then
!
      VisualGrid%ng = (na+1)*(nb+1)*(nc+1)
      if ( VisualGrid%ng <= 1 ) then 
         call ErrorHandler('setVisualGrid','Number of grid points',   &
                        VisualGrid%ng )
      endif
!
      length = getVecLength(3,vb3(1:3) - vb0(1:3))
      if ( length <= TEN2m12 ) then 
         call ErrorHandler('setVisualGrid','length =< 0', length )
      endif
!
      VisualGrid%ngc = nc+1
      VisualGrid%step_c=length/real(nc,kind=RealKind)
      VisualGrid%vec_c(1:3) = (vb3(1:3) - vb0(1:3))/length
!
   endif
!
   if ( present(radius) ) then
      DefRad = .true.
   else
      DefRad = .false.
   endif
   VGAInfo_size = VisualGrid%ng
   allocate( VGAInfo(1:VisualGrid%ng) )
   allocate( nAtomsGP(1:VisualGrid%ng) )
   do ic = 1, VisualGrid%ng
      nullify(VGAInfo(ic)%nextAtom)
      VGAInfo(ic)%r = ZERO
      VGAInfo(ic)%indexAtom = 0
      VGAInfo(ic)%VP_point = -1
   enddo
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
      do kc = 0, VisualGrid%ngc-1
         do jc = 0, VisualGrid%ngb-1
            do ic = 0, VisualGrid%nga-1
               gCounter = gCounter + 1
               rg(1:3) = VisualGrid%step_a * VisualGrid%vec_a * ic + &
                         VisualGrid%step_b * VisualGrid%vec_b * jc + &
                         VisualGrid%step_c * VisualGrid%vec_c * kc + &
                         VisualGrid%vec_o
!           ===================================
!           Loop over local atoms
!           ===================================
            currentVGAInfo => VGAInfo(gCounter)
            previousVGAInfo => currentVGAInfo
            do lp = -1, 1
               do mp = -1, 1
                  do np = -1, 1
!
                     do ia = 1, nLocal
                        r(1:3) = getLocalAtomPosition(ia)+lp*bravais(1:3,1)+ &
                                 mp*bravais(1:3,2)+np*bravais(1:3,3)
                        sqrtr2 = sqrt( (rg(1)-r(1))**2 + (rg(2)-r(2))**2 + &
                                 (rg(3)-r(3))**2 )
                        if ( radius >= sqrtr2 ) then
                           if (nAtomsGP(gCounter) == 0) then
                              VGAInfo(gCounter)%indexAtom = ia
                              VGAInfo(gCounter)%r = rg(1:3) - r(1:3)
                              allocate(VGAInfo(gCounter)%nextAtom)
                              previousVGAInfo => currentVGAInfo
                              currentVGAInfo => VGAInfo(gCounter)%nextAtom
                              currentVGAInfo%VP_point = -1
                              nullify(currentVGAInfo%nextAtom)
                           else
                              currentVGAInfo%indexAtom = ia
                              currentVGAInfo%r = rg(1:3) - r(1:3)
                              allocate(currentVGAInfo%nextAtom)
                              previousVGAInfo => currentVGAInfo
                              currentVGAInfo => currentVGAInfo%nextAtom
                              currentVGAInfo%VP_point = -1
                              nullify(currentVGAInfo%nextAtom)
                           endif
                           nAtomsGP(gCounter) = nAtomsGP(gCounter)+1
                           if ( isExternalPoint(ia, rg(1)-r(1), &
                                rg(2)-r(2), rg(3)-r(3)) ) then
                              previousVGAInfo%VP_point = -1
                           else if ( isSurfacePoint(ia, rg(1)-r(1), &
                                rg(2)-r(2), rg(3)-r(3)) ) then
                              previousVGAInfo%VP_point = 0         ! on the surface of VP
                           else
                              previousVGAInfo%VP_point = 1         ! inside VP
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
      do kc = 0, VisualGrid%ngc-1
         do jc = 0, VisualGrid%ngb-1
            do ic = 0, VisualGrid%nga-1
               gCounter = gCounter + 1
               rg(1:3) = VisualGrid%step_a * VisualGrid%vec_a * ic + &
                         VisualGrid%step_b * VisualGrid%vec_b * jc + &
                         VisualGrid%step_c * VisualGrid%vec_c * kc + &
                         VisualGrid%vec_o
!           ===================================
!           Loop over local atoms
!           ===================================
            currentVGAInfo => VGAInfo(gCounter)
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
                              VGAInfo(gCounter)%indexAtom = ia
                              VGAInfo(gCounter)%r = rg(1:3) - r(1:3)
                              allocate(VGAInfo(gCounter)%nextAtom)
                              previousVGAInfo => currentVGAInfo
                              currentVGAInfo => VGAInfo(gCounter)%nextAtom
                              currentVGAInfo%VP_point = -1
                              nullify(currentVGAInfo%nextAtom)
                           else
                              currentVGAInfo%indexAtom = ia
                              currentVGAInfo%r = rg(1:3) - r(1:3)
                              allocate(currentVGAInfo%nextAtom)
                              previousVGAInfo => currentVGAInfo
                              currentVGAInfo => currentVGAInfo%nextAtom
                              currentVGAInfo%VP_point = -1
                              nullify(currentVGAInfo%nextAtom)
                           endif
                           nAtomsGP(gCounter) = nAtomsGP(gCounter)+1
                           if ( isSurfacePoint(ia, rg(1)-r(1), &
                                rg(2)-r(2), rg(3)-r(3)) ) then
                              previousVGAInfo%VP_point = 0         ! on the surface of VP
                           else
                              previousVGAInfo%VP_point = 1         ! inside VP
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
   end subroutine setVisualGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printVisualGrid()
!  ===================================================================
   implicit none
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('printVisualGrid', &
                         'need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   write(6,'(/,80(''-''))')
   write(6,'(/,24x,a)')' *********************************'
   write(6,'(24x,a)')  ' *  Output from printVisualGrid  *'
   write(6,'(24x,a,/)')' *********************************'
   write(6,'(a)')'============================================================'
   write(6,'(  ''dimension   '',t40,''='',i3)') VisualGrid%dim
   write(6,'(  ''cell side a '',t40,''='',3d13.6)') VisualGrid%vec_a
   if (VisualGrid%dim > 1) then
      write(6,'(  ''cell side b '',t40,''='',3d13.6)')VisualGrid%vec_b
   endif
   if (VisualGrid%dim > 2) then
      write(6,'(  ''cell side c '',t40,''='',3d13.6)')VisualGrid%vec_c
   endif
   write(6,'(  ''nga        '',t40,''='',i5)') VisualGrid%nga
   if (VisualGrid%dim > 1) then
      write(6,'(  ''ngb        '',t40,''='',i5)') VisualGrid%ngb
   endif
   if (VisualGrid%dim > 2) then
      write(6,'(  ''ngc        '',t40,''='',i5)') VisualGrid%ngc
   endif
   write(6,'(  ''number of grid points'',t40,''='',i10)') VisualGrid%ng
   write(6,'(  ''step_a     '',t40,''='',d13.6)') VisualGrid%step_a
   if (VisualGrid%dim > 1) then
      write(6,'(  ''step_b     '',t40,''='',d13.6)') VisualGrid%step_b
   endif
   if (VisualGrid%dim > 2) then
      write(6,'(  ''step_c     '',t40,''='',d13.6)') VisualGrid%step_c
   endif
   write(6,'(a)')'============================================================'

!
   end subroutine printVisualGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printGPAtoms()
!  ===================================================================
   use MPPModule, only : MyPE
   implicit none
   integer (kind=IntKind) :: i, j, nLocal
   type (VGAtomInfoStruct), pointer :: currentVGAInfo
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('printGPAtoms','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   write(6,'(/,80(''-''))')
   write(6,'(/,24x,a)')' *********************************'
   write(6,'(24x,a)')  ' *  Output from printGPAtoms     *'
   write(6,'(24x,a,/)')' *********************************'
   write(6,'(a)')'============================================================'
   nLocal = getLocalNumAtoms()
   write(6, '(''**** MyPE = '', i5)') MyPE
   do j = 1, VisualGrid%ng
    currentVGAInfo => VGAInfo(j)
    write(6, '(''grid point, number of atoms = '', 3i5)') j, nAtomsGP(j)
    LoopNAGP: do i = 1, nAtomsGP(j)
     write(6, '(t5,''{i, atom index}'',t20, '' = {'',i5, i5, ''} '')') &
              i, currentVGAInfo%indexAtom
     if (associated(currentVGAInfo%nextAtom)) then
       currentVGAInfo => currentVGAInfo%nextAtom
     else
       exit LoopNAGP
     endif
    enddo LoopNAGP
   enddo
!
   end subroutine printGPAtoms
!  ===================================================================
!

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVisualDim() result(gridDimension)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: gridDimension
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVisualGrid','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   gridDimension = VisualGrid%dim
!
   end function getVisualDim
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVisualGrid() result(gp)
!  ===================================================================
   implicit none
!
   type (VisualGridStruct), pointer :: gp
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVisualGrid','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   gp => VisualGrid
!
   end function getVisualGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getInfoPoint( gCounter, ig )               result(visInfo)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: gCounter
   integer (kind=IntKind), intent(out) :: ig
!
   type (VGAtomInfoStruct), pointer :: visInfo
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGPAtomsI','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   ig = nAtomsGP(gCounter)
   visInfo => VGAInfo(gCounter)
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
   integer (kind=IntKind) :: ia
   type (VGAtomInfoStruct), pointer :: currentVGAInfo
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGPAtomsI','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   currentVGAInfo => VGAInfo(gCounter)
   la: do ia = 1, ig
      n = currentVGAInfo%indexAtom
      if (associated(currentVGAInfo%nextAtom)) then
         currentVGAInfo => currentVGAInfo%nextAtom
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
      call ErrorHandler('getGPAtomsN','need to call initVisualGrid first')
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
   integer (kind=IntKind) :: ia
   type (VGAtomInfoStruct), pointer :: currentVGAInfo
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGPAtomsR','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   currentVGAInfo => VGAInfo(gCounter)
   la:do ia = 1, ig
      r = currentVGAInfo%r
      if (associated(currentVGAInfo%nextAtom)) then
         currentVGAInfo => currentVGAInfo%nextAtom
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
   integer, pointer :: isVP
!
   integer (kind=IntKind) :: ia
   type (VGAtomInfoStruct), pointer :: currentVGAInfo
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGPAtomsR','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   currentVGAInfo => VGAInfo(gCounter)
   la:do ia = 1, ig
      isVP => currentVGAInfo%VP_point
      if (associated(currentVGAInfo%nextAtom)) then
         currentVGAInfo => currentVGAInfo%nextAtom
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
   integer (kind=IntKind) :: ia, ing, ng
   type (VGAtomInfoStruct), pointer :: currentVGAInfo
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGPAtomsR','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   currentVGAInfo => VGAInfo(gCounter)
   if ( DefRad ) then
      ng = nAtomsGP(gCounter)
   else
      ng = ig
   endif
   ing = 0
   la:do ia = 1, ng
      if ( currentVGAInfo%VP_point>=0 ) then
         ing = ing+1
         if (ing==ig) then
            r = currentVGAInfo%r
            exit la
         else
            currentVGAInfo => currentVGAInfo%nextAtom
            cycle la
         endif
      endif
      if (associated(currentVGAInfo%nextAtom)) then
         currentVGAInfo => currentVGAInfo%nextAtom
      else
         exit la
      endif
   enddo la
!
   end function getGPAtomsR_VP
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
      call ErrorHandler('getNumGridA','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   n = VisualGrid%nga
!
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
      call ErrorHandler('getNumGridB','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   n = VisualGrid%ngb
!
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
      call ErrorHandler('getNumGridC','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   n = VisualGrid%ngc
!
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
      call ErrorHandler('getNumGridPoints','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   n = VisualGrid%ng
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
      call ErrorHandler('getGridStepA','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   h = VisualGrid%step_a
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
      call ErrorHandler('getGridStepB','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   h = VisualGrid%step_b
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
      call ErrorHandler('getGridStepC','need to call initVisualGrid first')
!     ----------------------------------------------------------------
   endif
!
   h = VisualGrid%step_c
!
   end function getGridStepC
!  ===================================================================

end module VisualGridModule
