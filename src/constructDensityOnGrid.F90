!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine constructDensityOnGrid( density, den_type, lmax, spin,  &
                                 grid_type, InitModeFFT, InitModeVis )
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler
   use MathParamModule, only : ZERO
   use GroupCommModule, only : getGroupID, GlobalSumInGroup
!
   use Atom2ProcModule, only : getLocalNumAtoms, getGlobalIndex
!
   use DataServiceCenterModule, only : createDataStorage,     &
                                       getDataStorage,        &
                                       isDataStorageExisting, &
                                       RealType, IntegerType, &
                                       RealMark, IntegerMark
!
  use PublicTypeDefinitionsModule, only : VisualGridStruct, VGAtomInfoStruct, &
                                          FFTGridStruct, FFTGAtomInfoStruct
!
   use VisualGridModule, only : getVisualGrid
   use VisualGridModule, only : getVisAtN => getGPAtomsN
   use VisualGridModule, only : getVisAtR => getGPAtomsR, getVisAtI => getGPAtomsI
   use VisualGridModule, only : getVisAtR_VP => getGPAtomsR_VP
   use VisualGridModule, only : isVP_Point_Vis => isVP_Point
   use VisualGridModule, only : getInfoPointVis => getInfoPoint
!
   use FFTGridModule, only : getFFTGrid
   use FFTGridModule, only : getFFTAtN => getGPAtomsN
   use FFTGridModule, only : getFFTAtR => getGPAtomsR, getFFTAtI => getGPAtomsI
   use FFTGridModule, only : getFFTAtR_VP => getGPAtomsR_VP
   use FFTGridModule, only : isVP_Point_FFT => isVP_Point
   use FFTGridModule, only : getInfoPointFFT => getInfoPoint
!
   use AtomModule, only : getLocalnumSpecies, getLocalSpeciesContent
!
   use ChargeDensityModule, only : getChargeDensity, getMomentDensity
!
   use PotentialGenerationModule, only : getNewPotential => getPotential
!
!v2.0   use SingleSiteSolverModule, only : getWaveFunction
!
   implicit none
!
!  Variables
!
   integer (kind=IntKind), intent(in) :: lmax, spin
   integer (kind=IntKind), intent(inout) :: InitModeFFT, InitModeVis
!
   character(len=*), intent(in) :: density
   character(len=*), intent(in) :: den_type
   character(len=*), intent(in) :: grid_type
!
   logical :: isUniformGrid
   logical :: isCharge
   logical :: isPotential
   logical :: isMoment
   logical :: isWave
!
   character(len=3)  :: Lchar
   character(len=7)  :: gridKey
   character(len=20) :: keyDen
!
   character(len=50) :: keyAtomIndex, keyAtomNumber, keyAtomPosi
   character(len=50) :: keyDenData
   integer (kind=IntKind) :: lenKey = 50
   integer (kind=IntKind) :: lenKeyType
!
   integer (kind=IntKind) :: gCounter, ig, jmax
   integer (kind=IntKind) :: nga, ngb, ngc, ng
   integer (kind=IntKind) :: n, nc, ic, jc, kc, id, ia
   integer (kind=IntKind) :: atomsOnPoint, VP_point, n_mult
   integer (kind=IntKind) :: offset, is
   integer (kind=IntKind) :: GroupID
   integer (kind=IntKind), pointer :: indexAtom(:,:,:)
   integer (kind=IntKind), pointer :: nAtoms(:,:,:)
!
   real (kind=RealKind) :: r(3), sqrtr2, den_r
   real (kind=RealKind), pointer :: denOnGrid(:,:,:)
   type (FFTGridStruct), pointer :: fftgp
   type (FFTGAtomInfoStruct), pointer :: fft_point
   type (VisualGridStruct), pointer:: visualgp
   type (VGAtomInfoStruct), pointer :: vis_point
!
   GroupID = getGroupID('Unit Cell')
!
   isCharge = .false.
   isPotential=.false.
   isMoment = .false.
   isWave = .false.
   offset = 100
   jmax = (lmax+1)*(lmax+2)/2
   if ( spin==1 .or. spin==2 ) then
      is = spin
   else
      is = -1
   endif
! 
   gridKey(1:7) = " "
   if ( grid_type(1:10) == "VisualGrid" ) then
      gridKey(1:6) = "Visual"
      nc = 6
      visualgp => getVisualGrid()
      nga = visualgp%nga
      ngb = visualgp%ngb
      ngc = visualgp%ngc
      ng  = visualgp%ng
      isUniformGrid = .true.
   else if ( grid_type(1:11) == "UniformGrid" ) then
      gridKey(1:7) = "Uniform"
      nc = 7
      fftgp => getFFTGrid()
      nga = fftgp%nga
      ngb = fftgp%ngb
      ngc = fftgp%ngc
      ng  = fftgp%ng
      isUniformGrid = .false.
   else
      call ErrorHandler("constructDensityOnGrid","Unknown grid type:",grid_type)
   endif
!
   keyDen(1:20) = " "
   keyAtomIndex(1:lenKey)     = " "
   keyAtomNumber(1:lenKey)    = " "
   keyAtomPosi(1:lenKey)      = " "
   keyDenData(1:lenKey)       = " "
!
   keyAtomIndex  = "Atom Index on "//gridKey(1:nc)//" Mesh"
   if ( .not.isDataStorageExisting( trim(keyAtomIndex) ) ) then
      call createDataStorage( trim(keyAtomIndex), &
                             ng, IntegerType)
   endif
!
   keyAtomNumber = "Atom Number on "//gridKey(1:nc)//" Mesh"
   if ( .not.isDataStorageExisting( trim(keyAtomNumber) ) ) then
      call createDataStorage( trim(keyAtomNumber), &
                              ng, IntegerType)
   endif
!
   indexAtom => getDataStorage( trim(keyAtomIndex),                   &
                                nga, ngb, ngc, IntegerMark)
   nAtoms => getDataStorage( trim(keyAtomNumber),                     &
                             nga, ngb, ngc, IntegerMark)
!
!  Generate the atom degeneracy, index and relative position on the uniform grid points
!      InitMode == 0  => It is first time or the atom positions have changed
!
   if ( isUniformGrid .and. InitModeVis == 0 ) then
      n = getLocalNumAtoms()
      indexAtom = 0
      nAtoms = 0
      gCounter = 0
      id = -1
      do kc = 1, ngc
         do jc = 1, ngb
            do ic = 1, nga
               gCounter = gCounter + 1
               vis_point => getInfoPointVis( gCounter, atomsOnPoint )
               nAtoms(ic,jc,kc) = atomsOnPoint
!
               do ig = 1, atomsOnPoint
!
!                 =================================================
!                 loop over atoms on this grid point
!                 =================================================
                  r(1:3) = vis_point%r(1:3)
                  sqrtr2 = sqrt( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) )
                  id = vis_point%indexAtom
                  VP_point = vis_point%VP_Point
                  indexAtom(ic,jc,kc) = getGlobalIndex(id)
!                  rAtoms(ic,jc,kc)    = sqrtr2
                  if ( VP_point<0 ) then ! point inside or on the surface of Voronoi Polyhedron
                     nAtoms(ic,jc,kc) = nAtoms(ic,jc,kc) -1
                  endif
                  vis_point => vis_point%nextAtom
!
               enddo
               nullify( vis_point )
            enddo
         enddo
      enddo
      call GlobalSumInGroup(GroupID, nAtoms, nga, ngb, ngc)
      InitModeVis = 1
   endif
!
   if ( (.not.isUniformGrid) .and. InitModeFFT == 0 ) then
      n = getLocalNumAtoms()
      indexAtom = 0
      nAtoms = 0
      gCounter = 0
      id = -1
      do kc = 1, ngc
         do jc = 1, ngb
            do ic = 1, nga
               gCounter = gCounter + 1
               fft_point => getInfoPointFFT( gCounter, atomsOnPoint )
               nAtoms(ic,jc,kc) = atomsOnPoint
!
               do ig = 1, atomsOnPoint
!                 =================================================
!                 loop over atoms on this grid point
!                 =================================================
                  r(1:3) = fft_point%r(1:3)
                  sqrtr2 = sqrt( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) )
                  id = fft_point%indexAtom
                  VP_point = fft_point%VP_Point
                  indexAtom(ic,jc,kc) = getGlobalIndex(id)
!                  rAtoms(ic,jc,kc)    = sqrtr2
                  if ( VP_point<0 ) then
                     nAtoms(ic,jc,kc) = nAtoms(ic,jc,kc) -1
                  endif
                  fft_point => fft_point%nextAtom
!
               enddo
               nullify( fft_point )
            enddo
         enddo
      enddo
      call GlobalSumInGroup(GroupID, nAtoms, nga, ngb, ngc)
      InitModeFFT = 1
   endif
!
   lenKeyType = len(den_type)
   keyDen(1:lenKeyType) = den_type(1:lenKeyType)
   if ( density(1:6)=="Charge" ) then
      isCharge = .true.
!
      if ( jmax > 0 ) then
         write(Lchar,'(i3)') offset + lmax
         Lchar(1:1) = "L"
         keyDenData = &
         keyDen(1:lenKeyType)//" Electron Density on "//gridKey(1:nc)//" Mesh "//Lchar
         if (.not.isDataStorageExisting( trim(keyDenData) ) ) then
            call createDataStorage( trim(keyDenData), ng, RealType )
         endif
      else
         keyDenData = keyDen(1:lenKeyType)//" Electron Density on "//gridKey(1:nc)//" Mesh"
         if (.not.isDataStorageExisting( trim(keyDenData) ) ) then
            call createDataStorage( trim(keyDenData), ng, RealType )
         endif
      endif
!
      denOnGrid => getDataStorage( trim(keyDenData),                  &
                                   nga, ngb, ngc, RealMark )
   else if ( density(1:6)=="Moment" ) then
      isMoment = .true.
!
      if ( jmax > 0 ) then
         write(Lchar,'(i3)') offset + lmax
         Lchar(1:1) = "L"
         keyDenData = &
         keyDen(1:lenKeyType)//" Moment Density on "//gridKey(1:nc)//" Mesh "//Lchar
         if (.not.isDataStorageExisting( trim(keyDenData) ) ) then
            call createDataStorage( trim(keyDenData), ng, RealType )
         endif
      else
         keyDenData = keyDen(1:lenKeyType)//" Moment Density on "//gridKey(1:nc)//" Mesh"
         if (.not.isDataStorageExisting( trim(keyDenData) ) ) then
            call createDataStorage( trim(keyDenData), ng, RealType )
         endif
      endif
!
      denOnGrid => getDataStorage( trim(keyDenData),                  &
                                   nga, ngb, ngc, RealMark )
   else if ( density(1:9)=="Potential" ) then
      isPotential=.true.
!
      if ( jmax > 0 ) then
         write(Lchar,'(i3)') offset + lmax
         Lchar(1:1) = "L"
         keyDenData = & 
         keyDen(1:lenKeyType)//" Potential on "//gridKey(1:nc)//" Mesh "//Lchar
         if (.not.isDataStorageExisting( trim(keyDenData) ) ) then
            call createDataStorage( trim(keyDenData), ng, RealType )
         endif
      else
         keyDenData = keyDen(1:lenKeyType)//" Potential on "//gridKey(1:nc)//" Mesh"
         if (.not.isDataStorageExisting( trim(keyDenData) ) ) then
            call createDataStorage( trim(keyDenData), ng, RealType )
         endif
      endif
!
      denOnGrid => getDataStorage( trim(keyDenData),                  &
                                   nga, ngb, ngc, RealMark )
   else if ( density(1:9)=="WaveFunct" ) then
      isWave=.true.
!
      if ( jmax > 0 ) then
         write(Lchar,'(i3)') offset + lmax
         Lchar(1:1) = "L"
         keyDenData = & 
         keyDen(1:lenKeyType)//" WaveFunct on "//gridKey(1:nc)//" Mesh "//Lchar
         if (.not.isDataStorageExisting( trim(keyDenData) ) ) then
            call createDataStorage( trim(keyDenData), ng, RealType )
         endif
      else
         keyDenData = keyDen(1:lenKeyType)//" WaveFunct on "//gridKey(1:nc)//" Mesh"
         if (.not.isDataStorageExisting( trim(keyDenData) ) ) then
            call createDataStorage( trim(keyDenData), ng, RealType )
         endif
      endif
!
      denOnGrid => getDataStorage( trim(keyDenData),                  &
                                   nga, ngb, ngc, RealMark )
   else
      call ErrorHandler("constructDensityOnGrid", &
                        "Unknown density type",density)
   endif
!
   n = getLocalNumAtoms()
   denOnGrid = ZERO
   gCounter = 0
   id = -1
!  ===================================================================
!  loop over grid points
!  ===================================================================
   if ( isUniformGrid ) then
      if ( isCharge ) then
         do kc = 1, ngc
            do jc = 1, ngb
               do ic = 1, nga
                  gCounter = gCounter + 1
                  vis_point => getInfoPointVis( gCounter, atomsOnPoint )
                  n_mult = nAtoms(ic,jc,kc)
!
                  do ig = 1, atomsOnPoint
!                    ==============================================
!                    loop over atoms on this grid point
!                    ==============================================
                     r(1:3) = vis_point%r(1:3)
                     id = vis_point%indexAtom
                     VP_point = vis_point%VP_Point
!
                     den_r = ZERO
                     do ia = 1, getLocalNumSpecies(id)
                        den_r = den_r + getChargeDensity( trim(keyDen), id, ia, r, jmax, &
                                                          VP_point, n_mult )*getLocalSpeciesContent(id,ia)
                     enddo
                     denOnGrid(ic,jc,kc) = denOnGrid(ic,jc,kc) + den_r/n_mult
!
                     vis_point => vis_point%nextAtom
!
                  enddo
               enddo
            enddo
         enddo
      else if ( isMoment ) then
         do kc = 1, ngc
            do jc = 1, ngb
               do ic = 1, nga
                  gCounter = gCounter + 1
                  vis_point => getInfoPointVis( gCounter, atomsOnPoint )
                  n_mult = nAtoms(ic,jc,kc)
!
                  do ig = 1, atomsOnPoint
!                    ==============================================
!                    loop over atoms on this grid point
!                    ==============================================
                     r(1:3) = vis_point%r(1:3)
                     id = vis_point%indexAtom
                     VP_point = vis_point%VP_Point
!
                     den_r = ZERO
                     do ia = 1, getLocalNumSpecies(id)
                        den_r = den_r + getMomentDensity( trim(keyDen), id, ia, r, jmax, &
                                                          VP_point, n_mult )*getLocalSpeciesContent(id,ia)
                     enddo
                     denOnGrid(ic,jc,kc) = denOnGrid(ic,jc,kc) + den_r/n_mult
!
                     vis_point => vis_point%nextAtom
!
                  enddo
               enddo
            enddo
         enddo
      else if ( isPotential ) then
         do kc = 1, ngc
            do jc = 1, ngb
               do ic = 1, nga
                  gCounter = gCounter + 1
                  vis_point => getInfoPointVis( gCounter, atomsOnPoint )
                  n_mult = nAtoms(ic,jc,kc)
!
                  do ig = 1, atomsOnPoint
!                    ==============================================
!                    loop over atoms on this grid point
!                    ==============================================
                     r(1:3) = vis_point%r(1:3)
                     id = vis_point%indexAtom
                     VP_point = vis_point%VP_Point
!
                     den_r = ZERO
                     do ia = 1, getLocalNumSpecies(id)
                        den_r = den_r + getNewPotential(id, ia, jmax, r, trim(keyDen), &
                                                        is, VP_point )*getLocalSpeciesContent(id,ia)
                     enddo
                     denOnGrid(ic,jc,kc) = denOnGrid(ic,jc,kc) +       &
                                           den_r/n_mult
!
                     vis_point => vis_point%nextAtom
!
                  enddo
               enddo
            enddo
         enddo
      else if ( isWave ) then
         do kc = 1, ngc
            do jc = 1, ngb
               do ic = 1, nga
                  gCounter = gCounter + 1
                  vis_point => getInfoPointVis( gCounter, atomsOnPoint )
                  n_mult = nAtoms(ic,jc,kc)
!
                  do ig = 1, atomsOnPoint
!                    ==============================================
!                    loop over atoms on this grid point
!                    ==============================================
                     r(1:3) = vis_point%r(1:3)
                     id = vis_point%indexAtom
                     VP_point = vis_point%VP_Point
!
                     den_r = ZERO
                     do ia = 1, getLocalNumSpecies(id)
!v2.0                   den_r = den_r + getWaveFunction(id, trim(keyDen), 7, r, &
!v2.0                                                   VP_point )*getLocalSpeciesContent(id,ia)
stop 'Under construction ...'
                     enddo
                     denOnGrid(ic,jc,kc) = denOnGrid(ic,jc,kc) +       &
                                           den_r/n_mult
!
                     vis_point => vis_point%nextAtom
!
                  enddo
               enddo
            enddo
         enddo
      endif
   else
      if ( isCharge ) then
         do kc = 1, ngc
            do jc = 1, ngb
               do ic = 1, nga
                  gCounter = gCounter + 1
                  fft_point => getInfoPointFFT( gCounter, atomsOnPoint )
                  n_mult = nAtoms(ic,jc,kc)
!
                  do ig = 1, atomsOnPoint
!                    ==============================================
!                    loop over atoms on this grid point
!                    ==============================================
                     r(1:3) = fft_point%r(1:3)
                     id = fft_point%indexAtom
                     VP_point = fft_point%VP_Point
!
                     den_r = ZERO
                     do ia = 1, getLocalNumSpecies(id)
                        den_r = den_r + getChargeDensity( trim(keyDen), id, ia, r, jmax, &
                                                          VP_point, n_mult )*getLocalSpeciesContent(id,ia)
                     enddo
                     denOnGrid(ic,jc,kc) = denOnGrid(ic,jc,kc) + den_r/n_mult
!
                     fft_point => fft_point%nextAtom
!
                  enddo
               enddo
            enddo
         enddo
      else if ( isMoment ) then
         do kc = 1, ngc
            do jc = 1, ngb
               do ic = 1, nga
                  gCounter = gCounter + 1
                  fft_point => getInfoPointFFT( gCounter, atomsOnPoint )
                  n_mult = nAtoms(ic,jc,kc)
!
                  do ig = 1, atomsOnPoint
!                    ==============================================
!                    loop over atoms on this grid point
!                    ==============================================
                     r(1:3) = fft_point%r(1:3)
                     id = fft_point%indexAtom
                     VP_point = fft_point%VP_Point
!
                     den_r = ZERO
                     do ia = 1, getLocalNumSpecies(id)
                        den_r = den_r + getMomentDensity( trim(keyDen), id, ia, r, jmax, &
                                                          VP_point, n_mult )*getLocalSpeciesContent(id,ia)
                     enddo
                     denOnGrid(ic,jc,kc) = denOnGrid(ic,jc,kc) + den_r/n_mult
!
                     fft_point => fft_point%nextAtom
!
                  enddo
               enddo
            enddo
         enddo
      else if (isPotential) then
         do kc = 1, ngc
            do jc = 1, ngb
               do ic = 1, nga
                  gCounter = gCounter + 1
                  fft_point => getInfoPointFFT( gCounter, atomsOnPoint )
                  n_mult = nAtoms(ic,jc,kc)
!
                  do ig = 1, atomsOnPoint
!                    ==============================================
!                    loop over atoms on this grid point
!                    ==============================================
                     r(1:3) = fft_point%r(1:3)
                     id = fft_point%indexAtom
                     VP_point = fft_point%VP_Point
!
                     den_r = ZERO
                     do ia = 1, getLocalNumSpecies(id)
                        den_r = den_r + getNewPotential(id, ia, jmax, r, trim(keyDen), &
                                                        is, VP_point )*getLocalSpeciesContent(id,ia)
                     enddo
                     denOnGrid(ic,jc,kc) = denOnGrid(ic,jc,kc) +       &
                                           den_r/n_mult
!
                     fft_point => fft_point%nextAtom
!
                  enddo
               enddo
            enddo
         enddo
      else if (isWave) then
         do kc = 1, ngc
            do jc = 1, ngb
               do ic = 1, nga
                  gCounter = gCounter + 1
                  fft_point => getInfoPointFFT( gCounter, atomsOnPoint )
                  n_mult = nAtoms(ic,jc,kc)
!
                  do ig = 1, atomsOnPoint
!                    ==============================================
!                    loop over atoms on this grid point
!                    ==============================================
                     r(1:3) = fft_point%r(1:3)
                     id = fft_point%indexAtom
                     VP_point = fft_point%VP_Point
!
                     den_r = ZERO
                     do ia = 1, getLocalNumSpecies(id)
!v2.0                   den_r = den_r + getWaveFunction(id, trim(keyDen), 7, r, &
!v2.0                                                   VP_point )*getLocalSpeciesContent(id,ia)
stop 'Under construction ...'
                     enddo
                     denOnGrid(ic,jc,kc) = denOnGrid(ic,jc,kc) +       &
                                           den_r/n_mult
!
                     fft_point => fft_point%nextAtom
!
                  enddo
               enddo
            enddo
         enddo
      endif
   endif
!
!  ==================================================================
!  Gather the density on Grid
!  ==================================================================
!
   call GlobalSumInGroup(GroupID, denOnGrid, nga, ngb, ngc)
!  ===================================================================
   end subroutine constructDensityOnGrid
