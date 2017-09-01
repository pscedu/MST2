!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printDensityOnGrid( density, den_type, lmax )
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler
!
   use MPPModule, only : MyPE
!
   use DataServiceCenterModule, only : createDataStorage,     &
                                       getDataStorage,        &
                                       isDataStorageExisting, &
                                       RealType, IntegerType, &
                                       RealMark, IntegerMark
!
   use PublicTypeDefinitionsModule, only : VisualGridStruct
!
   use VisualGridModule, only : getVisualGrid
   use OutputModule, only : getDensityPrintFlag,     &
                            getDensityPrintFile,     &
                            getDensityPrintFormat,   &
                            getPotentialPrintFlag,   &
                            getPotentialPrintFile,   &
                            getPotentialPrintFormat, &
                            getWavePrintFlag,        &
                            getWavePrintFile,        &
                            getWavePrintFormat

!
   implicit none
!
!  ==================================================================
!  Variables
!  ==================================================================
!
   character(len=*), intent(in) :: density
   character(len=*), intent(in) :: den_type
!
   integer (kind=IntKind), intent(in) :: lmax
!
   logical :: isCharge
   logical :: isPotential
   logical :: isWave
!
   character (len=90) :: fname
   character (len=3) :: Lchar
   character(len=20) :: keyDen
!
   character(len=50) :: keyAtomIndex, keyAtomNumber
   character(len=50) :: keyDenData
   integer (kind=IntKind) :: lenKey = 50
   integer (kind=IntKind) :: lenKeyDen
!
   integer (kind=IntKind) :: offset, jmax
   integer (kind=IntKind) :: nga, ngb, ngc, ng, i, j, k
   integer (kind=IntKind) :: printFormat, funit
   integer (kind=IntKind), pointer :: indexAtom(:,:,:), nAtoms(:,:,:)
!
   real (kind=RealKind) :: rg(3), sqrt_rg
   real (kind=RealKind), pointer :: denOnGrid(:,:,:)
!
   type (VisualGridStruct), pointer:: visualgp
!
   isCharge    = .false.
   isPotential = .false.
   isWave      = .false.
!
   offset  = 100
   keyDen(1:20) = " "
   keyAtomIndex(1:lenKey)     = " "
   keyAtomNumber(1:lenKey)    = " "
   keyDenData(1:lenKey)       = " "
!  ==================================================================
!  Body
!  ==================================================================
   visualgp => getVisualGrid()
   nga = visualgp%nga
   ngb = visualgp%ngb
   ngc = visualgp%ngc
   ng  = visualgp%ng
   jmax = (lmax+1)*(lmax+2)/2
!
   keyAtomIndex  = "Atom Index on Visual Mesh"
   if ( .not.isDataStorageExisting( trim(keyAtomIndex) ) ) then
      call ErrorHandler("printDensityOnGrid", &
                        "Visual Grid Atom Index not defined",density)
   endif
!
   keyAtomNumber = "Atom Number on Visual Mesh"
   if ( .not.isDataStorageExisting( trim(keyAtomNumber) ) ) then
      call ErrorHandler("printDensityOnGrid", &
                        "Visual Grid Atom Number not defined",density)
   endif
!
   if ( density(1:6)=="Charge" ) then
      isCharge = .true.
   else if ( density(1:9)=="Potential" ) then
      isPotential=.true.
   else if ( density(1:9)=="WaveFunct" ) then
      isWave=.true.
   else
      call ErrorHandler("printDensityOnGrid", &
                        "Unknown density type",density)
   endif
!
   indexAtom => getDataStorage( trim(keyAtomIndex),                   &
                                nga, ngb, ngc, IntegerMark)
   nAtoms => getDataStorage( trim(keyAtomNumber),                     &
                             nga, ngb, ngc, IntegerMark)
!
   lenkeyDen = len(den_type)
   keyDen(1:lenkeyDen) = den_type(1:lenkeyDen)
   if ( isCharge ) then
      printFormat = getDensityPrintFormat()
      fname = getDensityPrintFile(lmax,trim(keyDen))
!
      if ( jmax > 0 ) then
         write(Lchar,'(i3)') offset + lmax
         Lchar(1:1) = "L"
         keyDenData = keyDen(1:lenkeyDen)//" Electron Density on Visual Mesh "//Lchar
      else
         keyDenData = keyDen(1:lenkeyDen)//" Electron Density on Visual Mesh"
      endif
!
      denOnGrid => getDataStorage( trim(keyDenData),                  &
                                   nga, ngb, ngc, RealMark )

   else if ( isPotential ) then
      printFormat = getPotentialPrintFormat()
      fname = getPotentialPrintFile(lmax,trim(keyDen))
!
      if ( jmax > 0 ) then
         write(Lchar,'(i3)') offset + lmax
         Lchar(1:1) = "L"
         keyDenData = trim(keyDen(1:lenkeyDen))//" Potential on Visual Mesh "//Lchar
      else
         keyDenData = trim(keyDen(1:lenkeyDen))//" Potential on Visual Mesh"
      endif
!
      denOnGrid => getDataStorage( trim(keyDenData),                  &
                                   nga, ngb, ngc, RealMark )
   else if ( isWave ) then
      printFormat = getWavePrintFormat()
      fname = getWavePrintFile(lmax,trim(keyDen))
!
      if ( jmax > 0 ) then
         write(Lchar,'(i3)') offset + lmax
         Lchar(1:1) = "L"
         keyDenData = trim(keyDen(1:lenkeyDen))//" WaveFunct on Visual Mesh "//Lchar
      else
         keyDenData = trim(keyDen(1:lenkeyDen))//" WaveFunct on Visual Mesh"
      endif
!
      denOnGrid => getDataStorage( trim(keyDenData),                  &
                                   nga, ngb, ngc, RealMark )
   endif
!
   if ( MyPE==0 ) then
      funit = 11
      open(unit=funit, file=trim(fname), status='unknown')
      if ( printFormat == 0 ) then
!
         write(funit, *) "# xyz sqrt(r^2) format"
         do k = 1, visualgp%ngc
            do j = 1, visualgp%ngb
               do i = 1, visualgp%nga
                  rg = visualgp%step_a * visualgp%vec_a * (i-1) + &
                       visualgp%step_b * visualgp%vec_b * (j-1) + &
                       visualgp%step_c * visualgp%vec_c * (k-1) + &
                       visualgp%vec_o
                  sqrt_rg = sqrt(rg(1)**2+rg(2)**2+rg(3)**2)
                  WRITE(funit,'(5f15.8,2i5)') rg, sqrt_rg, &
                        denOnGrid(i,j,k), indexAtom(i,j,k), nAtoms(i,j,k)
               enddo
            enddo
         enddo
!
      else if (printFormat == 1) then
!
         write(funit, '(''# vtk DataFile Version 1.0'')')
         write(funit, '(''density on the grid example'')')
         write(funit, '(''ASCII'')')
         write(funit, *)
         write(funit, '(''DATASET STRUCTURED_GRID'')')
         write(funit, '(''DIMENSIONS '' , 3i5)') &
                          visualgp%nga, visualgp%ngb, visualgp%ngc
         write(funit, '(''POINTS '', i10, '' double'')') visualgp%ng
!
         do k = 1, visualgp%ngc
            do j = 1, visualgp%ngb
               do i = 1, visualgp%nga
                  rg = visualgp%step_a * visualgp%vec_a * (i-1) + &
                       visualgp%step_b * visualgp%vec_b * (j-1) + &
                       visualgp%step_c * visualgp%vec_c * (k-1) + &
                       visualgp%vec_o
                  WRITE(funit,'(3f10.6)') rg
               enddo
            enddo
         enddo
!
         write(funit, '(a, i10)') "POINT_DATA ", visualgp%ng
         write(funit, '(a)')      "SCALARS density double 1"
         write(funit, '(a)')      "LOOKUP_TABLE default"
!
         do k = 1, visualgp%ngc
            do j = 1, visualgp%ngb
               do i = 1, visualgp%nga
                  write (funit, '(f15.8)') denOnGrid(i,j,k)
               enddo
            enddo
         enddo
!
      endif
      call FlushFile(funit)
      close(funit)
   endif
!
!  ==================================================================
   end subroutine printDensityOnGrid
