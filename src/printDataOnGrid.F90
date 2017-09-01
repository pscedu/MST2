 
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printDataOnGrid(data_name, na, nb, nc, data )
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MPPModule, only : MyPE
!
   use PublicTypeDefinitionsModule, only : VisualGridStruct
!
   use VisualGridModule, only : getVisualGrid
!
   use OutputModule, only : getDensityPrintFlag, &
                            getDensityPrintFile, &
                            getDensityPrintFormat
!
   implicit none
!  ==================================================================
!  Variables
!  ==================================================================
!
   character (len=*), intent(in) :: data_name
   integer (kind=IntKind), intent(in) :: na, nb, nc
   real (kind=RealKind), target :: data(na,nb,nc)
   character (len=90) :: fname
!
   integer (kind=IntKind) :: printFormat, funit
   integer (kind=IntKind) :: i, j, k
!
   real (kind=RealKind) :: rg(3), sqrt_rg
!
   type (VisualGridStruct), pointer:: visualgp
!
!  ==================================================================
!  Body
!  ==================================================================
   visualgp => getVisualGrid()
   funit = 11
   fname = getDensityPrintFile(-1,data_name)
!  ==================================================================
!  Print the density
!  ==================================================================
   printFormat = getDensityPrintFormat()
!
   if ( MyPE==0 ) then
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
                  WRITE(funit,'(5f15.8)') rg, sqrt_rg, data(i,j,k)
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
                  write (funit, '(f19.12)') data(i,j,k)
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
   end subroutine printDataOnGrid