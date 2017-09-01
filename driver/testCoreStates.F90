program testCoreStates
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use SystemModule, only : getNumAtoms
!
   use ScfDataModule, only : n_spin_pola, n_spin_cant
   use ScfDataModule, only : istop, EvBottom, isNonRelativisticCore
!
   use AtomModule, only : getPhiLmax, getStepFuncLmax, getPotLmax
!
   use MathParamModule, only : ZERO, TEN2m6, ONE, CZERO, SQRTm1, TEN2m8, PI4, PI
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : endSphericalHarmonics
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use CoreStatesModule, only : initCoreStates, calCoreStates, endCoreStates
   use CoreStatesModule, only : readCoreStates, printCoreStates
!
   use PolyhedraModule, only : initPolyhedra, endPolyhedra
   use PolyhedraModule, only : getVolume, getInscrSphVolume
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use OutputModule, only : getStandardOutputLevel
!
   use Atom2ProcModule, only : getLocalNumAtoms
!
   use PotentialModule, only : initPotential, endPotential
   use PotentialModule, only : readPotential
!
   implicit   none
!
   integer (kind=IntKind) :: iprint = 0
   integer (kind=IntKind) :: NumAtoms, LocalNumAtoms
   integer (kind=IntKind) :: lmax_max
   integer (kind=IntKind) :: id, ig
!
   integer (kind=IntKind), allocatable :: atom_print_level(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:), lmax_step(:)
!
!  -------------------------------------------------------------------
   call startProcess()
   NumAtoms = getNumAtoms()
   LocalNumAtoms=getLocalNumAtoms()
!  -------------------------------------------------------------------
!
   allocate(atom_print_level(1:LocalNumAtoms))
   allocate(lmax_step(LocalNumAtoms),lmax_pot(LocalNumAtoms))
!
   lmax_max = 0
   do id = 1, LocalNumAtoms
      lmax_pot(id) = getPotLmax(id)
      lmax_step(id)  = getStepFuncLmax(id)
      lmax_max = max(lmax_max,2*getPhiLmax(id),lmax_step(id))
      atom_print_level(id) = getStandardOutputLevel(id)
   enddo
!
!  -------------------------------------------------------------------
   call initSphericalHarmonics(2*lmax_max)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initilize GauntFactors:
!  ===================================================================
   call initGauntFactors(lmax_max,istop,iprint)
!  -------------------------------------------------------------------
!
!  -------------------------------------------------------------------
   call setupRadGridAndCell(NumAtoms,lmax_max)
!  -------------------------------------------------------------------
   call initPotential(LocalNumAtoms,lmax_pot,lmax_step,               &
                      n_spin_pola,n_spin_cant,istop,atom_print_level)
!  -------------------------------------------------------------------
   call initCoreStates(LocalNumAtoms,EvBottom,n_spin_pola,         &
                       isNonRelativisticCore(),istop,iprint)
!  ===================================================================
!  read potential data
!  -------------------------------------------------------------------
   call readPotential()
!  -------------------------------------------------------------------
!
!  ===================================================================
!  read core states and density data of local atoms from file
!  ===================================================================
   call readCoreStates()
!  -------------------------------------------------------------------
   call calCoreStates()
!  -------------------------------------------------------------------
   do id = 1,LocalNumAtoms
      if (atom_print_level(id) >= 0) then
!        -------------------------------------------------------------
         call printCoreStates(id)
!        -------------------------------------------------------------
      endif
   enddo
!
   deallocate(atom_print_level,lmax_step,lmax_pot)
!
!  -------------------------------------------------------------------
   call endCoreStates()
   call endPotential()
   call endGauntFactors()
   call endSphericalHarmonics()
   call finishProcess()
!  -------------------------------------------------------------------
!
end program testCoreStates
