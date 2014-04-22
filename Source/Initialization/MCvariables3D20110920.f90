!****************************************************************************
!   HEADING: MC3D VARIABLES MODULE FOR UNSTRUCTURED GRID
!   AUTHOR: 
!   PURPOSE: ALL THE GLOBAL CONSTANTS AND VARIABLES ARE DEFINED.
!   DATE : 
!****************************************************************************
MODULE mod_VARIABLES
IMPLICIT NONE

!----fixed PARAMETERs
REAL*8,PARAMETER::M_PI=3.14159265358979323846D0
REAL*8,PARAMETER::M_PI_2=6.28318530717958647692D0
REAL*8,PARAMETER::M_PI_half=1.57079632679489661923D0
REAL*8,PARAMETER::M_PI_3=9.4247779607693797153879301498385D0
REAL*8,PARAMETER::kB=8.617d-2,barh=6.5822d-1 ! meV/K, meV.ps
REAL*8,PARAMETER::zero_tol=0.00000000001
INTEGER*4,PARAMETER::LRSi=110,LRGe=120,LR1=130,LR2=140,LW1=150,LW2=160,LW3=170,ER=100

!----FUNCTIONS
INTERFACE
    FUNCTION CROSS_PRODUCT(A,B)
    IMPLICIT NONE
    REAL*8::A(3),B(3)
    REAL*8::CROSS_PRODUCT(3)
    END FUNCTION
END INTERFACE

!****MATERIAL TABLE INFORMATION
INTEGER*4::N_Ge1,N_Ge2,N_Si1,N_Si2
REAL*8::dU_Ge,dU_Si,Ge_start,Si_start,rho(2)=(/5.323d3,2.329d3/) !rho_Ge=5.323d3,rho_Si=2.329d3 ! kg/m^3
REAL*8,ALLOCATABLE::Ge_table(:,:),Si_table(:,:)

!****GRID INFORMATION
INTEGER*4:: Nnodes, Nelement, Ngrains
INTEGER*4:: Propelem
INTEGER*4,ALLOCATABLE:: Element(:,:),GrainMt(:)
INTEGER*4,ALLOCATABLE:: Nnumcell(:),Nbgcell(:)
REAL*8,ALLOCATABLE::xyNodes(:,:)
REAL*8,ALLOCATABLE::Volume(:),center(:,:),normVec(:,:,:),vecCoords(:,:,:),invers_vec(:,:,:)
REAL*8,ALLOCATABLE::dEcell(:),dEunit(:),dVunit(:),MFP(:),dTemp(:),dEdiff(:)
REAL*8::dPP,dPPB,TL0,TR0
REAL*8:: dLdomain(3)

TYPE::cell
    INTEGER*4::N
    INTEGER*4,ALLOCATABLE::array(:)
END TYPE cell
TYPE(cell),ALLOCATABLE:: ElemGroup(:)
REAL*8,ALLOCATABLE::Vgroup(:)
!-------------------------------------------------------------------------------
!Nnodes, Nelement : number of nodes and elements in computational domain
!Popelem : number of properties of each element, =9 by default
!-------------------------------------------------------------------------------
!Element( 1st~9th properties , cell label)
!       1st : material of this element
!       2nd~5th: global labels of its four neighbors:
!                2 is its neighbor of face that consisted of node 6¡B7¡B8
!                3 is its neighbor of face that consisted of node 6¡B7¡B9
!                4 is its neighbor of face that consisted of node 6¡B8¡B9
!                5 is its neighbor of face that consisted of node 7¡B8¡B9
!                   = Nelement+1 means a periodic b.c. in the first direction
!                   = Nelement+2 means a periodic b.c. in the second direction
!                   = Nelement+3 means a periodic b.c. in the third direction
!                   = Nelement+4 means an adiabatic b.c. required
!                   = Nelement+5 means there is a prescribed heat flux 
!       6th~9th: global labels of its four nodes
!
!ATTENTION!! If the element is located on boundary, then node Element(6,element label) must on boundary surface.
!
!-------------------------------------------------------------------------------
!Nnumcell : number of phonons in each element
!Nbgcell : beging label of phonon of each element in matrix address
!xyNodes : xyz coordinates of each node
!Volume, center : volume and center of each element
!-------------------------------------------------------------------------------
!normVec( xyz coords. , 1st~4th vector , element label )
!       1 is normal unit vector of face 2 of the cell
!       2 is normal unit vector of face 3 of the cell
!       3 is normal unit vector of face 4 of the cell
!       4 is normal unit vector of face 5 of the cell
!-------------------------------------------------------------------------------
!vecCoords( xyz coords. ,1st~3rd vector , element label )
!       1 is vector from node Element(6,element label) to Element(7,element label)
!       2 is vector from node Element(6,element label) to Element(8,element label)
!       3 is vector from node Element(6,element label) to Element(9,element label)
!-------------------------------------------------------------------------------
!invers_vec(:,:,element label) : inverse matrix of vecCoords(:,:,element label), 
!                                for solving linear eqs. in subroutine Find_Periodic_Neighbor
!-------------------------------------------------------------------------------
!dEcell,dEunit,dVunit,MFP,dTemp,dEdiff:
!energy, phonon density, phonon velocity, mean free path,temperature,energy difference of each element
!-------------------------------------------------------------------------------
!
!****BOUNDARY CELL INFORMATION
INTEGER*4:: NBCL,NBCR,NBCyP,NBCyN,NBCzP,NBCzN 
INTEGER*4,ALLOCATABLE:: BCelementL(:),BCelementR(:),BCelementyP(:),BCelementyN(:),BCelementzP(:),BCelementzN(:)
REAL*8,ALLOCATABLE:: areaBCL(:),areaBCR(:)
!-------------------------------------------------------------------------------
!NBCL,NBCR,NBCyP,NBCyN,NBCzP,NBCzN:
!       the number of cells which located on each boundary
!-------------------------------------------------------------------------------
!BCelementL,BCelementR,BCelementyP,BCelementyN,BCelementzP,BCelementzN 
!       global label of boundary cell
!
!****heatcontrol information
INTEGER*4,ALLOCATABLE::mlostL(:),mlostR(:),NemitL(:),NemitR(:)
INTEGER*4::WAY_DIR,Npool,Nmakeup,truePoolL,truePoolR
!-------------------------------------------------------------------------------
! Npool=7 for cell-to-cell periodically assign
!      =5 for material-to material periodically assign
!       1 is remaining time
!       2~4 is moving direction
!       5 is material of the phonon
!       6~7 is the y¡Bz location of the phonon (only cell-to-cell)
!-------------------------------------------------------------------------------
REAL*8,ALLOCATABLE::dEinjectL(:),dVinjectL(:),dEheatfluxL(:),dElostL(:),dtheat(:)
REAL*8,ALLOCATABLE::dPoolL(:,:,:),dPoolR(:,:,:),dEinjectR(:),dVinjectR(:),dEheatfluxR(:),dElostR(:)

!****phonon information
INTEGER*4::Nph0,Nph,heatNph,Nprop
INTEGER*4, ALLOCATABLE::address(:)
REAL*8, ALLOCATABLE::phn(:,:),phnheat(:,:)


!****simulation parameters
INTEGER*4::nsteady,iter,iter0,cellinfoMethod
REAL*8:: dt,bundle(2),time0,time

!****output data
REAL*8::qcenter
REAL*8, ALLOCATABLE:: Tz(:),qbdyL(:),qbdyR(:)
CHARACTER*40::inputfilename,outputfilename,restartfilename,gridfilename,outputfilename2
!===================================================================================
! 3D simulation:: iNprop=10
! 1: x, 2: y, 3: z, 4: vx, 5: vy, 6:vz, 7: energy, 8: velocity, 9: energy_material, 10: element 
!-----------------------------------------------------------------------------------
END MODULE mod_VARIABLES