!****************************************************************************
!   HEADING: MC3D I/O MODULE FOR UNSTRUCTURED GRID
!   AUTHOR: 
!   PURPOSE: THIS MODULE HANDLES THE INPUT/OUTPUT PROCESSES.
!   DATE: 
!****************************************************************************
MODULE mod_IO
USE mod_VARIABLES
USE mod_ROUTINES
USE RNG
IMPLICIT NONE
CONTAINS
!=====================================================================!
!=====================================================================!
SUBROUTINE readtable
IMPLICIT NONE

    WRITE(*,*) '========================================================================'
    WRITE(*,*) 'READING MATERIAL PROPERTIES...'

    OPEN(UNIT=LRGe,FILE="Ge_real_table.txt")
    OPEN(UNIT=LRSi,FILE="Si_real_table.txt")

    READ(LRGe,*) N_Ge1,N_Ge2
    READ(LRSi,*) N_Si1,N_Si2

    ALLOCATE( Ge_table(N_Ge1,N_Ge2) )
    ALLOCATE( Si_table(N_Si1,N_Si2) )
    ALLOCATE( Elimit(2,2) )
    READ(LRGe,*) Ge_table
    READ(LRSi,*) Si_table

    Elimit(1,1)=Ge_table(3,1)
    Elimit(2,1)=Ge_table(3,N_Ge2)
    dU_Ge=Ge_table(3,2)-Ge_table(3,1)
    Elimit(1,2)=Si_table(3,1)
    Elimit(2,2)=Si_table(3,N_Si2)
    dU_Si=Si_table(3,2)-Si_table(3,1)

    WRITE(*,*) '    # of Ge data = ',N_Ge2
    WRITE(*,*) '    # of Si data = ',N_Si2
    WRITE(*,*) '    Ge min energy (meV) = ',Elimit(1,1)
    WRITE(*,*) '    Ge MAX energy (meV) = ',Elimit(2,1)
    WRITE(*,*) '    dU_Ge    (meV) = ',dU_Ge
    WRITE(*,*) '    Si min energy (meV) = ',Elimit(1,2)
    WRITE(*,*) '    Si MAX energy (meV) = ',Elimit(2,2)
    WRITE(*,*) '    Si_Ge    (meV) = ',dU_Si

    CLOSE(LRGe)
    CLOSE(LRSi)
    
    WRITE(*,*) 'FINISH...'
    WRITE(*,*) '========================================================================'
    
END SUBROUTINE readtable
!============================================================================
SUBROUTINE gridINFO
IMPLICIT NONE
    
    WRITE(*,*) '========================================================================'
    WRITE(*,*) 'READING GRID INFORMATIONS...'
    
    OPEN(UNIT=LR1,FILE=gridfilename)
    
    READ(LR1,*) Propelem,Nelement,Nnodes,Ngrains
    ALLOCATE( Element(Propelem,Nelement),xyNodes(3,Nnodes),center(3,Nelement),Volume(Nelement) )
    ALLOCATE( normVec(3,4,Nelement),vecCoords(3,3,Nelement),invers_vec(3,3,Nelement) )
    ALLOCATE( GrainMt(Ngrains) )
    READ(LR1,*) GrainMt
    READ(LR1,*) Element
    READ(LR1,*) xyNodes
    READ(LR1,*) center
    READ(LR1,*) Volume
    READ(LR1,*) normVec
    READ(LR1,*) vecCoords
    READ(LR1,*) invers_vec
    READ(LR1,*) NBCL,NBCR,NBCyP,NBCyN,NBCzP,NBCzN
    ALLOCATE( BCelementL(NBCL),BCelementR(NBCR) )
    ALLOCATE( BCelementyP(NBCyP),BCelementyN(NBCyN) )
    ALLOCATE( BCelementzP(NBCzP),BCelementzN(NBCzN) )
    ALLOCATE( areaBCL(NBCL),areaBCR(NBCR) )
    READ(LR1,*) BCelementL,BCelementR
    READ(LR1,*) BCelementyP,BCelementyN
    READ(LR1,*) BCelementzP,BCelementzN
    READ(LR1,*) areaBCL,areaBCR

    dLdomain(1)=MAXVAL(xyNodes(1,:))-MINVAL(xyNodes(1,:))
    dLdomain(2)=MAXVAL(xyNodes(2,:))-MINVAL(xyNodes(2,:))
    dLdomain(3)=MAXVAL(xyNodes(3,:))-MINVAL(xyNodes(3,:))

    CLOSE(LR1)
    
    WRITE(*,*) 'FINISH...'
    WRITE(*,*) '========================================================================'
    
END SUBROUTINE gridINFO
!============================================================================
SUBROUTINE initialize
IMPLICIT NONE
INTEGER*4::i
REAL*8::rannum

    WRITE(*,*) '========================================================================'
    WRITE(*,*) 'INITIALIZING...'
    !-----------------------------------------------------------------!
    !-----------------------------------------------------------------!
    OPEN(LR2,FILE=inputfilename)
    
    READ(LR2,*) bundle,tStep,time0,iter0
    READ(LR2,*) DPP,DPPB
    READ(LR2,*) Nprop,Nph0,Nph

    WRITE(*,*) '    # of phonon bundles = ',Nph
    WRITE(*,*) '    bundle = ',bundle
    WRITE(*,*) '    DPP = ',DPP
    WRITE(*,*) '    DPPB = ',DPPB
    WRITE(*,*) '    Time Step (ps) = ',tStep
    WRITE(*,*) '    iter0 = ',iter0
    WRITE(*,*) '    time (ps) = ',time0

    ALLOCATE( Nnumcell(Nelement),Nbgcell(Nelement),ElemGroup(Nelement),Vgroup(Nelement) )
    ALLOCATE( dEcell(Nelement),dTemp(Nelement),dEunit(Nelement) )
    ALLOCATE( dEdiff(Nelement),dVunit(Nelement),MFP(Nelement) )
    ALLOCATE( phn(Nprop,Nph0),address(Nph0) )
    
    IF (cellInfoMethod.eq.2) ALLOCATE( tempEcell(Nelement) )
    
    address=0
    phn=0
    dEdiff=0
    READ(LR2,*) phn
    READ(LR2,*) dEdiff
    
    IF (WAY_DIR.gt.0) THEN
        
        ALLOCATE( dEheatfluxL(NBCL),dVinjectL(2),dEinjectL(2),NemitL(NBCL),dElostL(NBCL),qbdyL(NBCL) )
        ALLOCATE( dEheatfluxR(NBCR),dVinjectR(2),dEinjectR(2),NemitR(NBCR),dElostR(NBCR),qbdyR(NBCR) )
        
        dElostL=0
        dElostR=0
        qbdyL=0
        qbdyR=0
        qcenter=0
        
        CALL proc_BC( TL0,TR0 )
        
        Nmakeup=0
        
        IF (WAY_DIR.lt.3) THEN
        
            READ(LR2,*) Npool,Nmakeup
            
            IF (WAY_DIR.eq.1) THEN
                ALLOCATE( mlostL(NBCL),mlostR(NBCR),dPoolL(Npool,Nmakeup,NBCL),dPoolR(Npool,Nmakeup,NBCR) )
                ALLOCATE( truePoolL(NBCL),truePoolR(NBCR) )
            ELSE IF (WAY_DIR.eq.2) THEN
                ALLOCATE( mlostL(2),mlostR(2),dPoolL(Npool,Nmakeup,2),dPoolR(Npool,Nmakeup,2) )
                ALLOCATE( truePoolL(2),truePoolR(2) )
            ENDIF
            
            mlostL=0
            mlostR=0
            truePoolL=0
            truePoolR=0
            
	        READ(LR2,*) dPoolL,dPoolR
	        WRITE(*,*) '    Nmakeup = ',Nmakeup
	    ELSE
	        ALLOCATE( truePoolL(1),truePoolR(1) )
	        truePoolL=0
            truePoolR=0
        ENDIF
    ENDIF
    
    CLOSE(LR2)
    !-----------------------------------------------------------------!
    !-----------------------------------------------------------------!
    IF (WAY_DIR.eq.1) THEN
        FORALL(i=1:NBCL) mlostL(i)=COUNT(dPoolL(5,:,i).gt.0)
        FORALL(i=1:NBCR) mlostR(i)=COUNT(dPoolR(5,:,i).gt.0)   
        FORALL(i=1:NBCL,mlostL(i).eq.Nmakeup) truePoolL(i)=1
        FORALL(i=1:NBCR,mlostR(i).eq.Nmakeup) truePoolR(i)=1
    ELSE IF (WAY_DIR.eq.2) THEN
        FORALL(i=1:2) mlostL(i)=COUNT(dPoolL(5,:,i).gt.0)
        FORALL(i=1:2) mlostR(i)=COUNT(dPoolR(5,:,i).gt.0)   
        FORALL(i=1:2,mlostL(i).eq.Nmakeup) truePoolL(i)=1
        FORALL(i=1:2,mlostR(i).eq.Nmakeup) truePoolR(i)=1
    ENDIF
    !-----------------------------------------------------------------!
    !-----------------------------------------------------------------!
    OPEN(LR2,FILE=groupfilename)
    READ(LR2,*) Vgroup
    READ(LR2,*) ElemGroup(:)%N
    DO i=1,Nelement
        ALLOCATE( ElemGroup(i)%array(ElemGroup(i)%N) )
        IF (ElemGroup(i)%N.eq.1) THEN
            ElemGroup(i)%array=i
        ELSE
            READ(LR2,*) ElemGroup(i)%array
        ENDIF
    ENDDO
    CLOSE(LR2)
    !-----------------------------------------------------------------!
    !-----------------------------------------------------------------!
    CALL RANDOM_SEED()
    ALLOCATE( RNseedMP(0:Ncore-1) )
    DO i=0,Ncore-1
        CALL RANDOM_NUMBER(rannum)
        rannum=rannum*1d8
        CALL RNG_SEED(RNseedMP(i),INT(rannum))
    ENDDO
    CALL RANDOM_NUMBER(rannum)
    rannum=rannum*1d8
    CALL RNG_SEED(RNseed0,INT(rannum))
    !-----------------------------------------------------------------!
    !-----------------------------------------------------------------!
    CALL proc_reorder_and_Ecell
    CALL cellinfo

    WRITE(*,*) '    MAX & min # of phonon bundles in a cell :', MAXVAL(Nnumcell),MINVAL(Nnumcell)
    WRITE(*,*) '    MAX T & min T :',MAXVAL(dTemp),MINVAL(dTemp)
    !-----------------------------------------------------------------!
    !-----------------------------------------------------------------!
    ALLOCATE( Tz(Nelement))
    Tz=0
    !-----------------------------------------------------------------!
    !-----------------------------------------------------------------!
    WRITE(*,*) 'FINISH...'
    
END SUBROUTINE initialize
!============================================================================
SUBROUTINE restart
IMPLICIT NONE
INTEGER*4::i

    OPEN(LW3,file=restartfilename)    
    WRITE(LW3,*) bundle,tStep,time,iter
    WRITE(LW3,*) DPP,DPPB
    WRITE(LW3,*) Nprop,Nph0,Nph
    WRITE(LW3,*) phn
    WRITE(LW3,*) dEdiff
    IF (WAY_DIR.gt.0) THEN
        WRITE(LW3,*) Npool,Nmakeup
        WRITE(LW3,*) dPoolL,dPoolR
    ENDIF
    CLOSE(LW3)
    
END SUBROUTINE restart
!============================================================================
END MODULE mod_IO
