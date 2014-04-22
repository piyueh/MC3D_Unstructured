!****************************************************************************
!   HEADING: MC3D MAIN PROGRAM FOR UNSTRUCTURED MESH
!   LIMINATION: VERTICAL BOUNDARIES AT BOTH ENDS REQUIRED
!   AUTHOR: 
!   PURPOSE: TO SIMULATE THERMAL CONDUCTIVITY OF NANOSTRUCTURED MATERIAL
!   DATE: 
!****************************************************************************
PROGRAM MC3D
USE mod_VARIABLES
USE mod_heatcontrol
USE mod_ROUTINES
USE mod_ADVANCE
USE mod_CreateDelete
USE mod_IO
USE omp_lib
IMPLICIT NONE
INTEGER*4::iterations,i
REAL*8::bg,ed,cputime
REAL*8::bg2,ed2,t(5)
LOGICAL::true

    WRITE(*,*) 'INPUT CASE NAME (NO .''TXT'',NO ''_GRIDFILE.TXT'',NO ''_INITIAL.TXT'') = '
    READ(*,"(A72)") casename
    
    inputfilename=casename(1:LEN_TRIM(casename))//'_initial.txt'
    gridfilename=casename(1:LEN_TRIM(casename))//'_gridfile.txt'
    groupfilename=casename(1:LEN_TRIM(casename))//'_group.txt'
    restartfilename=casename(1:LEN_TRIM(casename))//'_restart.txt'
    
    INQUIRE(FILE=restartfilename,EXIST=true)
    IF (true) THEN
        WRITE(*,*) 'THE RESTART FILE IS EXIST, USE RESTART FILE ?(1:YES  2:NO)'
        READ(*,*) i
        IF (i.eq.1) inputfilename=casename(1:LEN_TRIM(casename))//'_restart.txt'
    ENDIF
    

    WRITE(*,*) 'ENTER HOW MANY ITERATIONS WILL OUTPUT ONCE:'
    READ(*,*) nsteady

    iterations=10000000
    backup=10000

    WRITE(*,*) 'ENTER CELL INFORMATION CALCULATION METHOD:'
    WRITE(*,*) '    1:ORIGINAL'
    WRITE(*,*) '    2:LOCAL SMOOTHING'
    READ(*,*) cellinfoMethod

    WRITE(*,*) 'WAY TO ASSIGN THE DIRECTIONS OF INCIDENT PHONONS'
    WRITE(*,*) '    0: NO INJECTION OF PHONONS FROM BOUNDARIES'
    WRITE(*,*) '    1: CELL-TO-CELL PERIODICALLY ASSIGNED '
    WRITE(*,*) '    2: MATERIAL-TO-MATERIAL PERIODICALLY ASSIGNED '
    WRITE(*,*) '    3: RANDOMLY ASSIGNED'
    READ(*,*) WAY_DIR

    IF (WAY_DIR.gt.0) THEN
        WRITE(*,*) 'THE SPECIFIED TEMPERATURES (K): TL,TR = '
        READ*, TL0,TR0
    ENDIF

    WRITE(*,*) 'ENTER THE NUMBER OF CPU CORES TO BE USED :'
    READ(*,*) Ncore

    SELECT CASE(Ncore)
        CASE(1)
            CALL OMP_SET_NUM_THREADS(1)
        CASE(2)
            CALL OMP_SET_NUM_THREADS(2)
        CASE(3)
            CALL OMP_SET_NUM_THREADS(3)
        CASE(4)
            CALL OMP_SET_NUM_THREADS(4)
        CASE(5)
            CALL OMP_SET_NUM_THREADS(5)
        CASE(6)
            CALL OMP_SET_NUM_THREADS(6)
        CASE(7)
            CALL OMP_SET_NUM_THREADS(7)
        CASE(8)
            CALL OMP_SET_NUM_THREADS(8)
    ENDSELECT
    !-----------------------------------------------------------------!
    !-----------------------------------------------------------------!

    CALL readtable
    CALL gridINFO
    CALL initialize

    WRITE(*,*) 'PREPROCESSING DONE'
    WRITE(*,*) '========================================================================'
    
    !-----------------------------------------------------------------!
    !-----------------------------------------------------------------!
    
    time = time0
    cputime=0
    t=0
    
    DO iter=iter0+1,iterations
    
        bg=OMP_GET_WTIME()
        
        !*****************************************************************!
        bg2=OMP_GET_WTIME()
        !$OMP PARALLEL DEFAULT(SHARED),PRIVATE(i,CPUID)
        CPUID=OMP_GET_THREAD_NUM()
        !$OMP DO SCHEDULE(DYNAMIC)
        DO i=1,Nph0
            IF (phn(7,i).gt.0) THEN
                CALL proc_advection(phn(:,i),tStep,1,RNseedMP(CPUID))
            ENDIF
        ENDDO
        !$OMP END DO
        !$OMP END PARALLEL
        ed2=OMP_GET_WTIME()
        t(1)=ed2-bg2+t(1)
        
        bg2=OMP_GET_WTIME()
        IF (WAY_DIR.gt.0) CALL proc_heatcontrol
        ed2=OMP_GET_WTIME()
        t(2)=ed2-bg2+t(2)
        
        bg2=OMP_GET_WTIME()
        CALL proc_reorder_and_Ecell
        ed2=OMP_GET_WTIME()
        t(3)=ed2-bg2+t(3)
        
        bg2=OMP_GET_WTIME()
        CALL proc_createdelete
        ed2=OMP_GET_WTIME()
        t(4)=ed2-bg2+t(4)
        
        bg2=OMP_GET_WTIME()
        CALL cellinfo
        ed2=OMP_GET_WTIME()
        t(5)=ed2-bg2+t(5)
        !*****************************************************************!
    
        Tz=Tz+dTemp
        time=time+tStep
    
        IF (MOD(iter,nsteady).eq.0) THEN
            WRITE(outputfilename,"((I7.7),'.txt')") iter
            OPEN(LW1,file=outputfilename)
            WRITE(LW1,*) time,tStep
            WRITE(LW1,*) qcenter/DBLE(nsteady)
            WRITE(LW1,*) SUM(dEheatfluxL)-SUM(qbdyL)/DBLE(nsteady),SUM(qbdyR)/DBLE(nsteady)-SUM(dEheatfluxR)
            WRITE(LW1,*) dEheatfluxL-qbdyL/DBLE(nsteady),qbdyR/DBLE(nsteady)-dEheatfluxR
            WRITE(LW1,*) Tz/DBLE(nsteady)
            CLOSE(LW1) 
            Tz=0
            qbdyL=0
   	        qbdyR=0
   	        qcenter=0
        ENDIF
    
        IF (iter.eq.iterations.or.MOD(iter,backup).eq.0) CALL restart
        
        ed=OMP_GET_WTIME()
        cputime=cputime+ed-bg
    
        WRITE(*,*) 'iter',iter,'time',time,'CPU_TIME',cputime/DBLE(iter-iter0)
        WRITE(*,"(5(F13.11,1X))") t/DBLE(iter-iter0)
    ENDDO
    !-----------------------------------------------------------------!
    !-----------------------------------------------------------------!
    DEALLOCATE( RNseedMP )
    DEALLOCATE( Tz,qbdyL,qbdyR )
    DEALLOCATE( Ge_table,Si_table )
    DEALLOCATE( Element,GrainMt,Volume,xyNodes,BCelementL,BCelementR,BCelementyN,BCelementyP,BCelementzN,BCelementzP )
    DEALLOCATE( Nnumcell,Nbgcell,dEcell,dEdiff,dEunit,dVunit,MFP,dTemp,phn,address )
    IF (WAY_DIR.gt.0) THEN
        DEALLOCATE( dEheatfluxL,dEinjectL,dVinjectL,NemitL,dElostL )
        DEALLOCATE( dEheatfluxR,dEinjectR,dVinjectR,NemitR,dElostR )
    ENDIF
    IF (WAY_DIR.eq.1) DEALLOCATE( dPoolL,dPoolR,mlostL,mlostR )

END PROGRAM MC3D
!============================================================================