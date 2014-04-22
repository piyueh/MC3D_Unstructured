!****************************************************************************
!   HEADING: MC3D HEAT CONTROL MODULE FOR UNSTRUCTURED GRID
!   AUTHOR: 
!   PURPOSE: THIS MODULE MAINTAIN A CONSTANT HEAT FLOW RATE ALONG ONE DIRECTION.
!   DATE : 
!****************************************************************************
MODULE mod_HEATCONTROL
USE mod_VARIABLES
USE mod_ROUTINES
USE mod_ADVANCE
USE omp_lib
IMPLICIT NONE
!****************************************************************************
CONTAINS
!============================================================================
!左右邊界必須為垂直面邊界
SUBROUTINE proc_heatcontrol
IMPLICIT NONE
INTEGER*4:: s,k,m,j,bg,ed,cc
REAL*8:: tmp,rannum1
REAL*8,ALLOCATABLE:: rannum(:,:)
INTEGER*4,ALLOCATABLE:: mL(:,:),mR(:,:)

    dElostL=dElostL+dEheatfluxL !dElost:energy that should be injected
    dElostR=dElostR+dEheatfluxR
    NemitL=0
    NemitR=0

    k=2*MAXVAL(dElostL)/MINVAL(dEinjectL)
    s=2*MAXVAL(dElostR)/MINVAL(dEinjectR)
    ALLOCATE( mL(k,NBCL),mR(s,NBCR) )

    IF (WAY_DIR.eq.1) THEN
        DO j=1,NBCL
            IF (truePoolL(j).eq.0) THEN
                IF (dPoolL(5,Nmakeup,j).gt.0) truePoolL(j)=1
            ENDIF
        ENDDO
        DO j=1,NBCR
            IF (truePoolR(j).eq.0) THEN
                IF (dPoolR(5,Nmakeup,j).gt.0) truePoolR(j)=1
            ENDIF
        ENDDO
    ELSE IF (WAY_DIR.eq.2) THEN
        DO j=1,2
            IF (truePoolL(j).eq.0) THEN
                IF (dPoolL(5,Nmakeup,j).gt.0) truePoolL(j)=1
            ENDIF
            IF (truePoolR(j).eq.0) THEN
                IF (dPoolR(5,Nmakeup,j).gt.0) truePoolR(j)=1
            ENDIF
        ENDDO
    ELSE
        truePoolL=0
        truePoolR=0
    ENDIF
!------------------------------------------------------------------------------
    DO j=1,NBCL
        IF (WAY_DIR.eq.1) THEN
            s=j
        ELSE IF (WAY_DIR.eq.2) THEN
            s=GrainMt(Element(1,BCelementL(j)))
        ELSE
            s=1
        ENDIF
        
        IF (truePoolL(s).eq.1) THEN
            tmp=dEinjectL(GrainMt(Element(1,BCelementL(j))))/2d0 !half phonon energy
            cc=0
            DO WHILE(dElostL(j).ge.tmp)
                CALL RAN_NUM(RNseed0,rannum1)
	            m=MIN( INT(Nmakeup*rannum1)+1,Nmakeup )
                NemitL(j)=NemitL(j)+1
		        mL(NemitL(j),j)=m
		        dElostL(j)=dElostL(j)-dEinjectL(dPoolL(5,m,s))
		        cc=cc+1
		        IF (cc.gt.1000000) THEN
	                WRITE(*,*) 'WRONG IN SUBROUTINE proc_heatcontrol 1!!'
	                WRITE(*,*) 'j,dElostL(j),tmp = ',j,dElostL(j),tmp
	                PAUSE '無限迴圈'
	            ENDIF
            ENDDO
        ELSE
            NemitL(j)=MAX( 0, INT(dElostL(j)/dEinjectL(GrainMt(Element(1,BCelementL(j))))+0.5d0) )
            dElostL(j)=dElostL(j)-DBLE(NemitL(j))*dEinjectL(GrainMt(Element(1,BCelementL(j))))
        ENDIF
    ENDDO

    DO j=1,NBCR
        IF (WAY_DIR.eq.1) THEN
            s=j
        ELSE IF (WAY_DIR.eq.2) THEN
            s=GrainMt(Element(1,BCelementR(j)))
        ELSE
            s=1
        ENDIF
        
        IF (truePoolR(s).eq.1) THEN
            tmp=dEinjectR(GrainMt(Element(1,BCelementR(j))))/2d0
            cc=0
            DO WHILE(dElostR(j).ge.tmp)
                CALL RAN_NUM(RNseed0,rannum1)
	            m=MIN( INT(Nmakeup*rannum1)+1,Nmakeup )
                NemitR(j)=NemitR(j)+1
	            mR(NemitR(j),j)=m
	            dElostR(j)=dElostR(j)-dEinjectR(dPoolR(5,m,s))
	            cc=cc+1
	            IF (cc.gt.1000000) THEN
	                WRITE(*,*) 'WRONG IN SUBROUTINE proc_heatcontrol 2!!'
	                WRITE(*,*) 'j,dElostR(j),tmp = ',j,dElostR(j),tmp
	                PAUSE '無限迴圈'
	            ENDIF
            ENDDO
        ELSE
            NemitR(j)=MAX( 0, INT(dElostR(j)/dEinjectR(GrainMt(Element(1,BCelementR(j))))+0.5d0) )
            dElostR(j)=dElostR(j)-DBLE(NemitR(j))*dEinjectR(GrainMt(Element(1,BCelementR(j))))
        ENDIF
    ENDDO
!------------------------------------------------------------------------------
    heatNph = SUM(NemitL)+SUM(NemitR) !heatcontrol 要增加至DOMAIN的總聲子數
    ALLOCATE( phnheat(Nprop,heatNph),dtheat(heatNph) )
    
    DO j=1,NBCL+NBCR
        IF (j.le.NBCL) THEN
            IF (NemitL(j).gt.0) THEN
                bg=SUM(NemitL(1:j-1))+1
                ed=SUM(NemitL(1:j))
                s=GrainMt(Element(1,BCelementL(j)))
	            IF (truePoolL(j).eq.1) THEN ! First way to determine the properties of emitted phonons (periodically)
	                IF (WAY_DIR.eq.1) THEN
	                    k=j
                        phnheat(2,bg:ed)=dPoolL(6,mL(1:NemitL(j),j),k)
                        phnheat(3,bg:ed)=dPoolL(7,mL(1:NemitL(j),j),k)	                    
	                ELSE IF (WAY_DIR.eq.2) THEN
	                    k=s
	                    ALLOCATE( rannum(NemitL(j),2) )
                        CALL RAN_NUM(RNseed0,rannum)
                        rannum(:,2)=rannum(:,2)*(1d0-rannum(:,1))
	                    phnheat(2,bg:ed)=xynodes(2,Element(6,BCelementL(j)))+rannum(:,1)*vecCoords(2,1,BCelementL(j))+rannum(:,2)*vecCoords(2,2,BCelementL(j))
                        phnheat(3,bg:ed)=xynodes(3,Element(6,BCelementL(j)))+rannum(:,1)*vecCoords(3,1,BCelementL(j))+rannum(:,2)*vecCoords(3,2,BCelementL(j))
                        DEALLOCATE( rannum )
	                ENDIF
                    phnheat(1,bg:ed)=0d0
                    phnheat(4,bg:ed)=dPoolL(2,mL(1:NemitL(j),j),k)*dVinjectL(s)
                    phnheat(5,bg:ed)=dPoolL(3,mL(1:NemitL(j),j),k)*dVinjectL(s)
                    phnheat(6,bg:ed)=dPoolL(4,mL(1:NemitL(j),j),k)*dVinjectL(s)
	                phnheat(7,bg:ed)=dEinjectL(dPoolL(5,mL(1:NemitL(j),j),k))
                    phnheat(8,bg:ed)=dVinjectL(s)
	                phnheat(9,bg:ed)=dPoolL(5,mL(1:NemitL(j),j),k)
                    phnheat(10,bg:ed)=BCelementL(j)
                    dtheat(bg:ed)=dPoolL(1,mL(1:NemitL(j),j),k)
	                
                ELSE	   ! Second way to determine the properties of emitted phonons (randomly)
	                ALLOCATE( rannum(NemitL(j),5) )
                    CALL RAN_NUM(RNseed0,rannum)
                    rannum(:,2)=rannum(:,2)*(1d0-rannum(:,1))
	                phnheat(1,bg:ed)=0
                    phnheat(2,bg:ed)=xynodes(2,Element(6,BCelementL(j)))+rannum(:,1)*vecCoords(2,1,BCelementL(j))+rannum(:,2)*vecCoords(2,2,BCelementL(j))
                    phnheat(3,bg:ed)=xynodes(3,Element(6,BCelementL(j)))+rannum(:,1)*vecCoords(3,1,BCelementL(j))+rannum(:,2)*vecCoords(3,2,BCelementL(j))
                    phnheat(4,bg:ed)=DSQRT(rannum(:,3))*dVinjectL(s)  ! 0 < rannum(:,3)=cos(theta)^2 < 1
	                phnheat(5,bg:ed)=DSQRT(1d0-rannum(:,3))*DCOS(M_PI_2*rannum(:,4))*dVinjectL(s)
	                phnheat(6,bg:ed)=DSQRT(1d0-rannum(:,3))*DSIN(M_PI_2*rannum(:,4))*dVinjectL(s)	            
                    phnheat(7,bg:ed)=dEinjectL(s)
	                phnheat(8,bg:ed)=dVinjectL(s)
	                phnheat(9,bg:ed)=s
                    phnheat(10,bg:ed)=BCelementL(j)
	                dtheat(bg:ed)=tStep*rannum(:,5)
	                DEALLOCATE( rannum )
	            ENDIF
            ENDIF
!------------------------------------------------------------------------------
        ELSE IF (j.gt.NBCL) THEN
            m=j-NBCL
            IF (NemitR(m).gt.0) THEN
                bg=SUM(NemitL)+SUM(NemitR(1:m-1))+1
                ed=SUM(NemitL)+SUM(NemitR(1:m))        
                s=GrainMt(Element(1,BCelementR(m)))
	            IF (truePoolR(m).eq.1) THEN
	                IF (WAY_DIR.eq.1) THEN
	                    k=m
                        phnheat(2,bg:ed)=dPoolR(6,mR(1:NemitR(m),m),k)
                        phnheat(3,bg:ed)=dPoolR(7,mR(1:NemitR(m),m),k)	                    
	                ELSE IF (WAY_DIR.eq.2) THEN
	                    k=s
	                    ALLOCATE( rannum(NemitR(m),2) )
                        CALL RAN_NUM(RNseed0,rannum)
                        rannum(:,2)=rannum(:,2)*(1d0-rannum(:,1))
	                    phnheat(2,bg:ed)=xynodes(2,Element(6,BCelementR(m)))+rannum(:,1)*vecCoords(2,1,BCelementR(m))+rannum(:,2)*vecCoords(2,2,BCelementR(m))
                        phnheat(3,bg:ed)=xynodes(3,Element(6,BCelementR(m)))+rannum(:,1)*vecCoords(3,1,BCelementR(m))+rannum(:,2)*vecCoords(3,2,BCelementR(m))
                        DEALLOCATE( rannum )
	                ENDIF
                    phnheat(1,bg:ed)=dLdomain(1)
                    phnheat(4,bg:ed)=dPoolR(2,mR(1:NemitR(m),m),k)*dVinjectR(s)
                    phnheat(5,bg:ed)=dPoolR(3,mR(1:NemitR(m),m),k)*dVinjectR(s)
                    phnheat(6,bg:ed)=dPoolR(4,mR(1:NemitR(m),m),k)*dVinjectR(s)
	                phnheat(7,bg:ed)=dEinjectR(dPoolR(5,mR(1:NemitR(m),m),k))
                    phnheat(8,bg:ed)=dVinjectR(s)
	                phnheat(9,bg:ed)=dPoolR(5,mR(1:NemitR(m),m),k)
                    phnheat(10,bg:ed)=BCelementR(m)
                    dtheat(bg:ed)=dPoolR(1,mR(1:NemitR(m),m),k)
                ELSE	   ! Second way to determine the properties of emitted phonons
	                ALLOCATE( rannum(NemitR(m),5) )
                    CALL RAN_NUM(RNseed0,rannum)
                    rannum(:,2)=rannum(:,2)*(1d0-rannum(:,1))
	                phnheat(1,bg:ed)=dLdomain(1)
	                phnheat(2,bg:ed)=xynodes(2,Element(6,BCelementR(m)))+rannum(:,1)*vecCoords(2,1,BCelementR(m))+rannum(:,2)*vecCoords(2,2,BCelementR(m))
                    phnheat(3,bg:ed)=xynodes(3,Element(6,BCelementR(m)))+rannum(:,1)*vecCoords(3,1,BCelementR(m))+rannum(:,2)*vecCoords(3,2,BCelementR(m))
                    phnheat(4,bg:ed)=-DSQRT(rannum(:,3))*dVinjectR(s)
	                phnheat(5,bg:ed)=DSQRT(1d0-rannum(:,3))*DCOS(M_PI_2*rannum(:,4))*dVinjectR(s)
	                phnheat(6,bg:ed)=DSQRT(1d0-rannum(:,3))*DSIN(M_PI_2*rannum(:,4))*dVinjectR(s)
                    phnheat(7,bg:ed)=dEinjectR(s)
	                phnheat(8,bg:ed)=dVinjectR(s)
	                phnheat(9,bg:ed)=s
                    phnheat(10,bg:ed)=BCelementR(m)
	                dtheat(bg:ed)=tStep*rannum(:,5)
	                DEALLOCATE( rannum )
	            ENDIF
            ENDIF
        ENDIF
    ENDDO
    
    !$OMP PARALLEL DEFAULT(SHARED),PRIVATE(j,CPUID)
    CPUID=OMP_GET_THREAD_NUM()
    !$OMP DO SCHEDULE(DYNAMIC)
    DO j=1,heatNph
        CALL proc_advection(phnheat(:,j),dtheat(j),2,RNseedMP(CPUID))
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
    
    DEALLOCATE( dtheat )
    DEALLOCATE( mL,mR )
 
END SUBROUTINE proc_heatcontrol

!============================================================================
END MODULE mod_heatcontrol