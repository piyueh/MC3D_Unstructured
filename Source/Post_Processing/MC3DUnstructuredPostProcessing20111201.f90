
MODULE variables
IMPLICIT NONE
!----MATERIAL TABLE INFORMATION
INTEGER*4::N_Ge1,N_Ge2,N_Si1,N_Si2
REAL*8::Ge_start,Si_start,dU_Ge,dU_Si,rho(2)=(/5.323d3,2.329d3/) !rho_Ge=5.323d3,rho_Si=2.329d3 ! kg/m^3
REAL*8,ALLOCATABLE::Ge_table(:,:),Si_table(:,:)
REAL*8,PARAMETER::zero_tol=0.00000000001
INTEGER*4,PARAMETER::LRSi=110,LRGe=120,LR1=130
!-----FOR MESH INFORMATION-----!
INTEGER*4:: Nnodes,Nelement,Propelem,Ngrains
INTEGER*4:: NBCL,NBCR,NBCyP,NBCyN,NBCzP,NBCzN
INTEGER*4,ALLOCATABLE:: BCelementL(:),BCelementR(:),BCelementyP(:),BCelementyN(:),BCelementzP(:),BCelementzN(:),Element(:,:),GrainMt(:)
REAL*8,ALLOCATABLE::xyNodes(:,:),Volume(:),center(:,:),normVec(:,:,:),vecCoords(:,:,:),invers_vec(:,:,:), areaBCL(:),areaBCR(:)
!----SOLVER OUTPUT VARIABLES---!
REAL*8::heatflowL,heatflowR,heatflowC,dt
REAL*8,ALLOCATABLE::qbdyL(:),qbdyR(:),cellTemp(:)
!-----FOR TEMPORARILY USE------!
INTEGER*4::i,j,k,m,s,temp_INT
REAL*8::DomainL(3)
LOGICAL::true(5)
CHARACTER*72::casename,gridfile
!----INTERPOLATING VARIABLES---!
INTEGER*4::Nx,Ny,Nz,Nelement2
REAL*8::dx,dy,dz,CuboidEnergy
REAL*8,ALLOCATABLE::xxyyzz(:,:,:,:),CuboidTemp(:,:,:) !,weight(:,:,:,:)
!INTEGER*4,ALLOCATABLE::NonZeroV(:,:,:),Vol_Comp(:)
REAL*8::temp_REAL
REAL*8,ALLOCATABLE::Exam_Weight(:),Exam_Volume(:)

TYPE CompositeCuboid
    INTEGER*4::Ncomp
    INTEGER*4,ALLOCATABLE::CompLabel(:)
    REAL*8,ALLOCATABLE::CompVolume(:)
END TYPE CompositeCuboid
TYPE(CompositeCuboid),ALLOCATABLE::CuboidComp(:,:,:)

END MODULE variables



!========================================================================================================
PROGRAM main
USE variables
IMPLICIT NONE
CHARACTER*72::YN
CHARACTER*72::outputname


    !************************READ GRID INFORMATION*****************************
    
    WRITE(*,*) 'ENTER CASE NAME (no ''.txt'') :'
    READ(*,"(A72)") casename
    WRITE(*,*) 'ENTER THE TITLE WHILE PLOTING IN MATLAB :'
    READ(*,"(A72)") outputname
    WRITE(*,*) outputname
    
    gridfile=casename(1:LEN_TRIM(casename))//'_gridfile.txt'
    OPEN(UNIT=LR1,FILE=gridfile)
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
    CLOSE(LR1)
    
    IF ((DABS(MINVAL(xyNodes(1,:))).gt.0.0000000001).OR.(DABS(MINVAL(xyNodes(2,:))).gt.0.0000000001).OR.(DABS(MINVAL(xyNodes(3,:))).gt.0.0000000001)) THEN
        PAUSE '幾何圖形沒定位在原點'
    ENDIF
    
    DomainL(1)=MAXVAL(xyNodes(1,:))
    DomainL(2)=MAXVAL(xyNodes(2,:))
    DomainL(3)=MAXVAL(xyNodes(3,:))

    !****************READ SOLVER OUPUT AVERAGE STEADY PROPERTIES*************** 
    ALLOCATE( qbdyL(NBCL),qbdyR(NBCR),cellTemp(Nelement) )
    CALL elem_ave_prop(Nelement,NBCL,NBCR,qbdyL,qbdyR,cellTemp,heatflowL,heatflowR,heatflowC,dt,Volume)


    CALL interpolate(outputname)
    
    PAUSE 'PROGRAM FINISH'
   
END PROGRAM main
!==================================================================================================
!==================================================================================================
SUBROUTINE elem_ave_prop(Nelement,NBCL,NBCR,qbdyL,qbdyR,cellTemp,heatflowL,heatflowR,heatflowC,dt,Volume)
IMPLICIT NONE
!---------------------------!
INTEGER*4::Nelement,NBCL,NBCR
REAL*8::qbdyL(NBCL),qbdyR(NBCR),cellTemp(Nelement)
REAL*8::heatflowL,heatflowR,heatflowC,dt,Volume(Nelement)
!---------------------------!
REAL*8::err,time,tempR1
INTEGER*4::i,j,k,s,N,tempI1,Nnosteady
INTEGER*4::bg0,ed0,bg,ed,step,Nfiles
LOGICAL::true
CHARACTER*72::readfile
REAL*8::tempheatflowC,tempheatflowL,tempheatflowR,tempR
REAL*8,ALLOCATABLE::tempCellTemp(:),tempqbdyL(:),tempqbdyR(:)
REAL*8,ALLOCATABLE::aveheatflowC(:),aveheatflowL(:),aveheatflowR(:)
REAL*8,ALLOCATABLE::aveCellTemp(:,:),aveqbdyL(:,:),aveqbdyR(:,:)

    
    
    WRITE(*,*) '輸入起始ITER'
    READ(*,*) bg0
    WRITE(*,*) '輸入結束ITER'
    READ(*,*) ed0
    WRITE(*,*) '輸入STEP'
    READ(*,*) step
    WRITE(*,*) '輸入N'
    READ(*,*) N
    
    Nfiles=(ed0-bg0)/step+1
    
    WRITE(*,*) 'THERE ARE',Nfiles,'FILES.'
    
    tempI1=Nfiles/N
    
    ALLOCATE( AveHeatFlowC(N),AveHeatFlowL(N),AveHeatFlowR(N) )
    ALLOCATE( AveCellTemp(Nelement,N),AveQbdyL(NBCL,N),AveQbdyR(NBCR,N) )
    ALLOCATE( tempCellTemp(Nelement),tempQbdyL(NBCL),tempQbdyR(NBCR) )
    
    AveHeatFlowC=0
    AveHeatFlowL=0
    AveHeatFlowR=0
    AveCellTemp=0
    AveQbdyL=0
    AveQbdyR=0
    
    DO i=1,N
    
        ed=ed0-(i-1)*tempI1*step
        bg=ed0-(i*tempI1-1)*step
        
        IF (i.ne.1) THEN
            AveHeatFlowC(i)=AveHeatFlowC(i-1)
            AveHeatFlowL(i)=AveHeatFlowL(i-1)
            AveHeatFlowR(i)=AveHeatFlowR(i-1)
            AveQbdyL(:,i)=AveQbdyL(:,i-1)
            AveQbdyR(:,i)=AveQbdyR(:,i-1)
            AveCellTemp(:,i)=AveCellTemp(:,i-1)
        ENDIF
        
        DO j=ed,bg,-step
        
            WRITE(readfile,"((I7.7),'.txt')") j
            OPEN(UNIT=110,FILE=readfile)
            READ(110,*) time,dt
            READ(110,*) tempHeatFlowC
            READ(110,*) tempHeatFlowL,tempHeatFlowR
            READ(110,*) tempQbdyL,tempQbdyR
            READ(110,*) tempCellTemp
            CLOSE(110)
        
            AveHeatFlowC(i)=AveHeatFlowC(i)+tempHeatFlowC
            AveHeatFlowL(i)=AveHeatFlowL(i)+tempHeatFlowL
            AveHeatFlowR(i)=AveHeatFlowR(i)+tempHeatFlowR
            AveQbdyL(:,i)=AveQbdyL(:,i)+tempQbdyL
            AveQbdyR(:,i)=AveQbdyR(:,i)+tempQbdyR
            AveCellTemp(:,i)=AveCellTemp(:,i)+tempCellTemp
        
        ENDDO
        
    ENDDO
    
    tempHeatFlowC=10000
    tempHeatFlowL=10000
    tempHeatFlowR=10000
    tempQbdyL=10000
    tempQbdyR=10000
    tempCellTemp=10000
    
    AveHeatFlowC(1)=AveHeatFlowC(1)/DBLE(tempI1)
    AveHeatFlowL(1)=AveHeatFlowL(1)/DBLE(tempI1)
    AveHeatFlowR(1)=AveHeatFlowR(1)/DBLE(tempI1)
    AveQbdyL(:,1)=AveQbdyL(:,1)/DBLE(tempI1)
    AveQbdyR(:,1)=AveQbdyR(:,1)/DBLE(tempI1)
    AveCellTemp(:,1)=AveCellTemp(:,1)/DBLE(tempI1)
    
    DO i=2,N
    
        AveHeatFlowC(i)=AveHeatFlowC(i)/DBLE(i*tempI1)
        tempR=DABS((AveHeatFlowC(i)-AveHeatFlowC(i-1))/AveHeatFlowC(i-1))
        IF (tempR.lt.tempHeatFlowC) THEN
            tempHeatFlowC=tempR
            HeatFlowC=AveHeatFlowC(i)
        ENDIF
        
        AveHeatFlowL(i)=AveHeatFlowL(i)/DBLE(i*tempI1)
        tempR=DABS((AveHeatFlowL(i)-AveHeatFlowL(i-1))/AveHeatFlowL(i-1))
        IF (tempR.lt.tempHeatFlowL) THEN
            tempHeatFlowL=tempR
            HeatFlowL=AveHeatFlowL(i)
        ENDIF
        
        AveHeatFlowR(i)=AveHeatFlowR(i)/DBLE(i*tempI1)
        tempR=DABS((AveHeatFlowR(i)-AveHeatFlowR(i-1))/AveHeatFlowR(i-1))
        IF (tempR.lt.tempHeatFlowR) THEN
            tempHeatFlowR=tempR
            HeatFlowR=AveHeatFlowR(i)
        ENDIF
        
        AveQbdyL(:,i)=AveQbdyL(:,i)/DBLE(i*tempI1)
        DO j=1,NBCL
            tempR=DABS((AveQbdyL(j,i)-AveQbdyL(j,i-1))/AveQbdyL(j,i-1))
            IF (tempR.lt.tempQbdyL(j)) THEN
                tempQbdyL(j)=tempR
                qbdyL(j)=AveQbdyL(j,i)
            ENDIF
        ENDDO
        
        AveQbdyR(:,i)=AveQbdyR(:,i)/DBLE(i*tempI1)
        DO j=1,NBCR
            tempR=DABS((AveQbdyR(j,i)-AveQbdyR(j,i-1))/AveQbdyR(j,i-1))
            IF (tempR.lt.tempQbdyR(j)) THEN
                tempQbdyR(j)=tempR
                qbdyR(j)=AveQbdyR(j,i)
            ENDIF
        ENDDO
        
        AveCellTemp(:,i)=AveCellTemp(:,i)/DBLE(i*tempI1)
        DO j=1,Nelement
            tempR=DABS((AveCellTemp(j,i)-AveCellTemp(j,i-1))/AveCellTemp(j,i-1))
            IF (tempR.lt.tempCellTemp(j)) THEN
                tempCellTemp(j)=tempR
                CellTemp(j)=AveCellTemp(j,i)
            ENDIF
        ENDDO
            
        
    ENDDO
    WRITE(*,"('AVERAGE HEAT FLUX(CENTER): ',E15.8,', ERROR: ',E15.8)") heatflowC,tempHeatFlowC
    WRITE(*,"('AVERAGE HEAT FLUX(LEFT BOUNDARY): ',E15.8,', ERROR: ',E15.8)") heatflowL,tempHeatFlowL
    WRITE(*,"('AVERAGE HEAT FLUX(RIGHT BOUNDARY): ',E15.8,', ERROR: ',E15.8)") heatflowR,tempHeatFlowR
    WRITE(*,*)
    WRITE(*,*) '| FLUXL-FLUXR |/AVERAGE FLUX:',DABS(heatflowR-heatflowL)*2d0/(heatflowR+heatflowL)
    WRITE(*,*)
    WRITE(*,"('SUM OF qbdyL: ',E15.8,', AVERAGE OF ERROR: ',E15.8)") SUM(qbdyL),SUM(tempQbdyL)/DBLE(NBCL)
    WRITE(*,"('SUM OF qbdyR: ',E15.8,', AVERAGE OF ERROR: ',E15.8)") SUM(qbdyR),SUM(tempQbdyR)/DBLE(NBCR)
    WRITE(*,*)
    WRITE(*,"('THE AVERAGE TEMPERATURE BY VOLUME WEIGHTING IS: ',F15.11)") SUM(CellTemp*Volume)/SUM(Volume)
    WRITE(*,"('THE AVERAGE ERROR OF ELEMENT TEMPERATURE: ',E15.8)") SUM(tempCellTemp)/DBLE(Nelement)
    
        
    OPEN(UNIT=50,FILE='AveHeatFlowC.txt')
    OPEN(UNIT=51,FILE='AveHeatFlowL.txt')
    OPEN(UNIT=52,FILE='AveHeatFlowR.txt')
    WRITE(50,"(<N>(E15.9,/))") AveHeatFlowC
    WRITE(51,"(<N>(E15.9,/))") AveHeatFlowL
    WRITE(52,"(<N>(E15.9,/))") AveHeatFlowR
    CLOSE(50)
    CLOSE(51)
    CLOSE(52)
   
    OPEN(UNIT=50,FILE='Tet_Ave_Temp.txt')
    WRITE(50,"(<Nelement>(F15.11,/))") CellTemp
    CLOSE(50)
    
    
    OPEN(UNIT=50,FILE='Ave_HeatFlux.txt')
    WRITE(50,"('AVERAGE HEAT FLUX(CENTER): ',E15.8,', ERROR: ',E15.8)") HeatFlowC,tempHeatFlowC
    WRITE(50,"('AVERAGE HEAT FLUX(LEFT BOUNDARY): ',E15.8,', ERROR: ',E15.8)") HeatFlowL,tempHeatFlowL
    WRITE(50,"('AVERAGE HEAT FLUX(RIGHT BOUNDARY): ',E15.8,', ERROR: ',E15.8)") HeatFlowR,tempHeatFlowR
    WRITE(50,"('AVERAGE HEAT FLUX((LEFT+RIGHT)/2): ',E15.8)") (HeatFlowR+HeatFlowL)*0.5
    WRITE(50,*) '| FLUXL-FLUXR |/AVERAGE FLUX:',DABS(HeatFlowR-HeatFlowL)*2d0/(HeatFlowR+HeatFlowL)
    WRITE(50,*)
    WRITE(50,"('SUM OF qbdyL: ',E15.8,', AVERAGE OF ERROR: ',E15.8)") SUM(qbdyL),SUM(tempQbdyL)/DBLE(NBCL)
    WRITE(50,"('SUM OF qbdyR: ',E15.8,', AVERAGE OF ERROR: ',E15.8)") SUM(qbdyR),SUM(tempQbdyR)/DBLE(NBCR)
    CLOSE(50)
    
    DEALLOCATE( aveheatflowC,aveheatflowL,aveheatflowR )
    DEALLOCATE( aveCellTemp,aveqbdyL,aveqbdyR )
    
END SUBROUTINE elem_ave_prop
!==================================================================================================
!==================================================================================================
SUBROUTINE interpolate(outputname)
USE variables
IMPLICIT NONE
REAL*8::examDomainL(3),dV
CHARACTER::YN
CHARACTER*72::outputname

    !************************READ INTERPOLATING GRID INFORMATION***************
    
    gridfile=casename(1:LEN_TRIM(casename))//'_VolumeWeight.txt'
    
    OPEN(UNIT=LR1,FILE=gridfile)
    READ(LR1,*) examDomainL
    WRITE(*,*) 'The Domain Domension of Interpolating File is :'
    WRITE(*,*) examDomainL
    READ(LR1,*) Nx,Ny,Nz,Nelement2
    WRITE(*,*) 'The Element Number of Interpolating File is :'
    WRITE(*,*) Nelement2
    
    IF ((Nelement2.ne.Nelement).OR.ANY(DABS(examDomainL-DomainL).gt.0.00000001)) THEN
        WRITE(*,*) 'THIS INTERPOLATING FILE DOES MATCH THE GRID FILE!!'
        WRITE(*,*) 'Nelement of Original Case :',Nelement
        WRITE(*,*) 'Nelement of Interpolating File :',Nelement2
        WRITE(*,*) 'Domain Size of Original Case :'
        WRITE(*,*) DomainL
        WRITE(*,*) 'Domain Size of Interpolating File :'
        WRITE(*,*) examDomainL
        WRITE(*,*) Nelement2.ne.Nelement,ANY(DABS(examDomainL-DomainL).gt.0.00000001)
        WRITE(*,*) 'CONTINUE AND CHANGE DOMAIN TO CURRENT GRID FILE ? (Y/N)'
        READ(*,*) YN
        IF ((YN.eq.'y').OR.(YN.eq.'Y')) THEN
            examDomainL=DomainL
            WRITE(*,*) examDomainL
        ELSE
            PAUSE 'PRESS ANY KEY TO EXIT...'
            STOP
        ENDIF
    ENDIF
    
    dx=DomainL(1)/DBLE(Nx)
    dy=DomainL(2)/DBLE(Ny)
    dz=DomainL(3)/DBLE(Nz)
    dV=dx*dy*dz
    
    READ(LR1,*) Nelement2
    
    ALLOCATE( CuboidComp(Nx,Ny,Nz),xxyyzz(3,Nx,Ny,Nz),Exam_Weight(Nx*Ny*Nz),Exam_Volume(Nelement),CuboidTemp(Nx,Ny,Nz) )

    READ(LR1,*) CuboidComp(:,:,:)%Ncomp
    
    DO k=1,Nz
        DO j=1,Ny
            DO i=1,Nx
                temp_INT=CuboidComp(i,j,k)%Ncomp
                IF (temp_INT.eq.0) THEN
                    WRITE(*,*) i,j,k,temp_INT
                    PAUSE
                ENDIF
                ALLOCATE( CuboidComp(i,j,k)%CompLabel(temp_INT),CuboidComp(i,j,k)%CompVolume(temp_INT) )
                CuboidComp(i,j,k)%CompVolume=0
                CuboidComp(i,j,k)%CompLabel=-999
                CuboidComp(i,j,k)%Ncomp=0
            ENDDO
        ENDDO
    ENDDO
    
    DO m=1,Nelement2
        READ(LR1,*) i,j,k,temp_INT,temp_REAL
        CuboidComp(i,j,k)%Ncomp=CuboidComp(i,j,k)%Ncomp+1
        s=CuboidComp(i,j,k)%Ncomp
        CuboidComp(i,j,k)%CompLabel(s)=temp_INT
        CuboidComp(i,j,k)%CompVolume(s)=temp_REAL*dV
    ENDDO
    
    WRITE(*,*) Nx,Ny,Nz
    
    CLOSE(LR1)
    
    m=0
    Exam_Volume=0
    DO k=1,Nz
        DO j=1,Ny
            DO i=1,Nx
                IF (CuboidComp(i,j,k)%Ncomp.ne.SIZE(CuboidComp(i,j,k)%CompLabel)) THEN
                    WRITE(*,*) i,j,k
                    WRITE(*,*) CuboidComp(i,j,k)%Ncomp,SIZE(CuboidComp(i,j,k)%CompLabel)
                    PAUSE 'WRONG IN INTERPOLATING!!'
                ENDIF
                
                m=m+1
                Exam_Weight(m)=DABS(SUM(CuboidComp(i,j,k)%CompVolume)-dV)/dV
                xxyyzz(1,i,j,k)=(DBLE(i)-0.5)*dx
                xxyyzz(2,i,j,k)=(DBLE(j)-0.5)*dy
                xxyyzz(3,i,j,k)=(DBLE(k)-0.5)*dz
                
                Exam_Volume(CuboidComp(i,j,k)%CompLabel)=Exam_Volume(CuboidComp(i,j,k)%CompLabel)+CuboidComp(i,j,k)%CompVolume
            ENDDO
        ENDDO
    ENDDO
    
    Exam_Volume=DABS(Exam_Volume-Volume)/Volume
    
    WRITE(*,*) 'THE MAX AND MIN RELATIVE ERROR OF WEIGHT :'
    WRITE(*,*) MAXVAL(Exam_Weight),MINVAL(Exam_Weight)
    WRITE(*,*) 'THE MAX AND MIN RELATIVE ERROR OF ELEMENT VOLUME :'
    WRITE(*,*) MAXVAL(Exam_Volume),MINVAL(Exam_Volume)
    
    
    DEALLOCATE( Exam_Weight,Exam_Volume )
    
    
    !*********************INTERPOLATING****************************************
    
    CALL readtable
    
    CuboidTemp=0
    
    OPEN(UNIT=110,FILE='InterpolatingTemperature.txt')
    WRITE(110,*) DomainL(1),DomainL(2),DomainL(3)
    WRITE(110,*) Nx,Ny,Nz
    
    OPEN(UNIT=111,FILE='MatlabInterpTemperature.m')
    WRITE(111,"('% Domain Dimension:',F10.2,1X,F10.2,1X,F10.2)") DomainL(1),DomainL(2),DomainL(3)
    WRITE(111,"('% Element#:',I10,1X,I10,1X,I10)") Nx,Ny,Nz
    DO k=1,Nz
        DO j=1,Ny
            DO i=1,Nx
            
                CuboidEnergy=0
                temp_INT=CuboidComp(i,j,k)%Ncomp
                DO m=1,temp_INT
                    s=CuboidComp(i,j,k)%CompLabel(m)
                    CALL proc_energy(GrainMt(Element(1,s)),cellTemp(s),temp_REAL)
                    CuboidEnergy=CuboidEnergy+CuboidComp(i,j,k)%CompVolume(m)*temp_REAL
                ENDDO    
                
                CALL find_equivalent_temperature(temp_INT,CuboidComp(i,j,k)%CompVolume,CuboidEnergy,GrainMt(Element(1,CuboidComp(i,j,k)%CompLabel)),CuboidTemp(i,j,k))
                temp_REAL=DOT_PRODUCT(CuboidComp(i,j,k)%CompVolume,cellTemp(CuboidComp(i,j,k)%CompLabel))/dV ! 直接使用體積加權算平均溫度
                
                WRITE(110,*) CuboidTemp(i,j,k),temp_REAL
                WRITE(111,"('CuboidT(',I4,',',I4,',',I4,')=',F10.5,';')") j,i,k,CuboidTemp(i,j,k)
            ENDDO
        ENDDO
    ENDDO
    CLOSE(110)
    
    
    WRITE(111,*) 'casename=',"' "//outputname(1:LEN_TRIM(outputname))//"'",';'
    WRITE(111,*) 'DomainL(1)=',DomainL(1),';'
    WRITE(111,*) 'DomainL(2)=',DomainL(2),';'
    WRITE(111,*) 'DomainL(3)=',DomainL(3),';'
    WRITE(111,*) 'Ncuboid(1)=',Nx,';'
    WRITE(111,*) 'Ncuboid(2)=',Ny,';'
    WRITE(111,*) 'Ncuboid(3)=',Nz,';'
    WRITE(111,*) 'dt=',dt,';'
    WRITE(111,*) 'Tx=sum(sum(CuboidT,1),3)./(size(CuboidT,1)*size(CuboidT,3));'
    WRITE(111,"('a=',F10.4,':',F10.4,':',F10.4,';')") 0.5d0*dx,dx,(DBLE(Nx)-0.5d0)*dx
    WRITE(111,"('b=',F10.4,':',F10.4,':',F10.4,';')") 0.5d0*dy,dy,(DBLE(Ny)-0.5d0)*dy
    WRITE(111,"('c=',F10.4,':',F10.4,':',F10.4,';')") 0.5d0*dz,dz,(DBLE(Nz)-0.5d0)*dz
    WRITE(111,*) '[x y z]=meshgrid(a,b,c);'
    WRITE(111,*) 'clear a b c'
    WRITE(111,*) 'qfluxC=',(SUM(qbdyL)+SUM(qbdyR))*0.5d0,';'
    WRITE(111,*) 'qfluxL=',heatflowL,';'
    WRITE(111,*) 'qfluxR=',heatflowR,';'
    WRITE(111,*) 'fluxArea=',sum(areaBCL),';'
    CLOSE(111)
    
END SUBROUTINE interpolate
!==================================================================================================

SUBROUTINE find_equivalent_temperature(Ncomponent,ComponentVolume,Etotal,ComponentMaterial,Tout)
USE variables
IMPLICIT NONE
INTEGER*4::p,Ncomponent,ComponentMaterial(Ncomponent)
REAL*8 :: Tout,Eout,Ea,Eb,Ta,Tb,Etotal,LastT
REAL*8::ComponentVolume(Ncomponent),ComponentTemp(Ncomponent)
REAL*8,EXTERNAL::f
    
    Ta=200d0
    Tb=500d0
    Tout=0.5*(Ta+Tb)
    
    Ea=f(Ncomponent,ComponentVolume,ComponentMaterial,Etotal,Ta)
    Eb=f(Ncomponent,ComponentVolume,ComponentMaterial,Etotal,Tb)
    Eout=f(Ncomponent,ComponentVolume,ComponentMaterial,Etotal,Tout)
    
    IF (Ea*Eb.ge.0) THEN
        PAUSE 'WRONG IN BISECTION METHOD!'
    ELSE
        p=0
        DO WHILE ((DABS(Eout)/Etotal).gt.1d-12)
            p=p+1
            IF (p.gt.100000000) THEN
                WRITE(*,*) Tout,DABS(Eout)
                !PAUSE
                EXIT
            ENDIF
            IF (Eout.gt.0) THEN
                Tb=Tout
            ELSE IF (Eout.lt.0) THEN
                Ta=Tout
            ELSE IF (Eout.eq.0) THEN
                EXIT
            ENDIF
            Tout=0.5*(Ta+Tb)
            Eout=f(Ncomponent,ComponentVolume,ComponentMaterial,Etotal,Tout)
        ENDDO
    ENDIF

END SUBROUTINE find_equivalent_temperature
!==================================================================================================
FUNCTION f(Ncomponent,ComponentVolume,ComponentMaterial,Etotal,T)
IMPLICIT NONE
INTEGER*4::Ncomponent,p,ComponentMaterial(Ncomponent)
REAL*8::f,T,Etotal
REAL*8::ComponentVolume(Ncomponent),U(Ncomponent)
    
    DO p=1,Ncomponent
        CALL proc_energy(ComponentMaterial(p),T,U(p))
    ENDDO
    
    f=DOT_PRODUCT(ComponentVolume,U)-Etotal

END FUNCTION f
!==================================================================================================
SUBROUTINE proc_energy(nc,T0,Eout) ! given T0, find Eout
USE variables
IMPLICIT NONE
INTEGER*4 :: nc
REAL*8 :: T0,Eout,Ea,Eb,Ta,Tb,Tout

IF (nc.eq.1) THEN
   Ea = Ge_table(3,1)
   Eb = Ge_table(3,N_Ge2)
   Ta = Ge_table(1,1)-T0
   Tb = Ge_table(1,N_Ge2)-T0
ELSE IF (nc.eq.2) THEN
   Ea = Si_table(3,1)
   Eb = Si_table(3,N_Si2)
   Ta = Si_table(1,1)-T0
   Tb = Si_table(1,N_Si2)-T0
ENDIF
Eout=(Ea+Eb)/2d0
CALL ETable(nc,1,Eout,Tout)
Tout=Tout-T0
DO WHILE(ABS(Tout).gt.zero_tol)
   IF (Tout*Ta.lt.0d0) THEN
       Eb=Eout
	   Tb=Tout
   ELSE
	   Ea=Eout
	   Ta=Tout
   ENDIF
   Eout=(Ea+Eb)/2d0
   CALL ETable(nc,1,Eout,Tout)
   Tout=Tout-T0
ENDDO

END SUBROUTINE proc_energy
!============================================================================
SUBROUTINE Etable(mat,idx,E,out) !mat=1:Ge, mat=2:Si, idx=1:temperature, 2:number density, 4:velocity, 5:MFP, 6:specific heat
USE variables
IMPLICIT NONE
INTEGER*4::mat,idx,iU
REAL*8::E,out,out_a,out_b

IF (mat.eq.1) THEN
	iU = INT( (E-Ge_start)/dU_Ge )+1
	IF (iU.ge.N_Ge2.or.iU.le.0) PAUSE 'out of table!'
	out_a = Ge_table(idx,iU)
	out_b = Ge_table(idx,iU+1)
	out = out_a + (out_b-out_a)*(E-Ge_table(3,iU))/dU_Ge
ELSE IF (mat.eq.2) THEN
	iU = int( (E-Si_start)/dU_Si )+1
	IF (iU.ge.N_Si2.or.iU.le.0) PAUSE 'out of table'
	out_a = Si_table(idx,iU)
	out_b = Si_table(idx,iU+1)
	out = out_a + (out_b-out_a)*(E-Si_table(3,iU))/dU_Si
END IF

END SUBROUTINE Etable
!============================================================================
SUBROUTINE readtable
USE variables
IMPLICIT NONE

OPEN(unit=LRGe,file="Ge_real_table.txt")
OPEN(unit=LRSi,file="Si_real_table.txt")

READ(LRGe,*) N_Ge1,N_Ge2
READ(LRSi,*) N_Si1,N_Si2

ALLOCATE( Ge_table(N_Ge1,N_Ge2) )
ALLOCATE( Si_table(N_Si1,N_Si2) )
READ(LRGe,*) Ge_table
READ(LRSi,*) Si_table

Ge_start=Ge_table(3,1)
dU_Ge=Ge_table(3,2)-Ge_table(3,1)
Si_start=Si_table(3,1)
dU_Si=Si_table(3,2)-Si_table(3,1)

PRINT*,'# of Ge data = ',N_Ge2
PRINT*,'# of Si data = ',N_Si2
PRINT*,'Ge_start (meV) = ',Ge_start
PRINT*,'dU_Ge    (meV) = ',dU_Ge
PRINT*,'Si_start (meV) = ',Si_start
PRINT*,'Si_Ge    (meV) = ',dU_Si

CLOSE(LRGe)
CLOSE(LRSi)
END SUBROUTINE readtable
!============================================================================