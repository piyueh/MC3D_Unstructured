!============================================================================
PROGRAM initial2D
USE mod_VARIABLES
IMPLICIT NONE
INTEGER*4 :: NperCell,Ndomain
CHARACTER*72 :: casename,gridfile
CHARACTER::YN

    PRINT*,'INPUT THE CASE NAME(NO .''TXT'', NO ''_GRIDFILE''): '
    READ(*,"(A72)") casename

    gridfile=casename(1:LEN_TRIM(casename))//'_gridfile.txt'
    PRINT*, 'CALCULATE INFORMATIONS OF ELEMENT GROUPS ? (Y/N)'
    READ(*,*) YN
    
    OPEN(LR1,file=gridfile)
    CALL readtable
    CALL gridINFO
    CLOSE(LR1)

    ALLOCATE( dEdiff(Nelement),dVunit(Nelement),dEunit(Nelement),dTemp(Nelement) )
    ALLOCATE( Nnumcell(Nelement),Nbgcell(Nelement) )
    CALL random_seed()

    CALL domain(NperCell,Ndomain)
    CALL properties(NperCell,Ndomain)
    IF ((YN.eq.'Y').OR.(YN.eq.'y')) THEN
        CALL elem_group(Npercell)
    ENDIF
    CALL output(casename,YN)

    DEALLOCATE( Element,Volume,xyNodes,BCelementL,BCelementR )
    DEALLOCATE( dEdiff,dTemp,dVunit,dEunit,Ge_table,Si_Table,phn,dPoolL,dPoolR )
    
END PROGRAM initial2D
!============================================================================
SUBROUTINE gridINFO
USE mod_VARIABLES
IMPLICIT NONE
INTEGER*4:: k
    
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
    
    WRITE(*,*) 'dLdomain:',dLdomain
    
    WRITE(*,*) MAXVAL(volume),MINVAL(volume)
    
END SUBROUTINE gridINFO
!============================================================================
SUBROUTINE output(casename,YN)
USE mod_VARIABLES
IMPLICIT NONE
INTEGER*4::i
CHARACTER*72::casename
CHARACTER::YN

    OPEN(140,file=casename(1:LEN_TRIM(casename))//'_initial.txt')
    
    
    WRITE(140,*) bundle,dt,time0,iter0
    WRITE(140,*) DPP,DPPB
    WRITE(140,*) Nprop,Nph0,Nph
    WRITE(140,*) phn
    WRITE(140,*) dEdiff
    WRITE(140,*) Npool,Nmakeup
    WRITE(140,*) dPoolL,dPoolR
    
    CLOSE(140)
    
    IF ((YN.eq.'Y').OR.(YN.eq.'y')) THEN
        OPEN(141,file=casename(1:LEN_TRIM(casename))//'_group.txt')
        WRITE(141,*) Vgroup
        WRITE(141,*) ElemGroup(:)%N
        DO i=1,Nelement
            IF (ElemGroup(i)%N.ne.1) THEN
                WRITE(141,*) ElemGroup(i)%array
            ENDIF
        ENDDO
        CLOSE(141)
    ENDIF
    
END SUBROUTINE output
!============================================================================
SUBROUTINE domain(NperCell,Ndomain)
USE mod_VARIABLES
IMPLICIT NONE
INTEGER*4 :: NperCell,Ndomain
INTEGER*4::k,i

    WRITE(*,*) '輸入介面粗糙度'
    READ(*,*) DPP
    WRITE(*,*) '輸入邊界粗糙度'
    READ(*,*) DPPB

    Nprop=10 !FOR 3-D
    
    WRITE(*,*) '輸入NperCell:'
    READ(*,*) NperCell
    WRITE(*,*) '輸入Ndomain:'
    READ(*,*) Ndomain

    WRITE(*,*) '選擇聲子入射方法：'
    WRITE(*,*) '1:Cell-to-Cell週期入射法'
    WRITE(*,*) '2:Material-to-Material週期入射法'
    READ(*,*) WAY_DIR
    
    IF (WAY_DIR.eq.1) THEN
        Npool=7
    ELSE IF (WAY_DIR.eq.2) THEN
        Npool=5
    ENDIF
    
    time0=0
    iter0=0
    
    dTemp=330d0

END SUBROUTINE domain
!============================================================================
SUBROUTINE properties(NperCell,Ndomain)
USE mod_VARIABLES
IMPLICIT NONE
REAL*8 :: checkE,checkT,dt_min,sumT,dt1,dt2,dt3,dt4,dt5,dt6,dt7
INTEGER*4 :: NperCell,Ndomain,s,k,true(2),i,tempI1
INTEGER*8::m
REAL*8,ALLOCATABLE::rannum(:,:)
CHARACTER*72::YN

    true=0
    dt_min=100000
    tempI1=0
    m=0
    
    DO k=1,Nelement
        s=GrainMt(Element(1,k))
        true(s)=1
        CALL energy(s,dTemp(k),dEdiff(k))            ! U
        CALL ETable(s,4,dEdiff(k),dVunit(k))         ! V
        CALL ETable(s,2,dEdiff(k),dEunit(k))         ! N
        CALL ETable(s,5,dEdiff(k),dTemp(k))          ! MFP
   
        m=m+INT(dEunit(k)*Volume(k)+0.5d0)
    ENDDO
    
    dt1=MINVAL((Volume(:)**(1d0/3d0))*0.5/dVunit)
    dt2=MINVAL(Volume(:)**(1d0/3d0))*0.5/MAXVAL(dVunit)
    CALL energy(1,330d0,dt4)
    CALL Etable(1,4,dt4,dt_min)
    CALL Etable(1,5,dt4,sumT)
    dt3=MINVAL(Volume(:)**(1d0/3d0))*0.5/dt_min
    
    dt4=100000d0
    dt5=100000d0
    DO i=1,Nelement
        dt6=MAXVAL(xyNodes(1,Element(6:9,i)))-MINVAL(xyNodes(1,Element(6:9,i)))
        dt7=dt6*0.5/dVunit(i)
        IF (dt7.lt.dt4) dt4=dt7
        IF (dt6.lt.dt5) dt5=dt6
    ENDDO
    dt5=dt5*0.5/dt_min
    
    PRINT*,'dt1=',dt1
    PRINT*,'dt2=',dt2
    PRINT*,'dt3=',dt3
    PRINT*,'dt4=',dt4
    PRINT*,'dt5=',dt5

    dt=MIN(dt1,dt2,dt3,dt4,dt5)
    PRINT*,'dt = ',dt
    WRITE(*,*)
    WRITE(*,*) 'CHANGE dt? KEY IN NEW dt, OR KEY IN THE ORIGINAL dt.'
    READ(*,*) dt
    PRINT*,'dt = ',dt
!-----------------------------------------------------
    PRINT*,'min Volume = ',MINVAL(Volume)
    PRINT*,'MFP = ',MINVAL(dTemp)
    PRINT*,'tau = ',MINVAL(dTemp)/MAXVAL(dVunit)
    PRINT*,'dz/vel = ',MINVAL(Volume**(1d0/3d0))/MAXVAL(dVunit)
    PRINT*,'vel = ',MAXVAL(dVunit)

    WRITE(*,*) m
    bundle=DBLE(m)/DBLE(Ndomain)
    
    WRITE(*,*) 'BUNDLE IS :',bundle
    WRITE(*,*) 'CHANGE BUNDLE?(Y/N)'
    READ(*,*) YN
    IF ((YN(1:1).eq.'y').OR.(YN(1:1).eq.'Y')) THEN
        WRITE(*,*) 'INPUT NEW BUNDLE:'
        READ(*,*) bundle(1)
        bundle(2)=bundle(1)
    ENDIF
    PRINT*,'min # of phonon bundles of each cell = ', INT(MINVAL(dEunit*Volume)/bundle(s)+0.5d0)
    PRINT*,'MAX # of phonon bundles of each cell = ', INT(MAXVAL(dEunit*Volume)/bundle(s)+0.5d0)

    Nmakeup =INT(MAXVAL(dEunit)*MAXVAL(Volume)/bundle(1)+0.5d0)*3
    WRITE(*,*) 'Nmakeup=',Nmakeup
    
    
    IF (WAY_DIR.eq.1) THEN
        ALLOCATE( dPoolL(Npool,Nmakeup,NBCL),dPoolR(Npool,Nmakeup,NBCR) )  
    ELSE IF (WAY_DIR.eq.2) THEN
        ALLOCATE( dPoolL(Npool,Nmakeup,2),dPoolR(Npool,Nmakeup,2) ) ! for two materials and two boundaries
    ENDIF

    dPoolL=0
    dPoolR=0

    PRINT*,'Nmakeup ',Nmakeup
    
    PAUSE '按任意鍵繼續'

    dTemp=dEdiff ! U now
    dEdiff=0
    Nph=0

    DO k=1,Nelement
        s=GrainMt(Element(1,k))
        m=INT(dEunit(k)*Volume(k)/bundle(s)+0.5d0)
        !IF (m.lt.1) PAUSE
        dEdiff(k)=dTemp(k)*Volume(k)-DBLE(m)*bundle(s)*dTemp(k)/dEunit(k)
        Nbgcell(k)=Nph
        Nnumcell(k)=m
        Nph=Nph+m
    ENDDO
    
    IF (sum(Nnumcell).ne.Nph) PAUSE 'sum(Nnumcell).ne.Nph'
    IF ((Nbgcell(Nelement)+Nnumcell(Nelement)).ne.Nph) PAUSE '(Nbgcell(Nelement)+Nnumcell(Nelement)).ne.Nph'
    PRINT*,'total # of phonons = ',Nph
    PRINT*,'Nprop = ',Nprop

!------------------------------------
    Nph0=Nph
    ALLOCATE( phn(Nprop,Nph0) )
    phn=0

    WRITE(*,*) 'start to give phonons properties'
    DO k=1,Nelement
        ALLOCATE( rannum(Nnumcell(k),5) )
	    CALL random_number(rannum)
	    rannum(:,3)=(1d0-rannum(:,1))*(1d0-rannum(:,2))*rannum(:,3)
	    rannum(:,2)=(1d0-rannum(:,1))*rannum(:,2)
	    phn(1,Nbgcell(k)+1:Nbgcell(k)+Nnumcell(k))=xyNodes(1,Element(6,k))+rannum(:,1)*vecCoords(1,1,k)+rannum(:,2)*vecCoords(1,2,k)+rannum(:,3)*vecCoords(1,3,k)
        phn(2,Nbgcell(k)+1:Nbgcell(k)+Nnumcell(k))=xyNodes(2,Element(6,k))+rannum(:,1)*vecCoords(2,1,k)+rannum(:,2)*vecCoords(2,2,k)+rannum(:,3)*vecCoords(2,3,k)
	    phn(3,Nbgcell(k)+1:Nbgcell(k)+Nnumcell(k))=xyNodes(3,Element(6,k))+rannum(:,1)*vecCoords(3,1,k)+rannum(:,2)*vecCoords(3,2,k)+rannum(:,3)*vecCoords(3,3,k)
	    phn(4,Nbgcell(k)+1:Nbgcell(k)+Nnumcell(k))=2d0*rannum(:,4)-1d0
	    phn(5,Nbgcell(k)+1:Nbgcell(k)+Nnumcell(k))=DSQRT(1d0-phn(4,Nbgcell(k)+1:Nbgcell(k)+Nnumcell(k))**2)*DCOS(M_PI_2*rannum(:,5))*dVunit(k)
	    phn(6,Nbgcell(k)+1:Nbgcell(k)+Nnumcell(k))=DSQRT(1d0-phn(4,Nbgcell(k)+1:Nbgcell(k)+Nnumcell(k))**2)*DSIN(M_PI_2*rannum(:,5))*dVunit(k)
	    phn(4,Nbgcell(k)+1:Nbgcell(k)+Nnumcell(k))=phn(4,Nbgcell(k)+1:Nbgcell(k)+Nnumcell(k))*dVunit(k)
	    phn(7,Nbgcell(k)+1:Nbgcell(k)+Nnumcell(k))=dTemp(k)/dEunit(k)*bundle(GrainMt(Element(1,k)))
	    phn(8,Nbgcell(k)+1:Nbgcell(k)+Nnumcell(k))=dVunit(k)
	    phn(9,Nbgcell(k)+1:Nbgcell(k)+Nnumcell(k))=GrainMt(Element(1,k))
	    phn(10,Nbgcell(k)+1:Nbgcell(k)+Nnumcell(k))=k 
	    DEALLOCATE( rannum )
	    
	    DO m=Nbgcell(k)+1,Nbgcell(k)+Nnumcell(k)
	        CALL check_element(k,m)
	    ENDDO
    ENDDO

    OPEN(UNIT=300,FILE='check initial T.txt')
    sumT=0
    DO k=1,Nelement
        checkE=(SUM(phn(7,Nbgcell(k)+1:Nbgcell(k)+Nnumcell(k)))+dEdiff(k))/Volume(k)
        CALL Etable(GrainMt(Element(1,k)),1,checkE,checkT)
        WRITE(300,*) center(1,k),checkT
        sumT=sumT+checkT*Volume(k)
    ENDDO
    CLOSE(300)

    WRITE(*,*) 'The initial average temperary is:',sumT/SUM(Volume)
    PAUSE 'Press to continue'

END SUBROUTINE properties
!============================================================================
SUBROUTINE elem_group(Npercell)
USE mod_VARIABLES
IMPLICIT NONE
INTEGER*4::Npercell,i,j,k,s,temp(200)
REAL*8::Vlim,R,E,Rnode(4),R2

    ALLOCATE( ElemGroup(Nelement),Vgroup(Nelement) )
    Vgroup=Volume
    CALL energy(2,330d0,E)  !Si 的單位體積聲子數比Ge少
    CALL Etable(2,2,E,Vlim) !Vlim = 單位體積聲子數
    Vlim=Npercell*bundle(1)/Vlim !now Vlim = group volume
    R=(Vlim*3*0.25/M_PI)**(1d0/3d0)
    WRITE(*,*) "THE GROUP VOLUME AND R ARE:",Vlim,R
    
    DO i=1,Nelement
        IF ( (Vgroup(i)-Vlim)/Vlim.lt.-0.25 ) THEN
            temp(1)=i
            temp(2:)=Nelement+10
            Vgroup(i)=Volume(i)
            R2=R
            s=1
            DO WHILE(Vgroup(i).lt.Vlim)
                DO j=1,Nelement
                    IF ((Element(1,j).eq.Element(1,i)).AND.(ALL(temp.ne.j))) THEN
                        FORALL(k=1:4) Rnode(k)=DSQRT(SUM((xyNodes(:,Element(k+5,j))-center(:,i))**2))
                        IF (ANY(Rnode.le.R2)) THEN
                            s=s+1
                            temp(s)=j
                            Vgroup(i)=Vgroup(i)+Volume(j)
                        ENDIF
                    ENDIF
                    IF (Vgroup(i).gt.Vlim) EXIT
                ENDDO
                R2=R2+1.06*R
            ENDDO
            ElemGroup(i)%N=s
            ALLOCATE( ElemGroup(i)%array(s) )
            ElemGroup(i)%array=temp(1:s)        
        ELSE
            ElemGroup(i)%N=1
            ALLOCATE( ElemGroup(i)%array(1) )
            ElemGroup(i)%array(1)=i
            Vgroup(i)=Volume(i)
        ENDIF
        
        IF ( ALL(ElemGroup(i)%array.ne.i) ) PAUSE 'WRONG 1 IN ELEM_GROUP'
        IF ( (Vgroup(i)-Vlim)/Vlim.lt.-0.25 ) THEN
            WRITE(*,*) Vgroup(i),Volume(i)
            PAUSE 'WRONG 2 IN ELEM_GROUP'
        ENDIF
        
        WRITE(*,*) i,ElemGroup(i)%N
        
    ENDDO
    
    WRITE(*,*) 'MAX(ElemGroup(:)%N)=',MAXVAL(ElemGroup(:)%N)

END SUBROUTINE elem_group
!============================================================================
SUBROUTINE check_element(k,m)
USE mod_VARIABLES
IMPLICIT NONE
REAL*8::A(3),r(3)
INTEGER*4:: k,m

    r(:)=phn(1:3,m)-xyNodes(:,Element(6,k))
    A(:)=MATMUL(invers_vec(:,:,k),r)
    IF ((A(1).lt.0).OR.(A(2).lt.0).OR.(A(3).lt.0).OR.(SUM(A).gt.1)) THEN
        WRITE(*,*) k,m
        WRITE(*,*) phn(:,m)
        WRITE(*,*) Element(:,k)
        WRITE(*,*) xyNodes(:,Element(6,k))
        WRITE(*,*) 'vecCoords:'
        WRITE(*,*) vecCoords(1,:,k)
        WRITE(*,*) vecCoords(2,:,k)
        WRITE(*,*) vecCoords(3,:,k)
        WRITE(*,*) 'invers_vec:'
        WRITE(*,*) invers_vec(1,:,k)
        WRITE(*,*) invers_vec(2,:,k)
        WRITE(*,*) invers_vec(3,:,k)
        WRITE(*,*) 'r:'
        WRITE(*,*) r
        WRITE(*,*) 'A:'
        WRITE(*,*) A
        PAUSE 'Phonon is not in the element.'
    ENDIF

END SUBROUTINE check_element
!============================================================================
SUBROUTINE energy(nc,T0,Eout)
USE mod_VARIABLES
IMPLICIT NONE
INTEGER*4 :: nc
REAL*8 :: T0,Eout,E1,E2,T1,T2,Tout

    IF (nc.eq.1) THEN
        E1 = Ge_table(3,1)
        E2 = Ge_table(3,N_Ge2)
        T1 = Ge_table(1,1)-T0
        T2 = Ge_table(1,N_Ge2)-T0
    ELSE IF (nc.eq.2) THEN
        E1 = Si_table(3,1)
        E2 = Si_table(3,N_Si2)
        T1 = Si_table(1,1)-T0
        T2 = Si_table(1,N_Si2)-T0
    ENDIF
    Eout=(E1+E2)/2d0
    CALL ETable(nc,1,Eout,Tout)
    Tout=Tout-T0
    DO WHILE(ABS(Tout).gt.zero_tol)
        IF (Tout*T1.lt.0d0) THEN
            E2=Eout
	        T2=Tout
        ELSE
	        E1=Eout
	        T1=Tout
        ENDIF
        Eout=(E1+E2)/2d0
        CALL ETable(nc,1,Eout,Tout)
        Tout=Tout-T0
    ENDDO

END SUBROUTINE energy
!============================================================================
SUBROUTINE readtable
USE mod_VARIABLES
IMPLICIT NONE

    OPEN(unit=110,file="Ge_real_table.txt")
    OPEN(unit=120,file="Si_real_table.txt")

    READ(110,*) N_Ge1,N_Ge2
    READ(120,*) N_Si1,N_Si2

    ALLOCATE( Ge_table(N_Ge1,N_Ge2) )
    ALLOCATE( Si_table(N_Si1,N_Si2) )
    READ(110,*) Ge_table
    READ(120,*) Si_table

    Ge_start=Ge_table(3,1)
    dU_Ge=Ge_table(3,2)-Ge_table(3,1)
    Si_start=Si_table(3,1)
    dU_Si=Si_table(3,2)-Si_table(3,1)

END SUBROUTINE readtable
!============================================================================
!mat=1:Ge, mat=2:Si, idx=1:temperature, 2:number density, 4:velocity, 5:MFP, 6:specific heat
SUBROUTINE Etable(mat,idx,E,out) 
USE mod_VARIABLES
IMPLICIT NONE
INTEGER*4::mat,idx,iU
REAL*8::E,out,out_a,out_b

    IF (mat.eq.1) THEN
	    iU = INT( (E-Ge_start)/dU_Ge )+1
	    out_a = Ge_table(idx,iU)
	    out_b = Ge_table(idx,iU+1)
	    out = out_a + (out_b-out_a)*(E-Ge_table(3,iU))/dU_Ge
    ELSE IF (mat.eq.2) THEN
    	iU = int( (E-Si_start)/dU_Si )+1
	    out_a = Si_table(idx,iU)
	    out_b = Si_table(idx,iU+1)
	    out = out_a + (out_b-out_a)*(E-Si_table(3,iU))/dU_Si
    END IF

END SUBROUTINE Etable
!============================================================================