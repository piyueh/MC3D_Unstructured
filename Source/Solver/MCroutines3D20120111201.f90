!****************************************************************************
!   HEADING: MC3D ROUTINES MODULE for unstructured grid
!   AUTHOR: 
!   PURPOSE: This module contains all the routines for transmission and reflection.
!   DATE: 
!****************************************************************************
MODULE mod_ROUTINES
USE mod_VARIABLES
USE omp_lib
IMPLICIT NONE
!****************************************************************************
CONTAINS
!============================================================================
! Inelastic Acoustic Mismatch Model (IAMM)
SUBROUTINE Snells( k,neighbor,phm,dcosth1,dsinth1,dcosth2,dsinth2 )
IMPLICIT NONE
INTEGER*4::k,neighbor
REAL*8:: phm(Nprop)
REAL*8:: dcosth1,dsinth1,dsinth2,dcosth2,t(3)

! Vin = |Vin|*cos(theta1)*n+|Vin|*sin(theta1)*t1
! Vtr = |Vtr|*cos(theta2)*n+|Vin|*sin(theta2)*t2
! Vin is inject velocity vector, and Vtr is transmitted velocity vector.
! n is normal vector. t1 and t2 is the unit vector that perpendicular to n.
! Beacuse Velocities before and after transimisstion lie on the same projected plane,
! i.e. they lie on the same line on the z-t plane.
! => unit vector t1 = unit vector t2 = t

    t(:)=(phm(4:6)/phm(8)-dcosth1*normVec(:,neighbor-1,k))/dsinth1
    phm(4:6)=dcosth2*normVec(:,neighbor-1,k)+dsinth2*t(:)

!this is only the unit direction vector!!

END SUBROUTINE Snells
!============================================================================
SUBROUTINE diffuseB( k,neighbor,phm,nc,rannum )
! nc=1 : transmission
! nc=-1 : reflection
IMPLICIT NONE
INTEGER*4:: nc,k,neighbor,Nrand
REAL*8:: phm(Nprop),dSinthSinphi,dSinthCosphi,dCosth,t1(3),t2(3),rannum(2),temp

    IF (nc.eq.1) THEN
        dcosth=DSQRT(rannum(1))
    ELSE
        dcosth=-DSQRT(rannum(1))
    ENDIF
    
    dSinthCosphi=DSQRT(1d0-rannum(1))*DCOS(M_PI_2*rannum(2))
    dSinthSinphi=DSQRT(1d0-rannum(1))*DSIN(M_PI_2*rannum(2))
    
    !計算t1、t2
    IF ((neighbor.eq.2).OR.(neighbor.eq.3)) THEN
        t1(:)=vecCoords(:,1,k)  !vecCoords(:,1,k) is the vector from xynodes(:,Element(6,k)) to xynodes(:,Element(7,k))
    ELSE IF ((neighbor.eq.4).OR.(neighbor.eq.5)) THEN
        t1(:)=xynodes(:,Element(9,k))-xynodes(:,Element(8,k))
    ELSE
        PAUSE 'Neighbor is wrong in subroutine diffuseB!'
    ENDIF
    
    t1=t1/DSQRT(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3)) !normalized
    t2(:)=CROSS_PRODUCT(normVec(:,neighbor-1,k),t1)
        
    ! velocity = cos(theta)*normal vector + sin(theta)*cos(phi)*tangential vector1 
    !                    + sin(theta)*sin(phi)*tangential vector2
    
    phm(4:6)=dcosth*normVec(:,neighbor-1,k)+dSinthCosphi*t1(:)+dSinthSinphi*t2(:)
    !This is only unit direction vector!!!
    
END SUBROUTINE diffuseB
!============================================================================
SUBROUTINE proc_reorder_and_Ecell
IMPLICIT NONE
INTEGER*4:: k,tot,m,temp(Nelement)
REAL*8::tempR(Nelement)

    temp=0
    tempR=0
    !$OMP PARALLEL DEFAULT(SHARED),FIRSTPRIVATE(phn),PRIVATE(m,k)
    !$OMP DO REDUCTION(+:temp,tempR)
    DO m=1,Nph0
        IF (phn(7,m).gt.0) THEN
            k=INT(phn(10,m))
            temp(k)=temp(k)+1
            tempR(k)=tempR(k)+phn(7,m)
        ENDIF
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    Nnumcell=temp
    tempEcell=tempR+dEdiff
    
    DO m=1,heatNph
        IF (phnheat(7,m).gt.0) THEN
            k=INT(phnheat(10,m))
            tempEcell(k)=tempEcell(k)+phnheat(7,m)
        ENDIF
    ENDDO
    
    tot=0
    DO k=1,Nelement
        Nbgcell(k)=tot
        tot=tot+Nnumcell(k)
    ENDDO
    Nph=tot
    
    Nnumcell=0
    DO m=1,Nph0
        IF (phn(7,m).gt.0) THEN
            k=phn(10,m)
            Nnumcell(k)=Nnumcell(k)+1
            address(Nbgcell(k)+Nnumcell(k))=m
        ENDIF
    ENDDO

END SUBROUTINE proc_reorder_and_Ecell
!============================================================================
SUBROUTINE Etable(mat,idx,E,out,sit) !mat=1:Ge, mat=2:Si, idx=1:temperature, 2:number density, 4:velocity, 5:MFP, 6:specific heat
IMPLICIT NONE
INTEGER*4::mat,idx,iU,sit
REAL*8::E,out,out_a,out_b

    IF (mat.eq.1) THEN
	    iU = INT( (E-Elimit(1,1))/dU_Ge )+1
	    IF (iU.ge.N_Ge2.or.iU.le.0) THEN
	        WRITE(*,*) 'sit=',sit
	        WRITE(*,*) 'E=',E
	        WRITE(*,*) 'iU=',iU
	        WRITE(*,*) 'Ge_start=',Elimit(1,1)
	        WRITE(*,*) 'dU_Ge=',dU_Ge
	        WRITE(*,*) 'mat=',mat
	        WRITE(*,*) 'idx=',idx
	        PAUSE 'out of table!'
	    ENDIF
	    out_a = Ge_table(idx,iU)
	    out_b = Ge_table(idx,iU+1)
	    out = out_a + (out_b-out_a)*(E-Ge_table(3,iU))/dU_Ge
    ELSE IF (mat.eq.2) THEN
	    iU = int( (E-Elimit(1,2))/dU_Si )+1
	    IF (iU.ge.N_Si2.or.iU.le.0) THEN
	        WRITE(*,*) 'sit=',sit
	        WRITE(*,*) 'E=',E
	        WRITE(*,*) 'iU=',iU
	        WRITE(*,*) 'Si_start=',Elimit(1,2)
	        WRITE(*,*) 'dU_Si=',dU_Si
	        WRITE(*,*) 'mat=',mat
	        WRITE(*,*) 'idx=',idx
	        PAUSE 'out of table!'
	    ENDIF
	    out_a = Si_table(idx,iU)
	    out_b = Si_table(idx,iU+1)
	    out = out_a + (out_b-out_a)*(E-Si_table(3,iU))/dU_Si
    END IF

END SUBROUTINE Etable
!============================================================================
SUBROUTINE proc_energy(nc,T0,Eout)
IMPLICIT NONE
INTEGER*4 :: nc,i,ia,ib,imid,cc
REAL*8 :: T0,Eout,Ttemp,Tb,Tout
    
    IF (nc.eq.1) THEN
        IF ((T0.gt.Ge_table(1,N_Ge2)).OR.(T0.lt.Ge_table(1,1))) THEN
            WRITE(*,*) 'T0=',T0
            WRITE(*,*) 'material=',nc
            PAUSE 'WRONG IN proc_energy'
        ELSE
            ia=1
            ib=N_Ge2
            cc=0
            DO WHILE((ib-ia).ne.1)
                imid=(ia+ib)/2
                Ttemp=T0-Ge_table(1,imid)
                IF (Ttemp.gt.zero_tol) THEN !Ge_table(1,imid) < T0 < Ge_table(1,ib)
                    ia=imid
                ELSE IF (Ttemp.lt.-zero_tol) THEN !Ge_table(1,imid) > T0 > Ge_table(1,ia)
                    ib=imid
                ELSE IF (DABS(Ttemp).le.zero_tol) THEN !Ge_table(1,imid) ~ T0
                    ib=imid
                    ia=imid-1
                    EXIT
                ENDIF
                cc=cc+1
                IF (cc.gt.1000000) THEN
	                WRITE(*,*) 'WRONG IN SUBROUTINE proc_energy 1!!'
	                WRITE(*,*) 'ia,ib,ib-ia = ',ia,ib,ib-ia
	                PAUSE '無限迴圈'
	            ENDIF
            ENDDO

            Eout=(T0-Ge_table(1,ia))*dU_Ge/(Ge_table(1,ib)-Ge_table(1,ia))+Ge_table(3,ia)
        ENDIF
    ELSE IF (nc.eq.2) THEN
        IF ((T0.gt.Si_table(1,N_Si2)).OR.(T0.lt.Si_table(1,1))) THEN
            WRITE(*,*) 'T0=',T0
            WRITE(*,*) 'material=',nc
            PAUSE 'WRONG IN proc_energy'
        ELSE
            ia=1
            ib=N_Si2
            cc=0
            DO WHILE((ib-ia).ne.1)
                imid=(ia+ib)/2
                Ttemp=T0-Si_table(1,imid)
                IF (Ttemp.gt.zero_tol) THEN !Si_table(1,imid) < T0 < Si_table(1,ib)
                    ia=imid
                ELSE IF (Ttemp.lt.-zero_tol) THEN !Si_table(1,imid) > T0 > Si_table(1,ia)
                    ib=imid
                ELSE IF (DABS(Ttemp).le.zero_tol) THEN !Si_table(1,imid) ~ T0
                    ib=imid
                    ia=imid-1
                    EXIT
                ENDIF
                cc=cc+1
                IF (cc.gt.1000000) THEN
	                WRITE(*,*) 'WRONG IN SUBROUTINE proc_energy 2!!'
	                WRITE(*,*) 'ia,ib,ib-ia = ',ia,ib,ib-ia
	                PAUSE '無限迴圈'
	            ENDIF
            ENDDO
            Eout=(T0-Si_table(1,ia))*dU_Si/(Si_table(1,ib)-Si_table(1,ia))+Si_table(3,ia)
        ENDIF
    ENDIF
END SUBROUTINE proc_energy
!============================================================================
SUBROUTINE Find_Periodic_Neighbor(phm,NBC,BCelement,s,k)
IMPLICIT NONE 
INTEGER*4::s,i,k,NBC,BCelement(NBC)
REAL*8::A(3),r(3),phm(Nprop)

    s=-1
    DO i=1,NBC
        !If the phonon is in the i-th boundary cell
        ! => r(:)=A(1)*vecCoords(:,1,BCelement(i))+A(2)*vecCoords(:,2,BCelement(i))+A(3)*vecCoords(:,3,BCelement(i)), 
        ! Sove this equations, then:
        !     1. A(1)、A(2)、A(3) must be greater than 0 or equal to 0
        !     2. A(1)+A(2)+A(3) must be smaller than 1 or equal to 1
        r(:)=phm(1:3)-xynodes(:,Element(6,BCelement(i)))
        A=MATMUL(invers_vec(:,:,BCelement(i)),r)
        IF ( ALL(A.ge.-0.00001).AND.(SUM(A).le.1.00001)) THEN
            s=i
            EXIT
        ENDIF
    ENDDO
    
    IF (s.eq.-1) THEN
        WRITE(30,*) 'WRONG in subroutine Find_Periodic_Neighbor'
        WRITE(30,*) 'phm'
        WRITE(30,*) phm
        WRITE(30,*) '原本網格'
        WRITE(30,*) k
        WRITE(30,*) Element(:,k)
        WRITE(30,*) '網格節點'
        WRITE(30,"(4(3(F15.11),/))") xyNodes(:,Element(6:9,k))
        PAUSE 'WRONG in subroutine Find_Periodic_Neighbor'
    ENDIF

END SUBROUTINE Find_Periodic_Neighbor
!============================================================================
SUBROUTINE Find_BCelement(k,NBC,BCelement,s) !Find the local label in BCelement of a cell that its global lable is k
IMPLICIT NONE
INTEGER*4::k,s,NBC,BCelement(NBC)
INTEGER*4::i,temp(1)

    IF (COUNT(BCelement.eq.k).eq.1) THEN
        temp=PACK((/ (i,i=1,NBC) /),BCelement.eq.k)
        s=temp(1)
    ELSE
        PAUSE 'WRONG in subroutine Find_BCelement'
    ENDIF
    
END SUBROUTINE Find_BCelement
!=================================================================================
SUBROUTINE proc_BC(TL,TR)
IMPLICIT NONE
INTEGER*4,ALLOCATABLE:: ms(:)
INTEGER*4:: j,mt
REAL*8:: dUL(2),dUR(2)
REAL*8:: TL,TR,dNum

    DO mt=1,2
        CALL proc_energy(mt,TL,dUL(mt))
        CALL Etable(mt,4,dUL(mt),dVinjectL(mt),0)
        CALL Etable(mt,2,dUL(mt),dNum,0)
        dEinjectL(mt)=dUL(mt)/dNum*bundle(mt)
   
        CALL proc_energy(mt,TR,dUR(mt))
        CALL Etable(mt,4,dUR(mt),dVinjectR(mt),0)
        CALL Etable(mt,2,dUR(mt),dNum,0)
        dEinjectR(mt)=dUR(mt)/dNum*bundle(mt)
    ENDDO
    
    ALLOCATE( ms(NBCL) )
    ms=GrainMt(Element(1,BCelementL(:)))
    dEheatfluxL(:)=dUL(ms(:))*dVinjectL(ms(:))*tStep*areaBCL(:)/4d0   ! meV
    DEALLOCATE( ms )
    
    ALLOCATE( ms(NBCR) )
    ms=GrainMt(Element(1,BCelementR(:)))
    dEheatfluxR(:)=dUR(ms(:))*dVinjectR(ms(:))*tStep*areaBCR(:)/4d0   ! meV
    DEALLOCATE( ms )
    
END SUBROUTINE proc_BC
!============================================================================
!============================================================================
SUBROUTINE cross_center_section(oldx,newx,energy)
IMPLICIT NONE
REAL*8::oldx,newx,energy
REAL*8::x1,x2

    x1=oldx-dLdomain(1)/2d0
    x2=newx-dLdomain(1)/2d0

    IF ((x1.le.0).AND.(x2.gt.0)) THEN
        qcenter=qcenter+energy
    ELSE IF ((x1.ge.0).AND.(x2.lt.0)) THEN
        qcenter=qcenter-energy
    ENDIF

END SUBROUTINE cross_center_section
!============================================================================
!============================================================================
SUBROUTINE cellinfo
IMPLICIT NONE
INTEGER*4::k,mt,i,j
REAL*8::tempE

    SELECT CASE(cellinfoMethod)
    CASE(1)

        !$OMP PARALLEL DEFAULT(SHARED),PRIVATE(k,mt)
        !$OMP DO
        DO k=1,Nelement
            mt=GrainMt(Element(1,k))
            dEcell(k) = tempEcell(j)/Volume(k)  ! total energy per unit volume of element k
            CALL Etable(mt,2,dEcell(k),dEunit(k),4) ! compute the equilibrium number of phonons 
            CALL Etable(mt,4,dEcell(k),dVunit(k),5) ! compute the average phonon group velocity 
            CALL Etable(mt,5,dEcell(k),MFP(k),6) !calculate MFP in this cell
            CALL Etable(mt,1,dEcell(k),dTemp(k),7) ! calculate temperature
            dEunit(k)=dEcell(k)/dEunit(k)*bundle(mt)  ! average energy per phonon bundle
        ENDDO
        !$OMP END DO
        !$OMP END PARALLEL
    
    CASE(2)

        !$OMP PARALLEL DEFAULT(SHARED),PRIVATE(k,tempE,mt,i,j)
        !$OMP DO
        DO k=1,Nelement
            tempE=0
            DO i=1,ElemGroup(k)%N
                j=ElemGroup(k)%array(i)
                tempE=tempE+tempEcell(j)
            ENDDO
            mt=GrainMt(Element(1,k))
            dEcell(k)=tempE/Vgroup(k)
        
            CALL Etable(mt,2,dEcell(k),dEunit(k),4) ! compute the equilibrium number of phonons 
            CALL Etable(mt,4,dEcell(k),dVunit(k),5) ! compute the average phonon group velocity 
            CALL Etable(mt,5,dEcell(k),MFP(k),6) !calculate MFP in this cell
            CALL Etable(mt,1,dEcell(k),dTemp(k),7) ! calculate temperature
            dEunit(k)=dEcell(k)/dEunit(k)*bundle(mt)  ! average energy per phonon bundle
        ENDDO
        !$OMP END DO
        !$OMP END PARALLEL
    
    ENDSELECT
    
END SUBROUTINE cellinfo
!==================================================================================
FUNCTION CROSS_PRODUCT(A,B)
IMPLICIT NONE
REAL*8::A(3),B(3)
REAL*8::CROSS_PRODUCT(3)
    
    CROSS_PRODUCT(1)=A(2)*B(3)-A(3)*B(2)
    CROSS_PRODUCT(2)=A(3)*B(1)-A(1)*B(3)
    CROSS_PRODUCT(3)=A(1)*B(2)-A(2)*B(1)

END FUNCTION CROSS_PRODUCT
!=====================================================================!
!=====================================================================!
END MODULE mod_ROUTINES