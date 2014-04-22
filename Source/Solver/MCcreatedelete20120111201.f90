MODULE mod_CreateDelete
USE mod_VARIABLES
USE RNG
USE omp_lib
IMPLICIT NONE
CONTAINS
!=====================================================================!
!=====================================================================!
SUBROUTINE proc_createdelete
IMPLICIT NONE
INTEGER*4::k,p,m,s,newNph,i,cc,ntake,Ndelete
REAL*8::random1,tmp
INTEGER*4,ALLOCATABLE::nadd(:),voids(:)
REAL*8,ALLOCATABLE::rannum(:,:),newphn(:,:)

    ALLOCATE( nadd(Nelement) )
    nadd=0
    Ndelete=0
    
    DO k=1,Nelement
        tmp=0.5d0*dEunit(k)
        IF (dEdiff(k).ge.tmp) THEN
            nadd(k) = INT(dEdiff(k)/dEunit(k)+0.5d0)
	        dEdiff(k)=dEdiff(k)-DBLE(nadd(k))*dEunit(k)
        ELSE IF (dEdiff(k).lt.-tmp) THEN
            cc=0
            ntake=0
            DO WHILE ((dEdiff(k).lt.-tmp).AND.(Nnumcell(k).gt.ntake))
	            CALL RAN_NUM(RNseed0,random1)
		        p=Nbgcell(k)+MIN(INT(random1*Nnumcell(k))+1,Nnumcell(k))
		        IF (phn(7,address(p)).gt.0) THEN
		            dEdiff(k)=dEdiff(k)+phn(7,address(p))
		            phn(7,address(p))=0
		            ntake=ntake+1
	            ENDIF
	            cc=cc+1
	            IF (cc.gt.1000000) THEN
	                WRITE(*,*) 'WRONG IN SUBROUTINE proc_createdelete 1!!'
	                WRITE(*,*) 'k,dEdiff(k),-tmp = ',k,dEdiff(k),-tmp
	                PAUSE '無限迴圈'
	            ENDIF
	        ENDDO
	        Ndelete=Ndelete+ntake
        ENDIF
    ENDDO

    p=SUM(nadd)+heatNph

    IF (p.gt.0) THEN  ! add phonons to the simulations

        newNph=(Nph-Ndelete)+p     ! real phonon numbers in the domain now
        
!-----------------------------------------newNph > Nph0  phn矩陣不夠用---------------------------------------------
        IF (newNph.gt.Nph0) THEN   ! must to re-allocate matrix phn, because real phonon numbers is greater than phn
            ALLOCATE( newphn(Nprop,newNph) )
            s=0
            DO m=1,Nph0
                IF (phn(7,m).gt.0) THEN
	                s=s+1
		            newphn(:,s)=phn(:,m)
	            ENDIF
            ENDDO
            DEALLOCATE( phn,address )
            
            IF (WAY_DIR.gt.0) THEN
                newphn(:,s+1:s+heatNph)=phnheat
                s=s+heatNph
                DEALLOCATE( phnheat )
            ENDIF
            
            DO k=1,Nelement
                IF (nadd(k).gt.0) THEN
                    ALLOCATE( rannum(nadd(k),5) )
	                CALL RAN_NUM(RNseed0,rannum)
	                rannum(:,3)=(1d0-rannum(:,1))*(1d0-rannum(:,2))*rannum(:,3)
		            rannum(:,2)=(1d0-rannum(:,1))*rannum(:,2)
	                newphn(1,s+1:s+nadd(k))=xyNodes(1,Element(6,k))+rannum(:,1)*vecCoords(1,1,k)+rannum(:,2)*vecCoords(1,2,k)+rannum(:,3)*vecCoords(1,3,k)
		            newphn(2,s+1:s+nadd(k))=xyNodes(2,Element(6,k))+rannum(:,1)*vecCoords(2,1,k)+rannum(:,2)*vecCoords(2,2,k)+rannum(:,3)*vecCoords(2,3,k)
		            newphn(3,s+1:s+nadd(k))=xyNodes(3,Element(6,k))+rannum(:,1)*vecCoords(3,1,k)+rannum(:,2)*vecCoords(3,2,k)+rannum(:,3)*vecCoords(3,3,k)
		            newphn(4,s+1:s+nadd(k))=2d0*rannum(:,4)-1d0
		            newphn(5,s+1:s+nadd(k))=DSQRT(1d0-newphn(4,s+1:s+nadd(k))**2)*DCOS(M_PI_2*rannum(:,5))*dVunit(k)
		            newphn(6,s+1:s+nadd(k))=DSQRT(1d0-newphn(4,s+1:s+nadd(k))**2)*DSIN(M_PI_2*rannum(:,5))*dVunit(k)
		            newphn(7,s+1:s+nadd(k))=dEunit(k)
		            newphn(8,s+1:s+nadd(k))=dVunit(k)
		            newphn(9,s+1:s+nadd(k))=GrainMt(Element(1,k))
		            newphn(10,s+1:s+nadd(k))=k
		            newphn(4,s+1:s+nadd(k))=newphn(4,s+1:s+nadd(k))*dVunit(k)
	                s=s+nadd(k)
		            DEALLOCATE( rannum )
	            ENDIF
            ENDDO
            
            Nph=newNph
            Nph0=DBLE(Nph)*1.002d0
            ALLOCATE( phn(Nprop,Nph0),address(Nph0) )
            phn=0
            phn(1:Nprop,1:Nph)=newphn
            DEALLOCATE( newphn )
!-----------------------------------------newNph < Nph0  phn矩陣夠用，要把新的聲子填入phn的空位---------------------------------------------
        ELSE
            ALLOCATE( voids(p) ) ! the number of phonons that will add into domain
            i=0
            m=1
            cc=0
            DO WHILE (i.lt.p)
                IF (phn(7,m).le.0) THEN
                    i=i+1
	                voids(i)=m
                ENDIF
	            m=m+1
	            cc=cc+1
	            IF (cc.gt.Nph0) THEN
	                WRITE(*,*) 'WRONG IN SUBROUTINE proc_createdelete 2!!'
	                WRITE(*,*) 'i,p = ',i,p
	                PAUSE '無限迴圈'
	            ENDIF
            ENDDO
            
            s=0
            
            IF (WAY_DIR.gt.0) THEN
                phn(:,voids(s+1:s+heatNph))=phnheat
	            s=s+heatNph
	            DEALLOCATE( phnheat )
            ENDIF
   
            DO k=1,Nelement
                IF (nadd(k).gt.0) THEN
                    ALLOCATE( rannum(nadd(k),5) )
	                CALL RAN_NUM(RNseed0,rannum)
	                rannum(:,3)=(1d0-rannum(:,1))*(1d0-rannum(:,2))*rannum(:,3)
		            rannum(:,2)=(1d0-rannum(:,1))*rannum(:,2)
	                phn(1,voids(s+1:s+nadd(k)))=xyNodes(1,Element(6,k))+rannum(:,1)*vecCoords(1,1,k)+rannum(:,2)*vecCoords(1,2,k)+rannum(:,3)*vecCoords(1,3,k)
		            phn(2,voids(s+1:s+nadd(k)))=xyNodes(2,Element(6,k))+rannum(:,1)*vecCoords(2,1,k)+rannum(:,2)*vecCoords(2,2,k)+rannum(:,3)*vecCoords(2,3,k)
		            phn(3,voids(s+1:s+nadd(k)))=xyNodes(3,Element(6,k))+rannum(:,1)*vecCoords(3,1,k)+rannum(:,2)*vecCoords(3,2,k)+rannum(:,3)*vecCoords(3,3,k)
		            phn(4,voids(s+1:s+nadd(k)))=2d0*rannum(:,4)-1d0
		            phn(5,voids(s+1:s+nadd(k)))=DSQRT(1d0-phn(4,voids(s+1:s+nadd(k)))**2)*DCOS(M_PI_2*rannum(:,5))*dVunit(k)
		            phn(6,voids(s+1:s+nadd(k)))=DSQRT(1d0-phn(4,voids(s+1:s+nadd(k)))**2)*DSIN(M_PI_2*rannum(:,5))*dVunit(k)
		            phn(7,voids(s+1:s+nadd(k)))=dEunit(k)
		            phn(8,voids(s+1:s+nadd(k)))=dVunit(k)
		            phn(9,voids(s+1:s+nadd(k)))=GrainMt(Element(1,k))
 	                phn(10,voids(s+1:s+nadd(k)))=k
		            phn(4,voids(s+1:s+nadd(k)))=phn(4,voids(s+1:s+nadd(k)))*dVunit(k)
	                s=s+nadd(k)
		            DEALLOCATE( rannum )
	            ENDIF
            ENDDO
            
            Nph=newNph
            
            DEALLOCATE( voids )   
        ENDIF
!--------------------------------------------------------------------------------------
    ENDIF
 
    IF (ALLOCATED(phnheat)) DEALLOCATE(phnheat)
    DEALLOCATE( nadd )
    
END SUBROUTINE proc_createdelete
!=====================================================================!
!=====================================================================!
END MODULE mod_CreateDelete