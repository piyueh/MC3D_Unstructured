!****************************************************************************
!   HEADING: MC3D ADVANCE MODULE FOR UNSTRUCTURED GRID
!   AUTHOR: 
!   PURPOSE: THIS MODULE CONTAINS THE SUBROUTINES ABOUT PHONON ADVANCE
!   DATE :  
!****************************************************************************
MODULE mod_ADVANCE
USE mod_VARIABLES
USE mod_ROUTINES
USE RNG
USE omp_lib
IMPLICIT NONE
CONTAINS
!=====================================================================!
!=====================================================================!
SUBROUTINE proc_advection(cphn,dt,nc,RNseed)
IMPLICIT NONE
!-------------------------------------!
TYPE(rng_t),INTENT(inout)::RNseed
INTEGER*4,INTENT(in)::nc
REAL*8,INTENT(inout)::cphn(Nprop)
REAL*8,INTENT(in)::dt
!-------------------------------------!
INTEGER*4::phcell,neighbor,cc,true
REAL*8::vecA(3),vecB(3),old_phn(Nprop),old_dt,oldx
REAL*8::dtremain,dtused
!-------------------------------------!

    dtremain=dt
    phcell=cphn(10)
    IF (nc.eq.2) THEN
        old_phn=cphn
        old_dt=dt
    ENDIF
    
    vecA=xyNodes(:,Element(6,phcell))-cphn(1:3)
    vecB=xyNodes(:,Element(7,phcell))-cphn(1:3)
    
    cc=0
    DO WHILE (DABS(dtremain).gt.1d-7)
    
        CALL hit_cell_boundary(phcell,cphn,vecA,vecB,neighbor,dtused)
        
        IF (dtused.gt.dtremain) THEN
            dtused=dtremain
            dtremain=0
            oldx=cphn(1)
            cphn(1:3)=cphn(1:3)+dtused*cphn(4:6)
            CALL cross_center_section(oldx,cphn(1),cphn(7))
            CALL proc_intrinsicscattering(phcell,cphn,dtused,RNseed)
        ELSE
            oldx=cphn(1)
            cphn(1:3)=cphn(1:3)+dtused*cphn(4:6)
            CALL cross_center_section(oldx,cphn(1),cphn(7))
            dtremain=dtremain-dtused
            CALL proc_intrinsicscattering(phcell,cphn,dtused,RNseed)
            
            true=0
            ! CHECK IF THE PHONON STILL BE LEAVING THE ELEMENT.
            IF (DOT_PRODUCT(cphn(4:6),normVec(:,neighbor-1,phcell)).gt.0) THEN 
                ! CHECK IF THE PHONON IS LEAVING THE COMPUTATIONAL DOMAIN
                CALL proc_outdomain(nc,phcell,neighbor,cphn,dtremain,true,RNseed)
                ! IF true=0, PHONON IS STILL INSIDE THE DOMAIN, AND IT'S GOING TO TRANSFER TO THE OTHER ELEMENT
                IF (true.eq.0) CALL proc_transmissivity(phcell,neighbor,cphn,RNseed)
            ENDIF
            
            IF ((true.eq.-2).AND.(nc.eq.2)) THEN
                cphn=old_phn
                phcell=cphn(10)
                dtremain=old_dt
            ENDIF
            
            IF ((true.ne.-2).OR.(nc.ne.1)) THEN
                vecA=xynodes(:,Element(6,phcell))-cphn(1:3)
                vecB=xynodes(:,Element(7,phcell))-cphn(1:3)
            ENDIF
            
        ENDIF
    
        cc=cc+1
        IF (cc.gt.1d6) THEN
            WRITE(*,*) 'WRONG IN SUBROUTINE proc_advection!!'
	        WRITE(*,*) 'dtremain = ',dtremain
	        PAUSE '無限迴圈'
        ENDIF
    ENDDO

END SUBROUTINE proc_advection
!=====================================================================!
!=====================================================================!
SUBROUTINE hit_cell_boundary(k,phm,r1,r2,neighbor,dtbdy)
IMPLICIT NONE
REAL*8:: r1(3),r2(3),phm(Nprop),dtbdy,dt(4),r,v
INTEGER*4:: k,neighbor,i,hittingface(1)

    DO i=1,4
        v=DOT_PRODUCT(phm(4:6),normVec(:,i,k))
        IF ( v.le.0 ) THEN !dt 小於、等於0表示不會穿越此面
            dt(i)=10000000
        ELSE
            SELECT CASE(i)
            CASE(1:3)
                r=DOT_PRODUCT(r1(:),normVec(:,i,k))            
            CASE(4)
                r=DOT_PRODUCT(r2(:),normVec(:,i,k))
            ENDSELECT
            IF (r.lt.0) r=0 !理論上r必>=0，但有時會有數值誤差，固強迫為0
            dt(i)=r/v
        ENDIF
    ENDDO
    
    hittingface=MINLOC( dt )
    neighbor=hittingface(1)+1
    dtbdy=dt(hittingface(1))
    
END SUBROUTINE hit_cell_boundary
!=====================================================================!
!=====================================================================!
SUBROUTINE proc_intrinsicscattering(k,phm,dt1,RNseed)
IMPLICIT NONE
TYPE(rng_t),INTENT(inout)::RNseed
INTEGER*4::k
REAL*8::phm(Nprop),dt1,prob
REAL*8::rannum(2)

    prob=1d0-DEXP(-dt1*phm(8)/MFP(k))
    CALL RAN_NUM(RNseed,rannum(1))

    IF (rannum(1).le.prob) THEN
        CALL RAN_NUM(RNseed,rannum)
        dEdiff(k)=dEdiff(k)+phm(7)-dEunit(k)
        phm(7)=dEunit(k)
        phm(8)=dVunit(k)
        phm(9)=GrainMt(Element(1,k))
        phm(4)=2d0*rannum(1)-1d0
        phm(5)=DSQRT(1d0-phm(4)**2)*DCOS(M_PI_2*rannum(2))*phm(8)
        phm(6)=DSQRT(1d0-phm(4)**2)*DSIN(M_PI_2*rannum(2))*phm(8)
        phm(4)=phm(4)*phm(8)
    ENDIF

END SUBROUTINE proc_intrinsicscattering
!=====================================================================!
!=====================================================================!
SUBROUTINE proc_outdomain(nc,k,neighbor,phm,dtremain,true,RNseed)
IMPLICIT NONE
!-------------------------------------!
TYPE(rng_t),INTENT(inout)::RNseed
INTEGER*4,INTENT(in)::nc,k,neighbor
INTEGER*4,INTENT(inout)::true
REAL*8,INTENT(inout)::phm(Nprop),dtremain
!-------------------------------------!
INTEGER*4::s,p,tempI
REAL*8::a,b,rannum(2)
!---------------------------------------------------------------------!
!!!!!! = Nelement+1 means a periodic b.c. in the first direction
!!!!!! = Nelement+2 means a periodic b.c. in the second direction
!!!!!! = Nelement+3 means a periodic b.c. in the third direction
!!!!!! = Nelement+4 means an adiabatic b.c. required
!!!!!! = Nelement+5 means there is a prescribed heat flux 
!---------------------------------------------------------------------!

    tempI=Element(neighbor,k)-Nelement
    SELECT CASE(tempI)
    CASE DEFAULT
        RETURN
    CASE(1)
        IF (phm(4).gt.0) THEN
            phm(1)=phm(1)-dLdomain(1)+zero_tol
        ELSE
            phm(1)=phm(1)+dLdomain(1)-zero_tol
        ENDIF
    CASE(2)
        IF (phm(5).gt.0) THEN
            phm(2)=phm(2)-dLdomain(2)+zero_tol
        ELSE
            phm(2)=phm(2)+dLdomain(2)-zero_tol
        ENDIF
    CASE(3)
        IF (phm(6).gt.0) THEN
            phm(3)=phm(3)-dLdomain(3)+zero_tol
        ELSE
            phm(3)=phm(3)+dLdomain(3)-zero_tol
        ENDIF
!-----------------------------------------------------------------------------
    CASE(4)
        CALL RAN_NUM(RNseed,rannum(1))
        IF (rannum(1).le.DPPB) THEN
            phm(4:6)=phm(4:6)-2*DOT_PRODUCT(phm(4:6),normVec(:,neighbor-1,k))*normVec(:,neighbor-1,k)
        ELSE
            CALL RAN_NUM(RNseed,rannum)
            CALL diffuseB( k,neighbor,phm,-1,rannum ) !ATTENTION: phm(4:6) returned from diffuseB are only directions,not velocity vector
	        dEdiff(k)=dEdiff(k)+phm(7)-dEunit(k) !phonon will scattered after reflection, i.e. property will be the same with the cell
	        phm(9)=GrainMt(Element(1,k))
	        phm(8)=dVunit(k)
	        phm(7)=dEunit(k)
            phm(4:6)=phm(4:6)*phm(8)
        ENDIF
        true=1
!-------------------------------------------------------------------------------
    CASE(5)
        IF (nc.eq.1) THEN !nc=1 means that now is advance procedure,and 2 is heatcontral procedure
	        IF (phm(4).gt.0) THEN !leave domain from right side
	            CALL Find_BCelement(k,NBCR,BCelementR,s) ! for cell-to-cell periodically asign
	            qbdyR(s)=qbdyR(s)+phm(7)
	            IF (WAY_DIR.eq.1) THEN !perodically assigned
	                phm(1)=phm(1)-dLdomain(1)+zero_tol
                    CALL Find_Periodic_Neighbor(phm,NBCL,BCelementL,s,k)
                    !$OMP CRITICAL
	                mlostL(s)=mlostL(s)+1
			        IF (mlostL(s).gt.Nmakeup) mlostL(s)=1
			        p=mlostL(s)
			        !$OMP END CRITICAL
		            dPoolL(1,p,s)=dtremain
                    dPoolL(2:4,p,s)=phm(4:6)/phm(8)
			        dPoolL(5,p,s)=phm(9)
			        dPoolL(6:7,p,s)=phm(2:3)
			    ELSE IF (WAY_DIR.eq.2) THEN ! for material-to-material periodically asign
			        s=GrainMt(Element(1,k))
			        !$OMP CRITICAL
			        mlostL(s)=mlostL(s)+1
			        IF (mlostL(s).gt.Nmakeup) mlostL(s)=1
			        p=mlostL(s)
			        !$OMP END CRITICAL
			        dPoolL(1,p,s)=dtremain
                    dPoolL(2:4,p,s)=phm(4:6)/phm(8)
			        dPoolL(5,p,s)=phm(9)
			    ENDIF
	        ELSE !phonon leave domain from left side
	            CALL Find_BCelement(k,NBCL,BCelementL,s)
	            qbdyL(s)=qbdyL(s)+phm(7)
	            IF (WAY_DIR.eq.1) THEN
	                phm(1)=phm(1)+dLdomain(1)-zero_tol
	                CALL Find_Periodic_Neighbor(phm,NBCR,BCelementR,s,k)
	                !$OMP CRITICAL
	                mlostR(s)=mlostR(s)+1
			        IF (mlostR(s).gt.Nmakeup) mlostR(s)=1
			        p=mlostR(s)
			        !$OMP END CRITICAL
		            dPoolR(1,p,s)=dtremain
                    dPoolR(2:4,p,s)=phm(4:6)/phm(8)
			        dPoolR(5,p,s)=phm(9)
			        dPoolR(6:7,p,s)=phm(2:3)
			    ELSE IF (WAY_DIR.eq.2) THEN
			        s=GrainMt(Element(1,k))
			        !$OMP CRITICAL
			        mlostR(s)=mlostR(s)+1
			        IF (mlostR(s).gt.Nmakeup) mlostR(s)=1
			        p=mlostR(s)
			        !$OMP END CRITICAL
		            dPoolR(1,p,s)=dtremain
                    dPoolR(2:4,p,s)=phm(4:6)/phm(8)
			        dPoolR(5,p,s)=phm(9)
			    ENDIF
	        ENDIF
	    ENDIF   
        phm(7)=0
	    dtremain=0
	    true=-2
    END SELECT

END SUBROUTINE proc_outdomain
!=====================================================================!
!=====================================================================!
SUBROUTINE proc_transmissivity(k,neighbor,phm,RNseed)
IMPLICIT NONE
!-------------------------------------!
TYPE(rng_t),INTENT(inout)::RNseed
INTEGER*4,INTENT(inout)::k
INTEGER*4,INTENT(in)::neighbor
REAL*8,INTENT(inout)::phm(Nprop)
!-------------------------------------!
INTEGER*4::k1,s
INTEGER*4::mtk1,mtk
REAL*8::tau12,ratio,rannum(2)
REAL*8::dcosth1,dsinth1,dcosth2,dsinth2
REAL*8::neighborE,neighborV,z1,z2
!-------------------------------------!

    k1=Element(neighbor,k)
    !-----------------------------------------------------------------!
    IF (k1.gt.Nelement) THEN
        k1=k1-Nelement
        SELECT CASE(k1)
        CASE(1)
            IF (phm(4).gt.0) THEN
                CALL Find_Periodic_Neighbor(phm,NBCL,BCelementL,s,k) !因為週期性BC，跳到另一邊的網格後，不知道新網格的編號
                k1=BCelementL(s)
            ELSE
                CALL Find_Periodic_Neighbor(phm,NBCR,BCelementR,s,k)
                k1=BCelementR(s)
            ENDIF
        CASE(2)
            IF (phm(5).gt.0) THEN
                CALL Find_Periodic_Neighbor(phm,NBCyN,BCelementyN,s,k) 
                k1=BCelementyN(s)
            ELSE
                CALL Find_Periodic_Neighbor(phm,NBCyP,BCelementyP,s,k)
                k1=BCelementyP(s)
            ENDIF
        CASE(3)
            IF (phm(6).gt.0) THEN
                CALL Find_Periodic_Neighbor(phm,NBCzN,BCelementzN,s,k)
                k1=BCelementzN(s)
            ELSE
                CALL Find_Periodic_Neighbor(phm,NBCzP,BCelementzP,s,k)
                k1=BCelementzP(s)
            ENDIF
        CASE DEFAULT
            PAUSE "WRONG"
        END SELECT
    ENDIF
    !-----------------------------------------------------------------!
    !-----------------------------------------------------------------!
    IF (Element(1,k1).eq.Element(1,k)) THEN ! NEIGHBOR IS IN THE SAME GRAIN, SO IT MUST BE THE SAME MATERIAL
        k=k1
        phm(10)=k1
        RETURN
    ELSE
        mtk1=GrainMt(Element(1,k1))
        mtk=GrainMt(Element(1,k))
        CALL proc_Energy(mtk1,dTemp(k),neighborE)
        CALL Etable(mtk1,4,neighborE,neighborV,3)
        
        CALL RAN_NUM(RNseed,rannum)
        
        IF (rannum(1).le.DPP) THEN ! SPECULAR RESPONSE
            ratio=(dEcell(k)*dVunit(k))/(neighborE*neighborV)
            dcosth1=DOT_PRODUCT(phm(4:6),normVec(:,neighbor-1,k))/phm(8)
            dsinth1=DSQRT(1d0-dcosth1**2)
	        dsinth2=dsinth1*DSQRT(ratio)
	        IF (dsinth2.lt.1d0) THEN
	            z1=rho(GrainMt(Element(1,k)))*dVunit(k)
		        z2=rho(GrainMt(Element(1,k1)))*neighborV
		        dcosth2=DSQRT(1d0-dsinth2**2)
		        tau12=(z2*dcosth2)/(z1*dcosth1)
		        tau12=1d0-((1d0-tau12)/(1d0+tau12))**2
		    ELSE
		        tau12=0
	        ENDIF
	        
	        IF (rannum(2).lt.tau12) THEN ! SPECULARLY TRANSMITTED (IAMM)
	        ! POSITION(1,2,3), ENERGY(7), AND ENERGY_MATERIAL(9) REMAIN UNCHANGED.
	            CALL Snells( k,neighbor,phm,dcosth1,dsinth1,dcosth2,dsinth2 )
	            phm(8)=dVunit(k1)
	            phm(4:6)=phm(4:6)*phm(8)
	            phm(10)=k1
	            k=k1
	        ELSE ! SPECULARLY REFLECTED
	            phm(4:6)=phm(4:6)-2*DOT_PRODUCT(phm(4:6),normVec(:,neighbor-1,k))*normVec(:,neighbor-1,k)
	            IF (k1.ne.Element(neighbor,k)) CALL back_to_original_side
	        ENDIF
        ELSE ! DIFFUSE RESPONSE
            tau12=(neighborE*neighborV)/(dEcell(k)*dVunit(k)+neighborE*neighborV)
            
            IF (rannum(2).lt.tau12) THEN ! DIFFUSELY TRANSMITTED (DAMM)
                CALL RAN_NUM(RNseed,rannum)
                CALL diffuseB( k,neighbor,phm,1,rannum )
                phm(8)=dVunit(k1)
                phm(4:6)=phm(4:6)*phm(8)
                k=k1
                phm(10)=k1
            ELSE ! DIFFUSELY REFLECTED (DAMM)
            ! FROM INELASTIC ACOUSTIC MISMATCH MODEL, THE PHONON ENERGY WON'T CHANGE.
                CALL RAN_NUM(RNseed,rannum)
                CALL diffuseB( k,neighbor,phm,-1,rannum )
                phm(8)=dVunit(k)
                phm(4:6)=phm(8)*phm(4:6)
                IF (k1.ne.Element(neighbor,k)) CALL back_to_original_side
            ENDIF
        ENDIF
    ENDIF
    !-----------------------------------------------------------------!
    !-----------------------------------------------------------------!
    CONTAINS
    SUBROUTINE back_to_original_side
    IMPLICIT NONE
        
        k1=Element(neighbor,k)-Nelement
        SELECT CASE(k1)
        CASE(1)
            IF (phm(4).gt.0) THEN
	            phm(1)=phm(1)-dLdomain(1)+zero_tol
            ELSE
                phm(1)=phm(1)+dLdomain(1)-zero_tol
            ENDIF
        CASE(2)
            IF (phm(5).gt.0) THEN
	            phm(2)=phm(2)-dLdomain(2)+zero_tol
            ELSE
                phm(2)=phm(2)+dLdomain(2)-zero_tol
            ENDIF
        CASE(3)
            IF (phm(6).gt.0) THEN
	            phm(3)=phm(3)-dLdomain(3)+zero_tol
            ELSE
                phm(3)=phm(3)+dLdomain(3)-zero_tol
            ENDIF
        CASE DEFAULT
            PAUSE "WRONG"
        END SELECT
    
    END SUBROUTINE back_to_original_side


END SUBROUTINE proc_transmissivity
!=====================================================================!
!=====================================================================!
END MODULE mod_ADVANCE