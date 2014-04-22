!*****************************************************************
!*****************************************************************
!*****************************************************************
MODULE variables
IMPLICIT NONE
!-----FOR MESH INFORMATION-----!
INTEGER*4:: Nnodes,Nelement,Propelem,Ngrains
INTEGER*4:: NBCL,NBCR,NBCyP,NBCyN,NBCzP,NBCzN
INTEGER*4,ALLOCATABLE:: BCelementL(:),BCelementR(:),BCelementyP(:),BCelementyN(:),BCelementzP(:),BCelementzN(:),Element(:,:),GrainMt(:)
REAL*8,ALLOCATABLE::xyNodes(:,:),Volume(:),center(:,:),normVec(:,:,:),vecCoords(:,:,:),invers_vec(:,:,:), areaBCL(:),areaBCR(:)
!-----FOR COMPOTING OTHER INFORMATIONS-----!
INTEGER*4::BCx,BCy,BCz
INTEGER*4::NBCLnode,NBCRnode,NBCyPnode,NBCyNnode,NBCzPnode,NBCzNnode
INTEGER*4,ALLOCATABLE::BCLnode(:),BCRnode(:),BCyPnode(:),BCyNnode(:),BCzPnode(:),BCzNnode(:)
REAL*8::dLdomain(2,3)
!================================================
CONTAINS
!------------------------------------------------
FUNCTION CROSS_PRODUCT(A,B)
IMPLICIT NONE
REAL*8::A(3),B(3)
REAL*8::CROSS_PRODUCT(3)
    
    CROSS_PRODUCT(1)=A(2)*B(3)-A(3)*B(2)
    CROSS_PRODUCT(2)=A(3)*B(1)-A(1)*B(3)
    CROSS_PRODUCT(3)=A(1)*B(2)-A(2)*B(1)

END FUNCTION CROSS_PRODUCT
!------------------------------------------------
SUBROUTINE read_tmp(LR,tmp)
IMPLICIT NONE
INTEGER*4 stat,i,LR
CHARACTER*72::tmp
    
    i=0
    stat=0
    tmp=' '

    DO WHILE(stat.eq.0)
        i=i+1
        READ(unit=LR,fmt="(A1)",iostat=stat,advance='NO') tmp(i:i)
    ENDDO

END SUBROUTINE read_tmp
!------------------------------------------------
END MODULE variables
!*****************************************************************
!*****************************************************************
!*****************************************************************

PROGRAM main
INCLUDE 'link_fnl_static_hpc.h'
USE IMSL_LIBRARIES
USE variables
IMPLICIT NONE
INTEGER*4:: LR1=100,i,j,k,NGe
INTEGER*4,ALLOCATABLE:: tempMatrix(:,:)
REAL*8::r,rhoGe
CHARACTER*72::temp,readfile,outputfile,GrainMtFile

!=======================================================================================!
    WRITE(*,*) 'ENTER THE INPUT FLIE NAME (NO ".txt") : '
    READ(*,"(A72)") readfile

    readfile=readfile(1:LEN_TRIM(readfile))//'.txt'
    WRITE(*,*) 'INPUT FILE :',readfile
    WRITE(*,*)
    WRITE(*,*) 'ENTER THE GRAIN MATERIAL DISTRIBUTION FILE NAME (NO ".txt") : '
    READ(*,"(A72)") GrainMtFile
    GrainMtFile=GrainMtFile(1:LEN_TRIM(GrainMtFile))//'.txt'
    WRITE(*,*) 'GRAIN MATERIAL FILE:', GrainMtFile
    WRITE(*,*)
    
    WRITE(*,*) 'ENTER BC OF x-DIRECTION (1.PERIODICAL 4.ADIABATIC 5.THERMAL FLOW):'
    READ(*,*) BCx
    WRITE(*,*) 'ENTER BC OF y-DIRECTION (2.PERIODICAL 4.ADIABATIC):'
    READ(*,*) BCy
    WRITE(*,*) 'ENTER BC OF z-DIRECTION (3.PERIODICAL 4.ADIABATIC):'
    READ(*,*) BCz
    
    Propelem=9
    
    outputfile=readfile(1:INDEX(readfile,'.txt')-1)//'_gridfile.txt'
    WRITE(*,*) 'THE OUTPUT FILE NAME WILL BE:',outputfile
!===============================READ ANSYS MESH FILE====================================!
    OPEN(UNIT=LR1,FILE=readfile)
    
    DO
        CALL read_tmp(LR1,temp)
        IF (temp(3:).eq.'NODES') THEN
            WRITE(*,*) '1.READING INFORMATION OF NODES......'
            CALL read_nodes(LR1)
        ELSE IF (temp(3:).eq.'ELEMENTS') THEN
            WRITE(*,*) '2.READING INFORMATION OF ELEMENTS......'
            CALL read_elements(LR1)
        ELSE IF (temp(3:).eq.'COMPONENTS') THEN
            WRITE(*,*) '3.READING BOUNDARY NODES......'
            CALL read_boundary_nodes(LR1)
        ELSE IF (temp.eq.'FINI') THEN
            EXIT
        ENDIF
    ENDDO
    
    CLOSE(LR1)
    
    BCx=BCx+Nelement
    BCy=BCy+Nelement
    BCz=BCz+Nelement
    
    !CALL delete_wrong_BC_nodes
!====================================指定晶粒材料=======================================!

    !-----For Different rhoGe Cases --------!
    WRITE(*,*) '4.ASSIGNNING MATERIAL OF GRAINS......'
    
    ALLOCATE( GrainMt(Ngrains) )
    OPEN(UNIT=999,FILE=GrainMtFile)
    READ(999,*) GrainMt
    CLOSE(999)
    
    WRITE(*,*) "Ngrains=",Ngrains
!====================FIND ELEMENTS NEIGHBOR AND BOUNDARY ELEMENTS=======================!
    
    WRITE(*,*) '5.CALCULATING OTHER INFORMATIONS......'
    CALL find_elem_nieghb_bound_elem
    
!============MAKE LOCAL NODE 6,7,8 LOCATED ON BOUNDARY FOR BCL/R ELEMENT===============!

    ALLOCATE( tempMatrix(Propelem,NBCL) )
    tempMatrix=Element(:,BCelementL(:))
    CALL make_face2_on_bound(NBCL,tempMatrix,BCx)
    Element(:,BCelementL(:))=tempMatrix
    DEALLOCATE( tempMatrix )
    
    ALLOCATE( tempMatrix(Propelem,NBCR) )
    tempMatrix=Element(:,BCelementR(:))
    CALL make_face2_on_bound(NBCR,tempMatrix,BCx)
    Element(:,BCelementR(:))=tempMatrix
    DEALLOCATE( tempMatrix )
    
!=======================CALCULATE THE VECTORS NEEDED IN SOLVER==========================!

    ALLOCATE( Volume(Nelement),center(3,Nelement),normVec(3,4,Nelement),vecCoords(3,3,Nelement),invers_vec(3,3,Nelement) )

    CALL compute_vec_vol_centor
    
    WRITE(*,*) '  THE TOTAL VOLUME IS:',SUM(Volume)
    WRITE(*,*) '  THE MAX VOLUME IS:',MAXVAL(Volume),MAXLOC(Volume)
    WRITE(*,*) '  THE MIN VOLUME IS:',MINVAL(Volume),MINLOC(Volume)
    WRITE(*,*) '  THE MAX/MIN VOLUME RATIO IS:',MAXVAL(Volume)/MINVAL(Volume)
    
    PAUSE 'PRESS ANY KEY TO CONTINUE...'
!========================CALCULATE THE AREA OF NBCL/R ELEMENT===========================!
    ALLOCATE( areaBCL(NBCL),areaBCR(NBCR) )
    CALL compute_BCelem_area
!======================================OUT PUT==========================================!
    WRITE(*,*) '5.OUTPUTTING......'
    CALL output(outputfile)
    
END PROGRAM main
!******************************************************************************
!******************************************************************************
SUBROUTINE delete_wrong_BC_nodes
USE variables
IMPLICIT NONE
INTEGER*4::i,N
INTEGER*4,ALLOCATABLE::tempI(:)
    
    N=0
    ALLOCATE( tempI(NBCLnode) )
    DO i=1,NBCLnode
        IF ( DABS(xyNodes(1,BCLnode(i))-dLdomain(1,1)).gt.0.0000001 ) THEN
            N=N+1
            tempI(N)=BCLnode(i)
        ENDIF
    ENDDO
    DEALLOCATE( BCLnode )
    NBCLnode=N
    ALLOCATE( BCLnode(N) )
    BCLnode=tempI(1:N)
    DEALLOCATE( tempI )
    
    N=0
    ALLOCATE( tempI(NBCRnode) )
    DO i=1,NBCRnode
        IF ( DABS(xyNodes(1,BCRnode(i))-dLdomain(2,1)).gt.0.0000001 ) THEN
            N=N+1
            tempI(N)=BCRnode(i)
        ENDIF
    ENDDO
    DEALLOCATE( BCRnode )
    NBCRnode=N
    ALLOCATE( BCRnode(N) )
    BCRnode=tempI(1:N)
    DEALLOCATE( tempI )
    
    N=0
    ALLOCATE( tempI(NBCYNnode) )
    DO i=1,NBCYNnode
        IF ( DABS(xyNodes(2,BCYNnode(i))-dLdomain(1,2)).gt.0.0000001 ) THEN
            N=N+1
            tempI(N)=BCYNnode(i)
        ENDIF
    ENDDO
    DEALLOCATE( BCYNnode )
    NBCYNnode=N
    ALLOCATE( BCYNnode(N) )
    BCYNnode=tempI(1:N)
    DEALLOCATE( tempI )
    
    N=0
    ALLOCATE( tempI(NBCYPnode) )
    DO i=1,NBCYPnode
        IF ( DABS(xyNodes(2,BCYPnode(i))-dLdomain(2,2)).gt.0.0000001 ) THEN
            N=N+1
            tempI(N)=BCYPnode(i)
        ENDIF
    ENDDO
    DEALLOCATE( BCYPnode )
    NBCYPnode=N
    ALLOCATE( BCYPnode(N) )
    BCYPnode=tempI(1:N)
    DEALLOCATE( tempI )
    
    N=0
    ALLOCATE( tempI(NBCZNnode) )
    DO i=1,NBCZNnode
        IF ( DABS(xyNodes(3,BCZNnode(i))-dLdomain(1,3)).gt.0.0000001 ) THEN
            N=N+1
            tempI(N)=BCZNnode(i)
        ENDIF
    ENDDO
    DEALLOCATE( BCZNnode )
    NBCZNnode=N
    ALLOCATE( BCZNnode(N) )
    BCZNnode=tempI(1:N)
    DEALLOCATE( tempI )
    
    N=0
    ALLOCATE( tempI(NBCZPnode) )
    DO i=1,NBCZPnode
        IF ( DABS(xyNodes(3,BCZPnode(i))-dLdomain(2,3)).gt.0.0000001 ) THEN
            N=N+1
            tempI(N)=BCZPnode(i)
        ENDIF
    ENDDO
    DEALLOCATE( BCZPnode )
    NBCZPnode=N
    ALLOCATE( BCZPnode(N) )
    BCZPnode=tempI(1:N)
    DEALLOCATE( tempI )

END SUBROUTINE delete_wrong_BC_nodes
!******************************************************************************
!******************************************************************************
SUBROUTINE read_nodes(LR1)
USE variables
IMPLICIT NONE
INTEGER*4::LR1,TRUE,i,j,k,readstat
REAL*8::unit_T
REAL*8,ALLOCATABLE::ReadTemp1(:,:),ReadTemp2(:,:)

        READ(LR1,*)
        READ(LR1,*)

        TRUE=0
        i=1
        j=500
        k=1
        readstat=0
        ALLOCATE( ReadTemp1(4,k*j) )
    
        DO WHILE(TRUE.ne.1)
            READ(unit=LR1,FMT=*,iostat=readstat) ReadTemp1(:,i)
            IF (readstat.le.0) THEN
                i=i+1
                IF (i.gt.k*j) THEN
                    ALLOCATE( ReadTemp2(4,k*j) )
                    ReadTemp2=ReadTemp1
                    DEALLOCATE( ReadTemp1 )
                    k=k+1
                    ALLOCATE( ReadTemp1(4,k*j) )
                    ReadTemp1(:,1:(k-1)*j)=ReadTemp2
                    DEALLOCATE( ReadTemp2 )
                ENDIF
            ELSE
                TRUE=1
            ENDIF
        ENDDO
            
        Nnodes=MAXVAL( ReadTemp1(1,:) )
        ALLOCATE( xyNodes(3,Nnodes) )
        xyNodes=ReadTemp1(2:4,1:Nnodes)
        dLdomain(:,1)=(/ MINVAL(xyNodes(1,:)),MAXVAL(xyNodes(1,:)) /)
        dLdomain(:,2)=(/ MINVAL(xyNodes(2,:)),MAXVAL(xyNodes(2,:)) /)
        dLdomain(:,3)=(/ MINVAL(xyNodes(3,:)),MAXVAL(xyNodes(3,:)) /)
        WRITE(*,*) '  THE DOMAIN DIMENSION IS:'
        WRITE(*,*) dLdomain(2,1),'*',dLdomain(2,2),'*',dLdomain(2,3)
        WRITE(*,*) '  PLEASE INPUT UNIT TRANSFORM CONSTANT:(IF THE UNTI IS WRITE,INPUT 1)'
        READ(*,*) unit_T
        xyNodes=xyNodes*unit_T
        dLdomain(:,1)=(/ MINVAL(xyNodes(1,:)),MAXVAL(xyNodes(1,:)) /)
        dLdomain(:,2)=(/ MINVAL(xyNodes(2,:)),MAXVAL(xyNodes(2,:)) /)
        dLdomain(:,3)=(/ MINVAL(xyNodes(3,:)),MAXVAL(xyNodes(3,:)) /)
        WRITE(*,*) '  NOW THE DOMAIN DIMENSION IS:'
        WRITE(*,*) dLdomain(2,1),'*',dLdomain(2,2),'*',dLdomain(2,3)
        DEALLOCATE( ReadTemp1 )
        WRITE(*,*) '  Nnodes=',Nnodes
        WRITE(*,*) '  END READING INFORMATION OF NODES......'
        
        
        DO i=1,Nnodes
            DO j=i+1,Nnodes
                IF (ALL(DABS(xynodes(:,j)-xynodes(:,i)).le.0.0000001)) THEN
                    PAUSE 'yyyy'
                ENDIF
            ENDDO
        ENDDO
        
END SUBROUTINE read_nodes

!*******************************************************************************
!*******************************************************************************

SUBROUTINE read_elements(LR1)
USE variables
IMPLICIT NONE
INTEGER*4::LR1,TRUE,i,j,k,Ntemp
REAL*8,ALLOCATABLE::ReadTemp1(:,:),ReadTemp2(:,:)
CHARACTER*72::temp

    TRUE=0
    k=0
    
    DO WHILE (TRUE.eq.0)
        CALL read_tmp(LR1,temp)
        IF (temp(1:14).eq.'/com, Elements') THEN
            k=k+1
            CALL read_tmp(LR1,temp)
            i=INDEX(temp,',',.TRUE.) !尋找這行最後一個逗號的位置，最後一個逗號後面有element的數量
            READ(temp(i+1:LEN_TRIM(temp)),*) Ntemp
            Nelement=Ntemp-1
            
            IF ( ALLOCATED(Element) ) THEN
                Ntemp=Ntemp-SIZE(Element,2)-1
            ELSE
                Ntemp=Ntemp-1
            ENDIF
            
            WRITE(*,*) '  THERE ARE',Ntemp,'ELEMENTS IN BLOCK',k
            ALLOCATE( ReadTemp1(15,Ntemp) )
            CALL read_tmp(LR1,temp)
            
            DO i=1,Ntemp
                READ(LR1,*) ReadTemp1(:,i)
            ENDDO
            
            IF ( ALLOCATED(Element) ) THEN
                ALLOCATE( ReadTemp2(Propelem,Nelement-Ntemp) )
                ReadTemp2=Element
                DEALLOCATE( Element )
                ALLOCATE( Element(Propelem,Nelement) )
                Element(:,1:Nelement-Ntemp)=ReadTemp2
                DEALLOCATE( ReadTemp2 )
                Element(6:9,Nelement-Ntemp+1:Nelement)=ReadTemp1(12:15,:)
                Element(1,Nelement-Ntemp+1:Nelement)=k
            ELSE
                ALLOCATE( Element(Propelem,Nelement) )
                Element(1,:)=k
                Element(6:9,:)=ReadTemp1(12:15,:)
            ENDIF
            
            DEALLOCATE( ReadTemp1 )
            CALL read_tmp(LR1,temp)
            CALL read_tmp(LR1,temp)
        ELSE
            TRUE=1
            BACKSPACE(LR1)
        ENDIF
    ENDDO
    Ngrains=k
    WRITE(*,*) '  TOTAL NUMBER OF BLOCKS(GRAINS)=',k
    WRITE(*,*) '  TOTAL NUMBER OF ELEMENTS=',Nelement
    WRITE(*,*) '  END READING INFORMATION OF ELEMENTS......'
    
END SUBROUTINE read_elements

!*******************************************************************************
!******************************************************************************

SUBROUTINE read_boundary_nodes(LR1)
USE variables
IMPLICIT NONE
INTEGER*4::LR1,i
CHARACTER*72::temp,BCname
INTERFACE
    SUBROUTINE read_NBCnode_BCname(LR1,temp,NBCnode,BCname,BCnode)
        IMPLICIT NONE
        INTEGER*4::i,LR1,NBCnode
        INTEGER*4,ALLOCATABLE::BCnode(:)
        CHARACTER*72::temp,BCname
    END SUBROUTINE read_NBCnode_BCname
END INTERFACE
        
    READ(LR1,*)
    DO
        CALL read_tmp(LR1,temp)
        
        IF (temp(1:7).eq.'CMBLOCK') THEN
        
            i=INDEX(temp,',')
            READ(temp(i+1:),*) BCname
            
            IF (BCname.eq.'BCL') THEN
                CALL read_NBCnode_BCname(LR1,temp,NBCLnode,BCname,BCLnode)
            ELSE IF (BCname.eq.'BCR') THEN
                CALL read_NBCnode_BCname(LR1,temp,NBCRnode,BCname,BCRnode)
            ELSE IF (BCname.eq.'BCYP') THEN
                CALL read_NBCnode_BCname(LR1,temp,NBCYPnode,BCname,BCYPnode)
            ELSE IF (BCname.eq.'BCYN') THEN
                CALL read_NBCnode_BCname(LR1,temp,NBCYNnode,BCname,BCYNnode)
            ELSE IF (BCname.eq.'BCZP') THEN
                CALL read_NBCnode_BCname(LR1,temp,NBCZPnode,BCname,BCZPnode)
            ELSE IF (BCname.eq.'BCZN') THEN
                CALL read_NBCnode_BCname(LR1,temp,NBCZNnode,BCname,BCZNnode)
            ELSE
                WRITE(*,*) 'WRONG!!'
                WRITE(*,*) 'BOUNDARY NAME ',BCname,' is not acceptable for MC3D solver.'
                WRITE(*,*) 'PLEASE CHANGE NAME AND TRY AGAIN.'
                PAUSE 'PRESS ANY KEY TO EXIT THE PROGRAM'
                STOP
            ENDIF  
        ELSE IF (temp(1:).eq.'!') THEN                 
            EXIT
        ENDIF
    ENDDO
    WRITE(*,*) '  END READING BOUNDARY NODES......'
    
END SUBROUTINE read_boundary_nodes

!----------------------------------------------------------------
SUBROUTINE read_NBCnode_BCname(LR1,temp,NBCnode,BCname,BCnode)
IMPLICIT NONE
INTEGER*4::i,LR1,NBCnode
INTEGER*4,ALLOCATABLE::BCnode(:)
CHARACTER*72::temp,BCname

    i=INDEX(temp,',',.TRUE.)
    READ(temp(i+1:),*) NBCnode
    WRITE(*,"('   THERE ARE ',I5,' NODES ON THE BOUNDARY ',A4,'.')") NBCnode,BCname
    ALLOCATE( BCnode(NBCnode) )
    READ(LR1,*)
    READ(LR1,*) BCnode

END SUBROUTINE read_NBCnode_BCname
!----------------------------------------------------------------

!******************************************************************************
!******************************************************************************

SUBROUTINE find_elem_nieghb_bound_elem
USE variables
IMPLICIT NONE
INTEGER*4::i,j,TRUE,m(4)
    
    NBCL=0
    NBCR=0
    NBCYP=0
    NBCYN=0
    NBCZP=0
    NBCZN=0
    Element(2:5,:)=Nelement*2
    
    DO i=51,56
        OPEN(UNIT=i,STATUS='SCRATCH')
    ENDDO
    
    DO i=1,Nelement
        m=0
        DO j=1,Nelement
            !IF (j.ne.i) CALL find_neighbor(Element(1:5,i),xyNodes(:,Element(6:9,i)),xyNodes(:,Element(6:9,j)),j,m,TRUE)
            IF (j.ne.i) CALL find_neighbor(Element(:,i),4,Element(6:9,j),j,m,TRUE)
        ENDDO
        !WRITE(*,*) '1',Element(2:5,i)
        
        IF ( ANY(Element(2:5,i).eq.2*Nelement) ) THEN
            CALL search_bound_elem(i,Element(:,i),NBCLnode,BCLnode,BCx,m,TRUE,NBCL,51)
            !WRITE(*,*) '2',Element(2:5,i)
            CALL search_bound_elem(i,Element(:,i),NBCRnode,BCRnode,BCx,m,TRUE,NBCR,52)
            !WRITE(*,*) '3',Element(2:5,i)
            CALL search_bound_elem(i,Element(:,i),NBCYPnode,BCYPnode,BCy,m,TRUE,NBCYP,53)
            !WRITE(*,*) '4',Element(2:5,i)
            CALL search_bound_elem(i,Element(:,i),NBCYNnode,BCYNnode,BCy,m,TRUE,NBCYN,54)
            !WRITE(*,*) '5',Element(2:5,i)
            CALL search_bound_elem(i,Element(:,i),NBCZPnode,BCZPnode,BCz,m,TRUE,NBCZP,55)
            !WRITE(*,*) '6',Element(2:5,i)
            CALL search_bound_elem(i,Element(:,i),NBCZNnode,BCZNnode,BCz,m,TRUE,NBCZN,56)
            !WRITE(*,*) '7',Element(2:5,i)
        ENDIF
        
        !IF ( ANY(Element(2:5,i).eq.2*Nelement) ) THEN
        !    CALL search_bound_elem2(i,Element(:,i),xyNodes(:,Element(6:9,i)),1,dLdomain(1,1),BCx,m,TRUE,NBCL,51)
        !    !WRITE(*,*) '2',Element(2:5,i)
        !    CALL search_bound_elem2(i,Element(:,i),xyNodes(:,Element(6:9,i)),1,dLdomain(2,1),BCx,m,TRUE,NBCR,52)
        !    !WRITE(*,*) '3',Element(2:5,i)
        !    CALL search_bound_elem2(i,Element(:,i),xyNodes(:,Element(6:9,i)),2,dLdomain(2,2),BCy,m,TRUE,NBCYP,53)
        !    !WRITE(*,*) '4',Element(2:5,i)
        !    CALL search_bound_elem2(i,Element(:,i),xyNodes(:,Element(6:9,i)),2,dLdomain(1,2),BCy,m,TRUE,NBCYN,54)
        !    !WRITE(*,*) '5',Element(2:5,i)
        !    CALL search_bound_elem2(i,Element(:,i),xyNodes(:,Element(6:9,i)),3,dLdomain(2,3),BCz,m,TRUE,NBCZP,55)
        !    !WRITE(*,*) '6',Element(2:5,i)
        !    CALL search_bound_elem2(i,Element(:,i),xyNodes(:,Element(6:9,i)),3,dLdomain(1,3),BCz,m,TRUE,NBCZN,56)
        !    !WRITE(*,*) '7',Element(2:5,i)
        !ENDIF
        
        IF (ANY(m.ne.1)) THEN
            WRITE(*,*) i,m
            WRITE(*,*) xyNodes(:,Element(6,i))
            WRITE(*,*) xyNodes(:,Element(7,i))
            WRITE(*,*) xyNodes(:,Element(8,i))
            WRITE(*,*) xyNodes(:,Element(9,i))
            
            WRITE(*,*) 'THE',i,'-th ELEMENT HAS SOME FACES THAT HAVE MORE THAN ONE'
            WRITE(*,*) 'CONSTRAINT RELATION WITH OTHER ELEMENTS OR BOUNDARYS.'
            PAUSE 'PRESS ANY KEY TO EXIT THE PROGRAM'
            STOP
        ENDIF
        IF (MOD(i,INT(DBLE(Nelement)*0.1)).eq.0) WRITE(*,*) i
    ENDDO

    ALLOCATE( BCelementL(NBCL),BCelementR(NBCR),BCelementyP(NBCYP),BCelementyN(NBCYN),BCelementzP(NBCZP),BCelementzN(NBCZN) )
    
    CALL read_BCelement(51,NBCL,BCelementL)
    CALL read_BCelement(52,NBCR,BCelementR)
    CALL read_BCelement(53,NBCYP,BCelementyP)
    CALL read_BCelement(54,NBCYN,BCelementyN)
    CALL read_BCelement(55,NBCZP,BCelementzP)
    CALL read_BCelement(56,NBCZN,BCelementzN)

END SUBROUTINE find_elem_nieghb_bound_elem
!--------------------------------------------
SUBROUTINE search_bound_elem2(i,Elem,ElemNodes,x,dL,BC,check,TRUE,NBC,writenum)
IMPLICIT NONE
INTEGER*4::i,Elem(9),BC,check(4),TRUE,NBC,writenum,x,k
REAL*8::dL,ElemNodes(3,4)
LOGICAL::TrueNode(4)
    
    TrueNode=.false.
    TRUE=0
    DO k=1,4
        IF (DABS(ElemNodes(x,k)-dL).lt.0.0000001) TrueNode(k)=.true.
    ENDDO
    
    IF (TrueNode(1).AND.TrueNode(2).AND.TrueNode(3).AND.(.not.TrueNode(4))) THEN
        Elem(2)=BC
        check(1)=check(1)+1
        TRUE=1
    ELSE IF (TrueNode(1).AND.TrueNode(2).AND.TrueNode(4).AND.(.not.TrueNode(3))) THEN
        Elem(3)=BC
        check(2)=check(2)+1
        TRUE=1
    ELSE IF (TrueNode(1).AND.TrueNode(3).AND.TrueNode(4).AND.(.not.TrueNode(2))) THEN
        Elem(4)=BC
        check(3)=check(3)+1
        TRUE=1
    ELSE IF (TrueNode(2).AND.TrueNode(3).AND.TrueNode(4).AND.(.not.TrueNode(1))) THEN
        Elem(5)=BC
        check(4)=check(4)+1
        TRUE=1
    ELSE IF ( ALL(TrueNode) ) THEN
        PAUSE 'WRONG IN search_bound_elem'
    ELSE
        TRUE=0
    ENDIF
    
    IF (TRUE.eq.1) THEN
        NBC=NBC+1
        WRITE(writenum,*) i
    ENDIF

END SUBROUTINE search_bound_elem2
!--------------------------------------------
!--------------------------------------------
SUBROUTINE search_bound_elem(i,Elem,NBCnode,BCnode,BC,check,TRUE,NBC,writenum)
IMPLICIT NONE
INTEGER*4::i,Elem(9),NBCnode,BCnode(NBCnode),BC,check(4),TRUE,NBC,writenum
    
    CALL find_neighbor(Elem,NBCnode,BCnode,BC,check,TRUE)
    IF (TRUE.eq.1) THEN
        NBC=NBC+1
        WRITE(writenum,*) i
    ENDIF

END SUBROUTINE search_bound_elem
!--------------------------------------------
SUBROUTINE read_BCelement(readnum,NBC,BCelement)
IMPLICIT NONE
INTEGER*4::i,NBC,readnum,BCelement(NBC)

    REWIND(readnum)
    DO i=1,NBC
        READ(readnum,*) BCelement(i)
    ENDDO
    CLOSE(readnum)

END SUBROUTINE read_BCelement
!--------------------------------------------

!******************************************************************************
!******************************************************************************

SUBROUTINE make_face2_on_bound(NBC,BCelem,BC)
IMPLICIT NONE
INTEGER*4::i,j,k,NBC,BCelem(9,NBC),BC
    
    DO i=1,NBC
        DO j=2,5
            IF (BCelem(j,i).eq.BC) EXIT
        ENDDO   
        IF (j.eq.3) THEN
            k=BCelem(9,i)
            BCelem(9,i)=BCelem(8,i)
            BCelem(8,i)=k
            k=BCelem(2,i)
            BCelem(2,i)=BCelem(3,i)
            BCelem(3,i)=k
        ELSE IF (j.eq.4) THEN
            k=BCelem(9,i)
            BCelem(9,i)=BCelem(7,i)
            BCelem(7,i)=k
            k=BCelem(2,i)
            BCelem(2,i)=BCelem(4,i)
            BCelem(4,i)=k
        ELSE IF (j.eq.5) THEN
            IF (BCelem(5,i).ne.BC) PAUSE 'WRONG 3'
            k=BCelem(9,i)
            BCelem(9,i)=BCelem(6,i)
            BCelem(6,i)=k
            k=BCelem(2,i)
            BCelem(2,i)=BCelem(5,i)
            BCelem(5,i)=k
        ENDIF
    ENDDO
    

END SUBROUTINE make_face2_on_bound

!******************************************************************************
!******************************************************************************

SUBROUTINE compute_vec_vol_centor
USE IMSL_LIBRARIES
USE variables
IMPLICIT NONE
INTEGER*4::i,j
REAL*8::r78(3),r79(3),temp(3,3)

    DO i=1,Nelement
        vecCoords(:,1,i)=xyNodes(:,Element(7,i))-xyNodes(:,Element(6,i))
        vecCoords(:,2,i)=xyNodes(:,Element(8,i))-xyNodes(:,Element(6,i))
        vecCoords(:,3,i)=xyNodes(:,Element(9,i))-xyNodes(:,Element(6,i))
        r78(:)=xyNodes(:,Element(8,i))-xyNodes(:,Element(7,i))
        r79(:)=xyNodes(:,Element(9,i))-xyNodes(:,Element(7,i))
        normVec(:,1,i)=CROSS_PRODUCT(vecCoords(:,1,i),vecCoords(:,2,i))
        normVec(:,2,i)=CROSS_PRODUCT(vecCoords(:,1,i),vecCoords(:,3,i))
        normVec(:,3,i)=CROSS_PRODUCT(vecCoords(:,2,i),vecCoords(:,3,i))
        normVec(:,4,i)=CROSS_PRODUCT(r78,r79)
        DO j=1,4
            normVec(:,j,i)=normVec(:,j,i)/DSQRT(normVec(1,j,i)**2+normVec(2,j,i)**2+normVec(3,j,i)**2)
        ENDDO
        IF (DOT_PRODUCT(normVec(:,1,i),r79).ge.0) normVec(:,1,i)=-normVec(:,1,i)
        IF (DOT_PRODUCT(normVec(:,2,i),r78).ge.0) normVec(:,2,i)=-normVec(:,2,i)
        IF (DOT_PRODUCT(normVec(:,3,i),r79).le.0) normVec(:,3,i)=-normVec(:,3,i)
        IF (DOT_PRODUCT(normVec(:,4,i),vecCoords(:,1,i)).le.0) normVec(:,4,i)=-normVec(:,4,i)
        
        center(:,i)=(xyNodes(:,Element(6,i))+xyNodes(:,Element(7,i))+xyNodes(:,Element(8,i))+xyNodes(:,Element(9,i)))/4
        Volume(i)=DABS(DOT_PRODUCT(vecCoords(:,1,i),CROSS_PRODUCT(vecCoords(:,2,i),vecCoords(:,3,i))))/6
        
        temp(:,1)=vecCoords(:,1,i)
        temp(:,2)=vecCoords(:,2,i)
        temp(:,3)=vecCoords(:,3,i)
        
        CALL D_LINRG(temp,invers_vec(:,:,i))
    ENDDO

END SUBROUTINE compute_vec_vol_centor


!******************************************************************************
!******************************************************************************

SUBROUTINE compute_BCelem_area
USE variables
IMPLICIT NONE
INTEGER*4::i
REAL*8::temp(3)

    DO i=1,NBCL
        IF (Element(2,BCelementL(i)).ne.BCx) THEN
            PAUSE 'WRONG IN CALCULATING BCL AREAS!!!'
            STOP
        ENDIF
        temp=CROSS_PRODUCT(vecCoords(:,1,BCelementL(i)),vecCoords(:,2,BCelementL(i)))
        areaBCL(i)=DSQRT(temp(1)**2+temp(2)**2+temp(3)**2)/2
    ENDDO
    
    DO i=1,NBCR
        IF (Element(2,BCelementR(i)).ne.BCx) THEN
            PAUSE 'WRONG IN CALCULATING BCR AREAS!!!'
            STOP
        ENDIF
        temp=CROSS_PRODUCT(vecCoords(:,1,BCelementR(i)),vecCoords(:,2,BCelementR(i)))
        areaBCR(i)=DSQRT(temp(1)**2+temp(2)**2+temp(3)**2)/2
    ENDDO
    
END SUBROUTINE compute_BCelem_area

!******************************************************************************
!******************************************************************************

SUBROUTINE output(outputfile)
USE variables
IMPLICIT NONE
INTEGER*4::i,j,k,LW1=110
CHARACTER*72::outputfile

    OPEN(UNIT=LW1,FILE=outputfile)
    WRITE(LW1,*) Propelem,Nelement,Nnodes,Ngrains
    WRITE(LW1,*) GrainMt
    WRITE(LW1,*) Element
    WRITE(LW1,*) xyNodes
    WRITE(LW1,*) center
    WRITE(LW1,*) Volume
    WRITE(LW1,*) normVec
    WRITE(LW1,*) vecCoords
    WRITE(LW1,*) invers_vec
    WRITE(LW1,*) NBCL,NBCR,NBCyP,NBCyN,NBCzP,NBCzN
    WRITE(LW1,*) BCelementL,BCelementR
    WRITE(LW1,*) BCelementyP,BCelementyN
    WRITE(LW1,*) BCelementzP,BCelementzN
    WRITE(LW1,*) areaBCL,areaBCR
    CLOSE(LW1)
    
    OPEN(UNIT=200,FILE='center.txt')
    DO i=1,Nelement
        WRITE(200,*) center(:,i)
    ENDDO
    CLOSE(200)
    
    OPEN(UNIT=200,FILE='neighber label.txt')
    DO i=1,Nelement
        WRITE(200,"(4(1X,I6))") Element(2:5,i)
    ENDDO
    CLOSE(200)
    
    OPEN(UNIT=200,FILE='Tet_Nodes_mt.m')
    DO i=1,Nelement
        DO j=1,3
            WRITE(200,*) 'Tet(1,1,',j,',',i,')=',xyNodes(j,Element(6,i)),';'
            WRITE(200,*) 'Tet(1,2,',j,',',i,')=',xyNodes(j,Element(6,i)),';'
            WRITE(200,*) 'Tet(1,3,',j,',',i,')=',xyNodes(j,Element(6,i)),';'
            WRITE(200,*) 'Tet(1,4,',j,',',i,')=',xyNodes(j,Element(7,i)),';'
            
            WRITE(200,*) 'Tet(2,1,',j,',',i,')=',xyNodes(j,Element(7,i)),';'
            WRITE(200,*) 'Tet(2,2,',j,',',i,')=',xyNodes(j,Element(7,i)),';'
            WRITE(200,*) 'Tet(2,3,',j,',',i,')=',xyNodes(j,Element(8,i)),';'
            WRITE(200,*) 'Tet(2,4,',j,',',i,')=',xyNodes(j,Element(8,i)),';'
            
            WRITE(200,*) 'Tet(3,1,',j,',',i,')=',xyNodes(j,Element(8,i)),';'
            WRITE(200,*) 'Tet(3,2,',j,',',i,')=',xyNodes(j,Element(9,i)),';'
            WRITE(200,*) 'Tet(3,3,',j,',',i,')=',xyNodes(j,Element(9,i)),';'
            WRITE(200,*) 'Tet(3,4,',j,',',i,')=',xyNodes(j,Element(9,i)),';'
        ENDDO
        WRITE(200,*) 'mt(',i,')=',GrainMt(Element(1,i)),';'
        WRITE(200,*) 'grain(',i,')=',Element(1,i),';'
    ENDDO
    CLOSE(200)
    
    OPEN(UNIT=200,FILE='Tet_Draw_Matlab.m')
    WRITE(200,*) 'for i=1:1:',Nelement,','
    WRITE(200,*) 'fill3(Tet(:,:,1,i),Tet(:,:,2,i),Tet(:,:,3,i),mt(i));'
    WRITE(200,*) 'hold on;'
    WRITE(200,*) 'end'
    WRITE(200,*) 'axis equal;'
    WRITE(200,"('axis ([0,',F9.4,',0,',F9.4,',0,',F9.4,']);')") MAXVAL(xyNodes(1,:)),MAXVAL(xyNodes(2,:)),MAXVAL(xyNodes(3,:))
    WRITE(*,"('axis ([0,',F9.4,',0,',F9.4,',0,',F9.4,']);')") MAXVAL(xyNodes(1,:)),MAXVAL(xyNodes(2,:)),MAXVAL(xyNodes(3,:))
    WRITE(200,*) 'axis vis3d;'
    CLOSE(200)
    
    OPEN(UNIT=200,FILE='Tet_norm__vecCoords_inv.m')
    DO i=1,Nelement
        DO j=1,4
            DO k=1,3
                WRITE(200,*) 'norm(',k,',',j,',',i,')=',normVec(k,j,i),';'
            ENDDO
        ENDDO
        DO j=1,3
            DO k=1,3
                WRITE(200,*) 'vecCoord(',k,',',j,',',i,')=',vecCoords(k,j,i),';'
                WRITE(200,*) 'invVec(',k,',',j,',',i,')=',invers_vec(k,j,i),';'
            ENDDO
        ENDDO
    ENDDO
    CLOSE(200)
END SUBROUTINE output

!*******************************************************************************
!******************************************************************************
!SUBROUTINE find_neighbor(ElemA,ElemANodes,ElemBNodes,neighbor,check,true)
!IMPLICIT NONE
!INTEGER*4::true,neighbor,check(4),ElemA(5),i,j
!REAL*8::ElemANodes(3,4),ElemBNodes(3,4)
!LOGICAL::TrueNode(4)
    
!    true=0
!    TrueNode=.false.
!    DO i=1,4
!        DO j=1,4
!            IF (ALL(DABS(ElemANodes(:,i)-ElemBNodes(:,j)).lt.0.0000001)) THEN
!                IF (.not.TrueNode(i)) THEN
!                    TrueNode(i)=.true.
!                ELSE
!                    WRITE(*,*) i,j,TrueNode(i)
!                    WRITE(*,"(4(3(E15.7),/))") ElemANodes
!                    WRITE(*,*)
!                    WRITE(*,"(4(3(E15.7),/))") ElemBNodes
!                    PAUSE 'WRONG IN find_neighbor'
!                    STOP
!                ENDIF
!            ENDIF
!        ENDDO
!    ENDDO
    
!    IF ( TrueNode(1).AND.TrueNode(2).AND.TrueNode(3) ) THEN
!        ElemA(2)=neighbor
!        check(1)=check(1)+1
!        true=1
!    ELSE IF ( TrueNode(1).AND.TrueNode(2).AND.TrueNode(4) ) THEN
!        ElemA(3)=neighbor
!        check(2)=check(2)+1
!        true=1
!    ELSE IF ( TrueNode(1).AND.TrueNode(3).AND.TrueNode(4) ) THEN
!        ElemA(4)=neighbor
!        check(3)=check(3)+1
!        true=1
!    ELSE IF ( TrueNode(2).AND.TrueNode(3).AND.TrueNode(4) ) THEN
!        ElemA(5)=neighbor
!        check(4)=check(4)+1
!        true=1
!    ELSE
!        true=0
!    ENDIF

!END SUBROUTINE find_neighbor
!*******************************************************************************
!******************************************************************************
SUBROUTINE find_neighbor(A,n,B,neighbor,check,true)
IMPLICIT NONE
INTEGER*4::n,true
INTEGER*4::A(9),B(n),neighbor,check(4)
    !A矩陣是一個Element的性質，B矩陣是一群node的集合(ex:另一個Element的4個node，或是都在同一個邊界面上的node)，
    !由A這個Element有幾個node是屬於B這群node 來判斷A的四個面是屬於與另一個Element的共面或是在邊界面上
    !而check(1~4)是用來判斷A這個Element的第1~4個面  是否已經知道與其他網格或DOMAIN邊界面的關係
    !若=0則表示仍然未知，=1則表示已知，>1則表示這個面有多重關係，這是不合理的，表示程式有問題
    
    true=0
    
    IF ( ANY(B.eq.A(6)).AND.ANY(B.eq.A(7)).AND.ANY(B.eq.A(8)).AND.ALL(B.ne.A(9)) ) THEN
        A(2)=neighbor
        check(1)=check(1)+1
        true=1
    ELSE IF ( ANY(B.eq.A(6)).AND.ANY(B.eq.A(7)).AND.ANY(B.eq.A(9)).AND.ALL(B.ne.A(8)) ) THEN
        A(3)=neighbor
        check(2)=check(2)+1
        true=1
    ELSE IF ( ANY(B.eq.A(6)).AND.ANY(B.eq.A(8)).AND.ANY(B.eq.A(9)).AND.ALL(B.ne.A(7)) ) THEN
        A(4)=neighbor
        check(3)=check(3)+1
        true=1
    ELSE IF ( ANY(B.eq.A(7)).AND.ANY(B.eq.A(8)).AND.ANY(B.eq.A(9)).AND.ALL(B.ne.A(6)) ) THEN
        A(5)=neighbor
        check(4)=check(4)+1
        true=1
    ELSE
        true=0
    ENDIF
END SUBROUTINE find_neighbor

!*******************************************************************************
!******************************************************************************
