!**************************************************************************************************************
!**************************************************************************************************************
!**************************************************************************************************************
!**************************************************************************************************************
!**************************************************************************************************************
PROGRAM main
INCLUDE 'link_fnl_shared_hpc.h'
USE omp_lib
IMPLICIT NONE
INTEGER*4::i,j,k,m,s,p,et(0:3),ct(0:3)
INTEGER*4::Nelement,N(3),CPUID,bg,ed
REAL*8::DiagVec(3),TetNodes(3,4),CubStartNode(3),dV,tempV
INTEGER*4,ALLOCATABLE::Element(:,:),ElementBound(:,:,:),NonZeroV(:,:,:)
REAL*8,ALLOCATABLE::xyNodes(:,:),volume(:),invers_vec(:,:,:),Vijk(:,:,:),Vm(:)
CHARACTER*72::casename
LOGICAL::true
INTERFACE
    SUBROUTINE read_gridfile(Nelement,Element,xyNodes,volume,invers_vec,casename)
        IMPLICIT NONE
        INTEGER*4::Propelem,Nelement,Nnodes
        INTEGER*4::NBCL,NBCR,NBCyP,NBCyN,NBCzP,NBCzN
        INTEGER*4,ALLOCATABLE::Element(:,:),BCelementL(:),BCelementR(:),BCelementyP(:),BCelementyN(:),BCelementzP(:),BCelementzN(:)
        REAL*8,ALLOCATABLE::xyNodes(:,:),center(:,:),Volume(:)
        REAL*8,ALLOCATABLE::normVec(:,:,:),vecCoords(:,:,:),invers_vec(:,:,:)
        CHARACTER*72::casename
    END SUBROUTINE read_gridfile
END INTERFACE
!------------------------------------------------
    WRITE(*,*) 'ENTER THE PROJECT NAME(no ''.txt'',no ''_gridfile''):'
    READ(*,"(A72)") casename

    CALL read_gridfile(Nelement,Element,xyNodes,volume,invers_vec,casename)
    
    
    WRITE(*,*) 'THE DIMENSIONS OF THE DOMAIN IS:'
    WRITE(*,*) MAXVAL(xyNodes,2)
    WRITE(*,*) 'ENTER THE NUMBER OF STRUCTURAL ELEMENT ON EACH DIRECTION (Nx,Ny,Nz):'
    READ(*,*) N(1),N(2),N(3)

    ALLOCATE( ElementBound(2,3,Nelement) )
    DiagVec=MAXVAL(xyNodes,2)/DBLE(N)
    dV=DiagVec(1)*DiagVec(2)*DiagVec(3)

    CALL OMP_SET_NUM_THREADS(4)
    !$OMP PARALLEL DEFAULT(SHARED),PRIVATE(i,CubStartNode,j)
    !$OMP DO SCHEDULE(DYNAMIC)
    DO i=1,Nelement
        CubStartNode=MINVAL(xyNodes(:,Element(6:9,i)),2) !temporarily used
        ElementBound(1,:,i)=INT(CubStartNode/DiagVec)
        DO j=1,3
            IF (ElementBound(1,j,i).eq.0) ElementBound(1,j,i)=ElementBound(1,j,i)+1
        ENDDO
        CubStartNode=MAXVAL(xyNodes(:,Element(6:9,i)),2) !temporarily used
        ElementBound(2,:,i)=INT(CubStartNode/DiagVec)
        DO j=1,3
            IF (ElementBound(2,j,i).ne.N(j)) ElementBound(2,j,i)=ElementBound(2,j,i)+1
        ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
    
    ct=0
    et=0
    CALL OMP_SET_NUM_THREADS(4)
    !$OMP PARALLEL DEFAULT(SHARED),PRIVATE(i,j,k,m,CubStartNode,TetNodes,CPUID,bg,ed,tempV)
    CPUID=OMP_GET_THREAD_NUM()
    OPEN(UNIT=91+CPUID,STATUS='SCRATCH')
    
    SELECT CASE(CPUID)
        CASE(0)
            bg=1
            ed=Nelement/4
        CASE(1)
            bg=Nelement/4+1
            ed=(Nelement/4)*2
        CASE(2)
            bg=(Nelement/4)*2+1
            ed=(Nelement/4)*3
        CASE(3)
            bg=(Nelement/4)*3+1
            ed=Nelement
    ENDSELECT
    !----------------------------------------------------------------
    DO m=bg,ed
        TetNodes=xyNodes(:,Element(6:9,m))
        DO i=ElementBound(1,1,m),ElementBound(2,1,m)
            CubStartNode(1)=(DBLE(i)-1d0)*DiagVec(1)
            DO j=ElementBound(1,2,m),ElementBound(2,2,m)
                CubStartNode(2)=(DBLE(j)-1d0)*DiagVec(2)
                DO  k=ElementBound(1,3,m),ElementBound(2,3,m)
                    CubStartNode(3)=(DBLE(k)-1d0)*DiagVec(3)
                    CALL cuboid_int_tet_3D(CubStartNode,DiagVec,TetNodes,invers_vec(:,:,m),tempV,CPUID)
                    IF (tempV.ne.0) THEN
                        WRITE(91+CPUID,*) i,j,k,m,tempV
                        ct(CPUID)=ct(CPUID)+1
                    ENDIF
                ENDDO
            ENDDO
        ENDDO        
        et(CPUID)=et(CPUID)+1
        WRITE(*,*) m
    ENDDO
    !$OMP END PARALLEL
    !----------------------------------------------------------------
    WRITE(*,*) et
    
    
    
    ALLOCATE( NonZeroV(N(1),N(2),N(3)),Vijk(N(1),N(2),N(3)),Vm(Nelement) )
    NonZeroV=0
    Vijk=0
    Vm=0
    
    OPEN(UNIT=89,FILE='TempFile.txt')
    DO s=0,3
        INQUIRE(UNIT=91+s,EXIST=true)
        IF (true) THEN
            REWIND(91+s)
            DO p=1,ct(s)
                READ(91+s,*) i,j,k,m,tempV
                WRITE(89,*) i,j,k,m,tempV
                NonZeroV(i,j,k)=NonZeroV(i,j,k)+1
                Vijk(i,j,k)=Vijk(i,j,k)+tempV
                Vm(m)=Vm(m)+tempV
            ENDDO
        ENDIF
        CLOSE(91+s)
    ENDDO
    REWIND(89)
    
    OPEN(UNIT=90,FILE=casename(1:LEN_TRIM(casename))//'_VolumeWeight.txt')
    WRITE(90,*) MAXVAL(xyNodes,2)
    WRITE(90,*) N(1),N(2),N(3),Nelement
    WRITE(90,*) SUM(ct)
    WRITE(90,*) NonZeroV
    
    DO p=1,SUM(ct)
        READ(89,*) i,j,k,m,tempV
        WRITE(90,*) i,j,k,m,tempV/dV
    ENDDO
    CLOSE(90)
    CLOSE(89,STATUS='DELETE')
    
    Vijk=Vijk/dV
    
    OPEN(UNIT=91,FILE='Weight Examination(Regular V).txt')
    OPEN(UNIT=92,FILE='Weight Examination(Tetrahedron V).txt')

    DO k=1,N(3)
        DO j=1,N(2)
            DO i=1,N(1)
                WRITE(91,*) Vijk(i,j,k),DABS(1-Vijk(i,j,k))
            ENDDO
        ENDDO
    ENDDO
    
    DO m=1,Nelement
        WRITE(92,*) Vm(m),DABS(Vm(m)-volume(m))/volume(m)
    ENDDO

    CLOSE(91)
    CLOSE(92)
    
    DEALLOCATE( NonZeroV,Vijk,Vm )
    DEALLOCATE( Element,xyNodes,volume,invers_vec,ElementBound )    

END PROGRAM main
!**************************************************************************************************************
!**************************************************************************************************************
SUBROUTINE read_gridfile(Nelement,Element,xyNodes,volume,invers_vec,casename)
IMPLICIT NONE
INTEGER*4::Propelem,Nelement,Nnodes,Ngrains
INTEGER*4::NBCL,NBCR,NBCyP,NBCyN,NBCzP,NBCzN
INTEGER*4,ALLOCATABLE::Element(:,:),BCelementL(:),BCelementR(:),BCelementyP(:),BCelementyN(:),BCelementzP(:),BCelementzN(:),GrainMt(:)
REAL*8,ALLOCATABLE::xyNodes(:,:),center(:,:),Volume(:)
REAL*8,ALLOCATABLE::normVec(:,:,:),vecCoords(:,:,:),invers_vec(:,:,:)
CHARACTER*72::casename

    OPEN(UNIT=100,FILE=casename(1:LEN_TRIM(casename))//'_gridfile.txt')
    READ(100,*) Propelem,Nelement,Nnodes,Ngrains
    ALLOCATE( Element(Propelem,Nelement),xyNodes(3,Nnodes),center(3,Nelement),Volume(Nelement) )
    ALLOCATE( normVec(3,4,Nelement),vecCoords(3,3,Nelement),invers_vec(3,3,Nelement) )
    ALLOCATE( GrainMt(Ngrains) )
    READ(100,*) GrainMt
    READ(100,*) Element
    READ(100,*) xyNodes
    READ(100,*) center
    READ(100,*) Volume
    READ(100,*) normVec
    READ(100,*) vecCoords
    READ(100,*) invers_vec
    READ(100,*) NBCL,NBCR,NBCyP,NBCyN,NBCzP,NBCzN
    ALLOCATE( BCelementL(NBCL),BCelementR(NBCR),BCelementyP(NBCyP),BCelementyN(NBCyN),BCelementzP(NBCzP),BCelementzN(NBCzN) )
    READ(100,*) BCelementL,BCelementR
    READ(100,*) BCelementyP,BCelementyN
    READ(100,*) BCelementzP,BCelementzN
    READ(100,*) 
    CLOSE(100)
    
    DEALLOCATE( center,normVec,vecCoords )
    DEALLOCATE( BCelementL,BCelementR,BCelementyP,BCelementyN,BCelementzP,BCelementzN )

END SUBROUTINE read_gridfile
!==============================================================================================================
SUBROUTINE cuboid_int_tet_3D(CuboidO,DiagVec,TetNodes,invTetVec,Volume,CPUID)
IMPLICIT NONE
INTEGER*4::i,j,k,NInt,CPUID
REAL*8::CuboidO(3),DiagVec(3),TetNodes(3,4)
REAL*8::invTetVec(3,3),Volume
REAL*8::tempNodes(3,100),CuboidNodes(3,8),tempNodes2(3,2),FaceNorm(3)
REAL*8::FaceNodes1(3,3,4),LineNodes1(3,2,6),FaceNodes2(3,4,6),LineNodes2(3,2,12)
LOGICAL::true,true2,true3
CHARACTER*72::string,file1,file2

    Volume=0
    NInt=0
    CuboidNodes(:,1)=CuboidO                                    !0 0 0
    CuboidNodes(:,2)=CuboidO+(/DiagVec(1),0d0,0d0/)             !1 0 0
    CuboidNodes(:,3)=CuboidO+(/0d0,DiagVec(2),0d0/)             !0 1 0
    CuboidNodes(:,4)=CuboidO+(/DiagVec(1),DiagVec(2),0d0/)      !1 1 0
    CuboidNodes(:,5)=CuboidO+(/0d0,0d0,DiagVec(3)/)             !0 0 1
    CuboidNodes(:,6)=CuboidO+(/DiagVec(1),0d0,DiagVec(3)/)      !1 0 1
    CuboidNodes(:,7)=CuboidO+(/0d0,DiagVec(2),DiagVec(3)/)      !0 1 1
    CuboidNodes(:,8)=CuboidO+DiagVec                            !1 1 1
    
    FaceNodes1(:,:,1)=TetNodes(:,(/1,2,3/))
    FaceNodes1(:,:,2)=TetNodes(:,(/1,2,4/))
    FaceNodes1(:,:,3)=TetNodes(:,(/1,3,4/))
    FaceNodes1(:,:,4)=TetNodes(:,(/2,3,4/))
    
    LineNodes1(:,:,1)=TetNodes(:,(/1,2/))
    LineNodes1(:,:,2)=TetNodes(:,(/1,3/))
    LineNodes1(:,:,3)=TetNodes(:,(/1,4/))
    LineNodes1(:,:,4)=TetNodes(:,(/2,3/))
    LineNodes1(:,:,5)=TetNodes(:,(/2,4/))
    LineNodes1(:,:,6)=TetNodes(:,(/3,4/))
    
    FaceNodes2(:,:,1)=CuboidNodes(:,(/1,2,3,4/))
    FaceNodes2(:,:,2)=CuboidNodes(:,(/5,6,7,8/))
    FaceNodes2(:,:,3)=CuboidNodes(:,(/1,2,5,6/))
    FaceNodes2(:,:,4)=CuboidNodes(:,(/3,4,7,8/))
    FaceNodes2(:,:,5)=CuboidNodes(:,(/1,3,5,7/))
    FaceNodes2(:,:,6)=CuboidNodes(:,(/2,4,6,8/))
    
    LineNodes2(:,:,1)=CuboidNodes(:,(/1,2/))
    LineNodes2(:,:,2)=CuboidNodes(:,(/2,4/))
    LineNodes2(:,:,3)=CuboidNodes(:,(/4,3/))
    LineNodes2(:,:,4)=CuboidNodes(:,(/3,1/))
    LineNodes2(:,:,5)=CuboidNodes(:,(/5,6/))
    LineNodes2(:,:,6)=CuboidNodes(:,(/6,8/))
    LineNodes2(:,:,7)=CuboidNodes(:,(/8,7/))
    LineNodes2(:,:,8)=CuboidNodes(:,(/7,5/))
    LineNodes2(:,:,9)=CuboidNodes(:,(/1,5/))
    LineNodes2(:,:,10)=CuboidNodes(:,(/2,6/))
    LineNodes2(:,:,11)=CuboidNodes(:,(/3,7/))
    LineNodes2(:,:,12)=CuboidNodes(:,(/4,8/))  
 
    DO i=1,4
        CALL vertex_in_cuboid_3D(TetNodes(:,i),CuboidO,DiagVec,true)
        IF (true) THEN
            NInt=Nint+1
            tempNodes(:,NInt)=TetNodes(:,i)
        ENDIF
    ENDDO

    DO i=1,8
        CALL vertex_in_tet_3D(CuboidNodes(:,i),TetNodes(:,1),invTetVec,true)
        IF (true) THEN
            NInt=Nint+1
            tempNodes(:,NInt)=CuboidNodes(:,i)
        ENDIF
    ENDDO

    DO i=1,6
        DO j=1,6
            true=.false.
            CALL sort_node_on_plane_3D(4,FaceNodes2(:,:,j),true,FaceNorm)
            CALL line_int_face_3D(4,FaceNodes2(:,:,j),LineNodes1(:,:,i),k,tempNodes2,true,FaceNorm)
            IF ((k.gt.0).AND.(k.le.2)) THEN
                tempNodes(:,NInt+1:NInt+k)=tempNodes2(:,1:k)
                NInt=NInt+k
            ELSE IF ((k.lt.0).OR.(k.gt.2)) THEN
                PAUSE "WRONG 1 IN SUBROUTINE cuboid_int_tet_3D"
            ENDIF
        ENDDO
    ENDDO

    DO i=1,12
        DO j=1,4
            true=.false.
            CALL sort_node_on_plane_3D(3,FaceNodes1(:,:,j),true,FaceNorm)
            CALL line_int_face_3D(3,FaceNodes1(:,:,j),LineNodes2(:,:,i),k,tempNodes2,true,FaceNorm)
            IF ((k.gt.0).AND.(k.le.2)) THEN
                tempNodes(:,NInt+1:NInt+k)=tempNodes2(:,1:k)
                NInt=NInt+k
            ELSE IF ((k.lt.0).OR.(k.gt.2)) THEN
                PAUSE "WRONG 1 IN SUBROUTINE cuboid_int_tet_3D"
            ENDIF
        ENDDO
    ENDDO
    
    true2=.false.
    IF (NInt.gt.0) THEN
    
        CALL del_reduplicate_nodes_3D(NInt,tempNodes(:,1:NInt))
        IF (NInt.eq.4) THEN
            CALL tet_volume(tempNodes(:,1:4),Volume)
        ELSE IF (NInt.gt.4) THEN
            CALL exam_coplane_3D(NInt,tempNodes(:,1:NInt),true,FaceNorm,true2)
            IF (true.eq..false.) THEN
                WRITE(file1,"('temp',I1,'.txt')") CPUID
                WRITE(file2,"('temp2',I1,'.txt')") CPUID
                OPEN(UNIT=50+CPUID,FILE=file1)
                WRITE(50+CPUID,*) 3
                WRITE(50+CPUID,*) NInt
                DO i=1,NInt
                    WRITE(50+CPUID,*) tempNodes(:,i)
                ENDDO
                WRITE(file1,"('qconvex',I1,'<temp',I1,'.txt')") CPUID,CPUID
                CALL SYSTEM(file1//' QJ Pp FA>'//file2)
                CLOSE(UNIT=50+CPUID,STATUS='DELETE')
                
                INQUIRE(FILE=file2,EXIST=true3)
                DO WHILE(true3.eq..false.)
                    INQUIRE(FILE=file2,EXIST=true3)
                ENDDO
                
                OPEN(UNIT=55+CPUID,FILE=file2)
                true=.true.
                DO WHILE(true)
                    CALL read_tmp(55+CPUID,string)
                    IF (INDEX(string,'volume').ne.0) THEN
                        READ(string(INDEX(string,':')+1:),*) Volume
                        true=.false.
                    ENDIF
                ENDDO
                CLOSE(UNIT=55+CPUID,STATUS='DELETE')
            ENDIF
        ENDIF
    ENDIF

    !WRITE(*,*) NInt
    !WRITE(*,"(<NInt>(3(F15.12,1X),/))") tempNodes(:,1:NInt)

END SUBROUTINE cuboid_int_tet_3D
!=============================================================================!
!=============================================================================!
SUBROUTINE line_int_line_3D(LineNodes1,LineNodes2,IntNode,IntTrue,NORM_TF,FaceNorm)
IMPLICIT NONE
REAL*8::LineNodes1(3,2),LineNodes2(3,2),IntNode(3),FaceNorm(3)
REAL*8::t1(3),t2(3),tempNodes(3,4),A(2,2),r1(3),x(2)
LOGICAL::IntTrue,CO_PLANE_TF,NORM_TF

    tempNodes(:,1:2)=LineNodes1
    tempNodes(:,3:)=LineNodes2
    CALL exam_coplane_3D(4,tempNodes,CO_PLANE_TF,FaceNorm,NORM_TF)
    
    IntTrue=.false.
    IF (CO_PLANE_TF) THEN
    
        t1=LineNodes1(:,2)-LineNodes1(:,1)
        t1=t1/DSQRT(SUM(t1**2))
        t2=LineNodes2(:,2)-LineNodes2(:,1)
        t2=t2/DSQRT(SUM(t2**2))
        
        r1=LineNodes2(:,1)-LineNodes1(:,1)
        
        IF (DABS(t1(1)*t2(2)-t2(1)*t1(2)).lt.0.0000001) RETURN
        !A=inverse of 2*2 matrix ( t1(1) t2(1) ; t1(2) t2(2) )
        A(1,1)=t2(2)
        A(1,2)=-t2(1)
        A(2,1)=-t1(2)
        A(2,2)=t1(1)
        A=A/(t1(1)*t2(2)-t1(2)*t2(1))
        x=MATMUL(A,r1(1:2))
        IntNode=x(1)*t1+LineNodes1(:,1)
        IF (( DOT_PRODUCT( (LineNodes1(:,1)-IntNode),(LineNodes1(:,2)-IntNode) ).le.0.0000001 ).AND. &
            ( DOT_PRODUCT( (LineNodes2(:,1)-IntNode),(LineNodes2(:,2)-IntNode) ).le.0.0000001 )) IntTrue=.true.

    ENDIF

END SUBROUTINE line_int_line_3D
!=============================================================================!
!=============================================================================!
SUBROUTINE line_int_face_3D(N,FaceNodes,LineNodes,NInt,IntNodes,NORM_TF,FaceNorm)
IMPLICIT NONE
INTEGER*4::N,Nint
INTEGER*4::i,j,k
REAL*8::FaceNodes(3,N),LineNodes(3,2)
REAL*8::FaceNorm(3),tempNodes(3,N+2),IntNodes(3,2),A,B
LOGICAL::NORM_TF,CO_PLANE_TF
!.Τ絬琿籔Τキ闽玒
!   1.絬籔常琌Τ絛瞅獶礚┑珿程ㄢユ翴程ぶ箂
!   2.ㄢ翴だㄢ凹程ユ翴程ぶ箂(ユ翴ぃΤキ絛瞅ず)
!   3.ㄢ翴常凹箂ユ翴

    tempNodes(:,1:N)=FaceNodes
    tempNodes(:,N+1:)=LineNodes
    CALL exam_coplane_3D(N+2,tempNodes,CO_PLANE_TF,FaceNorm,NORM_TF)
    NInt=0
    IF (CO_PLANE_TF) THEN !絬籔キ
        j=2
        DO i=1,N
            IF (i.eq.N) THEN
                tempNodes(:,1)=FaceNodes(:,N)
                tempNodes(:,2)=FaceNodes(:,1)
            ELSE
                tempNodes(:,1:2)=FaceNodes(:,i:i+1)
            ENDIF
            CALL line_int_line_3D(tempNodes(:,1:2),LineNodes,IntNodes(:,1),CO_PLANE_TF,NORM_TF,FaceNorm)
            IF (CO_PLANE_TF) THEN
                j=j+1
                IF (j.gt.N+2) PAUSE "WRONG 1 IN SUBROUTINE line_int_face_3D."
                tempNodes(:,j)=IntNodes(:,1)
            ENDIF
        ENDDO
        NInt=j-2
        IF (NInt.gt.0) THEN
            CALL del_reduplicate_nodes_3D(NInt,tempNodes(:,3:j))
            IF (NInt.gt.2) THEN
                WRITE(*,*) NInt
                WRITE(*,"(<NInt>(3(F15.8,1X),/))") tempNodes(:,3:2+NInt)
                PAUSE "WRONG 2 IN SUBROUTINE line_int_face_3D."
                STOP
            ENDIF
            IntNodes(:,1:Nint)=tempNodes(:,3:2+NInt)
        ENDIF
    ELSE
        A=DOT_PRODUCT(FaceNorm,(LineNodes(:,1)-FaceNodes(:,1)))
        B=DOT_PRODUCT(FaceNorm,(LineNodes(:,2)-FaceNodes(:,1)))
        IF (DABS(A).le.0.0000001) THEN !LineNodes(:,1)キ
            CALL vertex_in_polygon_3D(LineNodes(:,1),N,FaceNodes,CO_PLANE_TF)
            IF (CO_PLANE_TF) THEN
                NInt=1
                IntNodes(:,1)=LineNodes(:,1)
            ENDIF
        ELSE IF (DABS(B).le.0.0000001) THEN !LineNodes(:,2)キ
            CALL vertex_in_polygon_3D(LineNodes(:,2),N,FaceNodes,CO_PLANE_TF)
            IF (CO_PLANE_TF) THEN
                NInt=1
                IntNodes(:,1)=LineNodes(:,2)
            ENDIF
        ELSE IF ((A*B).lt.0) THEN !LineNodes(:,1)LineNodes(:,2)キㄢ凹
            B=B-A
            tempNodes(:,1)=-A*(LineNodes(:,2)-LineNodes(:,1))/B+LineNodes(:,1)
            CALL vertex_in_polygon_3D(tempNodes(:,1),N,FaceNodes,CO_PLANE_TF)
            IF (CO_PLANE_TF) THEN
                NInt=1
                IntNodes(:,1)=tempNodes(:,1)
            ENDIF
        ENDIF
    ENDIF
    

END SUBROUTINE line_int_face_3D
!=============================================================================!
!=============================================================================!
SUBROUTINE vertex_in_polygon_3D(Vertex,N,FaceNodes,True)
IMPLICIT NONE
INTEGER*4::N,Nt,i,j
REAL*8::Vertex(3),FaceNodes(3,N) !FaceNodes(3,N) must be sorted, and Vertex must be co-plane with polygon
REAL*8::r1(3),r2(3),r(3),A(2,2),B(2)
LOGICAL::True

    Nt=N-2
    r=Vertex-FaceNodes(:,1)
    True=.false.
    DO i=1,Nt
        r1=FaceNodes(:,i+1)-FaceNodes(:,1)
        r2=FaceNodes(:,i+2)-FaceNodes(:,1)
        IF ((r1(1)*r2(2)-r2(1)*r1(2)).ne.0) THEN
            A(1,1)=r2(2)
            A(1,2)=-r2(1)
            A(2,1)=-r1(2)
            A(2,2)=r1(1)
            A=A/(r1(1)*r2(2)-r2(1)*r1(2))
            B=MATMUL(A,r(1:2))
        ELSE
            IF ((r1(1)*r2(3)-r2(1)*r1(3)).ne.0) THEN
                A(1,1)=r2(3)
                A(1,2)=-r2(1)
                A(2,1)=-r1(3)
                A(2,2)=r1(1)
                A=A/(r1(1)*r2(3)-r2(1)*r1(3))
                B=MATMUL(A,r((/1,3/)))   
            ELSE
                IF ((r1(2)*r2(3)-r2(2)*r1(3)).ne.0) THEN
                    A(1,1)=r2(3)
                    A(1,2)=-r2(2)
                    A(2,1)=-r1(3)
                    A(2,2)=r1(2)
                    A=A/(r1(2)*r2(3)-r2(2)*r1(3))
                    B=MATMUL(A,r(2:3))
                ELSE
                    PAUSE "WRONG IN SUBROUTINE vertex_in_polygon_3D"
                ENDIF
            ENDIF
        ENDIF
        IF (ALL(B.ge.-0.0000001).AND.(SUM(B).le.1.0000001)) THEN
            True=.true.
            EXIT
        ENDIF
    ENDDO

END SUBROUTINE vertex_in_polygon_3D
!=============================================================================!
!=============================================================================!
SUBROUTINE vertex_in_tet_3D(Vertex,TetO,invTetVec,True)
IMPLICIT NONE
REAL*8::Vertex(3),invTetVec(3,3),TetO(3),r(3),A(3)
LOGICAL::True

    True=.false.
    r=Vertex-TetO
    A=MATMUL(invTetVec,r)
    IF (ALL(A.ge.-0.0000001).AND.(SUM(A).le.1.0000001)) True=.true.

END SUBROUTINE vertex_in_tet_3D
!=============================================================================!
!=============================================================================!
SUBROUTINE vertex_in_Cuboid_3D(Vertex,CuboidO,DiagVec,True)
IMPLICIT NONE
REAL*8::Vertex(3),CuboidO(3),DiagVec(3),r(3),A(3)
LOGICAL::True

    True=.false.
    r=Vertex-CuboidO
    A=r/DiagVec
    IF (ALL(A.ge.-0.0000001).AND.ALL(A.le.1.0000001)) True=.true.

END SUBROUTINE vertex_in_Cuboid_3D
!=============================================================================!
!=============================================================================!
SUBROUTINE sort_node_on_plane_3D(N,Nodes,NORM_TF,PlaneNorm)
!USE imsl_libraries
USE sort_real_int
IMPLICIT NONE
INTEGER*4::N,i,j,s
REAL*8::Nodes(3,N),tempNodes(3,N),PlaneNorm(3)
REAL*8::center(3),r(3,N),rsin(N),rcos(N)
LOGICAL::CO_PLANE_TF,NORM_TF
INTEGER*4,ALLOCATABLE::tempI1(:)
REAL*8,ALLOCATABLE::tempR1(:)
    
    CALL exam_coplane_3D(N,Nodes,CO_PLANE_TF,PlaneNorm,NORM_TF)
    IF (CO_PLANE_TF.eq..false.) THEN
        WRITE(*,*) "THIS NODES ARE NOT ON THE SAME PLANE!!!"
        PAUSE
        STOP
    ELSE
        center=SUM(Nodes,2)/DBLE(N)
        r(:,1)=Nodes(:,1)-center
        r(:,1)=r(:,1)/DSQRT(SUM(r(:,1)**2))
        rcos(1)=1d0
        rsin(1)=0d0
        s=COUNT(DABS(PlaneNorm).gt.0.0000001)
        ALLOCATE( tempI1(s),tempR1(3) )
        tempI1=PACK( (/(i,i=1,3)/),DABS(PlaneNorm).gt.0.0000001 )
        DO i=2,N
            r(:,i)=Nodes(:,i)-center
            r(:,i)=r(:,i)/DSQRT(SUM(r(:,i)**2))
            tempR1(1)=r(2,1)*r(3,i)-r(3,1)*r(2,i)
            tempR1(2)=r(3,1)*r(1,i)-r(1,1)*r(3,i)
            tempR1(3)=r(1,1)*r(2,i)-r(2,1)*r(1,i)
            tempR1(tempI1)=tempR1(tempI1)/PlaneNorm(tempI1)
            SELECT CASE(s)
                CASE(1)
                    rsin(i)=tempR1(tempI1(1))
                CASE(2)
                    IF (DABS(tempR1(tempI1(1))-tempR1(tempI1(2))).le.0.000001) THEN
                        rsin(i)=tempR1(tempI1(1))
                    ELSE
                        PAUSE "WRONG 1 IN SUBROUTINE sort_node_on_plane_3D"
                        STOP
                    ENDIF
                CASE(3)
                    IF ((DABS(tempR1(tempI1(1))-tempR1(tempI1(2))).le.0.000001).AND. &
                        (DABS(tempR1(tempI1(2))-tempR1(tempI1(3))).le.0.000001).AND. &
                        (DABS(tempR1(tempI1(3))-tempR1(tempI1(1))).le.0.000001)) THEN
                        rsin(i)=tempR1(tempI1(1))
                    ELSE
                        WRITE(*,*) PlaneNorm(tempI1(1:3))
                        WRITE(*,*) tempR1(tempI1(1:3))
                        PAUSE "WRONG 2 IN SUBROUTINE sort_node_on_plane_3D"
                        STOP
                    ENDIF
                CASE DEFAULT
                    PAUSE "WRONG 3 IN SUBROUTINE sort_node_on_plane_3D"
                    STOP
            ENDSELECT
            rcos(i)=DOT_PRODUCT(r(:,i),r(:,1))
        ENDDO
        DEALLOCATE( tempI1,tempR1 )
        
        tempNodes(:,1)=Nodes(:,1)
        
        s=COUNT(rsin(2:N).ge.0)
        ALLOCATE( tempI1(s),tempR1(s) )
        
        tempI1=PACK( (/(i,i=2,N)/),(rsin(2:N).ge.0) )
        CALL SORT_REAL(rcos(tempI1(:)),tempR1,IPERM=tempI1)
        j=1
        DO i=s,1,-1
            j=j+1
            tempNodes(:,j)=Nodes(:,tempI1(i))
        ENDDO
        DEALLOCATE( tempI1,tempR1 )
        
        s=COUNT(rsin(2:N).lt.0)
        ALLOCATE( tempI1(s),tempR1(s) )
        tempI1=PACK( (/(i,i=2,N)/),rsin(2:N).lt.0 )
        CALL SORT_REAL(rcos(tempI1(:)),tempR1,IPERM=tempI1)
        DO i=1,s
            j=j+1
            tempNodes(:,j)=Nodes(:,tempI1(i))
        ENDDO
        DEALLOCATE( tempI1,tempR1 )
        Nodes=tempNodes
    ENDIF
    

END SUBROUTINE sort_node_on_plane_3D
!=============================================================================!
!=============================================================================!
SUBROUTINE exam_coplane_3D(N,Nodes,CO_True,PlaneNorm,Norm_True)
IMPLICIT NONE
INTEGER*4::N,i
REAL*8::Nodes(3,N),r1(3),r2(3),r3(3),PlaneNorm(3)
LOGICAL::CO_True,Norm_True

    r1=Nodes(:,2)-Nodes(:,1)
    r2=Nodes(:,3)-Nodes(:,1)
    IF (Norm_True.eq..false.) THEN
        PlaneNorm(1)=r1(2)*r2(3)-r1(3)*r2(2)
        PlaneNorm(2)=r1(3)*r2(1)-r1(1)*r2(3)
        PlaneNorm(3)=r1(1)*r2(2)-r1(2)*r2(1)
        PlaneNorm=PlaneNorm/DSQRT(SUM(PlaneNorm**2))
        Norm_True=.true.
    ENDIF
    CO_True=.true.
    DO i=3,N
        r3=Nodes(:,i)-Nodes(:,1)
        IF (DABS(DOT_PRODUCT(PlaneNorm,r3)).gt.0.0000001) THEN
            CO_True=.false.
            RETURN
        ENDIF
    ENDDO

END SUBROUTINE exam_coplane_3D
!=============================================================================!
!=============================================================================!
SUBROUTINE del_reduplicate_nodes_3D(N,Nodes)
IMPLICIT NONE
INTEGER*4::N,tempN,i,j,s
REAL*8::Nodes(3,N),tempNodes(3,N)
LOGICAL::true(N)

    true=.false.
    tempN=0
    tempNodes=-10000
    DO i=1,N
        IF (true(i).eq..false.) THEN
            tempN=tempN+1
            tempNodes(:,tempN)=Nodes(:,i)
            true(i)=.true.
            DO j=i+1,N
                IF (ALL(DABS(Nodes(:,j)-Nodes(:,i)).lt.0.0000005)) THEN
                    true(j)=.true.
                ENDIF
            ENDDO
        ENDIF
    ENDDO
    N=tempN
    Nodes=tempNodes
    
END SUBROUTINE del_reduplicate_nodes_3D
!=============================================================================!
!=============================================================================!
SUBROUTINE tet_volume(Nodes,Volume)
IMPLICIT NONE
REAL*8::Nodes(3,4),Volume,r1(3),r2(3),r3(3)
    
    Volume=0

    r1=Nodes(:,2)-Nodes(:,1)
    r2=Nodes(:,3)-Nodes(:,1)
    r3=Nodes(:,4)-Nodes(:,1)
    
    Volume = DABS(r3(1)*(r1(2)*r2(3)-r1(3)*r2(2)) + &
                  r3(2)*(r1(3)*r2(1)-r1(1)*r2(3)) + &
                  r3(3)*(r1(1)*r2(2)-r1(2)*r2(1)))/6d0

END SUBROUTINE tet_volume
!=============================================================================!
!=============================================================================!
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