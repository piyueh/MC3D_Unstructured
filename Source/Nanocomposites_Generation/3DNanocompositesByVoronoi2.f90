MODULE variables
IMPLICIT NONE
!---------------------------!
    TYPE CELL
        INTEGER*4::arraysize
        INTEGER*4,ALLOCATABLE::array(:)
    END TYPE CELL
!---------------------------!

    REAL*8::DomainL(3)
    INTEGER*4::Ngrains,Nnodes,Nedges,Nfaces
    
    REAL*8,ALLOCATABLE::GrainCenter(:,:)
    REAL*8,ALLOCATABLE::Vertexes(:,:)
    INTEGER*4,ALLOCATABLE::Edges(:,:)
    
    TYPE(CELL),ALLOCATABLE::Faces(:)
    TYPE(CELL),ALLOCATABLE::Nanograins(:)

END MODULE variables
!===========================!
!===========================!
!===========================!
PROGRAM main
USE variables
USE IFQWIN
IMPLICIT NONE
REAL*8::dL(3)
INTEGER*4::i,j,k,s,m,c
REAL*8,ALLOCATABLE::tempR1(:,:)
INTEGER*4,ALLOCATABLE::transf(:),tempI2(:)
CHARACTER*72::OutPutTempFile
LOGICAL::true

    WRITE(*,*) "ENTER MODEL DIMENSION (Lx,Ly,Lz) : "
    READ(*,*) DomainL
    
    WRITE(*,*) "ENTER THE NUMBER OF TOTAL GRAINS : "
    READ(*,*) Ngrains

    true=.false.
    ALLOCATE( GrainCenter(3,Ngrains))
    c=0
    DO WHILE (true.eq..false.)
        c=c+1
        ALLOCATE( tempR1(3,Ngrains*7) )
    
        CALL initialize( Ngrains,DomainL,GrainCenter,tempR1 )
    
        OutPutTempFile='tempOUTPUT.txt'
        CALL construct_voronoi(3,SIZE(tempR1,2),tempR1,OutPutTempFile)
        DEALLOCATE( tempR1 )
    
        OPEN(UNIT=100,FILE=OutPutTempFile)
        READ(100,*) i
        READ(100,*) Nnodes,j,k
    
        ALLOCATE( tempR1(3,0:Nnodes-1),transf(0:Nnodes-1) )
        READ(100,*) tempR1

        CALL construct_transf(DomainL,Nnodes,tempR1,transf,j)
    
        ALLOCATE( Vertexes(3,j) )
        FORALL( j=0:Nnodes-1, transf(j).ne.-999 ) Vertexes(:,transf(j))=tempR1(:,j)
        Nnodes=j
        DEALLOCATE( tempR1 )
    
        OutPutTempFile='model.txt'
        CALL construct_model(100,3,Ngrains,SIZE(transf),Nnodes,transf,Vertexes,OutPutTempFile,true)
        CLOSE(UNIT=100,STATUS='DELETE')
        DEALLOCATE( transf,Vertexes )
        WRITE(*,*) c
    ENDDO
    
    OPEN(UNIT=991,FILE='GrainCenter.txt')
    WRITE(991,*) GrainCenter
    CLOSE(991)
    
    CALL export_gambit_jou_file(OutPutTempFile)
    PAUSE 'PRESS ANY KEY TO FINISH'
    
END PROGRAM main
!===========================!
!===========================!
SUBROUTINE initialize( N,L,xyzNodes,temp )
IMPLICIT NONE
INTEGER*4::N,i,j,k,s,m,z
REAL*8::L(3),xyzNodes(3,N),temp(3,N*7),Ri,Ro
LOGICAL::true
REAL*8,ALLOCATABLE::xxyyzz(:,:,:),R(:)

    ALLOCATE( xxyyzz(3,N,7) )

    Ri=(L(1)*L(2)*L(3)*3d0*0.25d0/3.1415926/DBLE(N))**(1d0/3d0)*1.4
    Ro=(L(1)*L(2)*L(3)*3d0*0.25d0/3.1415926/DBLE(N))**(1d0/3d0)*1.65
    
    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER(xxyyzz(:,1,1))
    
    xxyyzz(:,1,1)=xxyyzz(:,1,1)*L(:)
    DO i=2,N
        true=.false.
        m=0
        ALLOCATE( R(i-1) )
        DO WHILE (true.eq..false.)
            m=m+1
            CALL random_number(xxyyzz(:,i,1))
            xxyyzz(:,i,1)=xxyyzz(:,i,1)*L(:)
            FORALL(j=1:i-1) R(j)=DSQRT(SUM((xxyyzz(:,i,1)-xxyyzz(:,j,1))**2))
            IF (ALL(R.ge.Ri).AND.(MINVAL(R).le.Ro)) true=.true.
        ENDDO
        DEALLOCATE( R )
        !WRITE(*,*) i,m
    ENDDO
    
    xxyyzz(1,:,2)=-xxyyzz(1,:,1)
    xxyyzz(2:3,:,2)=xxyyzz(2:3,:,1)
    
    xxyyzz(1,:,3)=2*L(1)-xxyyzz(1,:,1)
    xxyyzz(2:3,:,3)=xxyyzz(2:3,:,1)
    
    xxyyzz(1,:,4)=xxyyzz(1,:,1)
    xxyyzz(2,:,4)=-xxyyzz(2,:,1)
    xxyyzz(3,:,4)=xxyyzz(3,:,1)
    
    xxyyzz(1,:,5)=xxyyzz(1,:,1)
    xxyyzz(2,:,5)=2*L(2)-xxyyzz(2,:,1)
    xxyyzz(3,:,5)=xxyyzz(3,:,1)
    
    xxyyzz(1,:,6)=xxyyzz(1,:,1)
    xxyyzz(2,:,6)=xxyyzz(2,:,1)
    xxyyzz(3,:,6)=-xxyyzz(3,:,1)
    
    xxyyzz(1,:,7)=xxyyzz(1,:,1)
    xxyyzz(2,:,7)=xxyyzz(2,:,1)
    xxyyzz(3,:,7)=2*L(3)-xxyyzz(3,:,1)
    
    s=0
    m=N
    temp=-10000d0
    xyzNodes=-10000d0
    DO j=1,7
        DO i=1,N
            IF (j.eq.1) THEN
                s=s+1
                temp(:,s)=xxyyzz(:,i,j)
                xyzNodes(:,s)=xxyyzz(:,i,j)
            ELSE
                m=m+1
                temp(:,m)=xxyyzz(:,i,j)
            ENDIF
        ENDDO
    ENDDO

    IF ( ANY( temp.eq.-10000d0 ).OR.ANY(xyzNodes.eq.-10000d0) ) PAUSE 'WRONG 1 IN INITIALIZE!'
    IF ( ANY( DABS(temp(:,1:s)-xyzNodes).ge.0.0000001 ) ) PAUSE 'WRONG 2 IN INITIALIZE!'
    
    DEALLOCATE( xxyyzz )
    
END SUBROUTINE initialize
!===========================!
!===========================!
SUBROUTINE construct_voronoi(ND,NS,Nodes,VoronoiFile)
IMPLICIT NONE
INTEGER*4::ND,NS,i
REAL*8::Nodes(ND,NS)
CHARACTER*72::VoronoiFile

    OPEN(UNIT=100,FILE='temp.txt')
    WRITE(100,*) ND
    WRITE(100,*) NS
    DO i=1,NS
        WRITE(100,*) Nodes(:,i)
    ENDDO

    FLUSH(UNIT=100) !此文字檔尚未關閉，使用此指令使得別的程式此時也可以使用這文字檔

    CALL SYSTEM('qvoronoi<temp.txt o TO '//VoronoiFile)
    CLOSE(UNIT=100,STATUS='DELETE')

END SUBROUTINE construct_voronoi
!===========================!
!===========================!
SUBROUTINE construct_transf(L,N1,xyzNodes,transf,N2)
IMPLICIT NONE
INTEGER*4::N1,N2,transf(0:N1-1),i,j,s,m
REAL*8::L(3),xyzNodes(3,0:N1-1)
LOGICAL,ALLOCATABLE:: true(:)

    ALLOCATE( true(0:N1-1) )
    true=.true.
    s=0
    m=0
    transf=-999
    
    FORALL(i=0:N1-1,(ALL(xyzNodes(:,i).ge.-0.0000001).AND.ALL(xyzNodes(:,i)-L(:).le.0.0000001))) true(i)=.false.
    
    DO i=0,N1-1
        IF (true(i).eq..false.) THEN
            s=s+1
            transf(i)=s
            true(i)=.true.
            DO j=i+1,N1-1
                IF ( ALL(DABS(xyzNodes(:,j)-xyzNodes(:,i)).le.0.0000001) ) THEN
                    m=m+1
                    true(j)=.true.
                    transf(j)=s
                ENDIF
            ENDDO
        ENDIF
    ENDDO
    N2=s
    DEALLOCATE( true )

END SUBROUTINE construct_transf
!===========================!
!===========================!
SUBROUTINE construct_model(RF1,ND,N1,N2,N3,transf,xyzNodes,OutPutFile,true3)
IMPLICIT NONE
!---------------------------!
TYPE CELL
    INTEGER*4::arraysize
    INTEGER*4,ALLOCATABLE::array(:)
END TYPE CELL
!---------------------------!
REAL*8::xyzNodes(3,N3),volume(N1)
INTEGER*4::i,j,k,s,ND,N1,N2,N3,RF1,RF2,NF,NE
INTEGER*4::tempI1,transf(0:N2-1)
CHARACTER*72::OutPutFile,string
LOGICAL::true3
TYPE(CELL)::Polyhedron(N1)
!---------------------------!
INTEGER*4,ALLOCATABLE::tempI2(:),tempI3(:),transfF(:),edge(:,:),tempEdge(:,:)
LOGICAL,ALLOCATABLE::true1(:),true2(:)
TYPE(CELL),ALLOCATABLE::Polygon(:),tempC1(:)
!---------------------------!

    RF2=RF1+1
    NF=0
    ALLOCATE( Polygon(20*N1) )
    DO i=1,N1
        
        READ(UNIT=RF1,FMT="(I2)",ADVANCE='NO') tempI1 !tempI1=第i個多面體的頂點數
        ALLOCATE( tempI2(tempI1),tempI3(tempI1),true1(tempI1) )
        READ(RF1,*) tempI2 
        
        true1=.false.
        tempI2=transf(tempI2) !tempI2=第i個多面體的頂點編號
        s=0
        DO j=1,tempI1
            IF (true1(j).eq..false.) THEN
                s=s+1
                tempI3(s)=tempI2(j)
                true1(j)=.true.
                DO k=j+1,tempI1
                    IF ((tempI2(k)-tempI2(j)).eq.0) THEN
                        true1(k)=.true.
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
        
        DEALLOCATE( tempI2,true1 )
        ALLOCATE( tempI2(s) )
        tempI2=tempI3(1:s)
        DEALLOCATE( tempI3 )
        
        OPEN(UNIT=RF2,FILE='temp.txt')
        WRITE(RF2,*) ND
        WRITE(RF2,*) s
        DO j=1,s
            WRITE(RF2,*) xyzNodes(:,tempI2(j))
        ENDDO
        CALL SYSTEM('qhull<temp.txt Fv i TO temp2.txt')
        CALL SYSTEM('qconvex<temp.txt FA>volume.txt')
        OPEN(UNIT=RF2+1,FILE='volume.txt')
        true3=.true.
        DO WHILE(true3)
            CALL read_tmp(RF2+1,string)
            IF (INDEX(string,'volume').ne.0) THEN
                READ(string(INDEX(string,':')+1:),*) volume(i)
                true3=.false.
            ENDIF
        ENDDO
        CLOSE(UNIT=RF2+1,STATUS='DELETE')
        CLOSE(UNIT=RF2,STATUS='DELETE')
        
        OPEN(UNIT=RF2,FILE='temp2.txt')
        READ(RF2,*) Polyhedron(i)%arraysize !第i個多面體的面數
        ALLOCATE( Polyhedron(i)%array(Polyhedron(i)%arraysize) )
        DO j=1,Polyhedron(i)%arraysize
            NF=NF+1
            Polyhedron(i)%array(j)=NF
            READ(RF2,*) Polygon(NF)%arraysize
            ALLOCATE( Polygon(NF)%array(Polygon(NF)%arraysize) )
        ENDDO
        
        READ(RF2,*) j
        
        IF (j.ne.Polyhedron(i)%arraysize) THEN
            PAUSE 'WRONG'
            STOP
        ELSE
            DO j=1,Polyhedron(i)%arraysize
                READ(RF2,*) Polygon(Polyhedron(i)%array(j))%array !第i個多面體的第j個面的頂點在tempI2裡的編號
                Polygon(Polyhedron(i)%array(j))%array = tempI2(Polygon(Polyhedron(i)%array(j))%array+1) !把頂點編號轉變為GLOBAL編號
            ENDDO
        ENDIF
        CLOSE(UNIT=RF2,STATUS='DELETE')
        DEALLOCATE( tempI2 )
    ENDDO
    
    CALL exam_volume_distribution(N1,volume,true3)
    IF (true3.eq..false.) RETURN
    
    ALLOCATE( tempC1(NF),true1(NF),transfF(NF) )
    s=0
    true1=.false.
    DO i=1,NF
        IF (true1(i).eq..false.) THEN
            s=s+1
            true1(i)=.true.
            tempC1(s)=Polygon(i)
            transfF(i)=s
            DO j=i+1,NF
                IF (Polygon(j)%arraysize.eq.Polygon(i)%arraysize) THEN
                    ALLOCATE( true2(Polygon(i)%arraysize) )
                    true2=.false.
                    DO k=1,Polygon(i)%arraysize
                        IF (ANY( Polygon(j)%array.eq.Polygon(i)%array(k) )) true2(k)=.true.
                    ENDDO
                    IF (ALL( true2 )) THEN
                        true1(j)=.true.
                        transfF(j)=s
                    ENDIF
                    DEALLOCATE(true2)
                ENDIF
            ENDDO
        ENDIF
    ENDDO
    
    DEALLOCATE( Polygon )
    NF=s
    ALLOCATE( Polygon(NF) )
    Polygon=tempC1(1:NF)
    DEALLOCATE( tempC1,true1 )
    
    FORALL(i=1:N1) Polyhedron(i)%array=transfF(Polyhedron(i)%array)
    DEALLOCATE(transfF)

    ALLOCATE( Edge(2,10*NF),tempC1(NF) )
    tempC1=Polygon
    NE=0
    DO i=1,NF
        DO j=1,Polygon(i)%arraysize
            NE=NE+1
            IF (j.ne.Polygon(i)%arraysize) THEN
                Edge(:,NE)=(/ Polygon(i)%array(j),Polygon(i)%array(j+1) /)
            ELSE
                Edge(:,NE)=(/ Polygon(i)%array(j),Polygon(i)%array(1) /)
            ENDIF
            tempC1(i)%array(j)=NE
        ENDDO
    ENDDO
    Polygon=tempC1
    DEALLOCATE( tempC1 )
    
    ALLOCATE( tempEdge(2,NE),true1(NE),transfF(NE) )
    s=0
    true1=.false.
    DO i=1,NE
        IF (true1(i).eq..false.) THEN
            s=s+1
            tempEdge(:,s)=Edge(:,i)
            true1(i)=.true.
            transfF(i)=s
            DO j=i+1,NE
                IF ( ANY(Edge(:,j).eq.Edge(1,i)) .AND. ANY(Edge(:,j).eq.Edge(2,i)) ) THEN
                    true1(j)=.true.
                    transfF(j)=s
                ENDIF
            ENDDO
        ENDIF
    ENDDO
    
    NE=s
    DEALLOCATE( Edge,true1 )
    ALLOCATE( Edge(2,NE) )
    Edge=tempEdge(:,1:NE)
    FORALL(i=1:NF) Polygon(i)%array=transfF(Polygon(i)%array)
    DEALLOCATE( transfF,tempEdge )
    
    OPEN(UNIT=RF2,FILE=OutPutFile)
    WRITE(RF2,*) N1,NF,NE,N3
    WRITE(RF2,*) Polyhedron(:)%arraysize
    DO i=1,N1
        WRITE(RF2,*) Polyhedron(i)%array
    ENDDO
    WRITE(RF2,*) Polygon(:)%arraysize
    DO i=1,NF
        WRITE(RF2,*) Polygon(i)%array
    ENDDO
    WRITE(RF2,*) Edge
    WRITE(RF2,*) xyzNodes
    CLOSE(RF2)
    
END SUBROUTINE construct_model
!===========================!
!===========================!
SUBROUTINE export_gambit_jou_file(InputTXT)
IMPLICIT NONE
!---------------------------!
TYPE CELL
    INTEGER*4::arraysize
    INTEGER*4,ALLOCATABLE::array(:)
END TYPE CELL
!---------------------------!
INTEGER*4::NV,NE,NF,NB,i,j,k,tempI
CHARACTER*72::InputTXT
!---------------------------!
REAL*8,ALLOCATABLE::Vertex(:,:)
INTEGER*4,ALLOCATABLE::Edge(:,:)
TYPE(CELL),ALLOCATABLE::Face(:),Volume(:)
!---------------------------!

    OPEN(UNIT=100,FILE=InputTXT)

    READ(100,*) NB,NF,NE,NV
    ALLOCATE( Volume(NB),Face(NF),Edge(2,NE),Vertex(3,NV) )
    READ(100,*) Volume(:)%arraysize
    DO i=1,NB
        ALLOCATE( Volume(i)%array(Volume(i)%arraysize) )
        READ(100,*) Volume(i)%array
    ENDDO
    READ(100,*) Face(:)%arraysize
    DO i=1,NF
        ALLOCATE( Face(i)%array(Face(i)%arraysize) )
        READ(100,*) Face(i)%array
    ENDDO
    READ(100,*) Edge
    READ(100,*) Vertex
    
    CLOSE(100)
    !----------------------------!

    OPEN(UNIT=100,FILE='ForGambit.jou')
    DO i=1,NV
        WRITE(100,500) Vertex(:,i)
    ENDDO
    
    DO i=1,NE
        WRITE(UNIT=100,FMT=501,ADVANCE='NO') 
        DO j=1,2
            tempI=INT( LOG10(DBLE(Edge(j,i)))+1 )
            WRITE(UNIT=100,FMT=504,ADVANCE='NO')  Edge(j,i)
        ENDDO
        WRITE(UNIT=100,FMT=507,ADVANCE='YES')
    ENDDO
    
    DO i=1,NF
        WRITE(UNIT=100,FMT=502,ADVANCE='NO') 
        DO j=1,Face(i)%arraysize
            tempI=INT( LOG10(DBLE(Face(i)%array(j)))+1 )
            WRITE(UNIT=100,FMT=505,ADVANCE='NO')  Face(i)%array(j)
        ENDDO
        WRITE(UNIT=100,FMT=507,ADVANCE='YES')
    ENDDO
    
    DO i=1,NB
        WRITE(UNIT=100,FMT=503,ADVANCE='NO') 
        DO j=1,Volume(i)%arraysize
            tempI=INT( LOG10(DBLE(Volume(i)%array(j)))+1 )
            WRITE(UNIT=100,FMT=506,ADVANCE='NO')  Volume(i)%array(j)
        ENDDO
        WRITE(UNIT=100,FMT=507,ADVANCE='YES')
    ENDDO
    
    CLOSE(100)
    500 FORMAT ("vertex create coordinates ",3(E15.8,1x))
    501 FORMAT ('edge create straight ')
    502 FORMAT ('face create wireframe ')
    503 FORMAT ('volume create stitch ')
    504 FORMAT ('"vertex.',I<tempI>,'" ')
    505 FORMAT ('"edge.',I<tempI>,'" ')
    506 FORMAT ('"face.',I<tempI>,'" ')
    507 FORMAT ('real')

END SUBROUTINE export_gambit_jou_file
!===========================!
!===========================!
SUBROUTINE exam_volume_distribution(N1,volume,true)
IMPLICIT NONE
INTEGER*4::N1,i,m
REAL*8::volume(N1),tot_V,ave_V,STD
LOGICAL::true

    tot_V=SUM(volume)
    ave_V=tot_V/DBLE(N1)
    STD=DSQRT(SUM((volume-ave_V)**2)/DBLE(N1))
    
    m=0
    DO i=1,N1
        IF (DABS(volume(i)-ave_V).le.STD) m=m+1
    ENDDO
    
    WRITE(*,*) DBLE(m)/DBLE(N1),ave_V,STD
    
    IF ( DBLE(m)/DBLE(N1).ge.0.65 ) THEN
        true=.true.
    ELSE
        true=.false.
    ENDIF


END SUBROUTINE exam_volume_distribution
!===========================!
!===========================!
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
!===========================!
!===========================!