MODULE rng
IMPLICIT NONE
PRIVATE
PUBLIC::rng_t,rng_seed,RAN_NUM
!=================================================================!
INTEGER*4,PARAMETER::ns=4
INTEGER*4,PARAMETER,DIMENSION(ns)::default_seed=(/521288629,362436069,16163801,1131199299/)
!=================================================================!
! A DATA TYPE FOR STORING THE STATE OF THE RNG
TYPE::rng_t
    INTEGER*4,DIMENSION(ns)::state
END TYPE rng_t
!=================================================================!
INTERFACE RAN_NUM
    MODULE PROCEDURE RAN_NUM_ONE_DOUBLE
    MODULE PROCEDURE RAN_NUM_DOUBLE_ARRAY
    MODULE PROCEDURE RAN_NUM_DOUBLE_ARRAY_2D
    MODULE PROCEDURE RAN_NUM_DOUBLE_ARRAY_3D
END INTERFACE RAN_NUM
!=================================================================!
CONTAINS

    !=================================================================!
    ! SEEDS THE RNG USING A DOUBLE INTEGER AND A DEFAULT SEED VECTOR
    SUBROUTINE rng_seed(self,seed)
    IMPLICIT NONE
    TYPE(rng_t),INTENT(inout)::self
    INTEGER*4,INTENT(in)::seed
        
        self%state=default_seed
        self%state(1)=seed
        self%state(2:ns)=default_seed(2:ns)
    
    END SUBROUTINE rng_seed
    !=================================================================!
    ! DRAWS A UNIFORM REAL NUMBER ON [0,1]
    SUBROUTINE RAN_NUM_ONE_DOUBLE(self,u)
    IMPLICIT NONE
    TYPE(rng_t),INTENT(inout)::self
    REAL*8::u
    INTEGER*4::imz
    
        imz=self%state(1)-self%state(3)
        
        IF (imz < 0) imz=imz+2147483579
        
        self%state(1)=self%state(2)
        self%state(2)=self%state(3)
        self%state(3)=imz
        self%state(4)=69069*self%state(4)+1013904243
        
        imz=imz+self%state(4)
        u=0.5d0 + 0.23283064d-9 * DBLE(imz)
    
    END SUBROUTINE RAN_NUM_ONE_DOUBLE
    !=================================================================!
    SUBROUTINE RAN_NUM_DOUBLE_ARRAY(self,A)
    IMPLICIT NONE
    TYPE(rng_t),INTENT(inout)::self
    INTEGER*4::i,N(1)
    REAL*8::A(:)
    
        N=SIZE(A)
        DO i=1,N(1)
            CALL RAN_NUM_ONE_DOUBLE(self,A(i))
        ENDDO
    
    END SUBROUTINE RAN_NUM_DOUBLE_ARRAY
    !=================================================================!
    SUBROUTINE RAN_NUM_DOUBLE_ARRAY_2D(self,A)
    IMPLICIT NONE
    TYPE(rng_t),INTENT(inout)::self
    INTEGER*4::i,j,N(2)
    REAL*8::A(:,:)
    
        N=SHAPE(A)
        DO j=1,N(2)
            DO i=1,N(1)
                CALL RAN_NUM_ONE_DOUBLE(self,A(i,j))
            ENDDO
        ENDDO
    
    END SUBROUTINE RAN_NUM_DOUBLE_ARRAY_2D
    !=================================================================!
    !=================================================================!
    SUBROUTINE RAN_NUM_DOUBLE_ARRAY_3D(self,A)
    IMPLICIT NONE
    TYPE(rng_t),INTENT(inout)::self
    INTEGER*4::i,j,k,N(3)
    REAL*8::A(:,:,:)
    
        N=SHAPE(A)
        DO k=1,N(3)
            DO j=1,N(2)
                DO i=1,N(1)
                    CALL RAN_NUM_ONE_DOUBLE(self,A(i,j,k))
                ENDDO
            ENDDO
        ENDDO
    
    END SUBROUTINE RAN_NUM_DOUBLE_ARRAY_3D
    !=================================================================!
END MODULE rng
