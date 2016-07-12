      SUBROUTINE STERN (JJA,NST,STG,TAUP)                               
!                                                                       
!     STERN GENERATES THE STAR OF REC LATTICE VECTOR KZZ(1,JJA)         
!     THE STAR VECTORS ARE STORED IN STG, THE STAR-SIZE IN  NST         
!     IZ CONTAINS THE SYMMETRY-MATRICES                                 
!                                                                       
!                                                                       
      USE param
      IMPLICIT REAL*8 (A-H,O-Z)
!
!
      COMPLEX*16          TAUP,IMAG                                     
      INTEGER          G,STG,INDEX(48)                                
      COMMON /SYMETR/  TAU(3,48),KZZ(3,10000),IZ(3,3,48),             &
                       INUM(48),IORD                                  
!                                                                       
      DIMENSION        STG(3,48),TAUP(48),G(3)                      
      DATA IMAG        /(0.D0,1.D0)/                                    
!-------------------------------------------------------------------    
!                                                                       
      TPI=2.D0*ACOS(-1.D0)                                              
      G(1)=KZZ(1,JJA)                                                   
      G(2)=KZZ(2,JJA)                                                   
      G(3)=KZZ(3,JJA)                                                   
      NST=0                                                             
!                                                                       
!.....START LOOP OVER ALL SYMMETRY OPERATIONS                           
      DO 1 I=1,IORD                                                     
         TK=0.                                                          
         DO 2 J=1,3                                                     
         TK=TK+TAU(J,I)*G(J)*TPI                                     
           K=0                                                         
           DO 3 L=1,3                                                  
           K=IZ(J,L,I)*G(L)+K                                          
3     continue
           STG(J,I)=K                                                     
 2      continue
         IF(NST.EQ.0) GOTO 7                                            
!                                                                       
!........PROOF, IF THE VECTOR STG(J,I) IS A NEW STARMEMBER OR NOT       
         DO 4 M=1,NST                                                   
            DO 5 J=1,3                                                  
               IF(STG(J,M).NE.STG(J,I)) GOTO 4                          
   5        CONTINUE                                                    
!           STG(J,I) IS NOT A NEW STARMEMBER, IT IS EQUIV TO STG(J,M).  
!           BUT THE TAUPHASE OF STG(J,I) CAN BE NEW.  THEREFORE THE     
!           ALREADY DEFINED PHASE TAUP(M) IS AVERAGED WITH THE PHASE    
!           OF STG(J,M)                                                 
            TAUP(M)=TAUP(M)+EXP(TK*IMAG)                                
            INDEX(M)=INDEX(M)+1                                         
            GOTO 1                                                      
   4     CONTINUE                                                       
!                                                                       
!........NEW VECTOR FOUND ]]]                                           
   7     NST=NST+1                                                      
         DO 6 J=1,3                                                     
   6     STG(J,NST)=STG(J,I)                                            
         TAUP(NST)=EXP(TK*IMAG)                                         
         INDEX(NST)=1                                                   
   1  CONTINUE                                                          
!                                                                       
      DO 10 I=1,NST           
      TAUP(I)=TAUP(I)/INDEX(I)                                          
  10  continue
      RETURN                                                            
      END                                                               
