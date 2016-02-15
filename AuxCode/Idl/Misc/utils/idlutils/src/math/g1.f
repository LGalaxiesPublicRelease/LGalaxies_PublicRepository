      SUBROUTINE G1 (A,B,C,S,SIG)
      ZERO=0.
      ONE=1.
      IF (ABS(A).LE.ABS(B)) GO TO 10
      XR=B/A
      YR=SQRT(ONE+XR**2)
      C=ONE/YR
      IF (A.GE.0) GO TO 5
      C=-C
  5   CONTINUE
      S=C*XR
      SIG=ABS(A)*YR
      RETURN
 10   IF (B) 20,30,20
 20   XR=A/B
      YR=SQRT(ONE+XR**2)
      S=ONE/YR
      IF (B.GE.0) GO TO 25
      S=-S
 25   CONTINUE
      C=S*XR
      SIG=ABS(B)*YR
      RETURN
 30   SIG=ZERO
      C=ZERO
      S=ONE
C      C=ONE
C      S=ZERO
      RETURN
      END
