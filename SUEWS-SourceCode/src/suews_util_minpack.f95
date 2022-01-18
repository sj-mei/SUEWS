SUBROUTINE chkder(m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode, err)

!*****************************************************************************80
!
!! CHKDER checks the gradients of M functions of N variables.
!
!  Discussion:
!
!    CHKDER checks the gradients of M nonlinear functions in N variables,
!    evaluated at a point X, for consistency with the functions themselves.
!
!    The user calls CHKDER twice, first with MODE = 1 and then with MODE = 2.
!
!    MODE = 1.
!      On input,
!        X contains the point of evaluation.
!      On output,
!        XP is set to a neighboring point.
!
!    Now the user must evaluate the function and gradients at X, and the
!    function at XP.  Then the subroutine is called again:
!
!    MODE = 2.
!      On input,
!        FVEC contains the function values at X,
!        FJAC contains the function gradients at X.
!        FVECP contains the functions evaluated at XP.
!      On output,
!        ERR contains measures of correctness of the respective gradients.
!
!    The subroutine does not perform reliably if cancellation or
!    rounding errors cause a severe loss of significance in the
!    evaluation of a function.  Therefore, none of the components
!    of X should be unusually small (in particular, zero) or any
!    other value which may cause loss of significance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, is the number of functions.
!
!    Input, integer ( kind = 4 ) N, is the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian is to be
!    evaluated.
!
!    Input, real ( kind = 8 ) FVEC(M), is used only when MODE = 2.
!    In that case, it should contain the function values at X.
!
!    Input, real ( kind = 8 ) FJAC(LDFJAC,N), an M by N array.  When MODE = 2,
!    FJAC(I,J) should contain the value of dF(I)/dX(J).
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of FJAC.
!    LDFJAC must be at least M.
!
!    Output, real ( kind = 8 ) XP(N), on output with MODE = 1, is a neighboring
!    point of X, at which the function is to be evaluated.
!
!    Input, real ( kind = 8 ) FVECP(M), on input with MODE = 2, is the function
!    value at XP.
!
!    Input, integer ( kind = 4 ) MODE, should be set to 1 on the first call, and
!    2 on the second.
!
!    Output, real ( kind = 8 ) ERR(M).  On output when MODE = 2, ERR contains
!    measures of correctness of the respective gradients.  If there is no
!    severe loss of significance, then if ERR(I):
!      = 1.0D+00, the I-th gradient is correct,
!      = 0.0D+00, the I-th gradient is incorrect.
!      > 0.5D+00, the I-th gradient is probably correct.
!      < 0.5D+00, the I-th gradient is probably incorrect.
!
   IMPLICIT NONE

   INTEGER(kind=4) ldfjac
   INTEGER(kind=4) m
   INTEGER(kind=4) n

   REAL(kind=8) eps
   REAL(kind=8) epsf
   REAL(kind=8) epslog
   REAL(kind=8) epsmch
   REAL(kind=8) err(m)
   REAL(kind=8) fjac(ldfjac, n)
   REAL(kind=8) fvec(m)
   REAL(kind=8) fvecp(m)
   INTEGER(kind=4) i
   INTEGER(kind=4) j
   INTEGER(kind=4) mode
   REAL(kind=8) temp
   REAL(kind=8) x(n)
   REAL(kind=8) xp(n)

   epsmch = EPSILON(epsmch)
   eps = SQRT(epsmch)
!
!  MODE = 1.
!
   IF (mode == 1) THEN

      DO j = 1, n
         temp = eps*ABS(x(j))
         IF (temp == 0.0D+00) THEN
            temp = eps
         END IF
         xp(j) = x(j) + temp
      END DO
!
!  MODE = 2.
!
   ELSE IF (mode == 2) THEN

      epsf = 100.0D+00*epsmch
      epslog = LOG10(eps)

      err = 0.0D+00

      DO j = 1, n
         temp = ABS(x(j))
         IF (temp == 0.0D+00) THEN
            temp = 1.0D+00
         END IF
         err(1:m) = err(1:m) + temp*fjac(1:m, j)
      END DO

      DO i = 1, m

         temp = 1.0D+00

         IF (fvec(i) /= 0.0D+00 .AND. fvecp(i) /= 0.0D+00 .AND. &
             ABS(fvecp(i) - fvec(i)) >= epsf*ABS(fvec(i))) THEN
            temp = eps*ABS((fvecp(i) - fvec(i))/eps - err(i)) &
                   /(ABS(fvec(i)) + ABS(fvecp(i)))
         END IF

         err(i) = 1.0D+00

         IF (epsmch < temp .AND. temp < eps) THEN
            err(i) = (LOG10(temp) - epslog)/epslog
         END IF

         IF (eps <= temp) THEN
            err(i) = 0.0D+00
         END IF

      END DO

   END IF

   RETURN
END SUBROUTINE chkder
SUBROUTINE dogleg(n, r, lr, diag, qtb, delta, x)

!*****************************************************************************80
!
!! DOGLEG finds the minimizing combination of Gauss-Newton and gradient steps.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N nonsingular diagonal
!    matrix D, an M-vector B, and a positive number DELTA, the
!    problem is to determine the convex combination X of the
!    Gauss-Newton and scaled gradient directions that minimizes
!    (A*X - B) in the least squares sense, subject to the
!    restriction that the euclidean norm of D*X be at most DELTA.
!
!    This subroutine completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization of A.  That is, if A = Q*R, where Q has
!    orthogonal columns and R is an upper triangular matrix,
!    then DOGLEG expects the full upper triangle of R and
!    the first N components of Q'*B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix R.
!
!    Input, real ( kind = 8 ) R(LR), the upper triangular matrix R stored
!    by rows.
!
!    Input, integer ( kind = 4 ) LR, the size of the R array, which must be
!    no less than (N*(N+1))/2.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real ( kind = 8 ) QTB(N), the first N elements of the vector Q'* B.
!
!    Input, real ( kind = 8 ) DELTA, is a positive upper bound on the
!    euclidean norm of D*X(1:N).
!
!    Output, real ( kind = 8 ) X(N), the desired convex combination of the
!    Gauss-Newton direction and the scaled gradient direction.
!
   IMPLICIT NONE

   INTEGER(kind=4) lr
   INTEGER(kind=4) n

   REAL(kind=8) alpha
   REAL(kind=8) bnorm
   REAL(kind=8) delta
   REAL(kind=8) diag(n)
   REAL(kind=8) enorm
   REAL(kind=8) epsmch
   REAL(kind=8) gnorm
   INTEGER(kind=4) i
   INTEGER(kind=4) j
   INTEGER(kind=4) jj
   INTEGER(kind=4) k
   INTEGER(kind=4) l
   REAL(kind=8) qnorm
   REAL(kind=8) qtb(n)
   REAL(kind=8) r(lr)
   REAL(kind=8) sgnorm
   REAL(kind=8) sum2
   REAL(kind=8) temp
   REAL(kind=8) wa1(n)
   REAL(kind=8) wa2(n)
   REAL(kind=8) x(n)

   epsmch = EPSILON(epsmch)
!
!  Calculate the Gauss-Newton direction.
!
   jj = (n*(n + 1))/2 + 1

   DO k = 1, n

      j = n - k + 1
      jj = jj - k
      l = jj + 1
      sum2 = 0.0D+00

      DO i = j + 1, n
         sum2 = sum2 + r(l)*x(i)
         l = l + 1
      END DO

      temp = r(jj)

      IF (temp == 0.0D+00) THEN

         l = j
         DO i = 1, j
            temp = MAX(temp, ABS(r(l)))
            l = l + n - i
         END DO

         IF (temp == 0.0D+00) THEN
            temp = epsmch
         ELSE
            temp = epsmch*temp
         END IF

      END IF

      x(j) = (qtb(j) - sum2)/temp

   END DO
!
!  Test whether the Gauss-Newton direction is acceptable.
!
   wa1(1:n) = 0.0D+00
   wa2(1:n) = diag(1:n)*x(1:n)
   qnorm = enorm(n, wa2)

   IF (qnorm <= delta) THEN
      RETURN
   END IF
!
!  The Gauss-Newton direction is not acceptable.
!  Calculate the scaled gradient direction.
!
   l = 1
   DO j = 1, n
      temp = qtb(j)
      DO i = j, n
         wa1(i) = wa1(i) + r(l)*temp
         l = l + 1
      END DO
      wa1(j) = wa1(j)/diag(j)
   END DO
!
!  Calculate the norm of the scaled gradient.
!  Test for the special case in which the scaled gradient is zero.
!
   gnorm = enorm(n, wa1)
   sgnorm = 0.0D+00
   alpha = delta/qnorm

   IF (gnorm /= 0.0D+00) THEN
!
!  Calculate the point along the scaled gradient which minimizes the quadratic.
!
      wa1(1:n) = (wa1(1:n)/gnorm)/diag(1:n)

      l = 1
      DO j = 1, n
         sum2 = 0.0D+00
         DO i = j, n
            sum2 = sum2 + r(l)*wa1(i)
            l = l + 1
         END DO
         wa2(j) = sum2
      END DO

      temp = enorm(n, wa2)
      sgnorm = (gnorm/temp)/temp
!
!  Test whether the scaled gradient direction is acceptable.
!
      alpha = 0.0D+00
!
!  The scaled gradient direction is not acceptable.
!  Calculate the point along the dogleg at which the quadratic is minimized.
!
      IF (sgnorm < delta) THEN

         bnorm = enorm(n, qtb)
         temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta)
         temp = temp - (delta/qnorm)*(sgnorm/delta)**2 &
                + SQRT((temp - (delta/qnorm))**2 &
                       + (1.0D+00 - (delta/qnorm)**2) &
                       *(1.0D+00 - (sgnorm/delta)**2))

         alpha = ((delta/qnorm)*(1.0D+00 - (sgnorm/delta)**2)) &
                 /temp

      END IF

   END IF
!
!  Form appropriate convex combination of the Gauss-Newton
!  direction and the scaled gradient direction.
!
   temp = (1.0D+00 - alpha)*MIN(sgnorm, delta)

   x(1:n) = temp*wa1(1:n) + alpha*x(1:n)

   RETURN
END SUBROUTINE dogleg
FUNCTION enorm(n, x)

!*****************************************************************************80
!
!! ENORM computes the Euclidean norm of a vector.
!
!  Discussion:
!
!    This is an extremely simplified version of the original ENORM
!    routine, which has been renamed to "ENORM2".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the length of the vector.
!
!    Input, real ( kind = 8 ) X(N), the vector whose norm is desired.
!
!    Output, real ( kind = 8 ) ENORM, the Euclidean norm of the vector.
!
   IMPLICIT NONE

   INTEGER(kind=4) n
   REAL(kind=8) x(n)
   REAL(kind=8) enorm

   enorm = SQRT(SUM(x(1:n)**2))

   RETURN
END FUNCTION enorm
FUNCTION enorm2(n, x)

!*****************************************************************************80
!
!! ENORM2 computes the Euclidean norm of a vector.
!
!  Discussion:
!
!    This routine was named ENORM.  It has been renamed "ENORM2",
!    and a simplified routine has been substituted.
!
!    The Euclidean norm is computed by accumulating the sum of
!    squares in three different sums.  The sums of squares for the
!    small and large components are scaled so that no overflows
!    occur.  Non-destructive underflows are permitted.  Underflows
!    and overflows do not occur in the computation of the unscaled
!    sum of squares for the intermediate components.
!
!    The definitions of small, intermediate and large components
!    depend on two constants, RDWARF and RGIANT.  The main
!    restrictions on these constants are that RDWARF^2 not
!    underflow and RGIANT^2 not overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the length of the vector.
!
!    Input, real ( kind = 8 ) X(N), the vector whose norm is desired.
!
!    Output, real ( kind = 8 ) ENORM2, the Euclidean norm of the vector.
!
   IMPLICIT NONE

   INTEGER(kind=4) n

   REAL(kind=8) agiant
   REAL(kind=8) enorm2
   INTEGER(kind=4) i
   REAL(kind=8) rdwarf
   REAL(kind=8) rgiant
   REAL(kind=8) s1
   REAL(kind=8) s2
   REAL(kind=8) s3
   REAL(kind=8) x(n)
   REAL(kind=8) xabs
   REAL(kind=8) x1max
   REAL(kind=8) x3max

   rdwarf = SQRT(TINY(rdwarf))
   rgiant = SQRT(HUGE(rgiant))

   s1 = 0.0D+00
   s2 = 0.0D+00
   s3 = 0.0D+00
   x1max = 0.0D+00
   x3max = 0.0D+00
   agiant = rgiant/REAL(n, kind=8)

   DO i = 1, n

      xabs = ABS(x(i))

      IF (xabs <= rdwarf) THEN

         IF (x3max < xabs) THEN
            s3 = 1.0D+00 + s3*(x3max/xabs)**2
            x3max = xabs
         ELSE IF (xabs /= 0.0D+00) THEN
            s3 = s3 + (xabs/x3max)**2
         END IF

      ELSE IF (agiant <= xabs) THEN

         IF (x1max < xabs) THEN
            s1 = 1.0D+00 + s1*(x1max/xabs)**2
            x1max = xabs
         ELSE
            s1 = s1 + (xabs/x1max)**2
         END IF

      ELSE

         s2 = s2 + xabs**2

      END IF

   END DO
!
!  Calculation of norm.
!
   IF (s1 /= 0.0D+00) THEN

      enorm2 = x1max*SQRT(s1 + (s2/x1max)/x1max)

   ELSE IF (s2 /= 0.0D+00) THEN

      IF (x3max <= s2) THEN
         enorm2 = SQRT(s2*(1.0D+00 + (x3max/s2)*(x3max*s3)))
      ELSE
         enorm2 = SQRT(x3max*((s2/x3max) + (x3max*s3)))
      END IF

   ELSE

      enorm2 = x3max*SQRT(s3)

   END IF

   RETURN
END FUNCTION enorm2
SUBROUTINE fdjac1(fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn, m, prms)

!*****************************************************************************80
!
!! FDJAC1 estimates an N by N jacobian matrix using forward differences.
!
!  Discussion:
!
!    This subroutine computes a forward-difference approximation
!    to the N by N jacobian matrix associated with a specified
!    problem of N functions in N variables. If the jacobian has
!    a banded form, then function evaluations are saved by only
!    approximating the nonzero terms.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!    30 Aug 2017: added `prms` to pass Parameters for constructing `FCN`
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!      subroutine fcn ( n, x, fvec, iflag )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fvec(n)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer ( kind = 4 ) N, the number of functions and variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the jacobian is evaluated.
!
!    Input, real ( kind = 8 ) FVEC(N), the functions evaluated at X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), the N by N approximate
!    jacobian matrix.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of FJAC, which
!    must not be less than N.
!
!    Output, integer ( kind = 4 ) IFLAG, is an error flag returned by FCN.
!    If FCN returns a nonzero value of IFLAG, then this routine returns
!    immediately to the calling program, with the value of IFLAG.
!
!    Input, integer ( kind = 4 ) ML, MU, specify the number of subdiagonals and
!    superdiagonals within the band of the jacobian matrix.  If the
!    jacobian is not banded, set ML and MU to N-1.
!
!    Input, real ( kind = 8 ) EPSFCN, is used in determining a suitable step
!    length for the forward-difference approximation.  This approximation
!    assumes that the relative errors in the functions are of the order of
!    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
!    the relative errors in the functions are of the order of the machine
!    precision.
!
   IMPLICIT NONE

   INTEGER(kind=4) ldfjac
   INTEGER(kind=4) n
   INTEGER(kind=4) m

   REAL(kind=8) eps
   REAL(kind=8) epsfcn
   REAL(kind=8) epsmch
   EXTERNAL fcn
   REAL(kind=8) fjac(ldfjac, n)
   REAL(kind=8) fvec(n)
   REAL(kind=8) h
   INTEGER(kind=4) i
   INTEGER(kind=4) iflag
   INTEGER(kind=4) j
   INTEGER(kind=4) k
   INTEGER(kind=4) ml
   INTEGER(kind=4) msum
   INTEGER(kind=4) mu
   REAL(kind=8) temp
   REAL(kind=8) wa1(n)
   REAL(kind=8) wa2(n)
   REAL(kind=8) x(n)
   REAL(kind=8) prms(m)

   epsmch = EPSILON(epsmch)

   eps = SQRT(MAX(epsfcn, epsmch))
   msum = ml + mu + 1
!
!  Computation of dense approximate jacobian.
!
   IF (n <= msum) THEN

      DO j = 1, n

         temp = x(j)
         h = eps*ABS(temp)
         IF (h == 0.0D+00) THEN
            h = eps
         END IF

         iflag = 1
         x(j) = temp + h
         CALL fcn(n, x, wa1, iflag, m, prms)

         IF (iflag < 0) THEN
            EXIT
         END IF

         x(j) = temp
         fjac(1:n, j) = (wa1(1:n) - fvec(1:n))/h

      END DO

   ELSE
!
!  Computation of banded approximate jacobian.
!
      DO k = 1, msum

         DO j = k, n, msum
            wa2(j) = x(j)
            h = eps*ABS(wa2(j))
            IF (h == 0.0D+00) THEN
               h = eps
            END IF
            x(j) = wa2(j) + h
         END DO

         iflag = 1
         CALL fcn(n, x, wa1, iflag, m, prms)

         IF (iflag < 0) THEN
            EXIT
         END IF

         DO j = k, n, msum

            x(j) = wa2(j)

            h = eps*ABS(wa2(j))
            IF (h == 0.0D+00) THEN
               h = eps
            END IF

            fjac(1:n, j) = 0.0D+00

            DO i = 1, n
               IF (j - mu <= i .AND. i <= j + ml) THEN
                  fjac(i, j) = (wa1(i) - fvec(i))/h
               END IF
            END DO

         END DO

      END DO

   END IF

   RETURN
END SUBROUTINE fdjac1
SUBROUTINE fdjac2(fcn, m, n, x, xdat, ydat, fvec, fjac, ldfjac, iflag, epsfcn)

!*****************************************************************************80
!
!! FDJAC2 estimates an M by N jacobian matrix using forward differences.
!
!  Discussion:
!
!    This subroutine computes a forward-difference approximation
!    to the M by N jacobian matrix associated with a specified
!    problem of M functions in N variables.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!      subroutine fcn ( m, n, x, xdat, ydat, fvec, iflag ) ! xdat/ydat added for AnOHM, TS 20 Jul 2017
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fvec(m)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer ( kind = 4 ) M, is the number of functions.
!
!    Input, integer ( kind = 4 ) N, is the number of variables.
!    N must not exceed M.
!
!    Input, real ( kind = 8 ) X(N), the point where the jacobian is evaluated.
!
!    Input, real ( kind = 8 ) FVEC(M), the functions evaluated at X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), the M by N approximate
!    jacobian matrix.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of FJAC,
!    which must not be less than M.
!
!    Output, integer ( kind = 4 ) IFLAG, is an error flag returned by FCN.
!    If FCN returns a nonzero value of IFLAG, then this routine returns
!    immediately to the calling program, with the value of IFLAG.
!
!    Input, real ( kind = 8 ) EPSFCN, is used in determining a suitable
!    step length for the forward-difference approximation.  This approximation
!    assumes that the relative errors in the functions are of the order of
!    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
!    the relative errors in the functions are of the order of the machine
!    precision.
!
   IMPLICIT NONE

   INTEGER(kind=4) ldfjac
   INTEGER(kind=4) m
   INTEGER(kind=4) n

   REAL(kind=8) eps
   REAL(kind=8) epsfcn
   REAL(kind=8) epsmch
   EXTERNAL fcn
   REAL(kind=8) fjac(ldfjac, n)
   REAL(kind=8) fvec(m), xdat(m), ydat(m)
   REAL(kind=8) h
   ! integer ( kind = 4 ) i
   INTEGER(kind=4) iflag
   INTEGER(kind=4) j
   REAL(kind=8) temp
   REAL(kind=8) wa(m)
   REAL(kind=8) x(n)

   epsmch = EPSILON(epsmch)

   eps = SQRT(MAX(epsfcn, epsmch))

   DO j = 1, n

      temp = x(j)
      h = eps*ABS(temp)
      IF (h == 0.0D+00) THEN
         h = eps
      END IF

      iflag = 1
      x(j) = temp + h
      CALL fcn(m, n, x, xdat, ydat, wa, iflag)

      IF (iflag < 0) THEN
         EXIT
      END IF

      x(j) = temp
      fjac(1:m, j) = (wa(1:m) - fvec(1:m))/h

   END DO

   RETURN
END SUBROUTINE fdjac2

SUBROUTINE hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, &
                 factor, nprint, info, nfev, fjac, ldfjac, r, lr, qtf, m, prms)

!*****************************************************************************80
!
!! HYBRD seeks a zero of N nonlinear equations in N variables.
!
!  Discussion:
!
!    HYBRD finds a zero of a system of N nonlinear functions in N variables
!    by a modification of the Powell hybrid method.  The user must provide a
!    subroutine which calculates the functions.  The jacobian is
!    then calculated by a forward-difference approximation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!      subroutine fcn ( n, x, fvec, iflag )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fvec(n)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer ( kind = 4 ) N, the number of functions and variables.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(N), the functions evaluated at the output X.
!
!    Input, real ( kind = 8 ) XTOL.  Termination occurs when the relative error
!    between two consecutive iterates is at most XTOL.  XTOL should be
!    nonnegative.
!
!    Input, integer ( kind = 4 ) MAXFEV.  Termination occurs when the number of
!    calls to FCN is at least MAXFEV by the end of an iteration.
!
!    Input, integer ( kind = 4 ) ML, MU, specify the number of subdiagonals and
!    superdiagonals within the band of the jacobian matrix.  If the jacobian
!    is not banded, set ML and MU to at least n - 1.
!
!    Input, real ( kind = 8 ) EPSFCN, is used in determining a suitable step
!    length for the forward-difference approximation.  This approximation
!    assumes that the relative errors in the functions are of the order of
!    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
!    the relative errors in the functions are of the order of the machine
!    precision.
!
!    Input/output, real ( kind = 8 ) DIAG(N).  If MODE = 1, then DIAG is set
!    internally.  If MODE = 2, then DIAG must contain positive entries that
!    serve as multiplicative scale factors for the variables.
!
!    Input, integer ( kind = 4 ) MODE, scaling option.
!    1, variables will be scaled internally.
!    2, scaling is specified by the input DIAG vector.
!
!    Input, real ( kind = 8 ) FACTOR, determines the initial step bound.  This
!    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
!    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
!    in the interval (0.1, 100) with 100 the recommended value.
!
!    Input, integer ( kind = 4 ) NPRINT, enables controlled printing of
!    iterates if it is positive.  In this case, FCN is called with IFLAG = 0 at
!    the beginning of the first iteration and every NPRINT iterations thereafter
!    and immediately prior to return, with X and FVEC available
!    for printing.  If NPRINT is not positive, no special calls
!    of FCN with IFLAG = 0 are made.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG.
!    See the description of FCN.
!    Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, relative error between two consecutive iterates is at most XTOL.
!    2, number of calls to FCN has reached or exceeded MAXFEV.
!    3, XTOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    4, iteration is not making good progress, as measured by the improvement
!       from the last five jacobian evaluations.
!    5, iteration is not making good progress, as measured by the improvement
!       from the last ten iterations.
!
!    Output, integer ( kind = 4 ) NFEV, the number of calls to FCN.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an N by N array which contains
!    the orthogonal matrix Q produced by the QR factorization of the final
!    approximate jacobian.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of FJAC.
!    LDFJAC must be at least N.
!
!    Output, real ( kind = 8 ) R(LR), the upper triangular matrix produced by
!    the QR factorization of the final approximate jacobian, stored rowwise.
!
!    Input, integer ( kind = 4 ) LR, the size of the R array, which must be no
!    less than (N*(N+1))/2.
!
!    Output, real ( kind = 8 ) QTF(N), contains the vector Q'*FVEC.
!
   IMPLICIT NONE

   INTEGER(kind=4) ldfjac
   INTEGER(kind=4) lr
   INTEGER(kind=4) n
   INTEGER(kind=4) m

   REAL(kind=8) actred
   REAL(kind=8) delta
   REAL(kind=8) diag(n)
   REAL(kind=8) enorm
   REAL(kind=8) epsfcn
   REAL(kind=8) epsmch
   REAL(kind=8) factor
   EXTERNAL fcn
   REAL(kind=8) fjac(ldfjac, n)
   REAL(kind=8) fnorm
   REAL(kind=8) fnorm1
   REAL(kind=8) fvec(n)
   INTEGER(kind=4) i
   INTEGER(kind=4) iflag
   INTEGER(kind=4) info
   INTEGER(kind=4) iter
   INTEGER(kind=4) iwa(1)
   INTEGER(kind=4) j
   LOGICAL jeval
   INTEGER(kind=4) l
   INTEGER(kind=4) maxfev
   INTEGER(kind=4) ml
   INTEGER(kind=4) mode
   INTEGER(kind=4) msum
   INTEGER(kind=4) mu
   INTEGER(kind=4) ncfail
   INTEGER(kind=4) nslow1
   INTEGER(kind=4) nslow2
   INTEGER(kind=4) ncsuc
   INTEGER(kind=4) nfev
   INTEGER(kind=4) nprint
   LOGICAL pivot
   REAL(kind=8) pnorm
   REAL(kind=8) prered
   REAL(kind=8) qtf(n)
   REAL(kind=8) r(lr)
   REAL(kind=8) ratio
   LOGICAL sing
   REAL(kind=8) sum2
   REAL(kind=8) temp
   REAL(kind=8) wa1(n)
   REAL(kind=8) wa2(n)
   REAL(kind=8) wa3(n)
   REAL(kind=8) wa4(n)
   REAL(kind=8) x(n)
   REAL(kind=8) xnorm
   REAL(kind=8) xtol
   REAL(kind=8) prms(m)

   epsmch = EPSILON(epsmch)

   info = 0
   iflag = 0
   nfev = 0
!
!  Check the input parameters for errors.
!
   IF (n <= 0) THEN
      RETURN
   ELSE IF (xtol < 0.0D+00) THEN
      RETURN
   ELSE IF (maxfev <= 0) THEN
      RETURN
   ELSE IF (ml < 0) THEN
      RETURN
   ELSE IF (mu < 0) THEN
      RETURN
   ELSE IF (factor <= 0.0D+00) THEN
      RETURN
   ELSE IF (ldfjac < n) THEN
      RETURN
   ELSE IF (lr < (n*(n + 1))/2) THEN
      RETURN
   END IF

   IF (mode == 2) THEN

      DO j = 1, n
         IF (diag(j) <= 0.0D+00) THEN
            go to 300
         END IF
      END DO

   END IF
!
!  Evaluate the function at the starting point
!  and calculate its norm.
!
   iflag = 1
   CALL fcn(n, x, fvec, iflag, m, prms)
   nfev = 1

   IF (iflag < 0) THEN
      go to 300
   END IF

   fnorm = enorm(n, fvec)
!
!  Determine the number of calls to FCN needed to compute the jacobian matrix.
!
   msum = MIN(ml + mu + 1, n)
!
!  Initialize iteration counter and monitors.
!
   iter = 1
   ncsuc = 0
   ncfail = 0
   nslow1 = 0
   nslow2 = 0
!
!  Beginning of the outer loop.
!
30 CONTINUE

   jeval = .TRUE.
!
!  Calculate the jacobian matrix.
!
   iflag = 2
   CALL fdjac1(fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn, m, prms)

   nfev = nfev + msum

   IF (iflag < 0) THEN
      go to 300
   END IF
!
!  Compute the QR factorization of the jacobian.
!
   pivot = .FALSE.
   CALL qrfac(n, n, fjac, ldfjac, pivot, iwa, 1, wa1, wa2)
!
!  On the first iteration, if MODE is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
   IF (iter == 1) THEN

      IF (mode /= 2) THEN

         diag(1:n) = wa2(1:n)
         DO j = 1, n
            IF (wa2(j) == 0.0D+00) THEN
               diag(j) = 1.0D+00
            END IF
         END DO

      END IF
!
!  On the first iteration, calculate the norm of the scaled X
!  and initialize the step bound DELTA.
!
      wa3(1:n) = diag(1:n)*x(1:n)
      xnorm = enorm(n, wa3)
      delta = factor*xnorm
      IF (delta == 0.0D+00) THEN
         delta = factor
      END IF

   END IF
!
!  Form Q' * FVEC and store in QTF.
!
   qtf(1:n) = fvec(1:n)

   DO j = 1, n

      IF (fjac(j, j) /= 0.0D+00) THEN
         temp = -DOT_PRODUCT(qtf(j:n), fjac(j:n, j))/fjac(j, j)
         qtf(j:n) = qtf(j:n) + fjac(j:n, j)*temp
      END IF

   END DO
!
!  Copy the triangular factor of the QR factorization into R.
!
   sing = .FALSE.

   DO j = 1, n
      l = j
      DO i = 1, j - 1
         r(l) = fjac(i, j)
         l = l + n - i
      END DO
      r(l) = wa1(j)
      IF (wa1(j) == 0.0D+00) THEN
         sing = .TRUE.
      END IF
   END DO
!
!  Accumulate the orthogonal factor in FJAC.
!
   CALL qform(n, n, fjac, ldfjac)
!
!  Rescale if necessary.
!
   IF (mode /= 2) THEN
      DO j = 1, n
         diag(j) = MAX(diag(j), wa2(j))
      END DO
   END IF
!
!  Beginning of the inner loop.
!
180 CONTINUE
!
!  If requested, call FCN to enable printing of iterates.
!
   IF (0 < nprint) THEN
      iflag = 0
      IF (MOD(iter - 1, nprint) == 0) THEN
         CALL fcn(n, x, fvec, iflag, m, prms)
      END IF
      IF (iflag < 0) THEN
         go to 300
      END IF
   END IF
!
!  Determine the direction P.
!
   CALL dogleg(n, r, lr, diag, qtf, delta, wa1)
!
!  Store the direction P and X + P.
!  Calculate the norm of P.
!
   wa1(1:n) = -wa1(1:n)
   wa2(1:n) = x(1:n) + wa1(1:n)
   wa3(1:n) = diag(1:n)*wa1(1:n)

   pnorm = enorm(n, wa3)
!
!  On the first iteration, adjust the initial step bound.
!
   IF (iter == 1) THEN
      delta = MIN(delta, pnorm)
   END IF
!
!  Evaluate the function at X + P and calculate its norm.
!
   iflag = 1
   CALL fcn(n, wa2, wa4, iflag, m, prms)
   nfev = nfev + 1

   IF (iflag < 0) THEN
      go to 300
   END IF

   fnorm1 = enorm(n, wa4)
!
!  Compute the scaled actual reduction.
!
   actred = -1.0D+00
   IF (fnorm1 < fnorm) THEN
      actred = 1.0D+00 - (fnorm1/fnorm)**2
   END IF
!
!  Compute the scaled predicted reduction.
!
   l = 1
   DO i = 1, n
      sum2 = 0.0D+00
      DO j = i, n
         sum2 = sum2 + r(l)*wa1(j)
         l = l + 1
      END DO
      wa3(i) = qtf(i) + sum2
   END DO

   temp = enorm(n, wa3)
   prered = 0.0D+00
   IF (temp < fnorm) THEN
      prered = 1.0D+00 - (temp/fnorm)**2
   END IF
!
!  Compute the ratio of the actual to the predicted reduction.
!
   ratio = 0.0D+00
   IF (0.0D+00 < prered) THEN
      ratio = actred/prered
   END IF
!
!  Update the step bound.
!
   IF (ratio < 0.1D+00) THEN

      ncsuc = 0
      ncfail = ncfail + 1
      delta = 0.5D+00*delta

   ELSE

      ncfail = 0
      ncsuc = ncsuc + 1

      IF (0.5D+00 <= ratio .OR. 1 < ncsuc) THEN
         delta = MAX(delta, pnorm/0.5D+00)
      END IF

      IF (ABS(ratio - 1.0D+00) <= 0.1D+00) THEN
         delta = pnorm/0.5D+00
      END IF

   END IF
!
!  Test for successful iteration.
!
!  Successful iteration.
!  Update X, FVEC, and their norms.
!
   IF (0.0001D+00 <= ratio) THEN
      x(1:n) = wa2(1:n)
      wa2(1:n) = diag(1:n)*x(1:n)
      fvec(1:n) = wa4(1:n)
      xnorm = enorm(n, wa2)
      fnorm = fnorm1
      iter = iter + 1
   END IF
!
!  Determine the progress of the iteration.
!
   nslow1 = nslow1 + 1
   IF (0.001D+00 <= actred) THEN
      nslow1 = 0
   END IF

   IF (jeval) THEN
      nslow2 = nslow2 + 1
   END IF

   IF (0.1D+00 <= actred) THEN
      nslow2 = 0
   END IF
!
!  Test for convergence.
!
   IF (delta <= xtol*xnorm .OR. fnorm == 0.0D+00) THEN
      info = 1
   END IF

   IF (info /= 0) THEN
      go to 300
   END IF
!
!  Tests for termination and stringent tolerances.
!
   IF (maxfev <= nfev) THEN
      info = 2
   END IF

   IF (0.1D+00*MAX(0.1D+00*delta, pnorm) <= epsmch*xnorm) THEN
      info = 3
   END IF

   IF (nslow2 == 5) THEN
      info = 4
   END IF

   IF (nslow1 == 10) THEN
      info = 5
   END IF

   IF (info /= 0) THEN
      go to 300
   END IF
!
!  Criterion for recalculating jacobian approximation
!  by forward differences.
!
   IF (ncfail == 2) THEN
      go to 290
   END IF
!
!  Calculate the rank one modification to the jacobian
!  and update QTF if necessary.
!
   DO j = 1, n
      sum2 = DOT_PRODUCT(wa4(1:n), fjac(1:n, j))
      wa2(j) = (sum2 - wa3(j))/pnorm
      wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
      IF (0.0001D+00 <= ratio) THEN
         qtf(j) = sum2
      END IF
   END DO
!
!  Compute the QR factorization of the updated jacobian.
!
   CALL r1updt(n, n, r, lr, wa1, wa2, wa3, sing)
   CALL r1mpyq(n, n, fjac, ldfjac, wa2, wa3)
   CALL r1mpyq(1, n, qtf, 1, wa2, wa3)
!
!  End of the inner loop.
!
   jeval = .FALSE.
   go to 180

290 CONTINUE
!
!  End of the outer loop.
!
   go to 30

300 CONTINUE
!
!  Termination, either normal or user imposed.
!
   IF (iflag < 0) THEN
      info = iflag
   END IF

   iflag = 0

   IF (0 < nprint) THEN
      CALL fcn(n, x, fvec, iflag, m, prms)
   END IF

   RETURN
END SUBROUTINE hybrd

SUBROUTINE hybrd1(fcn, n, x, fvec, tol, info, m, prms)

!*****************************************************************************80
!
!! HYBRD1 seeks a zero of N nonlinear equations in N variables.
!
!  Discussion:
!
!    HYBRD1 finds a zero of a system of N nonlinear functions in N variables
!    by a modification of the Powell hybrid method.  This is done by using the
!    more general nonlinear equation solver HYBRD.  The user must provide a
!    subroutine which calculates the functions.  The jacobian is then
!    calculated by a forward-difference approximation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2016
!    30 Aug 2017: added `prms` to pass paramters for constructing `fcn`
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!      subroutine fcn ( n, x, fvec, iflag )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fvec(n)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!      integer ( kind = 4 ) m
!      real ( kind = 8 ) prm(m)
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer ( kind = 4 ) N, the number of functions and variables.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Input, real ( kind = 8 ) FVEC(N), the functions evaluated at the output X.
!
!    Input, integer ( kind = 4 ) M, the number of parameters for constructing FCN.
!
!    Input, real ( kind = 8 ) PRMS(M), static paramters for constructing FCN.
!
!    Input, real ( kind = 8 ) TOL.  Termination occurs when the algorithm
!    estimates that the relative error between X and the solution is at
!    most TOL.  TOL should be nonnegative.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See the
!    description of FCN.
!    Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    2, number of calls to FCN has reached or exceeded 200*(N+1).
!    3, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    4, the iteration is not making good progress.
!
   IMPLICIT NONE

   ! integer ( kind = 4 ) lwa
   INTEGER(kind=4) n
   INTEGER(kind=4) m

   REAL(kind=8) diag(n)
   REAL(kind=8) epsfcn
   REAL(kind=8) factor
   EXTERNAL fcn
   REAL(kind=8) fjac(n, n)
   REAL(kind=8) fvec(n)
   INTEGER(kind=4) info
   ! integer ( kind = 4 ) j
   INTEGER(kind=4) ldfjac
   INTEGER(kind=4) lr
   INTEGER(kind=4) maxfev
   INTEGER(kind=4) ml
   INTEGER(kind=4) mode
   INTEGER(kind=4) mu
   INTEGER(kind=4) nfev
   INTEGER(kind=4) nprint
   REAL(kind=8) qtf(n)
   REAL(kind=8) r((n*(n + 1))/2)
   REAL(kind=8) tol
   REAL(kind=8) x(n)
   REAL(kind=8) xtol
   REAL(kind=8) prms(m)

   IF (n <= 0) THEN
      info = 0
      RETURN
   END IF

   IF (tol < 0.0D+00) THEN
      info = 0
      RETURN
   END IF

   xtol = tol
   maxfev = 200*(n + 1)
   ml = n - 1
   mu = n - 1
   epsfcn = 0.0D+00
   diag(1:n) = 1.0D+00
   mode = 2
   factor = 100.0D+00
   nprint = 0
   info = 0
   nfev = 0
   fjac(1:n, 1:n) = 0.0D+00
   ldfjac = n
   r(1:(n*(n + 1))/2) = 0.0D+00
   lr = (n*(n + 1))/2
   qtf(1:n) = 0.0D+00

   CALL hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, &
              factor, nprint, info, nfev, fjac, ldfjac, r, lr, qtf, m, prms)

   IF (info == 5) THEN
      info = 4
   END IF

   RETURN
END SUBROUTINE hybrd1
SUBROUTINE hybrj(fcn, n, x, fvec, fjac, ldfjac, xtol, maxfev, diag, mode, &
                 factor, nprint, info, nfev, njev, r, lr, qtf)

!*****************************************************************************80
!
!! HYBRJ seeks a zero of N nonlinear equations in N variables.
!
!  Discussion:
!
!    HYBRJ finds a zero of a system of N nonlinear functions in N variables
!    by a modification of the Powell hybrid method.  The user must provide a
!    subroutine which calculates the functions and the jacobian.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the jacobian.  FCN should have the form:
!
!      subroutine fcn ( n, x, fvec, fjac, ldfjac, iflag )
!      integer ( kind = 4 ) ldfjac
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fjac(ldfjac,n)
!      real ( kind = 8 ) fvec(n)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!    return this matrix in FJAC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer ( kind = 4 ) N, the number of functions and variables.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(N), the functions evaluated at the output X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an N by N matrix, containing
!    the orthogonal matrix Q produced by the QR factorization
!    of the final approximate jacobian.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of the
!    array FJAC.  LDFJAC must be at least N.
!
!    Input, real ( kind = 8 ) XTOL.  Termination occurs when the relative error
!    between two consecutive iterates is at most XTOL.  XTOL should be
!    nonnegative.
!
!    Input, integer ( kind = 4 ) MAXFEV.  Termination occurs when the number of
!    calls to FCN is at least MAXFEV by the end of an iteration.
!
!    Input/output, real ( kind = 8 ) DIAG(N).  If MODE = 1, then DIAG is set
!    internally.  If MODE = 2, then DIAG must contain positive entries that
!    serve as multiplicative scale factors for the variables.
!
!    Input, integer ( kind = 4 ) MODE, scaling option.
!    1, variables will be scaled internally.
!    2, scaling is specified by the input DIAG vector.
!
!    Input, real ( kind = 8 ) FACTOR, determines the initial step bound.  This
!    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
!    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
!    in the interval (0.1, 100) with 100 the recommended value.
!
!    Input, integer ( kind = 4 ) NPRINT, enables controlled printing of iterates
!    if it is positive.  In this case, FCN is called with IFLAG = 0 at the
!    beginning of the first iteration and every NPRINT iterations thereafter
!    and immediately prior to return, with X and FVEC available
!    for printing.  If NPRINT is not positive, no special calls
!    of FCN with IFLAG = 0 are made.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG.
!    See the description of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, relative error between two consecutive iterates is at most XTOL.
!    2, number of calls to FCN with IFLAG = 1 has reached MAXFEV.
!    3, XTOL is too small.  No further improvement in
!       the approximate solution X is possible.
!    4, iteration is not making good progress, as measured by the
!       improvement from the last five jacobian evaluations.
!    5, iteration is not making good progress, as measured by the
!       improvement from the last ten iterations.
!
!    Output, integer ( kind = 4 ) NFEV, the number of calls to FCN
!    with IFLAG = 1.
!
!    Output, integer ( kind = 4 ) NJEV, the number of calls to FCN
!    with IFLAG = 2.
!
!    Output, real ( kind = 8 ) R(LR), the upper triangular matrix produced
!    by the QR factorization of the final approximate jacobian, stored rowwise.
!
!    Input, integer ( kind = 4 ) LR, the size of the R array, which must
!    be no less than (N*(N+1))/2.
!
!    Output, real ( kind = 8 ) QTF(N), contains the vector Q'*FVEC.
!
   IMPLICIT NONE

   INTEGER(kind=4) ldfjac
   INTEGER(kind=4) lr
   INTEGER(kind=4) n

   REAL(kind=8) actred
   REAL(kind=8) delta
   REAL(kind=8) diag(n)
   REAL(kind=8) enorm
   REAL(kind=8) epsmch
   REAL(kind=8) factor
   EXTERNAL fcn
   REAL(kind=8) fjac(ldfjac, n)
   REAL(kind=8) fnorm
   REAL(kind=8) fnorm1
   REAL(kind=8) fvec(n)
   INTEGER(kind=4) i
   INTEGER(kind=4) iflag
   INTEGER(kind=4) info
   INTEGER(kind=4) iter
   INTEGER(kind=4) iwa(1)
   INTEGER(kind=4) j
   LOGICAL jeval
   INTEGER(kind=4) l
   INTEGER(kind=4) maxfev
   INTEGER(kind=4) mode
   INTEGER(kind=4) ncfail
   INTEGER(kind=4) nslow1
   INTEGER(kind=4) nslow2
   INTEGER(kind=4) ncsuc
   INTEGER(kind=4) nfev
   INTEGER(kind=4) njev
   INTEGER(kind=4) nprint
   LOGICAL pivot
   REAL(kind=8) pnorm
   REAL(kind=8) prered
   REAL(kind=8) qtf(n)
   REAL(kind=8) r(lr)
   REAL(kind=8) ratio
   LOGICAL sing
   REAL(kind=8) sum2
   REAL(kind=8) temp
   REAL(kind=8) wa1(n)
   REAL(kind=8) wa2(n)
   REAL(kind=8) wa3(n)
   REAL(kind=8) wa4(n)
   REAL(kind=8) x(n)
   REAL(kind=8) xnorm
   REAL(kind=8) xtol

   epsmch = EPSILON(epsmch)

   info = 0
   iflag = 0
   nfev = 0
   njev = 0
!
!  Check the input parameters for errors.
!
   IF (n <= 0) THEN
      IF (iflag < 0) THEN
         info = iflag
      END IF
      iflag = 0
      IF (0 < nprint) THEN
         CALL fcn(n, x, fvec, fjac, ldfjac, iflag)
      END IF
      RETURN
   END IF

   IF (ldfjac < n .OR. &
       xtol < 0.0D+00 .OR. &
       maxfev <= 0 .OR. &
       factor <= 0.0D+00 .OR. &
       lr < (n*(n + 1))/2) THEN
      IF (iflag < 0) THEN
         info = iflag
      END IF
      iflag = 0
      IF (0 < nprint) THEN
         CALL fcn(n, x, fvec, fjac, ldfjac, iflag)
      END IF
      RETURN
   END IF

   IF (mode == 2) THEN
      DO j = 1, n
         IF (diag(j) <= 0.0D+00) THEN
            IF (iflag < 0) THEN
               info = iflag
            END IF
            iflag = 0
            IF (0 < nprint) THEN
               CALL fcn(n, x, fvec, fjac, ldfjac, iflag)
            END IF
            RETURN
         END IF
      END DO
   END IF
!
!  Evaluate the function at the starting point and calculate its norm.
!
   iflag = 1
   CALL fcn(n, x, fvec, fjac, ldfjac, iflag)
   nfev = 1

   IF (iflag < 0) THEN
      IF (iflag < 0) THEN
         info = iflag
      END IF
      iflag = 0
      IF (0 < nprint) THEN
         CALL fcn(n, x, fvec, fjac, ldfjac, iflag)
      END IF
      RETURN
   END IF

   fnorm = enorm(n, fvec)
!
!  Initialize iteration counter and monitors.
!
   iter = 1
   ncsuc = 0
   ncfail = 0
   nslow1 = 0
   nslow2 = 0
!
!  Beginning of the outer loop.
!
   DO

      jeval = .TRUE.
!
!  Calculate the jacobian matrix.
!
      iflag = 2
      CALL fcn(n, x, fvec, fjac, ldfjac, iflag)
      njev = njev + 1

      IF (iflag < 0) THEN
         info = iflag
         iflag = 0
         IF (0 < nprint) THEN
            CALL fcn(n, x, fvec, fjac, ldfjac, iflag)
         END IF
         RETURN
      END IF
!
!  Compute the QR factorization of the jacobian.
!
      pivot = .FALSE.
      CALL qrfac(n, n, fjac, ldfjac, pivot, iwa, 1, wa1, wa2)
!
!  On the first iteration, if MODE is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
      IF (iter == 1) THEN

         IF (mode /= 2) THEN
            diag(1:n) = wa2(1:n)
            DO j = 1, n
               IF (wa2(j) == 0.0D+00) THEN
                  diag(j) = 1.0D+00
               END IF
            END DO
         END IF
!
!  On the first iteration, calculate the norm of the scaled X
!  and initialize the step bound DELTA.
!
         wa3(1:n) = diag(1:n)*x(1:n)
         xnorm = enorm(n, wa3)
         delta = factor*xnorm
         IF (delta == 0.0D+00) THEN
            delta = factor
         END IF

      END IF
!
!  Form Q'*FVEC and store in QTF.
!
      qtf(1:n) = fvec(1:n)

      DO j = 1, n
         IF (fjac(j, j) /= 0.0D+00) THEN
            sum2 = 0.0D+00
            DO i = j, n
               sum2 = sum2 + fjac(i, j)*qtf(i)
            END DO
            temp = -sum2/fjac(j, j)
            DO i = j, n
               qtf(i) = qtf(i) + fjac(i, j)*temp
            END DO
         END IF
      END DO
!
!  Copy the triangular factor of the QR factorization into R.
!
      sing = .FALSE.

      DO j = 1, n
         l = j
         DO i = 1, j - 1
            r(l) = fjac(i, j)
            l = l + n - i
         END DO
         r(l) = wa1(j)
         IF (wa1(j) == 0.0D+00) THEN
            sing = .TRUE.
         END IF
      END DO
!
!  Accumulate the orthogonal factor in FJAC.
!
      CALL qform(n, n, fjac, ldfjac)
!
!  Rescale if necessary.
!
      IF (mode /= 2) THEN
         DO j = 1, n
            diag(j) = MAX(diag(j), wa2(j))
         END DO
      END IF
!
!  Beginning of the inner loop.
!
      DO
!
!  If requested, call FCN to enable printing of iterates.
!
         IF (0 < nprint) THEN

            iflag = 0
            IF (MOD(iter - 1, nprint) == 0) THEN
               CALL fcn(n, x, fvec, fjac, ldfjac, iflag)
            END IF

            IF (iflag < 0) THEN
               info = iflag
               iflag = 0
               IF (0 < nprint) THEN
                  CALL fcn(n, x, fvec, fjac, ldfjac, iflag)
               END IF
               RETURN
            END IF

         END IF
!
!  Determine the direction P.
!
         CALL dogleg(n, r, lr, diag, qtf, delta, wa1)
!
!  Store the direction P and X + P.
!  Calculate the norm of P.
!
         wa1(1:n) = -wa1(1:n)
         wa2(1:n) = x(1:n) + wa1(1:n)
         wa3(1:n) = diag(1:n)*wa1(1:n)
         pnorm = enorm(n, wa3)
!
!  On the first iteration, adjust the initial step bound.
!
         IF (iter == 1) THEN
            delta = MIN(delta, pnorm)
         END IF
!
!  Evaluate the function at X + P and calculate its norm.
!
         iflag = 1
         CALL fcn(n, wa2, wa4, fjac, ldfjac, iflag)
         nfev = nfev + 1

         IF (iflag < 0) THEN
            info = iflag
            iflag = 0
            IF (0 < nprint) THEN
               CALL fcn(n, x, fvec, fjac, ldfjac, iflag)
            END IF
            RETURN
         END IF

         fnorm1 = enorm(n, wa4)
!
!  Compute the scaled actual reduction.
!
         actred = -1.0D+00
         IF (fnorm1 < fnorm) THEN
            actred = 1.0D+00 - (fnorm1/fnorm)**2
         END IF
!
!  Compute the scaled predicted reduction.
!
         l = 1
         DO i = 1, n
            sum2 = 0.0D+00
            DO j = i, n
               sum2 = sum2 + r(l)*wa1(j)
               l = l + 1
            END DO
            wa3(i) = qtf(i) + sum2
         END DO

         temp = enorm(n, wa3)
         prered = 0.0D+00
         IF (temp < fnorm) THEN
            prered = 1.0D+00 - (temp/fnorm)**2
         END IF
!
!  Compute the ratio of the actual to the predicted reduction.
!
         IF (0.0D+00 < prered) THEN
            ratio = actred/prered
         ELSE
            ratio = 0.0D+00
         END IF
!
!  Update the step bound.
!
         IF (ratio < 0.1D+00) THEN

            ncsuc = 0
            ncfail = ncfail + 1
            delta = 0.5D+00*delta

         ELSE

            ncfail = 0
            ncsuc = ncsuc + 1

            IF (0.5D+00 <= ratio .OR. 1 < ncsuc) THEN
               delta = MAX(delta, pnorm/0.5D+00)
            END IF

            IF (ABS(ratio - 1.0D+00) <= 0.1D+00) THEN
               delta = pnorm/0.5D+00
            END IF

         END IF
!
!  Test for successful iteration.
!

!
!  Successful iteration.
!  Update X, FVEC, and their norms.
!
         IF (0.0001D+00 <= ratio) THEN
            x(1:n) = wa2(1:n)
            wa2(1:n) = diag(1:n)*x(1:n)
            fvec(1:n) = wa4(1:n)
            xnorm = enorm(n, wa2)
            fnorm = fnorm1
            iter = iter + 1
         END IF
!
!  Determine the progress of the iteration.
!
         nslow1 = nslow1 + 1
         IF (0.001D+00 <= actred) THEN
            nslow1 = 0
         END IF

         IF (jeval) THEN
            nslow2 = nslow2 + 1
         END IF

         IF (0.1D+00 <= actred) THEN
            nslow2 = 0
         END IF
!
!  Test for convergence.
!
         IF (delta <= xtol*xnorm .OR. fnorm == 0.0D+00) THEN
            info = 1
         END IF

         IF (info /= 0) THEN
            iflag = 0
            IF (0 < nprint) THEN
               CALL fcn(n, x, fvec, fjac, ldfjac, iflag)
            END IF
            RETURN
         END IF
!
!  Tests for termination and stringent tolerances.
!
         IF (maxfev <= nfev) THEN
            info = 2
         END IF

         IF (0.1D+00*MAX(0.1D+00*delta, pnorm) <= epsmch*xnorm) THEN
            info = 3
         END IF

         IF (nslow2 == 5) THEN
            info = 4
         END IF

         IF (nslow1 == 10) THEN
            info = 5
         END IF

         IF (info /= 0) THEN
            iflag = 0
            IF (0 < nprint) THEN
               CALL fcn(n, x, fvec, fjac, ldfjac, iflag)
            END IF
            RETURN
         END IF
!
!  Criterion for recalculating jacobian.
!
         IF (ncfail == 2) THEN
            EXIT
         END IF
!
!  Calculate the rank one modification to the jacobian
!  and update QTF if necessary.
!
         DO j = 1, n
            sum2 = DOT_PRODUCT(wa4(1:n), fjac(1:n, j))
            wa2(j) = (sum2 - wa3(j))/pnorm
            wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
            IF (0.0001D+00 <= ratio) THEN
               qtf(j) = sum2
            END IF
         END DO
!
!  Compute the QR factorization of the updated jacobian.
!
         CALL r1updt(n, n, r, lr, wa1, wa2, wa3, sing)
         CALL r1mpyq(n, n, fjac, ldfjac, wa2, wa3)
         CALL r1mpyq(1, n, qtf, 1, wa2, wa3)
!
!  End of the inner loop.
!
         jeval = .FALSE.

      END DO
!
!  End of the outer loop.
!
   END DO

END SUBROUTINE hybrj
SUBROUTINE hybrj1(fcn, n, x, fvec, fjac, ldfjac, tol, info)

!*****************************************************************************80
!
!! HYBRJ1 seeks a zero of N equations in N variables by Powell's method.
!
!  Discussion:
!
!    HYBRJ1 finds a zero of a system of N nonlinear functions in N variables
!    by a modification of the Powell hybrid method.  This is done by using the
!    more general nonlinear equation solver HYBRJ.  The user
!    must provide a subroutine which calculates the functions
!    and the jacobian.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the jacobian.  FCN should have the form:
!      subroutine fcn ( n, x, fvec, fjac, ldfjac, iflag )
!      integer ( kind = 4 ) ldfjac
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fjac(ldfjac,n)
!      real ( kind = 8 ) fvec(n)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!    return this matrix in FJAC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer ( kind = 4 ) N, the number of functions and variables.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(N), the functions evaluated at the output X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an N by N array which contains
!    the orthogonal matrix Q produced by the QR factorization of the final
!    approximate jacobian.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of  FJAC.
!    LDFJAC must be at least N.
!
!    Input, real ( kind = 8 ) TOL.  Termination occurs when the algorithm
!    estimates that the relative error between X and the solution is at most
!    TOL.  TOL should be nonnegative.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    2, number of calls to FCN with IFLAG = 1 has reached 100*(N+1).
!    3, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    4, iteration is not making good progress.
!
   IMPLICIT NONE

   INTEGER(kind=4) ldfjac
   INTEGER(kind=4) n

   REAL(kind=8) diag(n)
   REAL(kind=8) factor
   EXTERNAL fcn
   REAL(kind=8) fjac(ldfjac, n)
   REAL(kind=8) fvec(n)
   INTEGER(kind=4) info
   ! integer ( kind = 4 ) j
   INTEGER(kind=4) lr
   INTEGER(kind=4) maxfev
   INTEGER(kind=4) mode
   INTEGER(kind=4) nfev
   INTEGER(kind=4) njev
   INTEGER(kind=4) nprint
   REAL(kind=8) qtf(n)
   REAL(kind=8) r((n*(n + 1))/2)
   REAL(kind=8) tol
   REAL(kind=8) x(n)
   REAL(kind=8) xtol

   info = 0

   IF (n <= 0) THEN
      RETURN
   ELSE IF (ldfjac < n) THEN
      RETURN
   ELSE IF (tol < 0.0D+00) THEN
      RETURN
   END IF

   maxfev = 100*(n + 1)
   xtol = tol
   mode = 2
   diag(1:n) = 1.0D+00
   factor = 100.0D+00
   nprint = 0
   lr = (n*(n + 1))/2

   CALL hybrj(fcn, n, x, fvec, fjac, ldfjac, xtol, maxfev, diag, mode, &
              factor, nprint, info, nfev, njev, r, lr, qtf)

   IF (info == 5) THEN
      info = 4
   END IF

   RETURN
END SUBROUTINE hybrj1
SUBROUTINE lmder(fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
                 diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf)

!*****************************************************************************80
!
!! LMDER minimizes M functions in N variables by the Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMDER minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    The user must provide a subroutine which calculates the functions
!    and the jacobian.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the jacobian.  FCN should have the form:
!      subroutine fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
!      integer ( kind = 4 ) ldfjac
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fjac(ldfjac,n)
!      real ( kind = 8 ) fvec(m)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!    return this matrix in FJAC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer ( kind = 4 ) M, is the number of functions.
!
!    Input, integer ( kind = 4 ) N, is the number of variables.
!    N must not exceed M.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(M), the functions evaluated at the output X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an M by N array.  The upper
!    N by N submatrix of FJAC contains an upper triangular matrix R with
!    diagonal elements of nonincreasing magnitude such that
!      P' * ( JAC' * JAC ) * P = R' * R,
!    where P is a permutation matrix and JAC is the final calculated jacobian.
!    Column J of P is column IPVT(J) of the identity matrix.  The lower
!    trapezoidal part of FJAC contains information generated during
!    the computation of R.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of FJAC.
!    LDFJAC must be at least M.
!
!    Input, real ( kind = 8 ) FTOL.  Termination occurs when both the actual
!    and predicted relative reductions in the sum of squares are at most FTOL.
!    Therefore, FTOL measures the relative error desired in the sum of
!    squares.  FTOL should be nonnegative.
!
!    Input, real ( kind = 8 ) XTOL.  Termination occurs when the relative error
!    between two consecutive iterates is at most XTOL.  XTOL should be
!    nonnegative.
!
!    Input, real ( kind = 8 ) GTOL.  Termination occurs when the cosine of the
!    angle between FVEC and any column of the jacobian is at most GTOL in
!    absolute value.  Therefore, GTOL measures the orthogonality desired
!    between the function vector and the columns of the jacobian.  GTOL should
!    be nonnegative.
!
!    Input, integer ( kind = 4 ) MAXFEV.  Termination occurs when the number of
!    calls to FCN with IFLAG = 1 is at least MAXFEV by the end of an iteration.
!
!    Input/output, real ( kind = 8 ) DIAG(N).  If MODE = 1, then DIAG is set
!    internally.  If MODE = 2, then DIAG must contain positive entries that
!    serve as multiplicative scale factors for the variables.
!
!    Input, integer ( kind = 4 ) MODE, scaling option.
!    1, variables will be scaled internally.
!    2, scaling is specified by the input DIAG vector.
!
!    Input, real ( kind = 8 ) FACTOR, determines the initial step bound.  This
!    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
!    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
!    in the interval (0.1, 100) with 100 the recommended value.
!
!    Input, integer ( kind = 4 ) NPRINT, enables controlled printing of iterates
!    if it is positive.  In this case, FCN is called with IFLAG = 0 at the
!    beginning of the first iteration and every NPRINT iterations thereafter
!    and immediately prior to return, with X and FVEC available
!    for printing.  If NPRINT is not positive, no special calls
!    of FCN with IFLAG = 0 are made.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, both actual and predicted relative reductions in the sum of
!       squares are at most FTOL.
!    2, relative error between two consecutive iterates is at most XTOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, the cosine of the angle between FVEC and any column of the jacobian
!       is at most GTOL in absolute value.
!    5, number of calls to FCN with IFLAG = 1 has reached MAXFEV.
!    6, FTOL is too small.  No further reduction in the sum of squares
!       is possible.
!    7, XTOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    8, GTOL is too small.  FVEC is orthogonal to the columns of the
!       jacobian to machine precision.
!
!    Output, integer ( kind = 4 ) NFEV, the number of calls to FCN with
!    IFLAG = 1.
!
!    Output, integer ( kind = 4 ) NJEV, the number of calls to FCN with
!    IFLAG = 2.
!
!    Output, integer ( kind = 4 ) IPVT(N), defines a permutation matrix P
!    such that JAC*P = Q*R, where JAC is the final calculated jacobian, Q is
!    orthogonal (not stored), and R is upper triangular with diagonal
!    elements of nonincreasing magnitude.  Column J of P is column
!    IPVT(J) of the identity matrix.
!
!    Output, real ( kind = 8 ) QTF(N), contains the first N elements of Q'*FVEC.
!
   IMPLICIT NONE

   INTEGER(kind=4) ldfjac
   INTEGER(kind=4) m
   INTEGER(kind=4) n

   REAL(kind=8) actred
   REAL(kind=8) delta
   REAL(kind=8) diag(n)
   REAL(kind=8) dirder
   REAL(kind=8) enorm
   REAL(kind=8) epsmch
   REAL(kind=8) factor
   EXTERNAL fcn
   REAL(kind=8) fjac(ldfjac, n)
   REAL(kind=8) fnorm
   REAL(kind=8) fnorm1
   REAL(kind=8) ftol
   REAL(kind=8) fvec(m)
   REAL(kind=8) gnorm
   REAL(kind=8) gtol
   ! integer ( kind = 4 ) i
   INTEGER(kind=4) iflag
   INTEGER(kind=4) info
   INTEGER(kind=4) ipvt(n)
   INTEGER(kind=4) iter
   INTEGER(kind=4) j
   INTEGER(kind=4) l
   INTEGER(kind=4) maxfev
   INTEGER(kind=4) mode
   INTEGER(kind=4) nfev
   INTEGER(kind=4) njev
   INTEGER(kind=4) nprint
   REAL(kind=8) par
   LOGICAL pivot
   REAL(kind=8) pnorm
   REAL(kind=8) prered
   REAL(kind=8) qtf(n)
   REAL(kind=8) ratio
   REAL(kind=8) sum2
   REAL(kind=8) temp
   REAL(kind=8) temp1
   REAL(kind=8) temp2
   REAL(kind=8) wa1(n)
   REAL(kind=8) wa2(n)
   REAL(kind=8) wa3(n)
   REAL(kind=8) wa4(m)
   REAL(kind=8) xnorm
   REAL(kind=8) x(n)
   REAL(kind=8) xtol

   epsmch = EPSILON(epsmch)

   info = 0
   iflag = 0
   nfev = 0
   njev = 0
!
!  Check the input parameters for errors.
!
   IF (n <= 0) THEN
      go to 300
   END IF

   IF (m < n) THEN
      go to 300
   END IF

   IF (ldfjac < m &
       .OR. ftol < 0.0D+00 .OR. xtol < 0.0D+00 .OR. gtol < 0.0D+00 &
       .OR. maxfev <= 0 .OR. factor <= 0.0D+00) THEN
      go to 300
   END IF

   IF (mode == 2) THEN
      DO j = 1, n
         IF (diag(j) <= 0.0D+00) THEN
            go to 300
         END IF
      END DO
   END IF
!
!  Evaluate the function at the starting point and calculate its norm.
!
   iflag = 1
   CALL fcn(m, n, x, fvec, fjac, ldfjac, iflag)
   nfev = 1
   IF (iflag < 0) THEN
      go to 300
   END IF

   fnorm = enorm(m, fvec)
!
!  Initialize Levenberg-Marquardt parameter and iteration counter.
!
   par = 0.0D+00
   iter = 1
!
!  Beginning of the outer loop.
!
30 CONTINUE
!
!  Calculate the jacobian matrix.
!
   iflag = 2
   CALL fcn(m, n, x, fvec, fjac, ldfjac, iflag)

   njev = njev + 1

   IF (iflag < 0) THEN
      go to 300
   END IF
!
!  If requested, call FCN to enable printing of iterates.
!
   IF (0 < nprint) THEN
      iflag = 0
      IF (MOD(iter - 1, nprint) == 0) THEN
         CALL fcn(m, n, x, fvec, fjac, ldfjac, iflag)
      END IF
      IF (iflag < 0) THEN
         go to 300
      END IF
   END IF
!
!  Compute the QR factorization of the jacobian.
!
   pivot = .TRUE.
   CALL qrfac(m, n, fjac, ldfjac, pivot, ipvt, n, wa1, wa2)
!
!  On the first iteration and if mode is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
   IF (iter == 1) THEN

      IF (mode /= 2) THEN
         diag(1:n) = wa2(1:n)
         DO j = 1, n
            IF (wa2(j) == 0.0D+00) THEN
               diag(j) = 1.0D+00
            END IF
         END DO
      END IF
!
!  On the first iteration, calculate the norm of the scaled X
!  and initialize the step bound DELTA.
!
      wa3(1:n) = diag(1:n)*x(1:n)

      xnorm = enorm(n, wa3)
      delta = factor*xnorm
      IF (delta == 0.0D+00) THEN
         delta = factor
      END IF

   END IF
!
!  Form Q'*FVEC and store the first N components in QTF.
!
   wa4(1:m) = fvec(1:m)

   DO j = 1, n

      IF (fjac(j, j) /= 0.0D+00) THEN
         sum2 = DOT_PRODUCT(wa4(j:m), fjac(j:m, j))
         temp = -sum2/fjac(j, j)
         wa4(j:m) = wa4(j:m) + fjac(j:m, j)*temp
      END IF

      fjac(j, j) = wa1(j)
      qtf(j) = wa4(j)

   END DO
!
!  Compute the norm of the scaled gradient.
!
   gnorm = 0.0D+00

   IF (fnorm /= 0.0D+00) THEN

      DO j = 1, n
         l = ipvt(j)
         IF (wa2(l) /= 0.0D+00) THEN
            sum2 = DOT_PRODUCT(qtf(1:j), fjac(1:j, j))/fnorm
            gnorm = MAX(gnorm, ABS(sum2/wa2(l)))
         END IF
      END DO

   END IF
!
!  Test for convergence of the gradient norm.
!
   IF (gnorm <= gtol) THEN
      info = 4
      go to 300
   END IF
!
!  Rescale if necessary.
!
   IF (mode /= 2) THEN
      DO j = 1, n
         diag(j) = MAX(diag(j), wa2(j))
      END DO
   END IF
!
!  Beginning of the inner loop.
!
200 CONTINUE
!
!  Determine the Levenberg-Marquardt parameter.
!
   CALL lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2)
!
!  Store the direction p and x + p. calculate the norm of p.
!
   wa1(1:n) = -wa1(1:n)
   wa2(1:n) = x(1:n) + wa1(1:n)
   wa3(1:n) = diag(1:n)*wa1(1:n)

   pnorm = enorm(n, wa3)
!
!  On the first iteration, adjust the initial step bound.
!
   IF (iter == 1) THEN
      delta = MIN(delta, pnorm)
   END IF
!
!  Evaluate the function at x + p and calculate its norm.
!
   iflag = 1
   CALL fcn(m, n, wa2, wa4, fjac, ldfjac, iflag)

   nfev = nfev + 1

   IF (iflag < 0) THEN
      go to 300
   END IF

   fnorm1 = enorm(m, wa4)
!
!  Compute the scaled actual reduction.
!
   actred = -1.0D+00
   IF (0.1D+00*fnorm1 < fnorm) THEN
      actred = 1.0D+00 - (fnorm1/fnorm)**2
   END IF
!
!  Compute the scaled predicted reduction and
!  the scaled directional derivative.
!
   DO j = 1, n
      wa3(j) = 0.0D+00
      l = ipvt(j)
      temp = wa1(l)
      wa3(1:j) = wa3(1:j) + fjac(1:j, j)*temp
   END DO

   temp1 = enorm(n, wa3)/fnorm
   temp2 = (SQRT(par)*pnorm)/fnorm
   prered = temp1**2 + temp2**2/0.5D+00
   dirder = -(temp1**2 + temp2**2)
!
!  Compute the ratio of the actual to the predicted reduction.
!
   IF (prered /= 0.0D+00) THEN
      ratio = actred/prered
   ELSE
      ratio = 0.0D+00
   END IF
!
!  Update the step bound.
!
   IF (ratio <= 0.25D+00) THEN

      IF (0.0D+00 <= actred) THEN
         temp = 0.5D+00
      END IF

      IF (actred < 0.0D+00) THEN
         temp = 0.5D+00*dirder/(dirder + 0.5D+00*actred)
      END IF

      IF (0.1D+00*fnorm1 >= fnorm .OR. temp < 0.1D+00) THEN
         temp = 0.1D+00
      END IF

      delta = temp*MIN(delta, pnorm/0.1D+00)
      par = par/temp

   ELSE

      IF (par == 0.0D+00 .OR. ratio >= 0.75D+00) THEN
         delta = 2.0D+00*pnorm
         par = 0.5D+00*par
      END IF

   END IF
!
!  Successful iteration.
!
!  Update X, FVEC, and their norms.
!
   IF (0.0001D+00 <= ratio) THEN
      x(1:n) = wa2(1:n)
      wa2(1:n) = diag(1:n)*x(1:n)
      fvec(1:m) = wa4(1:m)
      xnorm = enorm(n, wa2)
      fnorm = fnorm1
      iter = iter + 1
   END IF
!
!  Tests for convergence.
!
   IF (ABS(actred) <= ftol .AND. &
       prered <= ftol .AND. &
       0.5D+00*ratio <= 1.0D+00) THEN
      info = 1
   END IF

   IF (delta <= xtol*xnorm) THEN
      info = 2
   END IF

   IF (ABS(actred) <= ftol .AND. prered <= ftol &
       .AND. 0.5D+00*ratio <= 1.0D+00 .AND. info == 2) THEN
      info = 3
   END IF

   IF (info /= 0) THEN
      go to 300
   END IF
!
!  Tests for termination and stringent tolerances.
!
   IF (nfev >= maxfev) THEN
      info = 5
   END IF

   IF (ABS(actred) <= epsmch .AND. prered <= epsmch &
       .AND. 0.5D+00*ratio <= 1.0D+00) THEN
      info = 6
   END IF

   IF (delta <= epsmch*xnorm) THEN
      info = 7
   END IF

   IF (gnorm <= epsmch) THEN
      info = 8
   END IF

   IF (info /= 0) THEN
      go to 300
   END IF
!
!  End of the inner loop. repeat if iteration unsuccessful.
!
   IF (ratio < 0.0001D+00) THEN
      go to 200
   END IF
!
!  End of the outer loop.
!
   go to 30

300 CONTINUE
!
!  Termination, either normal or user imposed.
!
   IF (iflag < 0) THEN
      info = iflag
   END IF

   iflag = 0

   IF (0 < nprint) THEN
      CALL fcn(m, n, x, fvec, fjac, ldfjac, iflag)
   END IF

   RETURN
END SUBROUTINE lmder
SUBROUTINE lmder1(fcn, m, n, x, fvec, fjac, ldfjac, tol, info)

!*****************************************************************************80
!
!! LMDER1 minimizes M functions in N variables by Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMDER1 minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    This is done by using the more general least-squares solver LMDER.
!    The user must provide a subroutine which calculates the functions
!    and the jacobian.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the jacobian.  FCN should have the form:
!      subroutine fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
!      integer ( kind = 4 ) ldfjac
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fjac(ldfjac,n)
!      real ( kind = 8 ) fvec(m)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!    return this matrix in FJAC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer ( kind = 4 ) M, the number of functions.
!
!    Input, integer ( kind = 4 ) N, is the number of variables.
!    N must not exceed M.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(M), the functions evaluated at the output X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an M by N array.  The upper
!    N by N submatrix contains an upper triangular matrix R with
!    diagonal elements of nonincreasing magnitude such that
!      P' * ( JAC' * JAC ) * P = R' * R,
!    where P is a permutation matrix and JAC is the final calculated
!    jacobian.  Column J of P is column IPVT(J) of the identity matrix.
!    The lower trapezoidal part of FJAC contains information generated during
!    the computation of R.
!
!    Input, integer ( kind = 4 ) LDFJAC, is the leading dimension of FJAC,
!    which must be no less than M.
!
!    Input, real ( kind = 8 ) TOL.  Termination occurs when the algorithm
!    estimates either that the relative error in the sum of squares is at
!    most TOL or that the relative error between X and the solution is at
!    most TOL.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, algorithm estimates that the relative error in the sum of squares
!       is at most TOL.
!    2, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, FVEC is orthogonal to the columns of the jacobian to machine precision.
!    5, number of calls to FCN with IFLAG = 1 has reached 100*(N+1).
!    6, TOL is too small.  No further reduction in the sum of squares is
!       possible.
!    7, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!
   IMPLICIT NONE

   INTEGER(kind=4) ldfjac
   INTEGER(kind=4) m
   INTEGER(kind=4) n

   REAL(kind=8) diag(n)
   REAL(kind=8) factor
   EXTERNAL fcn
   REAL(kind=8) fjac(ldfjac, n)
   REAL(kind=8) ftol
   REAL(kind=8) fvec(m)
   REAL(kind=8) gtol
   INTEGER(kind=4) info
   INTEGER(kind=4) ipvt(n)
   INTEGER(kind=4) maxfev
   INTEGER(kind=4) mode
   INTEGER(kind=4) nfev
   INTEGER(kind=4) njev
   INTEGER(kind=4) nprint
   REAL(kind=8) qtf(n)
   REAL(kind=8) tol
   REAL(kind=8) x(n)
   REAL(kind=8) xtol

   info = 0

   IF (n <= 0) THEN
      RETURN
   ELSE IF (m < n) THEN
      RETURN
   ELSE IF (ldfjac < m) THEN
      RETURN
   ELSE IF (tol < 0.0D+00) THEN
      RETURN
   END IF

   factor = 100.0D+00
   maxfev = 100*(n + 1)
   ftol = tol
   xtol = tol
   gtol = 0.0D+00
   mode = 1
   nprint = 0

   CALL lmder(fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
              diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf)

   IF (info == 8) THEN
      info = 4
   END IF

   RETURN
END SUBROUTINE lmder1
SUBROUTINE lmdif(fcn, m, n, x, xdat, ydat, fvec, ftol, xtol, gtol, maxfev, epsfcn, &
                 diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf)

!*****************************************************************************80
!
!! LMDIF minimizes M functions in N variables by the Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMDIF minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    The user must provide a subroutine which calculates the functions.
!    The jacobian is then calculated by a forward-difference approximation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 Jul 2017
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!    adapted for AnOHM in SUEWS, Ting Sun, ting.sun@reading.ac.uk
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!      subroutine fcn ( m, n, x, xdat, ydat, fvec, iflag ) ! xdat, ydat added for AnOHM, TS 20 Jul 2017
!      integer ( kind = 4 ) m
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fvec(m),xdat(m),ydat(m)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer ( kind = 4 ) M, the number of functions.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!    N must not exceed M.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.

!    Input, real ( kind = 8 ) XDAT(M), YDAT(M).  On input, XDAT/YDAT must contain
!    observations to construct FVEC(M). ! added for AnOHM, TS 20 Jul 2017
!
!    Output, real ( kind = 8 ) FVEC(M), the functions evaluated at the output X.
!
!    Input, real ( kind = 8 ) FTOL.  Termination occurs when both the actual
!    and predicted relative reductions in the sum of squares are at most FTOL.
!    Therefore, FTOL measures the relative error desired in the sum of
!    squares.  FTOL should be nonnegative.
!
!    Input, real ( kind = 8 ) XTOL.  Termination occurs when the relative error
!    between two consecutive iterates is at most XTOL.  Therefore, XTOL
!    measures the relative error desired in the approximate solution.  XTOL
!    should be nonnegative.
!
!    Input, real ( kind = 8 ) GTOL. termination occurs when the cosine of the
!    angle between FVEC and any column of the jacobian is at most GTOL in
!    absolute value.  Therefore, GTOL measures the orthogonality desired
!    between the function vector and the columns of the jacobian.  GTOL should
!    be nonnegative.
!
!    Input, integer ( kind = 4 ) MAXFEV.  Termination occurs when the number of
!    calls to FCN is at least MAXFEV by the end of an iteration.
!
!    Input, real ( kind = 8 ) EPSFCN, is used in determining a suitable step
!    length for the forward-difference approximation.  This approximation
!    assumes that the relative errors in the functions are of the order of
!    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
!    the relative errors in the functions are of the order of the machine
!    precision.
!
!    Input/output, real ( kind = 8 ) DIAG(N).  If MODE = 1, then DIAG is set
!    internally.  If MODE = 2, then DIAG must contain positive entries that
!    serve as multiplicative scale factors for the variables.
!
!    Input, integer ( kind = 4 ) MODE, scaling option.
!    1, variables will be scaled internally.
!    2, scaling is specified by the input DIAG vector.
!
!    Input, real ( kind = 8 ) FACTOR, determines the initial step bound.
!    This bound is set to the product of FACTOR and the euclidean norm of
!    DIAG*X if nonzero, or else to FACTOR itself.  In most cases, FACTOR should
!    lie in the interval (0.1, 100) with 100 the recommended value.
!
!    Input, integer ( kind = 4 ) NPRINT, enables controlled printing of iterates
!    if it is positive.  In this case, FCN is called with IFLAG = 0 at the
!    beginning of the first iteration and every NPRINT iterations thereafter
!    and immediately prior to return, with X and FVEC available
!    for printing.  If NPRINT is not positive, no special calls
!    of FCN with IFLAG = 0 are made.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, both actual and predicted relative reductions in the sum of squares
!       are at most FTOL.
!    2, relative error between two consecutive iterates is at most XTOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, the cosine of the angle between FVEC and any column of the jacobian
!       is at most GTOL in absolute value.
!    5, number of calls to FCN has reached or exceeded MAXFEV.
!    6, FTOL is too small.  No further reduction in the sum of squares
!       is possible.
!    7, XTOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    8, GTOL is too small.  FVEC is orthogonal to the columns of the
!       jacobian to machine precision.
!
!    Output, integer ( kind = 4 ) NFEV, the number of calls to FCN.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an M by N array.  The upper
!    N by N submatrix of FJAC contains an upper triangular matrix R with
!    diagonal elements of nonincreasing magnitude such that
!
!      P' * ( JAC' * JAC ) * P = R' * R,
!
!    where P is a permutation matrix and JAC is the final calculated jacobian.
!    Column J of P is column IPVT(J) of the identity matrix.  The lower
!    trapezoidal part of FJAC contains information generated during
!    the computation of R.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of FJAC.
!    LDFJAC must be at least M.
!
!    Output, integer ( kind = 4 ) IPVT(N), defines a permutation matrix P such
!    that JAC * P = Q * R, where JAC is the final calculated jacobian, Q is
!    orthogonal (not stored), and R is upper triangular with diagonal
!    elements of nonincreasing magnitude.  Column J of P is column IPVT(J)
!    of the identity matrix.
!
!    Output, real ( kind = 8 ) QTF(N), the first N elements of Q'*FVEC.
!
   IMPLICIT NONE

   INTEGER(kind=4) ldfjac
   INTEGER(kind=4) m
   INTEGER(kind=4) n

   REAL(kind=8) actred
   REAL(kind=8) delta
   REAL(kind=8) diag(n)
   REAL(kind=8) dirder
   REAL(kind=8) enorm
   REAL(kind=8) epsfcn
   REAL(kind=8) epsmch
   REAL(kind=8) factor
   EXTERNAL fcn
   REAL(kind=8) fjac(ldfjac, n)
   REAL(kind=8) fnorm
   REAL(kind=8) fnorm1
   REAL(kind=8) ftol
   REAL(kind=8) :: fvec(m), xdat(m), ydat(m)
   REAL(kind=8) gnorm
   REAL(kind=8) gtol
   INTEGER(kind=4) i
   INTEGER(kind=4) iflag
   INTEGER(kind=4) iter
   INTEGER(kind=4) info
   INTEGER(kind=4) ipvt(n)
   INTEGER(kind=4) j
   INTEGER(kind=4) l
   INTEGER(kind=4) maxfev
   INTEGER(kind=4) mode
   INTEGER(kind=4) nfev
   INTEGER(kind=4) nprint
   REAL(kind=8) par
   LOGICAL pivot
   REAL(kind=8) pnorm
   REAL(kind=8) prered
   REAL(kind=8) qtf(n)
   REAL(kind=8) ratio
   REAL(kind=8) sum2
   REAL(kind=8) temp
   REAL(kind=8) temp1
   REAL(kind=8) temp2
   REAL(kind=8) wa1(n)
   REAL(kind=8) wa2(n)
   REAL(kind=8) wa3(n)
   REAL(kind=8) wa4(m)
   REAL(kind=8) x(n)
   REAL(kind=8) xnorm
   REAL(kind=8) xtol

   epsmch = EPSILON(epsmch)

   info = 0
   iflag = 0
   nfev = 0

   IF (n <= 0) THEN
      go to 300
   ELSE IF (m < n) THEN
      go to 300
   ELSE IF (ldfjac < m) THEN
      go to 300
   ELSE IF (ftol < 0.0D+00) THEN
      go to 300
   ELSE IF (xtol < 0.0D+00) THEN
      go to 300
   ELSE IF (gtol < 0.0D+00) THEN
      go to 300
   ELSE IF (maxfev <= 0) THEN
      go to 300
   ELSE IF (factor <= 0.0D+00) THEN
      go to 300
   END IF

   IF (mode == 2) THEN
      DO j = 1, n
         IF (diag(j) <= 0.0D+00) THEN
            go to 300
         END IF
      END DO
   END IF
!
!  Evaluate the function at the starting point and calculate its norm.
!
   iflag = 1
   CALL fcn(m, n, x, xdat, ydat, fvec, iflag)
   nfev = 1

   IF (iflag < 0) THEN
      go to 300
   END IF

   fnorm = enorm(m, fvec)
!
!  Initialize Levenberg-Marquardt parameter and iteration counter.
!
   par = 0.0D+00
   iter = 1
!
!  Beginning of the outer loop.
!
30 CONTINUE
!
!  Calculate the jacobian matrix.
!
   iflag = 2
   CALL fdjac2(fcn, m, n, x, xdat, ydat, fvec, fjac, ldfjac, iflag, epsfcn)
   nfev = nfev + n

   IF (iflag < 0) THEN
      go to 300
   END IF
!
!  If requested, call FCN to enable printing of iterates.
!
   IF (0 < nprint) THEN
      iflag = 0
      IF (MOD(iter - 1, nprint) == 0) THEN
         CALL fcn(m, n, x, xdat, ydat, fvec, iflag)
      END IF
      IF (iflag < 0) THEN
         go to 300
      END IF
   END IF
!
!  Compute the QR factorization of the jacobian.
!
   pivot = .TRUE.
   CALL qrfac(m, n, fjac, ldfjac, pivot, ipvt, n, wa1, wa2)
!
!  On the first iteration and if MODE is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
   IF (iter == 1) THEN

      IF (mode /= 2) THEN
         diag(1:n) = wa2(1:n)
         DO j = 1, n
            IF (wa2(j) == 0.0D+00) THEN
               diag(j) = 1.0D+00
            END IF
         END DO
      END IF
!
!  On the first iteration, calculate the norm of the scaled X
!  and initialize the step bound DELTA.
!
      wa3(1:n) = diag(1:n)*x(1:n)
      xnorm = enorm(n, wa3)
      delta = factor*xnorm
      IF (delta == 0.0D+00) THEN
         delta = factor
      END IF
   END IF
!
!  Form Q' * FVEC and store the first N components in QTF.
!
   wa4(1:m) = fvec(1:m)

   DO j = 1, n

      IF (fjac(j, j) /= 0.0D+00) THEN
         sum2 = DOT_PRODUCT(wa4(j:m), fjac(j:m, j))
         temp = -sum2/fjac(j, j)
         wa4(j:m) = wa4(j:m) + fjac(j:m, j)*temp
      END IF

      fjac(j, j) = wa1(j)
      qtf(j) = wa4(j)

   END DO
!
!  Compute the norm of the scaled gradient.
!
   gnorm = 0.0D+00

   IF (fnorm /= 0.0D+00) THEN

      DO j = 1, n

         l = ipvt(j)

         IF (wa2(l) /= 0.0D+00) THEN
            sum2 = 0.0D+00
            DO i = 1, j
               sum2 = sum2 + fjac(i, j)*(qtf(i)/fnorm)
            END DO
            gnorm = MAX(gnorm, ABS(sum2/wa2(l)))
         END IF

      END DO

   END IF
!
!  Test for convergence of the gradient norm.
!
   IF (gnorm <= gtol) THEN
      info = 4
      go to 300
   END IF
!
!  Rescale if necessary.
!
   IF (mode /= 2) THEN
      DO j = 1, n
         diag(j) = MAX(diag(j), wa2(j))
      END DO
   END IF
!
!  Beginning of the inner loop.
!
200 CONTINUE
!
!  Determine the Levenberg-Marquardt parameter.
!
   CALL lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2)
!
!  Store the direction P and X + P.
!  Calculate the norm of P.
!
   wa1(1:n) = -wa1(1:n)
   wa2(1:n) = x(1:n) + wa1(1:n)
   wa3(1:n) = diag(1:n)*wa1(1:n)

   pnorm = enorm(n, wa3)
!
!  On the first iteration, adjust the initial step bound.
!
   IF (iter == 1) THEN
      delta = MIN(delta, pnorm)
   END IF
!
!  Evaluate the function at X + P and calculate its norm.
!
   iflag = 1
   CALL fcn(m, n, wa2, xdat, ydat, wa4, iflag)
   nfev = nfev + 1
   IF (iflag < 0) THEN
      go to 300
   END IF
   fnorm1 = enorm(m, wa4)
!
!  Compute the scaled actual reduction.
!
   IF (0.1D+00*fnorm1 < fnorm) THEN
      actred = 1.0D+00 - (fnorm1/fnorm)**2
   ELSE
      actred = -1.0D+00
   END IF
!
!  Compute the scaled predicted reduction and the scaled directional derivative.
!
   DO j = 1, n
      wa3(j) = 0.0D+00
      l = ipvt(j)
      temp = wa1(l)
      wa3(1:j) = wa3(1:j) + fjac(1:j, j)*temp
   END DO

   temp1 = enorm(n, wa3)/fnorm
   temp2 = (SQRT(par)*pnorm)/fnorm
   prered = temp1**2 + temp2**2/0.5D+00
   dirder = -(temp1**2 + temp2**2)
!
!  Compute the ratio of the actual to the predicted reduction.
!
   ratio = 0.0D+00
   IF (prered /= 0.0D+00) THEN
      ratio = actred/prered
   END IF
!
!  Update the step bound.
!
   IF (ratio <= 0.25D+00) THEN

      IF (actred >= 0.0D+00) THEN
         temp = 0.5D+00
      END IF

      IF (actred < 0.0D+00) THEN
         temp = 0.5D+00*dirder/(dirder + 0.5D+00*actred)
      END IF

      IF (0.1D+00*fnorm1 >= fnorm .OR. temp < 0.1D+00) THEN
         temp = 0.1D+00
      END IF

      delta = temp*MIN(delta, pnorm/0.1D+00)
      par = par/temp

   ELSE

      IF (par == 0.0D+00 .OR. ratio >= 0.75D+00) THEN
         delta = 2.0D+00*pnorm
         par = 0.5D+00*par
      END IF

   END IF
!
!  Test for successful iteration.
!

!
!  Successful iteration. update X, FVEC, and their norms.
!
   IF (0.0001D+00 <= ratio) THEN
      x(1:n) = wa2(1:n)
      wa2(1:n) = diag(1:n)*x(1:n)
      fvec(1:m) = wa4(1:m)
      xnorm = enorm(n, wa2)
      fnorm = fnorm1
      iter = iter + 1
   END IF
!
!  Tests for convergence.
!
   IF (ABS(actred) <= ftol .AND. prered <= ftol &
       .AND. 0.5D+00*ratio <= 1.0D+00) THEN
      info = 1
   END IF

   IF (delta <= xtol*xnorm) THEN
      info = 2
   END IF

   IF (ABS(actred) <= ftol .AND. prered <= ftol &
       .AND. 0.5D+00*ratio <= 1.0D+00 .AND. info == 2) info = 3

   IF (info /= 0) THEN
      go to 300
   END IF
!
!  Tests for termination and stringent tolerances.
!
   IF (maxfev <= nfev) THEN
      info = 5
   END IF

   IF (ABS(actred) <= epsmch .AND. prered <= epsmch &
       .AND. 0.5D+00*ratio <= 1.0D+00) THEN
      info = 6
   END IF

   IF (delta <= epsmch*xnorm) THEN
      info = 7
   END IF

   IF (gnorm <= epsmch) THEN
      info = 8
   END IF

   IF (info /= 0) THEN
      go to 300
   END IF
!
!  End of the inner loop.  Repeat if iteration unsuccessful.
!
   IF (ratio < 0.0001D+00) THEN
      go to 200
   END IF
!
!  End of the outer loop.
!
   go to 30

300 CONTINUE
!
!  Termination, either normal or user imposed.
!
   IF (iflag < 0) THEN
      info = iflag
   END IF

   iflag = 0

   IF (0 < nprint) THEN
      CALL fcn(m, n, x, xdat, ydat, fvec, iflag)
   END IF

   RETURN
END SUBROUTINE lmdif
SUBROUTINE lmdif1(fcn, m, n, x, xdat, ydat, fvec, tol, info)

!*****************************************************************************80
!
!! LMDIF1 minimizes M functions in N variables using Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMDIF1 minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    This is done by using the more general least-squares solver LMDIF.
!    The user must provide a subroutine which calculates the functions.
!    The jacobian is then calculated by a forward-difference approximation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 Jul 2017
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!    adapted for AnOHM in SUEWS, Ting Sun, ting.sun@reading.ac.uk
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!      subroutine fcn ( m, n, x, xdat, ydat, fvec, iflag ) ! xdat, ydat added for AnOHM, TS 20 Jul 2017
!      integer ( kind = 4 ) m
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fvec(m),xdat(m),ydat(m)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer ( kind = 4 ) M, the number of functions.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!    N must not exceed M.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.

!    Input, real ( kind = 8 ) XDAT(M), YDAT(M).  On input, XDAT/YDAT must contain
!    observations to construct FVEC(M). ! added for AnOHM, TS 20 Jul 2017
!
!    Output, real ( kind = 8 ) FVEC(M), the functions evaluated at the output X.
!
!    Input, real ( kind = 8 ) TOL.  Termination occurs when the algorithm
!    estimates either that the relative error in the sum of squares is at
!    most TOL or that the relative error between X and the solution is at
!    most TOL.  TOL should be nonnegative.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, algorithm estimates that the relative error in the sum of squares
!       is at most TOL.
!    2, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, FVEC is orthogonal to the columns of the jacobian to machine precision.
!    5, number of calls to FCN has reached or exceeded 200*(N+1).
!    6, TOL is too small.  No further reduction in the sum of squares
!       is possible.
!    7, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!
   IMPLICIT NONE

   INTEGER(kind=4) m
   INTEGER(kind=4) n

   REAL(kind=8) diag(n)
   REAL(kind=8) epsfcn
   REAL(kind=8) factor
   EXTERNAL fcn
   REAL(kind=8) fjac(m, n)
   REAL(kind=8) ftol
   REAL(kind=8) ::fvec(m), xdat(m), ydat(m)
   REAL(kind=8) gtol
   INTEGER(kind=4) info
   INTEGER(kind=4) ipvt(n)
   INTEGER(kind=4) ldfjac
   INTEGER(kind=4) maxfev
   INTEGER(kind=4) mode
   INTEGER(kind=4) nfev
   INTEGER(kind=4) nprint
   REAL(kind=8) qtf(n)
   REAL(kind=8) tol
   REAL(kind=8) x(n)
   REAL(kind=8) xtol

   info = 0

   IF (n <= 0) THEN
      RETURN
   ELSE IF (m < n) THEN
      RETURN
   ELSE IF (tol < 0.0D+00) THEN
      RETURN
   END IF

   factor = 100.0D+00
   maxfev = 200*(n + 1)
   ftol = tol
   xtol = tol
   gtol = 0.0D+00
   epsfcn = 0.0D+00
   mode = 1
   nprint = 0
   ldfjac = m
! print*, 'x in limdif1',x
   CALL lmdif(fcn, m, n, x, xdat, ydat, fvec, ftol, xtol, gtol, maxfev, epsfcn, &
              diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf)

   IF (info == 8) THEN
      info = 4
   END IF

   RETURN
END SUBROUTINE lmdif1
SUBROUTINE lmpar(n, r, ldr, ipvt, diag, qtb, delta, par, x, sdiag)

!*****************************************************************************80
!
!! LMPAR computes a parameter for the Levenberg-Marquardt method.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N nonsingular diagonal
!    matrix D, an M-vector B, and a positive number DELTA,
!    the problem is to determine a value for the parameter
!    PAR such that if X solves the system
!
!      A*X = B,
!      sqrt ( PAR ) * D * X = 0,
!
!    in the least squares sense, and DXNORM is the euclidean
!    norm of D*X, then either PAR is zero and
!
!      ( DXNORM - DELTA ) <= 0.1 * DELTA,
!
!    or PAR is positive and
!
!      abs ( DXNORM - DELTA) <= 0.1 * DELTA.
!
!    This subroutine completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization, with column pivoting, of A.  That is, if
!    A*P = Q*R, where P is a permutation matrix, Q has orthogonal
!    columns, and R is an upper triangular matrix with diagonal
!    elements of nonincreasing magnitude, then LMPAR expects
!    the full upper triangle of R, the permutation matrix P,
!    and the first N components of Q'*B.  On output
!    LMPAR also provides an upper triangular matrix S such that
!
!      P' * ( A' * A + PAR * D * D ) * P = S'* S.
!
!    S is employed within LMPAR and may be of separate interest.
!
!    Only a few iterations are generally needed for convergence
!    of the algorithm.  If, however, the limit of 10 iterations
!    is reached, then the output PAR will contain the best
!    value obtained so far.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2014
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of R.
!
!    Input/output, real ( kind = 8 ) R(LDR,N),the N by N matrix.  The full
!    upper triangle must contain the full upper triangle of the matrix R.
!    On output the full upper triangle is unaltered, and the strict lower
!    triangle contains the strict upper triangle (transposed) of the upper
!    triangular matrix S.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of R.  LDR must be
!    no less than N.
!
!    Input, integer ( kind = 4 ) IPVT(N), defines the permutation matrix P
!    such that A*P = Q*R.  Column J of P is column IPVT(J) of the
!    identity matrix.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real ( kind = 8 ) QTB(N), the first N elements of the vector Q'*B.
!
!    Input, real ( kind = 8 ) DELTA, an upper bound on the euclidean norm
!    of D*X.  DELTA should be positive.
!
!    Input/output, real ( kind = 8 ) PAR.  On input an initial estimate of the
!    Levenberg-Marquardt parameter.  On output the final estimate.
!    PAR should be nonnegative.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution of the system
!    A*X = B, sqrt(PAR)*D*X = 0, for the output value of PAR.
!
!    Output, real ( kind = 8 ) SDIAG(N), the diagonal elements of the upper
!    triangular matrix S.
!
   IMPLICIT NONE

   INTEGER(kind=4) ldr
   INTEGER(kind=4) n

   REAL(kind=8) delta
   REAL(kind=8) diag(n)
   REAL(kind=8) dwarf
   REAL(kind=8) dxnorm
   REAL(kind=8) enorm
   REAL(kind=8) gnorm
   REAL(kind=8) fp
   ! integer ( kind = 4 ) i
   INTEGER(kind=4) ipvt(n)
   INTEGER(kind=4) iter
   INTEGER(kind=4) j
   INTEGER(kind=4) k
   INTEGER(kind=4) l
   INTEGER(kind=4) nsing
   REAL(kind=8) par
   REAL(kind=8) parc
   REAL(kind=8) parl
   REAL(kind=8) paru
   ! real ( kind = 8 ) qnorm
   REAL(kind=8) qtb(n)
   REAL(kind=8) r(ldr, n)
   REAL(kind=8) sdiag(n)
   REAL(kind=8) sum2
   REAL(kind=8) temp
   REAL(kind=8) wa1(n)
   REAL(kind=8) wa2(n)
   REAL(kind=8) x(n)
!
!  DWARF is the smallest positive magnitude.
!
   dwarf = TINY(dwarf)
!
!  Compute and store in X the Gauss-Newton direction.
!
!  If the jacobian is rank-deficient, obtain a least squares solution.
!
   nsing = n

   DO j = 1, n
      wa1(j) = qtb(j)
      IF (r(j, j) == 0.0D+00 .AND. nsing == n) THEN
         nsing = j - 1
      END IF
      IF (nsing < n) THEN
         wa1(j) = 0.0D+00
      END IF
   END DO

   DO k = 1, nsing
      j = nsing - k + 1
      wa1(j) = wa1(j)/r(j, j)
      temp = wa1(j)
      wa1(1:j - 1) = wa1(1:j - 1) - r(1:j - 1, j)*temp
   END DO

   DO j = 1, n
      l = ipvt(j)
      x(l) = wa1(j)
   END DO
!
!  Initialize the iteration counter.
!  Evaluate the function at the origin, and test
!  for acceptance of the Gauss-Newton direction.
!
   iter = 0
   wa2(1:n) = diag(1:n)*x(1:n)
   dxnorm = enorm(n, wa2)
   fp = dxnorm - delta

   IF (fp <= 0.1D+00*delta) THEN
      IF (iter == 0) THEN
         par = 0.0D+00
      END IF
      RETURN
   END IF
!
!  If the jacobian is not rank deficient, the Newton
!  step provides a lower bound, PARL, for the zero of
!  the function.
!
!  Otherwise set this bound to zero.
!
   parl = 0.0D+00

   IF (n <= nsing) THEN

      DO j = 1, n
         l = ipvt(j)
         wa1(j) = diag(l)*(wa2(l)/dxnorm)
      END DO

      DO j = 1, n
         sum2 = DOT_PRODUCT(wa1(1:j - 1), r(1:j - 1, j))
         wa1(j) = (wa1(j) - sum2)/r(j, j)
      END DO

      temp = enorm(n, wa1)
      parl = ((fp/delta)/temp)/temp

   END IF
!
!  Calculate an upper bound, PARU, for the zero of the function.
!
   DO j = 1, n
      sum2 = DOT_PRODUCT(qtb(1:j), r(1:j, j))
      l = ipvt(j)
      wa1(j) = sum2/diag(l)
   END DO

   gnorm = enorm(n, wa1)
   paru = gnorm/delta

   IF (paru == 0.0D+00) THEN
      paru = dwarf/MIN(delta, 0.1D+00)
   END IF
!
!  If the input PAR lies outside of the interval (PARL, PARU),
!  set PAR to the closer endpoint.
!
   par = MAX(par, parl)
   par = MIN(par, paru)
   IF (par == 0.0D+00) THEN
      par = gnorm/dxnorm
   END IF
!
!  Beginning of an iteration.
!
   DO

      iter = iter + 1
!
!  Evaluate the function at the current value of PAR.
!
      IF (par == 0.0D+00) THEN
         par = MAX(dwarf, 0.001D+00*paru)
      END IF

      wa1(1:n) = SQRT(par)*diag(1:n)

      CALL qrsolv(n, r, ldr, ipvt, wa1, qtb, x, sdiag)

      wa2(1:n) = diag(1:n)*x(1:n)
      dxnorm = enorm(n, wa2)
      temp = fp
      fp = dxnorm - delta
!
!  If the function is small enough, accept the current value of PAR.
!
      IF (ABS(fp) <= 0.1D+00*delta) THEN
         EXIT
      END IF
!
!  Test for the exceptional cases where PARL
!  is zero or the number of iterations has reached 10.
!
      IF (parl == 0.0D+00 .AND. fp <= temp .AND. temp < 0.0D+00) THEN
         EXIT
      ELSE IF (iter == 10) THEN
         EXIT
      END IF
!
!  Compute the Newton correction.
!
      DO j = 1, n
         l = ipvt(j)
         wa1(j) = diag(l)*(wa2(l)/dxnorm)
      END DO

      DO j = 1, n
         wa1(j) = wa1(j)/sdiag(j)
         temp = wa1(j)
         wa1(j + 1:n) = wa1(j + 1:n) - r(j + 1:n, j)*temp
      END DO

      temp = enorm(n, wa1)
      parc = ((fp/delta)/temp)/temp
!
!  Depending on the sign of the function, update PARL or PARU.
!
      IF (0.0D+00 < fp) THEN
         parl = MAX(parl, par)
      ELSE IF (fp < 0.0D+00) THEN
         paru = MIN(paru, par)
      END IF
!
!  Compute an improved estimate for PAR.
!
      par = MAX(parl, par + parc)
!
!  End of an iteration.
!
   END DO
!
!  Termination.
!
   IF (iter == 0) THEN
      par = 0.0D+00
   END IF

   RETURN
END SUBROUTINE lmpar
SUBROUTINE lmstr(fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
                 diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf)

!*****************************************************************************80
!
!! LMSTR minimizes M functions in N variables using Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMSTR minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm
!    which uses minimal storage.
!
!    The user must provide a subroutine which calculates the functions and
!    the rows of the jacobian.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the rows of the jacobian.
!    FCN should have the form:
!      subroutine fcn ( m, n, x, fvec, fjrow, iflag )
!      integer ( kind = 4 ) m
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fjrow(n)
!      real ( kind = 8 ) fvec(m)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If the input value of IFLAG is I > 1, calculate the (I-1)-st row of
!    the jacobian at X, and return this vector in FJROW.
!    To terminate the algorithm, set the output value of IFLAG negative.
!
!    Input, integer ( kind = 4 ) M, the number of functions.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!    N must not exceed M.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(M), the functions evaluated at the output X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an N by N array.  The upper
!    triangle of FJAC contains an upper triangular matrix R such that
!
!      P' * ( JAC' * JAC ) * P = R' * R,
!
!    where P is a permutation matrix and JAC is the final calculated jacobian.
!    Column J of P is column IPVT(J) of the identity matrix.  The lower
!    triangular part of FJAC contains information generated during
!    the computation of R.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of FJAC.
!    LDFJAC must be at least N.
!
!    Input, real ( kind = 8 ) FTOL.  Termination occurs when both the actual and
!    predicted relative reductions in the sum of squares are at most FTOL.
!    Therefore, FTOL measures the relative error desired in the sum of
!    squares.  FTOL should be nonnegative.
!
!    Input, real ( kind = 8 ) XTOL.  Termination occurs when the relative error
!    between two consecutive iterates is at most XTOL.  XTOL should be
!    nonnegative.
!
!    Input, real ( kind = 8 ) GTOL. termination occurs when the cosine of the
!    angle between FVEC and any column of the jacobian is at most GTOL in
!    absolute value.  Therefore, GTOL measures the orthogonality desired
!    between the function vector and the columns of the jacobian.  GTOL should
!    be nonnegative.
!
!    Input, integer ( kind = 4 ) MAXFEV.  Termination occurs when the number
!    of calls to FCN with IFLAG = 1 is at least MAXFEV by the end of
!    an iteration.
!
!    Input/output, real ( kind = 8 ) DIAG(N).  If MODE = 1, then DIAG is set
!    internally.  If MODE = 2, then DIAG must contain positive entries that
!    serve as multiplicative scale factors for the variables.
!
!    Input, integer ( kind = 4 ) MODE, scaling option.
!    1, variables will be scaled internally.
!    2, scaling is specified by the input DIAG vector.
!
!    Input, real ( kind = 8 ) FACTOR, determines the initial step bound.  This
!    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
!    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
!    in the interval (0.1, 100) with 100 the recommended value.
!
!    Input, integer ( kind = 4 ) NPRINT, enables controlled printing of iterates
!    if it is positive.  In this case, FCN is called with IFLAG = 0 at the
!    beginning of the first iteration and every NPRINT iterations thereafter
!    and immediately prior to return, with X and FVEC available
!    for printing.  If NPRINT is not positive, no special calls
!    of FCN with IFLAG = 0 are made.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See the
!    description of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, both actual and predicted relative reductions in the sum of squares
!       are at most FTOL.
!    2, relative error between two consecutive iterates is at most XTOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, the cosine of the angle between FVEC and any column of the jacobian
!       is at most GTOL in absolute value.
!    5, number of calls to FCN with IFLAG = 1 has reached MAXFEV.
!    6, FTOL is too small.  No further reduction in the sum of squares is
!       possible.
!    7, XTOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    8, GTOL is too small.  FVEC is orthogonal to the columns of the
!       jacobian to machine precision.
!
!    Output, integer ( kind = 4 ) NFEV, the number of calls to FCN
!    with IFLAG = 1.
!
!    Output, integer ( kind = 4 ) NJEV, the number of calls to FCN
!    with IFLAG = 2.
!
!    Output, integer ( kind = 4 ) IPVT(N), defines a permutation matrix P such
!    that JAC * P = Q * R, where JAC is the final calculated jacobian, Q is
!    orthogonal (not stored), and R is upper triangular.
!    Column J of P is column IPVT(J) of the identity matrix.
!
!    Output, real ( kind = 8 ) QTF(N), contains the first N elements of Q'*FVEC.
!
   IMPLICIT NONE

   INTEGER(kind=4) ldfjac
   INTEGER(kind=4) m
   INTEGER(kind=4) n

   REAL(kind=8) actred
   REAL(kind=8) delta
   REAL(kind=8) diag(n)
   REAL(kind=8) dirder
   REAL(kind=8) enorm
   REAL(kind=8) epsmch
   REAL(kind=8) factor
   EXTERNAL fcn
   REAL(kind=8) fjac(ldfjac, n)
   REAL(kind=8) fnorm
   REAL(kind=8) fnorm1
   REAL(kind=8) ftol
   REAL(kind=8) fvec(m)
   REAL(kind=8) gnorm
   REAL(kind=8) gtol
   INTEGER(kind=4) i
   INTEGER(kind=4) iflag
   INTEGER(kind=4) info
   INTEGER(kind=4) ipvt(n)
   INTEGER(kind=4) iter
   INTEGER(kind=4) j
   INTEGER(kind=4) l
   INTEGER(kind=4) maxfev
   INTEGER(kind=4) mode
   INTEGER(kind=4) nfev
   INTEGER(kind=4) njev
   INTEGER(kind=4) nprint
   REAL(kind=8) par
   LOGICAL pivot
   REAL(kind=8) pnorm
   REAL(kind=8) prered
   REAL(kind=8) qtf(n)
   REAL(kind=8) ratio
   LOGICAL sing
   REAL(kind=8) sum2
   REAL(kind=8) temp
   REAL(kind=8) temp1
   REAL(kind=8) temp2
   REAL(kind=8) wa1(n)
   REAL(kind=8) wa2(n)
   REAL(kind=8) wa3(n)
   REAL(kind=8) wa4(m)
   REAL(kind=8) x(n)
   REAL(kind=8) xnorm
   REAL(kind=8) xtol

   epsmch = EPSILON(epsmch)

   info = 0
   iflag = 0
   nfev = 0
   njev = 0
   xnorm = 0
!
!  Check the input parameters for errors.
!
   IF (n <= 0) THEN
      go to 340
   ELSE IF (m < n) THEN
      go to 340
   ELSE IF (ldfjac < n) THEN
      go to 340
   ELSE IF (ftol < 0.0D+00) THEN
      go to 340
   ELSE IF (xtol < 0.0D+00) THEN
      go to 340
   ELSE IF (gtol < 0.0D+00) THEN
      go to 340
   ELSE IF (maxfev <= 0) THEN
      go to 340
   ELSE IF (factor <= 0.0D+00) THEN
      go to 340
   END IF

   IF (mode == 2) THEN
      DO j = 1, n
         IF (diag(j) <= 0.0D+00) THEN
            go to 340
         END IF
      END DO
   END IF
!
!  Evaluate the function at the starting point and calculate its norm.
!
   iflag = 1
   CALL fcn(m, n, x, fvec, wa3, iflag)
   nfev = 1

   IF (iflag < 0) THEN
      go to 340
   END IF

   fnorm = enorm(m, fvec)
!
!  Initialize Levenberg-Marquardt parameter and iteration counter.
!
   par = 0.0D+00
   iter = 1
!
!  Beginning of the outer loop.
!
30 CONTINUE
!
!  If requested, call FCN to enable printing of iterates.
!
   IF (0 < nprint) THEN
      iflag = 0
      IF (MOD(iter - 1, nprint) == 0) THEN
         CALL fcn(m, n, x, fvec, wa3, iflag)
      END IF
      IF (iflag < 0) THEN
         go to 340
      END IF
   END IF
!
!  Compute the QR factorization of the jacobian matrix calculated one row
!  at a time, while simultaneously forming Q'* FVEC and storing
!  the first N components in QTF.
!
   qtf(1:n) = 0.0D+00
   fjac(1:n, 1:n) = 0.0D+00
   iflag = 2

   DO i = 1, m
      CALL fcn(m, n, x, fvec, wa3, iflag)
      IF (iflag < 0) THEN
         go to 340
      END IF
      temp = fvec(i)
      CALL rwupdt(n, fjac, ldfjac, wa3, qtf, temp, wa1, wa2)
      iflag = iflag + 1
   END DO

   njev = njev + 1
!
!  If the jacobian is rank deficient, call QRFAC to
!  reorder its columns and update the components of QTF.
!
   sing = .FALSE.
   DO j = 1, n
      IF (fjac(j, j) == 0.0D+00) THEN
         sing = .TRUE.
      END IF
      ipvt(j) = j
      wa2(j) = enorm(j, fjac(1, j))
   END DO

   IF (sing) THEN

      pivot = .TRUE.
      CALL qrfac(n, n, fjac, ldfjac, pivot, ipvt, n, wa1, wa2)

      DO j = 1, n

         IF (fjac(j, j) /= 0.0D+00) THEN

            sum2 = DOT_PRODUCT(qtf(j:n), fjac(j:n, j))
            temp = -sum2/fjac(j, j)
            qtf(j:n) = qtf(j:n) + fjac(j:n, j)*temp

         END IF

         fjac(j, j) = wa1(j)

      END DO

   END IF
!
!  On the first iteration
!    if mode is 1,
!      scale according to the norms of the columns of the initial jacobian.
!    calculate the norm of the scaled X,
!    initialize the step bound delta.
!
   IF (iter == 1) THEN

      IF (mode /= 2) THEN

         diag(1:n) = wa2(1:n)
         DO j = 1, n
            IF (wa2(j) == 0.0D+00) THEN
               diag(j) = 1.0D+00
            END IF
         END DO

      END IF

      wa3(1:n) = diag(1:n)*x(1:n)
      xnorm = enorm(n, wa3)
      delta = factor*xnorm
      IF (delta == 0.0D+00) THEN
         delta = factor
      END IF

   END IF
!
!  Compute the norm of the scaled gradient.
!
   gnorm = 0.0D+00

   IF (fnorm /= 0.0D+00) THEN

      DO j = 1, n
         l = ipvt(j)
         IF (wa2(l) /= 0.0D+00) THEN
            sum2 = DOT_PRODUCT(qtf(1:j), fjac(1:j, j))/fnorm
            gnorm = MAX(gnorm, ABS(sum2/wa2(l)))
         END IF
      END DO

   END IF
!
!  Test for convergence of the gradient norm.
!
   IF (gnorm <= gtol) THEN
      info = 4
      go to 340
   END IF
!
!  Rescale if necessary.
!
   IF (mode /= 2) THEN
      DO j = 1, n
         diag(j) = MAX(diag(j), wa2(j))
      END DO
   END IF
!
!  Beginning of the inner loop.
!
240 CONTINUE
!
!  Determine the Levenberg-Marquardt parameter.
!
   CALL lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2)
!
!  Store the direction P and X + P.
!  Calculate the norm of P.
!
   wa1(1:n) = -wa1(1:n)
   wa2(1:n) = x(1:n) + wa1(1:n)
   wa3(1:n) = diag(1:n)*wa1(1:n)
   pnorm = enorm(n, wa3)
!
!  On the first iteration, adjust the initial step bound.
!
   IF (iter == 1) THEN
      delta = MIN(delta, pnorm)
   END IF
!
!  Evaluate the function at X + P and calculate its norm.
!
   iflag = 1
   CALL fcn(m, n, wa2, wa4, wa3, iflag)
   nfev = nfev + 1
   IF (iflag < 0) THEN
      go to 340
   END IF
   fnorm1 = enorm(m, wa4)
!
!  Compute the scaled actual reduction.
!
   IF (0.1D+00*fnorm1 < fnorm) THEN
      actred = 1.0D+00 - (fnorm1/fnorm)**2
   ELSE
      actred = -1.0D+00
   END IF
!
!  Compute the scaled predicted reduction and
!  the scaled directional derivative.
!
   DO j = 1, n
      wa3(j) = 0.0D+00
      l = ipvt(j)
      temp = wa1(l)
      wa3(1:j) = wa3(1:j) + fjac(1:j, j)*temp
   END DO

   temp1 = enorm(n, wa3)/fnorm
   temp2 = (SQRT(par)*pnorm)/fnorm
   prered = temp1**2 + temp2**2/0.5D+00
   dirder = -(temp1**2 + temp2**2)
!
!  Compute the ratio of the actual to the predicted reduction.
!
   IF (prered /= 0.0D+00) THEN
      ratio = actred/prered
   ELSE
      ratio = 0.0D+00
   END IF
!
!  Update the step bound.
!
   IF (ratio <= 0.25D+00) THEN

      IF (actred >= 0.0D+00) THEN
         temp = 0.5D+00
      ELSE
         temp = 0.5D+00*dirder/(dirder + 0.5D+00*actred)
      END IF

      IF (0.1D+00*fnorm1 >= fnorm .OR. temp < 0.1D+00) THEN
         temp = 0.1D+00
      END IF

      delta = temp*MIN(delta, pnorm/0.1D+00)
      par = par/temp

   ELSE

      IF (par == 0.0D+00 .OR. ratio >= 0.75D+00) THEN
         delta = pnorm/0.5D+00
         par = 0.5D+00*par
      END IF

   END IF
!
!  Test for successful iteration.
!
   IF (ratio >= 0.0001D+00) THEN
      x(1:n) = wa2(1:n)
      wa2(1:n) = diag(1:n)*x(1:n)
      fvec(1:m) = wa4(1:m)
      xnorm = enorm(n, wa2)
      fnorm = fnorm1
      iter = iter + 1
   END IF
!
!  Tests for convergence, termination and stringent tolerances.
!
   IF (ABS(actred) <= ftol .AND. prered <= ftol &
       .AND. 0.5D+00*ratio <= 1.0D+00) THEN
      info = 1
   END IF

   IF (delta <= xtol*xnorm) THEN
      info = 2
   END IF

   IF (ABS(actred) <= ftol .AND. prered <= ftol &
       .AND. 0.5D+00*ratio <= 1.0D+00 .AND. info == 2) THEN
      info = 3
   END IF

   IF (info /= 0) THEN
      go to 340
   END IF

   IF (nfev >= maxfev) THEN
      info = 5
   END IF

   IF (ABS(actred) <= epsmch .AND. prered <= epsmch &
       .AND. 0.5D+00*ratio <= 1.0D+00) THEN
      info = 6
   END IF

   IF (delta <= epsmch*xnorm) THEN
      info = 7
   END IF

   IF (gnorm <= epsmch) THEN
      info = 8
   END IF

   IF (info /= 0) THEN
      go to 340
   END IF
!
!  End of the inner loop.  Repeat if iteration unsuccessful.
!
   IF (ratio < 0.0001D+00) THEN
      go to 240
   END IF
!
!  End of the outer loop.
!
   go to 30

340 CONTINUE
!
!  Termination, either normal or user imposed.
!
   IF (iflag < 0) THEN
      info = iflag
   END IF

   iflag = 0

   IF (0 < nprint) THEN
      CALL fcn(m, n, x, fvec, wa3, iflag)
   END IF

   RETURN
END SUBROUTINE lmstr
SUBROUTINE lmstr1(fcn, m, n, x, fvec, fjac, ldfjac, tol, info)

!*****************************************************************************80
!
!! LMSTR1 minimizes M functions in N variables using Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMSTR1 minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm
!    which uses minimal storage.
!
!    This is done by using the more general least-squares solver
!    LMSTR.  The user must provide a subroutine which calculates
!    the functions and the rows of the jacobian.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2016
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the rows of the jacobian.
!    FCN should have the form:
!      subroutine fcn ( m, n, x, fvec, fjrow, iflag )
!      integer ( kind = 4 ) m
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fjrow(n)
!      real ( kind = 8 ) fvec(m)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If the input value of IFLAG is I > 1, calculate the (I-1)-st row of
!    the jacobian at X, and return this vector in FJROW.
!    To terminate the algorithm, set the output value of IFLAG negative.
!
!    Input, integer ( kind = 4 ) M, the number of functions.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!    N must not exceed M.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(M), the functions evaluated at the output X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an N by N array.  The upper
!    triangle contains an upper triangular matrix R such that
!
!      P' * ( JAC' * JAC ) * P = R' * R,
!
!    where P is a permutation matrix and JAC is the final calculated
!    jacobian.  Column J of P is column IPVT(J) of the identity matrix.
!    The lower triangular part of FJAC contains information generated
!    during the computation of R.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of FJAC.
!    LDFJAC must be at least N.
!
!    Input, real ( kind = 8 ) TOL. Termination occurs when the algorithm
!    estimates either that the relative error in the sum of squares is at
!    most TOL or that the relative error between X and the solution is at
!    most TOL.  TOL should be nonnegative.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, algorithm estimates that the relative error in the sum of squares
!       is at most TOL.
!    2, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, FVEC is orthogonal to the columns of the jacobian to machine precision.
!    5, number of calls to FCN with IFLAG = 1 has reached 100*(N+1).
!    6, TOL is too small.  No further reduction in the sum of squares
!       is possible.
!    7, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!
   IMPLICIT NONE

   INTEGER(kind=4) ldfjac
   INTEGER(kind=4) m
   INTEGER(kind=4) n

   REAL(kind=8) diag(n)
   REAL(kind=8) factor
   EXTERNAL fcn
   REAL(kind=8) fjac(ldfjac, n)
   REAL(kind=8) ftol
   REAL(kind=8) fvec(m)
   REAL(kind=8) gtol
   INTEGER(kind=4) info
   INTEGER(kind=4) ipvt(n)
   INTEGER(kind=4) maxfev
   INTEGER(kind=4) mode
   INTEGER(kind=4) nfev
   INTEGER(kind=4) njev
   INTEGER(kind=4) nprint
   REAL(kind=8) qtf(n)
   REAL(kind=8) tol
   REAL(kind=8) x(n)
   REAL(kind=8) xtol

   IF (n <= 0) THEN
      info = 0
      RETURN
   END IF

   IF (m < n) THEN
      info = 0
      RETURN
   END IF

   IF (ldfjac < n) THEN
      info = 0
      RETURN
   END IF

   IF (tol < 0.0D+00) THEN
      info = 0
      RETURN
   END IF

   fvec(1:n) = 0.0D+00
   fjac(1:ldfjac, 1:n) = 0.0D+00
   ftol = tol
   xtol = tol
   gtol = 0.0D+00
   maxfev = 100*(n + 1)
   diag(1:n) = 0.0D+00
   mode = 1
   factor = 100.0D+00
   nprint = 0
   info = 0
   nfev = 0
   njev = 0
   ipvt(1:n) = 0
   qtf(1:n) = 0.0D+00

   CALL lmstr(fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
              diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf)

   IF (info == 8) THEN
      info = 4
   END IF

   RETURN
END SUBROUTINE lmstr1
SUBROUTINE qform(m, n, q, ldq)

!*****************************************************************************80
!
!! QFORM produces the explicit QR factorization of a matrix.
!
!  Discussion:
!
!    The QR factorization of a matrix is usually accumulated in implicit
!    form, that is, as a series of orthogonal transformations of the
!    original matrix.  This routine carries out those transformations,
!    to explicitly exhibit the factorization constructed by QRFAC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, is a positive integer input variable set
!    to the number of rows of A and the order of Q.
!
!    Input, integer ( kind = 4 ) N, is a positive integer input variable set
!    to the number of columns of A.
!
!    Input/output, real ( kind = 8 ) Q(LDQ,M).  Q is an M by M array.
!    On input the full lower trapezoid in the first min(M,N) columns of Q
!    contains the factored form.
!    On output, Q has been accumulated into a square matrix.
!
!    Input, integer ( kind = 4 ) LDQ, is a positive integer input variable
!    not less than M which specifies the leading dimension of the array Q.
!
   IMPLICIT NONE

   INTEGER(kind=4) ldq
   INTEGER(kind=4) m
   INTEGER(kind=4) n

   INTEGER(kind=4) j
   INTEGER(kind=4) k
   INTEGER(kind=4) l
   INTEGER(kind=4) minmn
   REAL(kind=8) q(ldq, m)
   REAL(kind=8) temp
   REAL(kind=8) wa(m)

   minmn = MIN(m, n)

   DO j = 2, minmn
      q(1:j - 1, j) = 0.0D+00
   END DO
!
!  Initialize remaining columns to those of the identity matrix.
!
   q(1:m, n + 1:m) = 0.0D+00

   DO j = n + 1, m
      q(j, j) = 1.0D+00
   END DO
!
!  Accumulate Q from its factored form.
!
   DO l = 1, minmn

      k = minmn - l + 1

      wa(k:m) = q(k:m, k)

      q(k:m, k) = 0.0D+00
      q(k, k) = 1.0D+00

      IF (wa(k) /= 0.0D+00) THEN

         DO j = k, m
            temp = DOT_PRODUCT(wa(k:m), q(k:m, j))/wa(k)
            q(k:m, j) = q(k:m, j) - temp*wa(k:m)
         END DO

      END IF

   END DO

   RETURN
END SUBROUTINE qform
SUBROUTINE qrfac(m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm)

!*****************************************************************************80
!
!! QRFAC computes a QR factorization using Householder transformations.
!
!  Discussion:
!
!    This subroutine uses Householder transformations with column
!    pivoting (optional) to compute a QR factorization of the
!    M by N matrix A.  That is, QRFAC determines an orthogonal
!    matrix Q, a permutation matrix P, and an upper trapezoidal
!    matrix R with diagonal elements of nonincreasing magnitude,
!    such that A*P = Q*R.  The Householder transformation for
!    column K, K = 1,2,...,min(M,N), is of the form
!
!      I - ( 1 / U(K) ) * U * U'
!
!    where U has zeros in the first K-1 positions.  The form of
!    this transformation and the method of pivoting first
!    appeared in the corresponding LINPACK subroutine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, real ( kind = 8 ) A(LDA,N), the M by N array.
!    On input, A contains the matrix for which the QR factorization is to
!    be computed.  On output, the strict upper trapezoidal part of A contains
!    the strict upper trapezoidal part of R, and the lower trapezoidal
!    part of A contains a factored form of Q (the non-trivial elements of
!    the U vectors described above).
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A, which must
!    be no less than M.
!
!    Input, logical PIVOT, is TRUE if column pivoting is to be carried out.
!
!    Output, integer ( kind = 4 ) IPVT(LIPVT), defines the permutation matrix P
!    such that A*P = Q*R.  Column J of P is column IPVT(J) of the identity
!    matrix.  If PIVOT is false, IPVT is not referenced.
!
!    Input, integer ( kind = 4 ) LIPVT, the dimension of IPVT, which should
!    be N if pivoting is used.
!
!    Output, real ( kind = 8 ) RDIAG(N), contains the diagonal elements of R.
!
!    Output, real ( kind = 8 ) ACNORM(N), the norms of the corresponding
!    columns of the input matrix A.  If this information is not needed,
!    then ACNORM can coincide with RDIAG.
!
   IMPLICIT NONE

   INTEGER(kind=4) lda
   INTEGER(kind=4) lipvt
   INTEGER(kind=4) m
   INTEGER(kind=4) n

   REAL(kind=8) a(lda, n)
   REAL(kind=8) acnorm(n)
   REAL(kind=8) ajnorm
   REAL(kind=8) enorm
   REAL(kind=8) epsmch
   ! integer ( kind = 4 ) i
   INTEGER(kind=4) i4_temp
   INTEGER(kind=4) ipvt(lipvt)
   INTEGER(kind=4) j
   INTEGER(kind=4) k
   INTEGER(kind=4) kmax
   INTEGER(kind=4) minmn
   LOGICAL pivot
   REAL(kind=8) r8_temp(m)
   REAL(kind=8) rdiag(n)
   REAL(kind=8) temp
   REAL(kind=8) wa(n)

   epsmch = EPSILON(epsmch)
!
!  Compute the initial column norms and initialize several arrays.
!
   DO j = 1, n
      acnorm(j) = enorm(m, a(1:m, j))
   END DO

   rdiag(1:n) = acnorm(1:n)
   wa(1:n) = acnorm(1:n)

   IF (pivot) THEN
      DO j = 1, n
         ipvt(j) = j
      END DO
   END IF
!
!  Reduce A to R with Householder transformations.
!
   minmn = MIN(m, n)

   DO j = 1, minmn
!
!  Bring the column of largest norm into the pivot position.
!
      IF (pivot) THEN

         kmax = j

         DO k = j, n
            IF (rdiag(kmax) < rdiag(k)) THEN
               kmax = k
            END IF
         END DO

         IF (kmax /= j) THEN

            r8_temp(1:m) = a(1:m, j)
            a(1:m, j) = a(1:m, kmax)
            a(1:m, kmax) = r8_temp(1:m)

            rdiag(kmax) = rdiag(j)
            wa(kmax) = wa(j)

            i4_temp = ipvt(j)
            ipvt(j) = ipvt(kmax)
            ipvt(kmax) = i4_temp

         END IF

      END IF
!
!  Compute the Householder transformation to reduce the
!  J-th column of A to a multiple of the J-th unit vector.
!
      ajnorm = enorm(m - j + 1, a(j, j))

      IF (ajnorm /= 0.0D+00) THEN

         IF (a(j, j) < 0.0D+00) THEN
            ajnorm = -ajnorm
         END IF

         a(j:m, j) = a(j:m, j)/ajnorm
         a(j, j) = a(j, j) + 1.0D+00
!
!  Apply the transformation to the remaining columns and update the norms.
!
         DO k = j + 1, n

            temp = DOT_PRODUCT(a(j:m, j), a(j:m, k))/a(j, j)

            a(j:m, k) = a(j:m, k) - temp*a(j:m, j)

            IF (pivot .AND. rdiag(k) /= 0.0D+00) THEN

               temp = a(j, k)/rdiag(k)
               rdiag(k) = rdiag(k)*SQRT(MAX(0.0D+00, 1.0D+00 - temp**2))

               IF (0.05D+00*(rdiag(k)/wa(k))**2 <= epsmch) THEN
                  rdiag(k) = enorm(m - j, a(j + 1, k))
                  wa(k) = rdiag(k)
               END IF

            END IF

         END DO

      END IF

      rdiag(j) = -ajnorm

   END DO

   RETURN
END SUBROUTINE qrfac
SUBROUTINE qrsolv(n, r, ldr, ipvt, diag, qtb, x, sdiag)

!*****************************************************************************80
!
!! QRSOLV solves a rectangular linear system A*x=b in the least squares sense.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N diagonal matrix D,
!    and an M-vector B, the problem is to determine an X which
!    solves the system
!
!      A*X = B
!      D*X = 0
!
!    in the least squares sense.
!
!    This subroutine completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization, with column pivoting, of A.  That is, if
!    Q*P = Q*R, where P is a permutation matrix, Q has orthogonal
!    columns, and R is an upper triangular matrix with diagonal
!    elements of nonincreasing magnitude, then QRSOLV expects
!    the full upper triangle of R, the permutation matrix p,
!    and the first N components of Q'*B.
!
!    The system is then equivalent to
!
!      R*Z = Q'*B
!      P'*D*P*Z = 0
!
!    where X = P*Z.  If this system does not have full rank,
!    then a least squares solution is obtained.  On output QRSOLV
!    also provides an upper triangular matrix S such that
!
!      P'*(A'*A + D*D)*P = S'*S.
!
!    S is computed within QRSOLV and may be of separate interest.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of R.
!
!    Input/output, real ( kind = 8 ) R(LDR,N), the N by N matrix.
!    On input the full upper triangle must contain the full upper triangle
!    of the matrix R.  On output the full upper triangle is unaltered, and
!    the strict lower triangle contains the strict upper triangle
!    (transposed) of the upper triangular matrix S.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of R, which must be
!    at least N.
!
!    Input, integer ( kind = 4 ) IPVT(N), defines the permutation matrix P such
!    that A*P = Q*R.  Column J of P is column IPVT(J) of the identity matrix.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real ( kind = 8 ) QTB(N), the first N elements of the vector Q'*B.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution.
!
!    Output, real ( kind = 8 ) SDIAG(N), the diagonal elements of the upper
!    triangular matrix S.
!
   IMPLICIT NONE

   INTEGER(kind=4) ldr
   INTEGER(kind=4) n

   REAL(kind=8) c
   REAL(kind=8) cotan
   REAL(kind=8) diag(n)
   INTEGER(kind=4) i
   INTEGER(kind=4) ipvt(n)
   INTEGER(kind=4) j
   INTEGER(kind=4) k
   INTEGER(kind=4) l
   INTEGER(kind=4) nsing
   REAL(kind=8) qtb(n)
   REAL(kind=8) qtbpj
   REAL(kind=8) r(ldr, n)
   REAL(kind=8) s
   REAL(kind=8) sdiag(n)
   REAL(kind=8) sum2
   REAL(kind=8) t
   REAL(kind=8) temp
   REAL(kind=8) wa(n)
   REAL(kind=8) x(n)
!
!  Copy R and Q'*B to preserve input and initialize S.
!
!  In particular, save the diagonal elements of R in X.
!
   DO j = 1, n
      r(j:n, j) = r(j, j:n)
      x(j) = r(j, j)
   END DO

   wa(1:n) = qtb(1:n)
!
!  Eliminate the diagonal matrix D using a Givens rotation.
!
   DO j = 1, n
!
!  Prepare the row of D to be eliminated, locating the
!  diagonal element using P from the QR factorization.
!
      l = ipvt(j)

      IF (diag(l) /= 0.0D+00) THEN

         sdiag(j:n) = 0.0D+00
         sdiag(j) = diag(l)
!
!  The transformations to eliminate the row of D
!  modify only a single element of Q'*B
!  beyond the first N, which is initially zero.
!
         qtbpj = 0.0D+00

         DO k = j, n
!
!  Determine a Givens rotation which eliminates the
!  appropriate element in the current row of D.
!
            IF (sdiag(k) /= 0.0D+00) THEN

               IF (ABS(r(k, k)) < ABS(sdiag(k))) THEN
                  cotan = r(k, k)/sdiag(k)
                  s = 0.5D+00/SQRT(0.25D+00 + 0.25D+00*cotan**2)
                  c = s*cotan
               ELSE
                  t = sdiag(k)/r(k, k)
                  c = 0.5D+00/SQRT(0.25D+00 + 0.25D+00*t**2)
                  s = c*t
               END IF
!
!  Compute the modified diagonal element of R and
!  the modified element of (Q'*B,0).
!
               r(k, k) = c*r(k, k) + s*sdiag(k)
               temp = c*wa(k) + s*qtbpj
               qtbpj = -s*wa(k) + c*qtbpj
               wa(k) = temp
!
!  Accumulate the tranformation in the row of S.
!
               DO i = k + 1, n
                  temp = c*r(i, k) + s*sdiag(i)
                  sdiag(i) = -s*r(i, k) + c*sdiag(i)
                  r(i, k) = temp
               END DO

            END IF

         END DO

      END IF
!
!  Store the diagonal element of S and restore
!  the corresponding diagonal element of R.
!
      sdiag(j) = r(j, j)
      r(j, j) = x(j)

   END DO
!
!  Solve the triangular system for Z.  If the system is
!  singular, then obtain a least squares solution.
!
   nsing = n

   DO j = 1, n

      IF (sdiag(j) == 0.0D+00 .AND. nsing == n) THEN
         nsing = j - 1
      END IF

      IF (nsing < n) THEN
         wa(j) = 0.0D+00
      END IF

   END DO

   DO j = nsing, 1, -1
      sum2 = DOT_PRODUCT(wa(j + 1:nsing), r(j + 1:nsing, j))
      wa(j) = (wa(j) - sum2)/sdiag(j)
   END DO
!
!  Permute the components of Z back to components of X.
!
   DO j = 1, n
      l = ipvt(j)
      x(l) = wa(j)
   END DO

   RETURN
END SUBROUTINE qrsolv
SUBROUTINE r1mpyq(m, n, a, lda, v, w)

!*****************************************************************************80
!
!! R1MPYQ computes A*Q, where Q is the product of Householder transformations.
!
!  Discussion:
!
!    Given an M by N matrix A, this subroutine computes A*Q where
!    Q is the product of 2*(N - 1) transformations
!
!      GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
!
!    and GV(I), GW(I) are Givens rotations in the (I,N) plane which
!    eliminate elements in the I-th and N-th planes, respectively.
!    Q itself is not given, rather the information to recover the
!    GV, GW rotations is supplied.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, real ( kind = 8 ) A(LDA,N), the M by N array.
!    On input, the matrix A to be postmultiplied by the orthogonal matrix Q.
!    On output, the value of A*Q.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A, which must not
!    be less than M.
!
!    Input, real ( kind = 8 ) V(N), W(N), contain the information necessary
!    to recover the Givens rotations GV and GW.
!
   IMPLICIT NONE

   INTEGER(kind=4) lda
   INTEGER(kind=4) m
   INTEGER(kind=4) n

   REAL(kind=8) a(lda, n)
   REAL(kind=8) c
   INTEGER(kind=4) i
   INTEGER(kind=4) j
   REAL(kind=8) s
   REAL(kind=8) temp
   REAL(kind=8) v(n)
   REAL(kind=8) w(n)
!
!  Apply the first set of Givens rotations to A.
!
   DO j = n - 1, 1, -1

      IF (1.0D+00 < ABS(v(j))) THEN
         c = 1.0D+00/v(j)
         s = SQRT(1.0D+00 - c**2)
      ELSE
         s = v(j)
         c = SQRT(1.0D+00 - s**2)
      END IF

      DO i = 1, m
         temp = c*a(i, j) - s*a(i, n)
         a(i, n) = s*a(i, j) + c*a(i, n)
         a(i, j) = temp
      END DO

   END DO
!
!  Apply the second set of Givens rotations to A.
!
   DO j = 1, n - 1

      IF (ABS(w(j)) > 1.0D+00) THEN
         c = 1.0D+00/w(j)
         s = SQRT(1.0D+00 - c**2)
      ELSE
         s = w(j)
         c = SQRT(1.0D+00 - s**2)
      END IF

      DO i = 1, m
         temp = c*a(i, j) + s*a(i, n)
         a(i, n) = -s*a(i, j) + c*a(i, n)
         a(i, j) = temp
      END DO

   END DO

   RETURN
END SUBROUTINE r1mpyq
SUBROUTINE r1updt(m, n, s, ls, u, v, w, sing)

!*****************************************************************************80
!
!! R1UPDT re-triangularizes a matrix after a rank one update.
!
!  Discussion:
!
!    Given an M by N lower trapezoidal matrix S, an M-vector U, and an
!    N-vector V, the problem is to determine an orthogonal matrix Q such that
!
!      (S + U * V' ) * Q
!
!    is again lower trapezoidal.
!
!    This subroutine determines Q as the product of 2 * (N - 1)
!    transformations
!
!      GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
!
!    where GV(I), GW(I) are Givens rotations in the (I,N) plane
!    which eliminate elements in the I-th and N-th planes,
!    respectively.  Q itself is not accumulated, rather the
!    information to recover the GV and GW rotations is returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of S.
!
!    Input, integer ( kind = 4 ) N, the number of columns of S.
!    N must not exceed M.
!
!    Input/output, real ( kind = 8 ) S(LS).  On input, the lower trapezoidal
!    matrix S stored by columns.  On output S contains the lower trapezoidal
!    matrix produced as described above.
!
!    Input, integer ( kind = 4 ) LS, the length of the S array.  LS must be at
!    least (N*(2*M-N+1))/2.
!
!    Input, real ( kind = 8 ) U(M), the U vector.
!
!    Input/output, real ( kind = 8 ) V(N).  On input, V must contain the
!    vector V.  On output V contains the information necessary to recover the
!    Givens rotations GV described above.
!
!    Output, real ( kind = 8 ) W(M), contains information necessary to
!    recover the Givens rotations GW described above.
!
!    Output, logical SING, is set to TRUE if any of the diagonal elements
!    of the output S are zero.  Otherwise SING is set FALSE.
!
   IMPLICIT NONE

   INTEGER(kind=4) ls
   INTEGER(kind=4) m
   INTEGER(kind=4) n

   REAL(kind=8) cos
   REAL(kind=8) cotan
   REAL(kind=8) giant
   INTEGER(kind=4) i
   INTEGER(kind=4) j
   INTEGER(kind=4) jj
   INTEGER(kind=4) l
   REAL(kind=8) s(ls)
   REAL(kind=8) sin
   LOGICAL sing
   REAL(kind=8) tan
   REAL(kind=8) tau
   REAL(kind=8) temp
   REAL(kind=8) u(m)
   REAL(kind=8) v(n)
   REAL(kind=8) w(m)
!
!  GIANT is the largest magnitude.
!
   giant = HUGE(giant)
!
!  Initialize the diagonal element pointer.
!
   jj = (n*(2*m - n + 1))/2 - (m - n)
!
!  Move the nontrivial part of the last column of S into W.
!
   l = jj
   DO i = n, m
      w(i) = s(l)
      l = l + 1
   END DO
!
!  Rotate the vector V into a multiple of the N-th unit vector
!  in such a way that a spike is introduced into W.
!
   DO j = n - 1, 1, -1

      jj = jj - (m - j + 1)
      w(j) = 0.0D+00

      IF (v(j) /= 0.0D+00) THEN
!
!  Determine a Givens rotation which eliminates the
!  J-th element of V.
!
         IF (ABS(v(n)) < ABS(v(j))) THEN
            cotan = v(n)/v(j)
            sin = 0.5D+00/SQRT(0.25D+00 + 0.25D+00*cotan**2)
            cos = sin*cotan
            tau = 1.0D+00
            IF (ABS(cos)*giant > 1.0D+00) THEN
               tau = 1.0D+00/cos
            END IF
         ELSE
            tan = v(j)/v(n)
            cos = 0.5D+00/SQRT(0.25D+00 + 0.25D+00*tan**2)
            sin = cos*tan
            tau = sin
         END IF
!
!  Apply the transformation to V and store the information
!  necessary to recover the Givens rotation.
!
         v(n) = sin*v(j) + cos*v(n)
         v(j) = tau
!
!  Apply the transformation to S and extend the spike in W.
!
         l = jj
         DO i = j, m
            temp = cos*s(l) - sin*w(i)
            w(i) = sin*s(l) + cos*w(i)
            s(l) = temp
            l = l + 1
         END DO

      END IF

   END DO
!
!  Add the spike from the rank 1 update to W.
!
   w(1:m) = w(1:m) + v(n)*u(1:m)
!
!  Eliminate the spike.
!
   sing = .FALSE.

   DO j = 1, n - 1

      IF (w(j) /= 0.0D+00) THEN
!
!  Determine a Givens rotation which eliminates the
!  J-th element of the spike.
!
         IF (ABS(s(jj)) < ABS(w(j))) THEN

            cotan = s(jj)/w(j)
            sin = 0.5D+00/SQRT(0.25D+00 + 0.25D+00*cotan**2)
            cos = sin*cotan

            IF (1.0D+00 < ABS(cos)*giant) THEN
               tau = 1.0D+00/cos
            ELSE
               tau = 1.0D+00
            END IF

         ELSE

            tan = w(j)/s(jj)
            cos = 0.5D+00/SQRT(0.25D+00 + 0.25D+00*tan**2)
            sin = cos*tan
            tau = sin

         END IF
!
!  Apply the transformation to S and reduce the spike in W.
!
         l = jj
         DO i = j, m
            temp = cos*s(l) + sin*w(i)
            w(i) = -sin*s(l) + cos*w(i)
            s(l) = temp
            l = l + 1
         END DO
!
!  Store the information necessary to recover the Givens rotation.
!
         w(j) = tau

      END IF
!
!  Test for zero diagonal elements in the output S.
!
      IF (s(jj) == 0.0D+00) THEN
         sing = .TRUE.
      END IF

      jj = jj + (m - j + 1)

   END DO
!
!  Move W back into the last column of the output S.
!
   l = jj
   DO i = n, m
      s(l) = w(i)
      l = l + 1
   END DO

   IF (s(jj) == 0.0D+00) THEN
      sing = .TRUE.
   END IF

   RETURN
END SUBROUTINE r1updt
SUBROUTINE r8vec_print(n, a, title)

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
   IMPLICIT NONE

   INTEGER(kind=4) n

   REAL(kind=8) a(n)
   INTEGER(kind=4) i
   CHARACTER(len=*) title

   WRITE (*, '(a)') ' '
   WRITE (*, '(a)') TRIM(title)
   WRITE (*, '(a)') ' '
   DO i = 1, n
      WRITE (*, '(2x,i8,2x,g16.8)') i, a(i)
   END DO

   RETURN
END SUBROUTINE r8vec_print
SUBROUTINE rwupdt(n, r, ldr, w, b, alpha, c, s)

!*****************************************************************************80
!
!! RWUPDT computes the decomposition of triangular matrix augmented by one row.
!
!  Discussion:
!
!    Given an N by N upper triangular matrix R, this subroutine
!    computes the QR decomposition of the matrix formed when a row
!    is added to R.  If the row is specified by the vector W, then
!    RWUPDT determines an orthogonal matrix Q such that when the
!    N+1 by N matrix composed of R augmented by W is premultiplied
!    by Q', the resulting matrix is upper trapezoidal.
!    The matrix Q' is the product of N transformations
!
!      G(N)*G(N-1)* ... *G(1)
!
!    where G(I) is a Givens rotation in the (I,N+1) plane which eliminates
!    elements in the (N+1)-st plane.  RWUPDT also computes the product
!    Q'*C where C is the (N+1)-vector (B,ALPHA).  Q itself is not
!    accumulated, rather the information to recover the G rotations is
!    supplied.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of R.
!
!    Input/output, real ( kind = 8 ) R(LDR,N).  On input the upper triangular
!    part of R must contain the matrix to be updated.  On output R contains the
!    updated triangular matrix.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of the array R.
!    LDR must not be less than N.
!
!    Input, real ( kind = 8 ) W(N), the row vector to be added to R.
!
!    Input/output, real ( kind = 8 ) B(N).  On input, the first N elements
!    of the vector C.  On output the first N elements of the vector Q'*C.
!
!    Input/output, real ( kind = 8 ) ALPHA.  On input, the (N+1)-st element
!    of the vector C.  On output the (N+1)-st element of the vector Q'*C.
!
!    Output, real ( kind = 8 ) C(N), S(N), the cosines and sines of the
!    transforming Givens rotations.
!
   IMPLICIT NONE

   INTEGER(kind=4) ldr
   INTEGER(kind=4) n

   REAL(kind=8) alpha
   REAL(kind=8) b(n)
   REAL(kind=8) c(n)
   REAL(kind=8) cotan
   INTEGER(kind=4) i
   INTEGER(kind=4) j
   REAL(kind=8) r(ldr, n)
   REAL(kind=8) rowj
   REAL(kind=8) s(n)
   REAL(kind=8) tan
   REAL(kind=8) temp
   REAL(kind=8) w(n)

   DO j = 1, n

      rowj = w(j)
!
!  Apply the previous transformations to R(I,J), I=1,2,...,J-1, and to W(J).
!
      DO i = 1, j - 1
         temp = c(i)*r(i, j) + s(i)*rowj
         rowj = -s(i)*r(i, j) + c(i)*rowj
         r(i, j) = temp
      END DO
!
!  Determine a Givens rotation which eliminates W(J).
!
      c(j) = 1.0D+00
      s(j) = 0.0D+00

      IF (rowj /= 0.0D+00) THEN

         IF (ABS(r(j, j)) < ABS(rowj)) THEN
            cotan = r(j, j)/rowj
            s(j) = 0.5D+00/SQRT(0.25D+00 + 0.25D+00*cotan**2)
            c(j) = s(j)*cotan
         ELSE
            tan = rowj/r(j, j)
            c(j) = 0.5D+00/SQRT(0.25D+00 + 0.25D+00*tan**2)
            s(j) = c(j)*tan
         END IF
!
!  Apply the current transformation to R(J,J), B(J), and ALPHA.
!
         r(j, j) = c(j)*r(j, j) + s(j)*rowj
         temp = c(j)*b(j) + s(j)*alpha
         alpha = -s(j)*b(j) + c(j)*alpha
         b(j) = temp

      END IF

   END DO

   RETURN
END SUBROUTINE rwupdt
SUBROUTINE timestamp()

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
   IMPLICIT NONE

   CHARACTER(len=8) ampm
   INTEGER(kind=4) d
   INTEGER(kind=4) h
   INTEGER(kind=4) m
   INTEGER(kind=4) mm
   CHARACTER(len=9), PARAMETER, DIMENSION(12) :: month = (/ &
                                                 'January  ', 'February ', 'March    ', 'April    ', &
                                                 'May      ', 'June     ', 'July     ', 'August   ', &
                                                 'September', 'October  ', 'November ', 'December '/)
   INTEGER(kind=4) n
   INTEGER(kind=4) s
   INTEGER(kind=4) values(8)
   INTEGER(kind=4) y

   CALL DATE_AND_TIME(values=values)

   y = values(1)
   m = values(2)
   d = values(3)
   h = values(5)
   n = values(6)
   s = values(7)
   mm = values(8)

   IF (h < 12) THEN
      ampm = 'AM'
   ELSE IF (h == 12) THEN
      IF (n == 0 .AND. s == 0) THEN
         ampm = 'Noon'
      ELSE
         ampm = 'PM'
      END IF
   ELSE
      h = h - 12
      IF (h < 12) THEN
         ampm = 'PM'
      ELSE IF (h == 12) THEN
         IF (n == 0 .AND. s == 0) THEN
            ampm = 'Midnight'
         ELSE
            ampm = 'AM'
         END IF
      END IF
   END IF

   WRITE (*, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)') &
      d, TRIM(month(m)), y, h, ':', n, ':', s, '.', mm, TRIM(ampm)

   RETURN
END SUBROUTINE timestamp
