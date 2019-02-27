*
* File: locate.f (from SR)
*
c routine taken out of Numerical Recipes
      Subroutine Locate(xx,n,x,j)
      Integer j,n
      Double Precision x,xx(n)
      Integer jl,ju,jm
      jl = 0
      ju = n+1
 10   If (ju - jl .gt. 1) then
         jm = (ju + jl)/2
c         If ((xx(n) .gt. xx(1)) .eqv. (x .gt. xx(jm))) then
         If ((xx(n) .gt. xx(1)) .and. (x .gt. xx(jm))) then
            jl = jm
         else
            ju = jm
         endif
         goto 10
      endif
      j = jl
      return
      End

         
