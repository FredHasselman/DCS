      subroutine jspd(m, n, k, Adata, jb, Bdata, d, tv, na, nb, nd)
C   Central part of spdiags for R
C     m and n are row and column sizes of A (underlying matrix)
C     jb will be number of returned diagonals
C     returns jb, Bdata, d
      integer m, n, na, nb, nd, jb, d(nd)
      integer i, j, k, kend, q, mn, kk1, js, je, qx
      double precision Adata(na), Bdata(nb), tv(k)
      LOGICAL not0
C      k = min(m, n)
C ??  check if k=1
C      Bdata<-NULL # start with nothing in B matrix (as vector)
      jb = 0 
      kend=0
C     column index of "last" column saved for B
C   d = NULL # index vector of diagonals from A 
C   # d contains 0 for the principal diagonal, -i for i'th lower
C   # diagonal (prefaced with zeros), +j for j'th upper diagonal
C   # (suffixed by zeros)
      q = (m-1)+n 
C  There are m-1 subdiagonals and n-1 superdiagonals + main diagonal
C   assume we have already built Adata for tall or fat matrix
C   # Augment the data with columns of zeros fore and aft
      mn = m*n
C      print *,"Original Adata"
C      print 1000, (Adata(i), i=1,mn)
      do 10 i=1,mn
        Bdata(i) = Adata(i)
 10   continue
      kk1 = (k-1)*k
      do 15 i=1,kk1
         Adata(i)=0.0
 15   continue
      js = kk1+1
      je = kk1+mn
      do 20 i=js,je
         Adata(i) = Bdata(i-js+1)
 20   continue
      js = je+1
      je = je+kk1
      do 25 i=js,je
         Adata(i) = 0.0
 25   continue
C      print *, "Augmented Adata with ", je,"  elements"
C      print 1000,(Adata(i), i=1,je)
 1000 FORMAT(1H ,25f4.0)
      do 100 i=1,q
          qx = k*(i-1)+1
C          print *,"Top element is no ",qx,"=",Adata(qx)
          tv(1) = Adata(qx)
          do 30 j=1,(k-1)
             tv(j+1) = 0.0
 30       continue
          qx = min(q-i, k-1) 
C          print *,"qx=",qx
          if (qx .gt. 0) then
C  qx will be 0 when we are at last superdiagonal,i.e., i == q
             do 35 j=1,qx
                tv(j+1) =  Adata(k*(i+j-1)+j+1)
 35          continue
          endif
          not0 = .FALSE.
          do 40 j=1,k
             if (tv(j) .ne. 0.0) not0 = .TRUE.
 40       continue
C          if (.NOT. not0) print *," Zero column for i=",i
          if (not0) then 
C             print *,"Nonzero column for i=",i
             jb = jb+1
             d(jb) = (i-k)
C             print *, " New jb=",jb,"   d element =",d(jb)
C             print 1000,(tv(j), j=1,k)
             do 45 j=1,k
                Bdata(kend+j) = tv(j)
 45          continue
             kend = kend + k
C  save the diagonal as column of B in vector form
          endif
 100  continue
C      print *,"Bdata ", jb, " columns"
C      qx = k*jb
C      print 1000,(Bdata(j), j=1,qx)
      return
      end
