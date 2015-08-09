c  Copyright (C) 1998
c  Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
c  $Id: achck.f,v 1.2 1998/07/23 05:05:58 bturlach Exp $
c 
c  This library is free software; you can redistribute it and/or
c  modify it under the terms of the GNU Library General Public
c  License as published by the Free Software Foundation; either
c  version 2 of the License, or (at your option) any later version.
c  
c  This library is distributed in the hope that it will be useful,
c  but WITHOUT ANY WARRANTY; without even the implied warranty of
c  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c  Library General Public License for more details.
c  
c  You should have received a copy of the GNU Library General Public
c  License along with this library; if not, write to the Free Software
c  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
c  MA 02111-1307 USA
c 
c  this routine checks whether all constraints are fulfilled.
c
      subroutine achck(sol,n,amat,aind,bvec,m,q,meq,prec,ok)
      implicit none
      integer m, aind(m+1,*), q, meq, n, i, j
      double precision sol(*), amat(m,*), bvec(*), sum, prec
      logical ok
      ok = .FALSE.
      do i=1,q
         sum = -bvec(i)
         do j=1,aind(1,i)
            sum = sum + amat(j,i)*sol(aind(j+1,i))
         enddo
         if( i .LE. meq ) sum = -abs(sum)
         if( sum .LT. -prec ) return
      enddo
      ok = .TRUE.
      return
      end
