c  Copyright (C) 1997, 1998
c  Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
c  $Id: aind.f,v 1.3 1998/07/23 05:06:32 bturlach Exp $
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
c  this routine checks whether Aind has valid entries, i.e., 
c    1) 1<= Aind(1,i) <= n for i=1,...,q (number of constraints)
c    2) 1<= Aind(j,i) <= n for j=2,...,Aind(1,i)+1, i=1,...,q
c
c  Aind is a m times q matrix
c
      subroutine aind(ind,m,q,n,ok)
      implicit none
      integer m, ind(m,*), q, n, i, j
      logical ok
      ok = .FALSE.
      do i=1,q
         if( ind(1,i) .LT. 1 .OR. ind(1,i) .GT. n ) return
         do j=2,ind(1,i)+1
            if( ind(j,i) .LT. 1 .OR. ind(j,i) .GT. n ) return
         enddo
      enddo
      ok = .TRUE.
      return
      end
