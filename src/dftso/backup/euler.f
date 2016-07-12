Subroutine Euler(Rot_new,a,b,c)
  Implicit REAL*8 (A-H,O-Z)
  DIMENSION Rot_old(3,3),Rot_new(3,3),z(3),zz(3)
  DIMENSION y(3),yy(3),yyy(3),pom(3),x(3),xx(3)

  zero=0d0
  pi=acos(-1d0)
  do i=1,3
     do j=1,3
	Rot_old(i,j)=0
        if (i.eq.j) Rot_old(i,i)=1
     end do
  end do
3 format(20X,3(F10.8,2X))

  do j=1,3
     y(j)=Rot_old(j,2)
     yyy(j)=Rot_new(j,2)
     z(j)=Rot_old(j,3)
     zz(j)=Rot_new(j,3)
  end do

  call vecprod(z,zz,yy)
  y_norm=dsqrt(dot(yy,yy))

  if (y_norm.lt.1d-10) then
     !	if (abs(dot(y,yyy)).gt.1d0) then
     !	aa=dot(y,yyy)/abs(dot(y,yyy))
     !	a=acos(aa)
     !	else
     !	a=acos(dot(y,yyy))
     !	end if
     !
     !......modified 12.2.2002
     a=ang(y,yyy)

     if (dot(z,zz).gt.zero) then
	c=zero
        b=zero
        if (yyy(1).gt.zero) a=2*pi-a
     else
	c=a
        a=zero
	b=pi
        if (yyy(1).lt.zero) c=2*pi-c
     end if
  else
     do j=1,3
	yy(j)=yy(j)/y_norm
     end do

     !	aa=dot(y,yy)
     !	bb=dot(z,zz)
     !	cc=dot(yy,yyy)
     !	if (abs(aa).gt.1d0) aa=aa/abs(aa)
     !       if (abs(bb).gt.1d0) bb=bb/abs(bb)
     !       if (abs(cc).gt.1d0) cc=cc/abs(cc)
     !	b=acos(bb)
     !	a=acos(aa)
     !	c=acos(cc)
     !
     !.....modifeied 12.2.2002
     a=ang(y,yy)
     b=ang(z,zz)
     c=ang(yy,yyy)
     if (yy(1).gt.zero) a=2*pi-a
     call vecprod(yy,yyy,pom)
     if (dot(pom,zz).lt.zero) c=2*pi-c
  end if
END Subroutine Euler


Subroutine vecprod(a,b,c)
  IMPLICIT REAL*8 (A-H)
  dimension a(3),b(3),c(3)
  
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=a(3)*b(1)-a(1)*b(3)
  c(3)=a(1)*b(2)-a(2)*b(1)
end Subroutine vecprod

REAL*8 Function dot(a,b)
  IMPLICIT REAL*8 (A-H)
  dimension a(3),b(3)
  dot=0
  do i=1,3
     dot=dot+a(i)*b(i)
  end do
end Function dot

REAL*8 function ang(a,b)
  IMPLICIT REAL*8 (A-H,O-Z)
  dimension a(3),b(3),c(3)
  
  pi=dacos(-1.d0)
  aa=dsqrt(dot(a,a))
  bb=dsqrt(dot(b,b))
  e=dot(a,b)/(aa*bb)
  if (abs(e).gt.0.8) then
     call vecprod(a,b,c)
     cc=dsqrt(dot(c,c))
     ang=dasin(cc/(aa*bb))
     if (e.lt.0.d0) ang=pi-ang
  else
     ang=dacos(e)
  end if
end function ang
