*
* File: dilog.f (from SR)
*
       double precision function dilog(x)
       implicit double precision  (a-z)
       dimension b(8)
       integer ncall
       data ncall/0/,pi6/1.644934066848226d+00/,een,vier/1.d+00,.25d+00/
       ncall = 0
       if(ncall.eq.0)go to 2
1      if(x.lt.0)go to 3
       if(x.gt.0.5)go to 4
       z=-dlog(1.-x)
7      z2=z*z
       dilog=z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een)-z2*vier
       if(x.gt.een)dilog=-dilog-.5*u*u+2.*pi6
       return
2      b(1)=een
       b(2)=een/36.
       b(3)=-een/3600.
       b(4)=een/211680.
       b(5)=-een/(30.*362880.d+00)
       b(6)=5./(66.*39916800.d+00)
       b(7)=-691./(2730.*39916800.d+00*156.)
       b(8)=een/(39916800.d+00*28080.)
       ncall=1
       go to 1
3      if(x.gt.-een)go to 5
       y=een/(een-x)
       z=-dlog(een-y)
       z2=z*z
       u=dlog(y)
       dilog=z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een)-z2*vier-u*(z+.5*u)-pi6
       return
4      if(x.ge.een)go to 10
       y=een-x
       z=-dlog(x)
6      u=dlog(y)
       z2=z*z
       dilog=-z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een-u)+z2*vier+pi6
       if(x.gt.een)dilog=-dilog-.5*z*z+pi6*2.
       return
5      y=een/(een-x)
       z=-dlog(y)
       z2=z*z
       dilog=-z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b(8)+b(7))+b(6))
     1 +b(5))+b(4))+b(3))+b(2))+een)-z2*vier
       return
10     if(x.eq.een)go to 20
       xx=1./x
       if(x.gt.2.)go to 11
       z=dlog(x)
       y=1.-xx
       go to 6
11     u=dlog(x)
       z=-dlog(1.-xx)
       go to 7
20     dilog=pi6
       return
       end
