B(x, t, g)
   g(1)=0.2D0*x(1)*x(2)+0.125142D1*cos(t)-0.911111D0*sin(t)-0.2D0*sin(t)*cos(t)
   g(2)=0.3D0*x(1)**2-0.1D0*x(2)**2-0.748585D0*sin(t)-0.288888D0*cos(t)-0.3D0*sin(t)**2+0.1D0*cos(t)**2

A(eig)
   cg1(1,1)=0.222222D-1+0.888889D0*eig
   cg1(1,2)=0.628539D-1-0.314269D0*eig
   cg1(2,1)=0.628536D-1-0.314269D0*eig
   cg1(2,2)=0.177777D0+0.111111D0*eig
   cg1=-cg1
   
