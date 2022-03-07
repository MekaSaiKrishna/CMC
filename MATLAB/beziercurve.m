% Bezier Curve Cubic Area

syms r1 r2 t real
P0=[0 1];   B0=(1-t)^3;
P1=[r1 r1]; B1=3*(1-t)^2*t;
P2=[r2 r2]; B2=3*(1-t)*t^2;
P3=[1 0];   B3=t^3;

curve = P0*B0 + P1*B1 + P2*B2 + P3*B3

x = curve(1,1)
y = curve(1,2)

dx = diff(x,t)
Area = int(y*dx,t,0,1)
