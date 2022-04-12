%Performs rotation of Stress Vector (2D)

function vdash=rot_stress(v,t)

c=cosd(t); s=sind(t);

R=[c^2,     s^2,    2*s*c;
   s^2,     c^2,   -2*s*c;
  -s*c,     s*c,  c^2-s^2];

vdash=R*v;

end
