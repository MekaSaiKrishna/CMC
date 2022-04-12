%Finds Qbar matrices

%NOTE: The inputs E1,E2,G12,nu12 must be manually entered in this file before proceeding to the Main File

function QA=Qbar(t)

c=cosd(t); s=sind(t);

E1A=158490e6;    E2A=9017e6; G12A=5703e6;    v12A=0.316;

v21A=E2A*v12A/E1A;

Q11=E1A/(1-(v12A*v21A)); 
Q12=(v12A*E2A)/(1-(v12A*v21A));
Q22=E2A/(1-(v12A*v21A));
Q66=G12A;

% [Q11,Q12,Q22,Q16,Q26,Q66]

QA(1,1) = (Q11*c^4) + 2*(Q12+2*Q66)*(s^2)*(c^2) + Q22*(s^4);
QA(1,2) = ((Q11+Q22-4*Q66)*(s^2)*(c^2)) + Q12*((s^4)+(c^4));
QA(2,2) = (Q11*s^4) + (2*(Q12 + 2*Q66)*(s^2*c^2)) + (Q22*c^4);
QA(1,3) = ((Q11-Q12-2*Q66)*s*c^3) + ((Q12 - Q22 + 2*Q66)*s^3*c);
QA(2,3) = ((Q11-Q12-2*Q66)*c*s^3) + ((Q12 - Q22 + 2*Q66)*c^3*s);
QA(3,3) = ((Q11 + Q22 - 2*Q12 - 2*Q66)*s^2*c^2) + Q66*(s^4 + c^4);
QA(2,1) = QA(1,2);
QA(3,2) = QA(2,3);
QA(3,1) = QA(1,3);

end
