% Finds the A,B,D Matrices

function [A,B,D] =ABDmat(Q,z)
N=numel(z)-1; %number of layers
A=zeros(3,3);
B=zeros(3,3);
D=zeros(3,3);
for i=1:N
    A_indu=Q(:,:,i)*(z(i+1)-z(i));
    A=A+A_indu;    
    B_indu=Q(:,:,i)*0.5*((z(i+1))^2-(z(i))^2);
    B=B+B_indu;
    D_indu=Q(:,:,i)*(1/3)*((z(i+1))^3-(z(i))^3);
    D=D+D_indu;
end
