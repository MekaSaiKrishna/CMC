%Problem Description:

% Material Constants and Geometric Constants
E1=158490e6;  E2=9017e6;  G12=5703e6; v12=0.316;
h=0.125e-3;           %thickness of each ply
L = 200e-3;           %side length of the lamina plate, here it is assumed to be a square panel

%Obtaining Qbar
Q1bar=Qbar(45); 
Q2bar=Qbar(-45);
z=[-2*h,-h,0,h,2*h]; %there are 4 plies, it is a symmetrical laminate

%Laminate Sequence: [+45/-45/-45/+45]

%Defining Q_G matrix
Q_G=zeros(3,3,4); %here '4' is the number of plies, this has to be changed based on our layup configuration

Q_G(:,:,1)=Q1bar; %  45deg
Q_G(:,:,2)=Q2bar; % -45deg
Q_G(:,:,3)=Q2bar; % -45deg
Q_G(:,:,4)=Q1bar; %  45deg

%Obtaining A matrix
[A,B,D]=ABDmat(Q_G,z);

ABD=[A,B;
     B,D];
ABDinv=inv(ABD);

Ainv=ABDinv(1:3,1:3);
Binv=ABDinv(4:6,1:3);

%Midplane Strains
%syms Ny
sigma0 = 10000e6;     %Applied stress, go to the "readme.MD" file for more information
area = (4*h)*L
Ny = sigma0*area/L;
N=[0;Ny;0];

eps=vpa(Ainv*N)
curv=vpa(Binv*N)

%(b) The [Global] stresses in the 45° and -45° plies 
%syms z;
z=0; %this is because 'curv' is very small we want to ignore its contribution
%NOTE: Before assuming z=0, check is the 'curv' is actually very small or not!
  
sigma1=Q1bar*(eps+z*curv); %in  45deg -h<z<0
sigma2=Q2bar*(eps+z*curv); %in -45deg  0<z<h

%% %(b) The [Local] stresses in the 45° and -45° plies
sigma1_local = rot_stress(sigma1,45);
sigma2_local = rot_stress(sigma2,45);

%% Results, printing results

fprintf("<strong>RESULTS</strong>:")
fprintf("\n")
fprintf("--------")
fprintf("\n")
fprintf("Stresses in the +45 Plies (Ply-1 and Ply-4):")
fprintf("\n")
fprintf("--------------------------------------------")
fprintf("\n")
fprintf("Global:")
fprintf("\n")
fprintf('sigmaX = %0.2f MPa',sigma1(1)*1e-6)
fprintf("\n")
fprintf('sigmaY = %0.2f MPa',sigma1(2)*1e-6)
fprintf("\n")
fprintf('sigmaXY = %0.2f MPa',sigma1(3)*1e-6)
fprintf("\n")
fprintf("********")
fprintf("\n")
fprintf("Local:")
fprintf("\n")
fprintf('sigma1 = %0.2f MPa',sigma1_local(1)*1e-6)
fprintf("\n")
fprintf('sigma2 = %0.2f MPa',sigma1_local(2)*1e-6)
fprintf("\n")
fprintf('sigma12 = %0.2f MPa',sigma1_local(3)*1e-6)
fprintf("\n")


fprintf("\n")
fprintf("Stresses in the -45 Plies (Ply-2 and Ply-3):")
fprintf("\n")
fprintf("--------------------------------------------")
fprintf("\n")
fprintf("Global:")
fprintf("\n")
fprintf('sigmaX = %0.2f MPa',sigma2(1)*1e-6)
fprintf("\n")
fprintf('sigmaY = %0.2f MPa',sigma2(2)*1e-6)
fprintf("\n")
fprintf('sigmaXY = %0.2f MPa',sigma2(3)*1e-6)
fprintf("\n")
fprintf("********")
fprintf("\n")
fprintf("Local:")
fprintf("\n")
fprintf('sigma1 = %0.2f MPa',sigma2_local(1)*1e-6)
fprintf("\n")
fprintf('sigma2 = %0.2f MPa',sigma2_local(2)*1e-6)
fprintf("\n")
fprintf('sigma12 = %0.2f MPa',sigma2_local(3)*1e-6)
fprintf("\n")




