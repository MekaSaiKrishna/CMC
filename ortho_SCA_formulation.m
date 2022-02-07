% Orthotropic SCA formulation:
%-----------------------------

% sigma   = [σ11, σ22, σ33, σ12, σ13, σ23]'
% epsilon = [є11, є22, є33, γ12, γ13, γ23]'

% aij - direction cosines 
syms a11 a12 a13 a21 a22 a23 a31 a32 a33 real

K1 = [a11^2, a12^2, a13^2;
      a21^2, a22^2, a23^2;
      a31^2, a32^2, a33^2];

K2 = [a11*a12, a11*a13, a12*a13;
      a21*a22, a21*a23, a22*a23;
      a31*a32, a31*a33, a32*a33];

K3 = [a11*a21, a12*a22, a13*a23;
      a11*a31, a12*a32, a13*a33;
      a21*a31, a22*a32, a23*a33];

K4 = [a11*a22+a12*a21, a11*a23+a13*a21, a12*a23+a13*a22;
      a11*a32+a12*a31, a11*a33+a13*a31, a12*a33+a13*a32;
      a21*a32+a22*a31, a21*a33+a23*a31, a22*a33+a23*a32];

% [σ'] = [R_sigma][σ]
R_sigma = [K1, 2*K2; K3, K4];

% [є'] = [R_eps][є]
R_eps = (inv(R_sigma))';

% Plane Strain Assumption: [є'22, є'33, γ'23] = 0
%------------------------------------------------
% [N]; where [є]=[N][є']_reduced & [σ']_reduced=[N]'[σ]
% [σ']_reduced = [σ'11, σ'12, σ'13]
% [є']_reduced = [є'11, γ'12, γ'13]
RT = R_sigma';
N = [RT(:,1), RT(:,4), RT(:,5)]; %because [є'22, є'33, γ'23] = 0

% Orthotropy Introduction:
%-------------------------

% Local co-ordinates == Global co-ordinates; Crack Normal = [1,0,0]
%------------------------------------------------------------------
N_new = subs(N,{a12,a13,a21,a31,a32,a23,a11,a22,a33},{0,0,0,0,0,0,1,1,1});

% Stress Deviator:
%------------------
syms s11 s12 s13 s21 s22 s23 s31 s32 s33 real
syms phi 
%phi = φ(J2): for J2-material

%s: stress deviator
s = phi*[s11,s12,s13;
         s21,s22,s23;
         s31,s32,s33];
     
% Pseudo Compliance Matrix for Plasticity
%-----------------------------------------
syms S11 S12 S13 S14 S15 S16 real
syms S22 S23 S24 S25 S26 real
syms S33 S34 S35 S36 real
syms S44 S45 S46 real
syms S55 S56 real
syms S66 real

S_plas = [S11,S12,S13,S14,S15,S16;
          S12,S22,S23,S24,S25,S26;
          S13,S23,S33,S34,S35,S36;
          S14,S24,S34,S44,S45,S46;
          S15,S25,S35,S45,S55,S56;
          S16,S26,S36,S46,S56,S66];
      
%Modified S_plas, of the form of S_ortho
S_plas_mod = subs(S_plas,{S14,S15,S16,S24,S25,S26,S34,S35,S36,S45,S46,S56},...
                 {0,0,0,0,0,0,0,0,0,0,0,0});

%% Orthotropic Constants:
%------------------------
syms E1 E2 E3 nu12 nu13 nu23 G23 G12 G13 real

nu21=(nu12/E1)*E2;
nu32=(nu23/E2)*E3;
nu31=(nu13/E1)*E3;

%Orthotropic complinace matrix
S_ortho = [1/E1,(-nu21/E2),(-nu31/E3), 0,    0,  0;
          (-nu12/E1), 1/E2,(-nu32/E3), 0,    0,  0;
          (-nu13/E1), (-nu23/E2), 1/E3,0,    0,  0;
                   0,          0,    0,0.5/G12, 0, 0;
                   0,          0,    0,0, 0.5/G13, 0;
                   0,          0,    0,0, 0, 0.5/G23];
%Orthotropic stiffness matrix
C_ortho = inv(S_ortho);

    %Case-1: Plane Stress
    S_ortho_reduced_plstress = [S_ortho(1,1), S_ortho(1,2), S_ortho(1,4);
                       S_ortho(2,1), S_ortho(2,2), S_ortho(2,4);
                       S_ortho(4,1), S_ortho(4,2), S_ortho(4,4)];

    C_ortho_redcued_plstress = inv(S_ortho_reduced_plstress);
    
    %Case-2: Plane Strain
    C_ortho_reduced_plstrain = [C_ortho(1,1), C_ortho(1,2), C_ortho(1,4);
                                C_ortho(2,1), C_ortho(2,2), C_ortho(2,4);
                                C_ortho(4,1), C_ortho(4,2), C_ortho(4,4)];
    
    S_ortho_reduced_plstrain = inv(C_ortho_reduced_plstrain);
    
    
%% Assumptions:
%-------------
% (1) Dda = [0], i.e. there is no damping to take care of the numerical
% instabilities
    
% Defining Dco and Dcr
%----------------------
% Assumption: Dco and Dcr are symmetric matrices
syms Dco11 Dco22 Dco33 Dco44 Dco55 Dco66 real
syms Dco12 Dco13 Dco14 Dco15 Dco16 real
syms Dco23 Dco24 Dco25 Dco26 real
syms Dco34 Dco35 Dco36 real
syms Dco45 Dco46 real
syms Dco56 real

syms Dcr11 Dcr12 Dcr13 Dcr22 Dcr23 Dcr33 real

Dco_def = [Dco11 Dco12 Dco13 Dco14 Dco15 Dco16;
       Dco12 Dco22 Dco23 Dco24 Dco25 Dco26;
       Dco13 Dco23 Dco33 Dco34 Dco35 Dco36;
       Dco14 Dco24 Dco34 Dco44 Dco45 Dco46;
       Dco15 Dco25 Dco35 Dco45 Dco55 Dco56;
       Dco16 Dco26 Dco36 Dco46 Dco56 Dco66];
   
Dcr_def = [Dcr11 Dcr12 Dcr13;
       Dcr12 Dcr22 Dcr23;
       Dcr13 Dcr23 Dcr33];
   
%% -------------------------------------------------------------------------
%   Values to Dco and Dcr
%  -------------------------------------------------------------------------   
 %S_plas:  a fictional compliance matrix btwn plastic strain and global stress 
 %S_ortho: generalized compliance matrix for orthotropic materials
 
 %Find S_plas components in terms of secant moduli and elastic moduli 
 
 Sco = S_plas + S_ortho; 
 Dco = inv(Sco);
 
 
% Global Stress-Strain Relationship (1):
%---------------------------------------
% inv(A)*b == A\b
% b*inv(A) == b/A
% C_eff_1 = (Dco - Dco*N_new*inv(N_new'*Dco*N_new + Dcr)*(N_new'*Dco))
temp1 = N_new'/Sco; %N_new'*Dco
temp2 = (temp1*N_new + Dcr)\(temp1);
C_eff_1 = (Dco - (Sco\N_new)*temp2);












