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

%%
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


%% Assumptions:
%-------------
% (1) Dda = [0], i.e. there is no damping to take care of the numerical
% instabilities
    
% Defining Dco and Dcr
%----------------------

%Linear Traction-Separation Law
syms sigf GIC epscr nu real

Ecr  = sigf*((1/epscr)-(sigf/(2*GIC)));
Gcr1 = 0.5*Ecr/(1+nu);
Gcr2 = Gcr1;

%Dcr matrix definition
Dcr = sym(zeros(3,3));

Dcr(1,1) = Ecr;
Dcr(2,2) = Gcr1;
Dcr(3,3) = Gcr2;

%Dco matrix definition
Dco = C_ortho;
 
 
%% ecr and epsilon relation
%---------------------------

eps = sym('eps',[1 6],'real'); 

ecr = ((N_new'*Dco*N_new + Dcr)\(N_new'*Dco))*eps';


epscr1 = ecr(1,1)
gammacr1 = ecr(2,1)
gammacr2 = ecr(3,1)


