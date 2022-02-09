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
syms sigf GIC e real
syms tau1f GIIC g1 real
syms tau2f GIIC g2 real

% sigf = y-intercept of the traction-separation law corresponding to GIC
% tau1f = y-intercept of the traction-separation law corresponding to GIIC
% tau2f = y-intercept of the traction-separation law corresponding to GIIC

%e = normal crack strain
%g1, g2 = sliding crack strains

Ecr  = sigf*((1/e)-(sigf/(2*GIC)));
Gcr1 = tau1f*((1/g1)-(tau1f/(2*GIIC))); 
Gcr2 = tau2f*((1/g2)-(tau2f/(2*GIIC)));


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


eq1 = e-simplify(ecr(1,1))==0;
eq2 = g1-simplify(ecr(2,1))==0;
eq3 = g2-simplify(ecr(3,1))==0;

sol_e = solve(eq1,e)
sol_g1 = solve(eq2,g1)
sol_g2 = solve(eq3,g2)

% Update Dcr with the solutions

Dcr = subs(Dcr,{e, g1, g2},{sol_e(2), sol_g1(2), sol_g2(2)})
fprintf("---------------------------------- \n")

%% ----------------------------------------------------------------------
% Uniaxial Tensile Loading Case

% Material Properties

% [E1 E2 E3 nu12 nu13 nu23 G23 G12 G13] = []
% [GIC GIIC] = []
% [sigf tau1f tau2f] = []

% Define Global Stress
    % Uniaxial Tensile Stress in 1-drxn i.e. ?11 !=0
    sigma_global = [1000 0 0 0 0 0]'; %SI Units

% Find Global Strain
%C_global = (Dco - Dco*N_new*((N_new'*Dco*N_new + Dcr)\(N_new'*Dco)));
strain_global = (Dco - Dco*N_new*((N_new'*Dco*N_new + Dcr)\(N_new'*Dco)))\sigma_global;

% Find Crack Strain Components

ecr_nn = subs(sol_e(2),{eps(1),eps(2),eps(3),eps(4),eps(5),eps(6)},{strain_global(1,1),...
    strain_global(2,1),strain_global(3,1),strain_global(4,1),strain_global(5,1)...
    strain_global(6,1)})

g1cr = subs(sol_g1(2),{eps(1),eps(2),eps(3),eps(4),eps(5),eps(6)},{strain_global(1,1),...
    strain_global(2,1),strain_global(3,1),strain_global(4,1),strain_global(5,1)...
    strain_global(6,1)})

g2cr = subs(sol_g2(2),{eps(1),eps(2),eps(3),eps(4),eps(5),eps(6)},{strain_global(1,1),...
    strain_global(2,1),strain_global(3,1),strain_global(4,1),strain_global(5,1)...
    strain_global(6,1)})
