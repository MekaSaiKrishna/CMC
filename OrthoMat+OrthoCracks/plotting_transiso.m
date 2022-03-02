% Orthotropic SCA formulation:
%-----------------------------

% sigma   = [σ11, σ22, σ33, τ12, τ13, τ23]'

% epsilon = [ε11, ε22, ε33, γ12, γ13, γ23]'


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


% [sigma'] = [R_sigma][sigma]
R_sigma = [K1, 2*K2; K3, K4];


% [eps'] = [R_eps][eps]
R_eps = (inv(R_sigma))';


%% Formulating the N matrix
%------------------------------------------------

% [N]; where [?]=[N][?']_reduced & [?']_reduced=[N]'[?]

% [?']_reduced = [?'11, ?'12, ?'13]

% [?']_reduced = [?'11, ?'12, ?'13]
RT = R_sigma';
N = [RT(:,1), RT(:,4), RT(:,5)]; %because [?'22, ?'33, ?'23] = 0


% Orthotropy Introduction:
%-------------------------


% Local co-ordinates == Global co-ordinates; Crack Normal = [1,0,0]
%------------------------------------------------------------------
N_new = subs(N,{a12,a13,a21,a31,a32,a23,a11,a22,a33},{0,0,0,0,0,0,1,1,1});



%% Orthotropic Constants:
%------------------------
syms E1 E2 nu12 nu23 G23 G12 real

G23  = E2*0.5/(1+nu23);

%Orthotropic complinace matrix
S_ortho = [1/E1,(-nu12/E1),(-nu12/E1), 0,    0,  0;
          (-nu12/E1), 1/E2,(-nu23/E2), 0,    0,  0;
          (-nu12/E1), (-nu23/E2), 1/E2,0,    0,  0;
                   0,          0,    0,1/G12, 0, 0;
                   0,          0,    0,0, 1/G12, 0;
                   0,          0,    0,0, 0, 1/G23];

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

% sigf  = y-intercept of the traction-separation law corresponding to GIC
% tau1f = y-intercept of the traction-separation law corresponding to GIIC
% tau2f = y-intercept of the traction-separation law corresponding to GIIC

%     e = normal crack strain
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
 
%% ecr and eps relation
%-----------------------
% ecr - crack strains [local]
% eps - global strains

eps = sym('eps',[1 6],'real'); 

ecr = ((N_new'*Dco*N_new + Dcr)\(N_new'*Dco))*eps';


eq1 = e-simplify(ecr(1,1))==0;
eq2 = g1-simplify(ecr(2,1))==0;
eq3 = g2-simplify(ecr(3,1))==0;

sol_e = solve(eq1,e);
sol_g1 = solve(eq2,g1);
sol_g2 = solve(eq3,g2);

% Update Dcr with the solutions
Dcr = subs(Dcr,{e, g1, g2},{sol_e(2), sol_g1(2), sol_g2(2)});

%% sigma(global) and strain(global) relation

%NOTE: eps(2) = eps(3) = -nu12*eps(1); because sigma1=!0 and all other
%stress components of global stress are zero

%ecr_updated = subs([sol_e(2), sol_g1(2), sol_g2(2)]',{eps(2),eps(3)},{-nu12*eps(1),-nu12*eps(1)})
ecr_updated = [sol_e(2),sol_g1(2),sol_g2(2)]';
%sigma = subs(Dco*(eps-N_new*ecr_updated),{eps(2),eps(3)},{-nu12*eps(1),-nu12*eps(1)});
sigma_updated = Dco*(eps-N_new*ecr_updated);
%sigma11_exp = sigma(1,1)
%sigma22_exp = sigma(2,2)


%% ----------------------------------------------------------------------
% Uniaxial Tensile Loading Case

% Material Properties

% Material Constants: [E1 E2 G12 nu12 nu23]
% Fracture Toughness: [GIC GIIC]
% Cr. Crack Stresses: [sigf tau1f tau2f]

%Material ConstantS
E1_val = 2.5e9;   %in Pa
E2_val = 2e9;     %in Pa
G12_val = 1.5e9;  %in Pa
nu12_val = 0.30;  %in Pa
nu23_val = 0.35;  %in Pa

%Fracture Toughness Values
GIC_val  = 1.5e3; %in N/m
GIIC_val = 1e3;   %in N/m

%Characteristic Length (h) ; gIC=GIC/h & gIIC=GIIC/h
%h = 0.01e-3; %h=0.01mm
h=1e-3;
gIC_val  = GIC_val/h;
gIIC_val = GIIC_val/h;

%crack stress intercepts from Traction-Separation Law
sigf_val  = 60e6; %in Pa
tau1f_val = 50e6; %in Pa
tau2f_val = 40e6; %in Pa


%% Plots of Stresses vs Strain before peak failure

array1 = (0:1000:sigf_val);
eps1_plot = zeros(1,length(array1));
sig11_plot = zeros(1,length(array1));
for i = 1:length(array1)
    sig11_plot(i) = array1(i);
    if eps1_plot < (sigf_val/E1_val)
        eps1_plot(i) = sig11_plot(i)/E1_val;
    end
end

figure('Name','sigma11 vs eps1')
plot(eps1_plot,sig11_plot)
grid on;
xlabel("\epsilon_1")
ylabel("\sigma_{11}")
title("Global Stress(\sigma_{11}) vs Global Strain(\epsilon_{1}) Pre-Peak")

eps2_plot = zeros(1,length(array1));
for i = 1:length(array1)
    sig11_plot(i) = array1(i);
    if abs(eps2_plot) < nu12_val*(sigf_val/E1_val)
        eps2_plot(i) = -nu12_val*sig11_plot(i)/E1_val;
    end
end
sig22_plot=zeros(1,length(eps2_plot));
figure('Name','sigma22 vs eps2')
plot(eps2_plot,sig22_plot)
grid on;
xlabel("\epsilon_2")
ylabel("\sigma_{22}")
title("Global Stress(\sigma_{22}) vs Global Strain(\epsilon_{2}) Pre-Peak")

eps3_plot = eps2_plot;
sig33_plot=zeros(1,length(eps3_plot));
figure('Name','sigma33 vs eps3')
plot(eps3_plot,sig33_plot)
grid on;
xlabel("\epsilon_3")
ylabel("\sigma_{33}")
title("Global Stress(\sigma_{33}) vs Global Strain(\epsilon_{3}) Pre-Peak")



%% Plotting of Crack Strains vs Global Strain Components

array1 = (0:0.0001:0.05);
% g1cr vs eps4
% Note: When eps4 < g1cr, g1cr = 0
% When eps4=tau1f/G12, then ecr_nn=sol_g1(2)

g1cr_exp = subs(sol_g1(2),{GIIC,G12,tau1f},{gIIC_val,G12_val,tau1f_val});
g1cr_plot = zeros(1,length(array1));

i=0;
for i = 1:length(array1)
    if array1(i) < (tau1f_val/G12_val)
        g1cr_plot(i) = 0; 
    else
        g1cr_plot(i) = subs(g1cr_exp,{eps(4)},{array1(i)});
    end
end

figure('Name','g1cr vs eps4')
plot(array1,g1cr_plot)
grid on;
xlabel("\epsilon_4")
ylabel("\gamma1_{cr}")
title("Crack Strain(\gamma_{1}^{cr}) vs Global Strain(\epsilon_{4})")

% g2cr vs eps5
% Note: When eps5 < g2cr, g2cr = 0
% When eps5=tau2f/G12, then g2cr=sol_g2(2)

g2cr_exp = subs(sol_g2(2),{GIIC,G12,tau2f},{gIIC_val,G12_val,tau2f_val});
g2cr_plot = zeros(1,length(array1));

i=0;
for i = 1:length(array1)
    if array1(i) < (tau2f_val/G12_val)
        g2cr_plot(i) = 0; 
    else
        g2cr_plot(i) = subs(g2cr_exp,{eps(5)},{array1(i)});
    end
end

figure('Name','g2cr vs eps5')
plot(array1,g2cr_plot)
grid on;
xlabel("\epsilon_5")
ylabel("\gamma2_{cr}")
title("Crack Strain(\gamma_{2}^{cr}) vs Global Strain(\epsilon_{5})")

% ecr_nn vs eps1
% Note: When eps1 < ecr_nn, ecr_nn = 0
% When eps1=sigf/E1, then ecr_nn=sol_e(2)

ecr_exp_temp = subs(sol_e(2),{eps(2),eps(3)},{-nu12*eps(1),-nu12*eps(1)});

ecr_exp = subs(ecr_exp_temp,{GIC,E1,E2,nu12,nu23,sigf},{gIC_val,...
    E1_val,E2_val,nu12_val,nu23_val,sigf_val});
ecr_plot = zeros(1,length(array1));

i=0;
for i = 1:length(array1)
    if array1(i) < (sigf_val/E1_val)
        ecr_plot(i) = 0; 
    else
        ecr_plot(i) = subs(ecr_exp,{eps(1)},{array1(i)});
    end
end

figure('Name','ecr vs eps1')
plot(array1,ecr_plot)
grid on;
xlabel("\epsilon_1")
ylabel("\epsilon1_{cr}")
title("Crack Strain(\epsilon^{cr}) vs Global Strain(\epsilon_{1})")

%% Plotting of Global Stresses vs Global Strains

strain11_global_1=zeros(1,60);
strain11_global_2=zeros(1,60);
sigma11_global=zeros(1,60);

for i=0:59
% Define Global Stress
%-----------------------
    % Uniaxial Tensile Stress in 1-drxn i.e. ?11 !=0
      sigma_global = 1e6*[(60-i) 0 0 0 0 0]'; %SI Units

% Find Global Strain
% C_global = (Dco - Dco*N_new*((N_new'*Dco*N_new + Dcr)\(N_new'*Dco)));
strain_global = (Dco - Dco*N_new*((N_new'*Dco*N_new + Dcr)\(N_new'*Dco)))\sigma_global;

% Solving for Global Strain Components
  EQ1 = eps(1) - strain_global(1)==0;
  EQ2 = eps(2) - strain_global(2)==0;
  EQ3 = eps(3) - strain_global(3)==0;
 %EQ4 = eps(4) - strain_global(4)==0;
 %EQ5 = eps(5) - strain_global(5)==0;
 %EQ6 = eps(6) - strain_global(6)==0;
 
%--------------------------------------------------------------------------
%NOTE: EQ4,EQ5,EQ6 are not necessary because we get eps(4)=eps(5)=eps(6)=0
%--------------------------------------------------------------------------

 eps_sol = solve([EQ1,EQ2,EQ3],[eps(1),eps(2),eps(3)]);
 
 %Extracting solutions for Global Strain
 global_sols = [eps_sol.eps1, eps_sol.eps2, eps_sol.eps3];
 
 %Solution Set-1
 global_sols_1 = global_sols(1,:);
 
 %Solution Set-2
 global_sols_2 = global_sols(2,:); 
 
 % Update Strain Global
   strain_global_1 = subs(([global_sols_1'; strain_global(4:6,1)]),...
       {E1,E2,GIC,sigf,nu12,nu23},{E1_val,E2_val,gIC_val,sigf_val,nu12_val,nu23_val});
   strain_global_2 = subs(([global_sols_2'; strain_global(4:6,1)]),...
       {E1,GIC,sigf,nu12},{E1_val,gIC_val,sigf_val,nu12_val});
   
   strain11_global_1(i+1)=strain_global_1(1);
   strain11_global_2(i+1)=strain_global_2(1);
   sigma11_global(i+1)=sigma_global(1);
   
end

%figure('Name','sigma11 vs eps1')
%plot(strain11_global_1,sigma11_global)
%grid on;
%xlabel("\epsilon_1 (Solution-1)")
%ylabel("\sigma_{11}")
%title("Global Stress(\sigma_{11}) vs Global Strain(\epsilon_{1}) Post-Peak")
%%
figure('Name','sigma11 vs eps1')
plot(eps1_plot,sig11_plot,'r')
hold on;
plot(strain11_global_2,sigma11_global,'b')
hold off;
grid on;
xlabel("\epsilon_1 (Solution-2)")
ylabel("\sigma_{11}")
legend('pre-peak','post-peak')
title("Global Stress(\sigma_{11}) vs Global Strain(\epsilon_{1})")

%% Importing ABAQUS Data
filename = 'S11vsE11fromABAQUS.xlsx';
A = readtable(filename);

E11_data = table2array(A(:,{'E11'}))';
S11_data = table2array(A(:,{'S11'}))';

%% Abaqus Plot of S11 vs E11

figure('Name','sigma11 vs eps1 (ABAQUS)')
plot(eps1_plot,sig11_plot*1e-6,'r')
hold on;
plot(strain11_global_2,sigma11_global*1e-6,'b')
hold on;
plot(E11_data,S11_data,'*')
grid on;
xlabel("\epsilon_1 (Solution-2)")
ylabel("\sigma_{11} (MPa)")
legend('pre-peak','post-peak','ABAQUS')
ax = gca
ax.FontSize = 18;
title("Global Stress(\sigma_{11}) vs Global Strain(\epsilon_{1})")
