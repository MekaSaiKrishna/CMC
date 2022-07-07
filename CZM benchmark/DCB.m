%% LEFM, VCCT, CZM for DCB [i.e. Mode-1]

% Corresponding explanatory figures can be found in the presentation 
% 'ExplanatoryFiguresMATLAB.pptx'

% Summary: This code contains the analytical LEFM, VCCT and CZM solutions

%% LEFM

%Note: Refer Figure-1 in slides

% Geometry:
a0= 30;                   %initial crack length
a = linspace(30,75,100);  %crack progression
b = 1;                    %in-plane width
h = 5;                    %thickness
I = b*(h^3)/12;           %Area moment of Inertia of beam

%note: All lengths are in mm

% Material:
E   = 70e3;                %Young's Modulus (MPa)
G1c = 7.48;                %Mode-1 Fracture Toughness (N/mm) 

% Force
F = (1./a)*sqrt(E*G1c*b*I);%Force

d = F.*(a.^3)/(3*E*I);     %"DELTA"=0.5*(Total Deflection)

% Plotting LEFM solution:
F1=[0,F]; %Force
a1=[0,a]; %Crack Length
d1=[0,d]; %0.5*(Total Deflection)

figure
plot(a1,F1,'-o','MarkerIndices',1:10:length(F1))
title("Force vs Crack Length")
ylabel("Force (F) [in N]")
xlabel("Crack Length (a) [in mm]")
legend({'LEFM'},'Location','northeast')
hold off
savefig('LEFM_FvsCrackLength.fig')

figure
plot(d1,F1,'-o','MarkerIndices',1:10:length(F1))
title("Force vs Deflection")
ylabel("Force (F) [in N]")
xlabel("Deflection (\Delta) [in mm]")
legend({'LEFM'},'Location','northeast')
hold off
savefig('LEFM_FvsDeflection.fig')

figure
tiledlayout(2,1)
%Top Plot
ax1 = nexttile;
plot(ax1,a1,F1)
title(ax1,'Force vs Crack Length')
ylabel(ax1,'Force (F) [in N]')
xlabel(ax1,'Crack Length (a) [in mm]')
legend(ax1,{'LEFM'},'Location','northeast')

%Bottom Plot
ax2 = nexttile;
plot(ax2,d1,F1)
title(ax2,'Force vs Deflection')
ylabel(ax2,'Force (F) [in N]')
xlabel(ax2,'Deflection (\Delta) [in mm]')
legend(ax2,{'LEFM'},'Location','northeast')

savefig('LEFM_DCB_AllPlots.fig')

%% CZM Analytical

