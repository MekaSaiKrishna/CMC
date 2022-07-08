%% LEFM, VCCT, CZM for ENF [i.e. Mode-2]

% Corresponding explanatory figures can be found in the presentation 
% 'ExplanatoryFiguresMATLAB.pptx'

% Summary: This code contains the analytical LEFM, VCCT and CZM solutions


%% ANALYTICAL SOLUTION

%Reference: (1) https://journals.sagepub.com/doi/pdf/10.1177/002199839803201401
%           (2) https://doi.org/10.1016/j.engfracmech.2006.03.006

E   = 70e3;      % in MPa
b   = 10;        % in mm
h   = 1.5;       % in mm
I   = b*(h^3)/12;
a   = 30;        % in mm; [Initial Crack Length]
G2c = 1.45;      % in N/mm
L   = 50;        % in mm [half length of beam]


% Intersection point of (OB) and (ABC)
Px1 = sqrt(64*G2c*b*E*I)/(a*(3*sqrt(3))^(1/3));

% Intersection point of (ABC) and (DE)
Px2= sqrt(64*G2c*b*E*I)/(L*(3*sqrt(3))^(1/3));

%-----------------------------------------------------
%Loading Line (OB)
P1=linspace(0,round(Px1),100);
d1 = P1.*(2*(L^3)+3*(a^3))/(96*E*I);

figure
p=plot(d1,P1,'LineWidth',1.5);
p.Color=[1 0.5 0];
hold on;

%-----------------------------------------------------
%Unloading Line (ABC) [a<L]
P2 = linspace(round(Px1),round(Px2),100);
d2 = (P2./(96*E*I)).*(2*L^3 + (((64*G2c*b*E*I)^1.5)./(sqrt(3)*P2.^3)));

p=plot(d2,P2,'--*b','MarkerIndices',1:15:length(d2),'LineWidth',2);
p.MarkerFaceColor = [0 0.7 0.7];
p.MarkerEdgeColor = [0 0.7 0.7];
p.MarkerSize = 2;

hold on;

%-----------------------------------------------------
% DE [a>L]
P3 = linspace(Px2,3*Px2,100);
d3 = (P3./(24*E*I)).*(2*L^3 - (((64*G2c*b*E*I)^1.5)./(4*sqrt(3)*P3.^3)));

p=plot(d3,P3,':*','MarkerIndices',1:20:length(d3),'LineWidth',2);
p.Color= [0.8 0.4 0];
p.MarkerFaceColor = [0 0.7 0.7];
p.MarkerEdgeColor = [0 0.7 0.7];
p.MarkerSize = 4;

%-----------------------------------------------------
% OE Completely Split
P4 = linspace(0,3*Px2,100);
d4 = (P4.*L^3)./(12*E*I);
p=plot(d4,P4,'--k','LineWidth',1);

legend({'OB','ABC','DE','OE'},'Location','northwest');
title("Load vs Deflection");
xlabel("Deflection in the middle of the specimen \Delta [mm]");
ylabel("Load P [N]");

savefig("Analytic_ENF.fig");





















