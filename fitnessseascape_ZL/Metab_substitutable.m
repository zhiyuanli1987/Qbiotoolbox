function metabfuns=Metab_substitutable(species_para)

metabfuns=[];

% load parameters
Ks=species_para.Ks;
Vs=species_para.Vs;
Es=species_para.Es;
Es=reshape(Es,length(Es),1);


% environment varibles: cs
% inner varibles: xs

% Define: Intake rates as a vector
Intakes=@(cs,xs)Es.*cs./(cs+Ks);
% Define: Growth  rate:
Growth=@(cs,xs)sum(Vs.*Intakes(cs,xs));
% Define: Within cell varibles:
dxdt=[];

%%%% If there are analytical solutions 
% Flux balance curve and Growth contour
%%%% for dimension 2
FB_cb_dim2=@(c_a,c_supplys) (c_a^2*Es(2) + (c_supplys(1)^2*Ks(1)^2*Es(2)^2 + 2*c_supplys(1)^2*Ks(1)*c_a*Es(2)^2 + c_supplys(1)^2*c_a^2*Es(2)^2 - 2*c_supplys(1)*c_supplys(2)*Ks(1)*c_a*Es(1)*Es(2) - 2*c_supplys(1)*c_supplys(2)*c_a^2*Es(1)*Es(2) - 2*c_supplys(1)*Ks(1)^2*c_a*Es(2)^2 + 2*c_supplys(1)*Ks(1)*Ks(2)*c_a*Es(1)*Es(2) - 4*c_supplys(1)*Ks(1)*c_a^2*Es(2)^2 + 2*c_supplys(1)*Ks(2)*c_a^2*Es(1)*Es(2) - 2*c_supplys(1)*c_a^3*Es(2)^2 + c_supplys(2)^2*c_a^2*Es(1)^2 + 2*c_supplys(2)*Ks(1)*c_a^2*Es(1)*Es(2) + 2*c_supplys(2)*Ks(2)*c_a^2*Es(1)^2 + 2*c_supplys(2)*c_a^3*Es(1)*Es(2) + Ks(1)^2*c_a^2*Es(2)^2 - 2*Ks(1)*Ks(2)*c_a^2*Es(1)*Es(2) + 2*Ks(1)*c_a^3*Es(2)^2 + Ks(2)^2*c_a^2*Es(1)^2 - 2*Ks(2)*c_a^3*Es(1)*Es(2) + c_a^4*Es(2)^2)^(1/2) - c_supplys(1)*Ks(1)*Es(2) - c_supplys(1)*c_a*Es(2) + c_supplys(2)*c_a*Es(1) + Ks(1)*c_a*Es(2) - Ks(2)*c_a*Es(1))/(2*c_a*Es(1));
GC_cb_dim2=@(c_a,g)-Ks(2)/(1-(Vs(2)*Es(2))/(g-Vs(1)*Es(1)*c_a/(c_a+Ks(1))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metabfuns.Intakes=Intakes;
metabfuns.Growth=Growth;
metabfuns.dxdt=dxdt;
metabfuns.intercell_dim=species_para.intercell_dim;

% this might be empty
metabfuns.analytical_FBGCincell_dim2={FB_cb_dim2;
    GC_cb_dim2;
    []
    };


