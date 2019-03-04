function metabfuns=Metab_RwImport(species_para)

metabfuns=[];

% load parameters
Ks=species_para.Ks;
env_dim=length(Ks);
Gamma=species_para.Gamma;
Es=reshape(species_para.Es,env_dim,1);

% environment varibles: cs
% inner varibles: xs

% Define: Intake rates as a vector
Intakes=@(cs,xs)Es.*cs./(cs+Ks);
% Define: Growth  rate:
Growth=@(cs,xs)Gamma/(sum(1./Intakes(cs,xs)));
% Define: Within cell varibles:
dxdt=[];

%%%% If there are analytical solutions 
% Flux balance curve and Growth contour
%%%%%%%%%%%%% dim 2
FB_cb_dim2=@(c_a,c_supplys) (c_a^2*Es(2) + (c_supplys(1)^2*Ks(1)^2*Es(2)^2 + 2*c_supplys(1)^2*Ks(1)*c_a*Es(2)^2 + c_supplys(1)^2*c_a^2*Es(2)^2 - 2*c_supplys(1)*c_supplys(2)*Ks(1)*c_a*Es(1)*Es(2) - 2*c_supplys(1)*c_supplys(2)*c_a^2*Es(1)*Es(2) - 2*c_supplys(1)*Ks(1)^2*c_a*Es(2)^2 + 2*c_supplys(1)*Ks(1)*Ks(2)*c_a*Es(1)*Es(2) - 4*c_supplys(1)*Ks(1)*c_a^2*Es(2)^2 + 2*c_supplys(1)*Ks(2)*c_a^2*Es(1)*Es(2) - 2*c_supplys(1)*c_a^3*Es(2)^2 + c_supplys(2)^2*c_a^2*Es(1)^2 + 2*c_supplys(2)*Ks(1)*c_a^2*Es(1)*Es(2) + 2*c_supplys(2)*Ks(2)*c_a^2*Es(1)^2 + 2*c_supplys(2)*c_a^3*Es(1)*Es(2) + Ks(1)^2*c_a^2*Es(2)^2 - 2*Ks(1)*Ks(2)*c_a^2*Es(1)*Es(2) + 2*Ks(1)*c_a^3*Es(2)^2 + Ks(2)^2*c_a^2*Es(1)^2 - 2*Ks(2)*c_a^3*Es(1)*Es(2) + c_a^4*Es(2)^2)^(1/2) - c_supplys(1)*Ks(1)*Es(2) - c_supplys(1)*c_a*Es(2) + c_supplys(2)*c_a*Es(1) + Ks(1)*c_a*Es(2) - Ks(2)*c_a*Es(1))/(2*c_a*Es(1));
GC_cb_dim2=@(c_a,g)            (Ks(2)/((1-Es(1) ) ))/(Gamma/g-1/((1-Es(1) ) )-1/(Es(1)*c_a/(c_a+Ks(1) )));

%%%% for dimension 3
FB_cbcc_dim3 = @(c_a,c_supplys)[(c_a^2*Es(2) + (c_supplys(1)^2*Ks(1)^2*Es(2)^2 + 2*c_supplys(1)^2*Ks(1)*c_a*Es(2)^2 + c_supplys(1)^2*c_a^2*Es(2)^2 - 2*c_supplys(1)*c_supplys(2)*Ks(1)*c_a*Es(1)*Es(2) - 2*c_supplys(1)*c_supplys(2)*c_a^2*Es(1)*Es(2) - 2*c_supplys(1)*Ks(1)^2*c_a*Es(2)^2 + 2*c_supplys(1)*Ks(1)*Ks(2)*c_a*Es(1)*Es(2) - 4*c_supplys(1)*Ks(1)*c_a^2*Es(2)^2 + 2*c_supplys(1)*Ks(2)*c_a^2*Es(1)*Es(2) - 2*c_supplys(1)*c_a^3*Es(2)^2 + c_supplys(2)^2*c_a^2*Es(1)^2 + 2*c_supplys(2)*Ks(1)*c_a^2*Es(1)*Es(2) + 2*c_supplys(2)*Ks(2)*c_a^2*Es(1)^2 + 2*c_supplys(2)*c_a^3*Es(1)*Es(2) + Ks(1)^2*c_a^2*Es(2)^2 - 2*Ks(1)*Ks(2)*c_a^2*Es(1)*Es(2) + 2*Ks(1)*c_a^3*Es(2)^2 + Ks(2)^2*c_a^2*Es(1)^2 - 2*Ks(2)*c_a^3*Es(1)*Es(2) + c_a^4*Es(2)^2)^(1/2) - c_supplys(1)*Ks(1)*Es(2) - c_supplys(1)*c_a*Es(2) + c_supplys(2)*c_a*Es(1) + Ks(1)*c_a*Es(2) - Ks(2)*c_a*Es(1))/(2*c_a*Es(1));
    (c_a^2*Es(3) + (c_supplys(1)^2*Ks(1)^2*Es(3)^2 + 2*c_supplys(1)^2*Ks(1)*c_a*Es(3)^2 + c_supplys(1)^2*c_a^2*Es(3)^2 - 2*c_supplys(1)*c_supplys(3)*Ks(1)*c_a*Es(1)*Es(3) - 2*c_supplys(1)*c_supplys(3)*c_a^2*Es(1)*Es(3) - 2*c_supplys(1)*Ks(1)^2*c_a*Es(3)^2 + 2*c_supplys(1)*Ks(1)*Ks(3)*c_a*Es(1)*Es(3) - 4*c_supplys(1)*Ks(1)*c_a^2*Es(3)^2 + 2*c_supplys(1)*Ks(3)*c_a^2*Es(1)*Es(3) - 2*c_supplys(1)*c_a^3*Es(3)^2 + c_supplys(3)^2*c_a^2*Es(1)^2 + 2*c_supplys(3)*Ks(1)*c_a^2*Es(1)*Es(3) + 2*c_supplys(3)*Ks(3)*c_a^2*Es(1)^2 + 2*c_supplys(3)*c_a^3*Es(1)*Es(3) + Ks(1)^2*c_a^2*Es(3)^2 - 2*Ks(1)*Ks(3)*c_a^2*Es(1)*Es(3) + 2*Ks(1)*c_a^3*Es(3)^2 + Ks(3)^2*c_a^2*Es(1)^2 - 2*Ks(3)*c_a^3*Es(1)*Es(3) + c_a^4*Es(3)^2)^(1/2) - c_supplys(1)*Ks(1)*Es(3) - c_supplys(1)*c_a*Es(3) + c_supplys(3)*c_a*Es(1) + Ks(1)*c_a*Es(3) - Ks(3)*c_a*Es(1))/(2*c_a*Es(1))];

GC_cc_dim3 = @(c_a,c_b,g) Ks(3)/(Es(3)*(Gamma/g-1/(Es(1)*c_a/(c_a+Ks(1)))-1/(Es(2)*c_b/(c_b+Ks(2))))-1);


%%%%%%%%%%%%%%%%%%%%%
metabfuns.Intakes=Intakes;
metabfuns.Growth=Growth;
metabfuns.dxdt=dxdt;
metabfuns.intercell_dim=species_para.intercell_dim;

metabfuns.analytical_FBGCincell_dim2={FB_cb_dim2;
    GC_cb_dim2;
    []
    };

metabfuns.analytical_FBGCincell_dim3={FB_cbcc_dim3;
    GC_cc_dim3;
    []
    };
