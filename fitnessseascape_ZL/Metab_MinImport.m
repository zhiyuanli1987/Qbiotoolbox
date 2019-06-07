function metabfuns=Metab_MinImport(species_para)

metabfuns=[];

% load parameters
metabfuns.r=species_para.r;
Ks=species_para.Ks;
Gamma=species_para.Gamma;
Es=species_para.Es;
Edim=length(Es);
Es=reshape(Es,Edim,1);

% environment varibles: cs
% inner varibles: xs

% Define: Intake rates as a vector
Intakes=@(cs,xs)Es.*cs./(cs+Ks);
% Define: Growth  rate:
Growth=@(cs,xs)Gamma*min(Intakes(cs,xs));
% Define: Within cell varibles:
dxdt=[];

metabfuns.Intakes=@(cs,xs)Intakes(cs,xs);
metabfuns.Growth=Growth;
metabfuns.dxdt=dxdt;
metabfuns.intercell_dim=species_para.intercell_dim;
metabfuns.Edim=Edim;

%%%% If there are analytical solutions 
% Flux balance curve and Growth contour
%%%% for dimension 2
FB_cb_dim2=@(c_a,c_supplys) (c_a^2*Es(2) + (c_supplys(1)^2*Ks(1)^2*Es(2)^2 + 2*c_supplys(1)^2*Ks(1)*c_a*Es(2)^2 + c_supplys(1)^2*c_a^2*Es(2)^2 - 2*c_supplys(1)*c_supplys(2)*Ks(1)*c_a*Es(1)*Es(2) - 2*c_supplys(1)*c_supplys(2)*c_a^2*Es(1)*Es(2) - 2*c_supplys(1)*Ks(1)^2*c_a*Es(2)^2 + 2*c_supplys(1)*Ks(1)*Ks(2)*c_a*Es(1)*Es(2) - 4*c_supplys(1)*Ks(1)*c_a^2*Es(2)^2 + 2*c_supplys(1)*Ks(2)*c_a^2*Es(1)*Es(2) - 2*c_supplys(1)*c_a^3*Es(2)^2 + c_supplys(2)^2*c_a^2*Es(1)^2 + 2*c_supplys(2)*Ks(1)*c_a^2*Es(1)*Es(2) + 2*c_supplys(2)*Ks(2)*c_a^2*Es(1)^2 + 2*c_supplys(2)*c_a^3*Es(1)*Es(2) + Ks(1)^2*c_a^2*Es(2)^2 - 2*Ks(1)*Ks(2)*c_a^2*Es(1)*Es(2) + 2*Ks(1)*c_a^3*Es(2)^2 + Ks(2)^2*c_a^2*Es(1)^2 - 2*Ks(2)*c_a^3*Es(1)*Es(2) + c_a^4*Es(2)^2)^(1/2) - c_supplys(1)*Ks(1)*Es(2) - c_supplys(1)*c_a*Es(2) + c_supplys(2)*c_a*Es(1) + Ks(1)*c_a*Es(2) - Ks(2)*c_a*Es(1))/(2*c_a*Es(1));
GC_cb_dim2=@(c_a,g)Ks(2)/((1-Es(1) )*Gamma/g-1)*( c_a>Ks(1)/(Gamma/g*Es(1)-1) ) + 100* ~(( c_a>Ks(1)/(Gamma/g*Es(1)-1) ));

%%%% for dimension 3
FB_cbcc_dim3 = @(c_a,c_supplys)[(c_a^2*Es(2) + (c_supplys(1)^2*Ks(1)^2*Es(2)^2 + 2*c_supplys(1)^2*Ks(1)*c_a*Es(2)^2 + c_supplys(1)^2*c_a^2*Es(2)^2 - 2*c_supplys(1)*c_supplys(2)*Ks(1)*c_a*Es(1)*Es(2) - 2*c_supplys(1)*c_supplys(2)*c_a^2*Es(1)*Es(2) - 2*c_supplys(1)*Ks(1)^2*c_a*Es(2)^2 + 2*c_supplys(1)*Ks(1)*Ks(2)*c_a*Es(1)*Es(2) - 4*c_supplys(1)*Ks(1)*c_a^2*Es(2)^2 + 2*c_supplys(1)*Ks(2)*c_a^2*Es(1)*Es(2) - 2*c_supplys(1)*c_a^3*Es(2)^2 + c_supplys(2)^2*c_a^2*Es(1)^2 + 2*c_supplys(2)*Ks(1)*c_a^2*Es(1)*Es(2) + 2*c_supplys(2)*Ks(2)*c_a^2*Es(1)^2 + 2*c_supplys(2)*c_a^3*Es(1)*Es(2) + Ks(1)^2*c_a^2*Es(2)^2 - 2*Ks(1)*Ks(2)*c_a^2*Es(1)*Es(2) + 2*Ks(1)*c_a^3*Es(2)^2 + Ks(2)^2*c_a^2*Es(1)^2 - 2*Ks(2)*c_a^3*Es(1)*Es(2) + c_a^4*Es(2)^2)^(1/2) - c_supplys(1)*Ks(1)*Es(2) - c_supplys(1)*c_a*Es(2) + c_supplys(2)*c_a*Es(1) + Ks(1)*c_a*Es(2) - Ks(2)*c_a*Es(1))/(2*c_a*Es(1));
    (c_a^2*Es(3) + (c_supplys(1)^2*Ks(1)^2*Es(3)^2 + 2*c_supplys(1)^2*Ks(1)*c_a*Es(3)^2 + c_supplys(1)^2*c_a^2*Es(3)^2 - 2*c_supplys(1)*c_supplys(3)*Ks(1)*c_a*Es(1)*Es(3) - 2*c_supplys(1)*c_supplys(3)*c_a^2*Es(1)*Es(3) - 2*c_supplys(1)*Ks(1)^2*c_a*Es(3)^2 + 2*c_supplys(1)*Ks(1)*Ks(3)*c_a*Es(1)*Es(3) - 4*c_supplys(1)*Ks(1)*c_a^2*Es(3)^2 + 2*c_supplys(1)*Ks(3)*c_a^2*Es(1)*Es(3) - 2*c_supplys(1)*c_a^3*Es(3)^2 + c_supplys(3)^2*c_a^2*Es(1)^2 + 2*c_supplys(3)*Ks(1)*c_a^2*Es(1)*Es(3) + 2*c_supplys(3)*Ks(3)*c_a^2*Es(1)^2 + 2*c_supplys(3)*c_a^3*Es(1)*Es(3) + Ks(1)^2*c_a^2*Es(3)^2 - 2*Ks(1)*Ks(3)*c_a^2*Es(1)*Es(3) + 2*Ks(1)*c_a^3*Es(3)^2 + Ks(3)^2*c_a^2*Es(1)^2 - 2*Ks(3)*c_a^3*Es(1)*Es(3) + c_a^4*Es(3)^2)^(1/2) - c_supplys(1)*Ks(1)*Es(3) - c_supplys(1)*c_a*Es(3) + c_supplys(3)*c_a*Es(1) + Ks(1)*c_a*Es(3) - Ks(3)*c_a*Es(1))/(2*c_a*Es(1))];

GC_cc_dim3 = @(c_a,c_b,g) Ks(3)/(Es(3)*Gamma/g-1)* ( c_a>Ks(1)/(Gamma/g*Es(1)-1)  & c_b>Ks(2)/(Gamma/g*Es(2)-1)) + 100 * ~( c_a>Ks(1)/(Gamma/g*Es(1)-1)  & c_b>Ks(2)/(Gamma/g*Es(2)-1));
        
metabfuns.analytical_FBGCincell_dim2={FB_cb_dim2;
    GC_cb_dim2;
    []
    };

metabfuns.analytical_FBGCincell_dim3={FB_cbcc_dim3;
    GC_cc_dim3;
    []
    };

%% optimization
% there is analytical solution: growth rate & the corresponding strategy
optEs_notnormalized=@(cs)1./(cs./(cs+Ks));
optg=@(env)Gamma/sum(optEs_notnormalized(env));
optEs=@(env)optEs_notnormalized(env)/sum(optEs_notnormalized(env));
analytical_opt_givenenv=@(env)[optg(env);
    optEs(env)];
metabfuns.analytical_opt_givenenv=analytical_opt_givenenv;

% optimal growth contours
opt_GC=[];
opt_GC.GC=@(c_a,D)Ks(2)/(Gamma/D-1-(c_a+Ks(1))/c_a);
opt_GC.unkowndim=2; % which value it presents;
opt_GC.Es=@(c_a,D)[(c_a+Ks(1))/c_a/(Gamma/D);
1-(c_a+Ks(1))/c_a/(Gamma/D)];
metabfuns.opt_GC=opt_GC;
