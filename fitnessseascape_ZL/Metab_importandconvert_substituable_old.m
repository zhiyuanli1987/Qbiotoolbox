function metabfuns=Metab_importandconvert_substituable(species_para)

metabfuns=[];

% load parameters
Ks=species_para.Ks;
V_import=species_para.V_import;
env_dim=length(Ks);

k=species_para.k;
Es=species_para.Es;

Es_import=reshape(Es(1:env_dim),env_dim,1);
Es_convert=reshape(Es((env_dim+1):end),env_dim,1);


% Define: Intake rates as a vector

Intakes=@(cs,xs)V_import*Es_import.*cs./(cs+Ks);
% Define: Growth  rate:
Growth=@(cs,xs)k*sum(Es_convert.*xs);
% Define: Within cell varibles
dxdt=@(cs,xs)Intakes(cs,xs)-Growth(cs,xs).*xs-k*Es_convert.*xs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metabfuns.Intakes=Intakes;
metabfuns.Growth=Growth;
metabfuns.withincell_dxdt=dxdt;
metabfuns.intercell_dim=species_para.intercell_dim;


%%%% If there are analytical solutions 
% Flux balance curve and Growth contour
%%%% for dimension 2
FB_cb_dim2=@(a,c_supplys) ((Es_import(1)^2*c_supplys(2)^2*a^2 + 2*Es_import(1)^2*c_supplys(2)*Ks(2)*a^2 + Es_import(1)^2*Ks(2)^2*a^2 - 2*Es_import(1)*Es_import(2)*c_supplys(1)*c_supplys(2)*Ks(1)*a - 2*Es_import(1)*Es_import(2)*c_supplys(1)*c_supplys(2)*a^2 + 2*Es_import(1)*Es_import(2)*c_supplys(1)*Ks(1)*Ks(2)*a + 2*Es_import(1)*Es_import(2)*c_supplys(1)*Ks(2)*a^2 + 2*Es_import(1)*Es_import(2)*c_supplys(2)*Ks(1)*a^2 + 2*Es_import(1)*Es_import(2)*c_supplys(2)*a^3 - 2*Es_import(1)*Es_import(2)*Ks(1)*Ks(2)*a^2 - 2*Es_import(1)*Es_import(2)*Ks(2)*a^3 + Es_import(2)^2*c_supplys(1)^2*Ks(1)^2 + 2*Es_import(2)^2*c_supplys(1)^2*Ks(1)*a + Es_import(2)^2*c_supplys(1)^2*a^2 - 2*Es_import(2)^2*c_supplys(1)*Ks(1)^2*a - 4*Es_import(2)^2*c_supplys(1)*Ks(1)*a^2 - 2*Es_import(2)^2*c_supplys(1)*a^3 + Es_import(2)^2*Ks(1)^2*a^2 + 2*Es_import(2)^2*Ks(1)*a^3 + Es_import(2)^2*a^4)^(1/2) + Es_import(2)*a^2 - Es_import(2)*c_supplys(1)*Ks(1) + Es_import(1)*c_supplys(2)*a - Es_import(2)*c_supplys(1)*a - Es_import(1)*Ks(2)*a + Es_import(2)*Ks(1)*a)/(2*Es_import(1)*a);
GC_cb_dim2=@(c_a,D)(Ks(2))/(1/((D/(k*V_import)-(Es_convert(1)*Es_import(1))/(D+Es_convert(1)*k)*c_a/(c_a+Ks(1)))/((Es_convert(2)*Es_import(2))/(D+Es_convert(2)*k)))-1);
metabfuns.analytical_FBGCincell_dim2={FB_cb_dim2;
    GC_cb_dim2;
    []
    };

FB_cb_dim3=@(a,c_supplys) [((Es_import(1)^2*c_supplys(2)^2*a^2 + 2*Es_import(1)^2*c_supplys(2)*Ks(2)*a^2 + Es_import(1)^2*Ks(2)^2*a^2 - 2*Es_import(1)*Es_import(2)*c_supplys(1)*c_supplys(2)*Ks(1)*a - 2*Es_import(1)*Es_import(2)*c_supplys(1)*c_supplys(2)*a^2 + 2*Es_import(1)*Es_import(2)*c_supplys(1)*Ks(1)*Ks(2)*a + 2*Es_import(1)*Es_import(2)*c_supplys(1)*Ks(2)*a^2 + 2*Es_import(1)*Es_import(2)*c_supplys(2)*Ks(1)*a^2 + 2*Es_import(1)*Es_import(2)*c_supplys(2)*a^3 - 2*Es_import(1)*Es_import(2)*Ks(1)*Ks(2)*a^2 - 2*Es_import(1)*Es_import(2)*Ks(2)*a^3 + Es_import(2)^2*c_supplys(1)^2*Ks(1)^2 + 2*Es_import(2)^2*c_supplys(1)^2*Ks(1)*a + Es_import(2)^2*c_supplys(1)^2*a^2 - 2*Es_import(2)^2*c_supplys(1)*Ks(1)^2*a - 4*Es_import(2)^2*c_supplys(1)*Ks(1)*a^2 - 2*Es_import(2)^2*c_supplys(1)*a^3 + Es_import(2)^2*Ks(1)^2*a^2 + 2*Es_import(2)^2*Ks(1)*a^3 + Es_import(2)^2*a^4)^(1/2) + Es_import(2)*a^2 - Es_import(2)*c_supplys(1)*Ks(1) + Es_import(1)*c_supplys(2)*a - Es_import(2)*c_supplys(1)*a - Es_import(1)*Ks(2)*a + Es_import(2)*Ks(1)*a)/(2*Es_import(1)*a);
((Es_import(1)^2*c_supplys(3)^2*a^2 + 2*Es_import(1)^2*c_supplys(3)*Ks(3)*a^2 + Es_import(1)^2*Ks(3)^2*a^2 - 2*Es_import(1)*Es_import(3)*c_supplys(1)*c_supplys(3)*Ks(1)*a - 2*Es_import(1)*Es_import(3)*c_supplys(1)*c_supplys(3)*a^2 + 2*Es_import(1)*Es_import(3)*c_supplys(1)*Ks(1)*Ks(3)*a + 2*Es_import(1)*Es_import(3)*c_supplys(1)*Ks(3)*a^2 + 2*Es_import(1)*Es_import(3)*c_supplys(3)*Ks(1)*a^2 + 2*Es_import(1)*Es_import(3)*c_supplys(3)*a^3 - 2*Es_import(1)*Es_import(3)*Ks(1)*Ks(3)*a^2 - 2*Es_import(1)*Es_import(3)*Ks(3)*a^3 + Es_import(3)^2*c_supplys(1)^2*Ks(1)^2 + 2*Es_import(3)^2*c_supplys(1)^2*Ks(1)*a + Es_import(3)^2*c_supplys(1)^2*a^2 - 2*Es_import(3)^2*c_supplys(1)*Ks(1)^2*a - 4*Es_import(3)^2*c_supplys(1)*Ks(1)*a^2 - 2*Es_import(3)^2*c_supplys(1)*a^3 + Es_import(3)^2*Ks(1)^2*a^2 + 2*Es_import(3)^2*Ks(1)*a^3 + Es_import(3)^2*a^4)^(1/2) + Es_import(3)*a^2 - Es_import(3)*c_supplys(1)*Ks(1) + Es_import(1)*c_supplys(3)*a - Es_import(3)*c_supplys(1)*a - Es_import(1)*Ks(3)*a + Es_import(3)*Ks(1)*a)/(2*Es_import(1)*a);];
GC_cb_dim3=@(c_a,c_b,d)Ks(3)/(((Es_convert(3)*Es_import(3))/(d+Es_convert(3)*k))/(d/(V_import*k)-(Es_convert(1)*Es_import(1))/(d+Es_convert(1)*k)*c_a/(c_a+Ks(1))-(Es_convert(2)*Es_import(2))/(d+Es_convert(2)*k)*c_b/(c_b+Ks(2)) )-1);
metabfuns.analytical_FBGCincell_dim3={FB_cb_dim3;
    GC_cb_dim3;
    []
    };

