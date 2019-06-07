function metabfuns=Metab_assimilation(species_para)

metabfuns=[];
r=species_para.r;
metabfuns.r=r;


% load parameters
Ks=species_para.Ks;
V_impot=species_para.V_import;
env_dim=length(Ks);

r=species_para.r;
k=species_para.k;
Es=species_para.Es;

Es_impot=reshape(Es(1:env_dim),env_dim,1);
Es_convet=reshape(Es((env_dim+1):end),env_dim,1);


% Define: Intake rates as a vector
J1=@(cs,xs)V_impot*Es_impot.*cs./(cs+Ks);
Intakes=@(cs,xs)J1(cs,xs);
% Define: Growth  rate:
Growth=@(cs,xs)sum(k*Es_convet.*xs)/r;
% Define: Within cell varibles
dxdt=@(cs,xs)J1(cs,xs)-Growth(cs,xs).*xs-k*Es_convet.*xs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metabfuns.Intakes=Intakes;
metabfuns.Growth=Growth;
metabfuns.withincell_dxdt=dxdt;
metabfuns.intercell_dim=species_para.intercell_dim;


FB_cb_dim3=@(a,c_supplys) [((Es_impot(1)^2*c_supplys(2)^2*a^2 + 2*Es_impot(1)^2*c_supplys(2)*Ks(2)*a^2 + Es_impot(1)^2*Ks(2)^2*a^2 - 2*Es_impot(1)*Es_impot(2)*c_supplys(1)*c_supplys(2)*Ks(1)*a - 2*Es_impot(1)*Es_impot(2)*c_supplys(1)*c_supplys(2)*a^2 + 2*Es_impot(1)*Es_impot(2)*c_supplys(1)*Ks(1)*Ks(2)*a + 2*Es_impot(1)*Es_impot(2)*c_supplys(1)*Ks(2)*a^2 + 2*Es_impot(1)*Es_impot(2)*c_supplys(2)*Ks(1)*a^2 + 2*Es_impot(1)*Es_impot(2)*c_supplys(2)*a^3 - 2*Es_impot(1)*Es_impot(2)*Ks(1)*Ks(2)*a^2 - 2*Es_impot(1)*Es_impot(2)*Ks(2)*a^3 + Es_impot(2)^2*c_supplys(1)^2*Ks(1)^2 + 2*Es_impot(2)^2*c_supplys(1)^2*Ks(1)*a + Es_impot(2)^2*c_supplys(1)^2*a^2 - 2*Es_impot(2)^2*c_supplys(1)*Ks(1)^2*a - 4*Es_impot(2)^2*c_supplys(1)*Ks(1)*a^2 - 2*Es_impot(2)^2*c_supplys(1)*a^3 + Es_impot(2)^2*Ks(1)^2*a^2 + 2*Es_impot(2)^2*Ks(1)*a^3 + Es_impot(2)^2*a^4)^(1/2) + Es_impot(2)*a^2 - Es_impot(2)*c_supplys(1)*Ks(1) + Es_impot(1)*c_supplys(2)*a - Es_impot(2)*c_supplys(1)*a - Es_impot(1)*Ks(2)*a + Es_impot(2)*Ks(1)*a)/(2*Es_impot(1)*a);
((Es_impot(1)^2*c_supplys(3)^2*a^2 + 2*Es_impot(1)^2*c_supplys(3)*Ks(3)*a^2 + Es_impot(1)^2*Ks(3)^2*a^2 - 2*Es_impot(1)*Es_impot(3)*c_supplys(1)*c_supplys(3)*Ks(1)*a - 2*Es_impot(1)*Es_impot(3)*c_supplys(1)*c_supplys(3)*a^2 + 2*Es_impot(1)*Es_impot(3)*c_supplys(1)*Ks(1)*Ks(3)*a + 2*Es_impot(1)*Es_impot(3)*c_supplys(1)*Ks(3)*a^2 + 2*Es_impot(1)*Es_impot(3)*c_supplys(3)*Ks(1)*a^2 + 2*Es_impot(1)*Es_impot(3)*c_supplys(3)*a^3 - 2*Es_impot(1)*Es_impot(3)*Ks(1)*Ks(3)*a^2 - 2*Es_impot(1)*Es_impot(3)*Ks(3)*a^3 + Es_impot(3)^2*c_supplys(1)^2*Ks(1)^2 + 2*Es_impot(3)^2*c_supplys(1)^2*Ks(1)*a + Es_impot(3)^2*c_supplys(1)^2*a^2 - 2*Es_impot(3)^2*c_supplys(1)*Ks(1)^2*a - 4*Es_impot(3)^2*c_supplys(1)*Ks(1)*a^2 - 2*Es_impot(3)^2*c_supplys(1)*a^3 + Es_impot(3)^2*Ks(1)^2*a^2 + 2*Es_impot(3)^2*Ks(1)*a^3 + Es_impot(3)^2*a^4)^(1/2) + Es_impot(3)*a^2 - Es_impot(3)*c_supplys(1)*Ks(1) + Es_impot(1)*c_supplys(3)*a - Es_impot(3)*c_supplys(1)*a - Es_impot(1)*Ks(3)*a + Es_impot(3)*Ks(1)*a)/(2*Es_impot(1)*a);];
GC_cb_dim3=@(c_a,c_b,d)Ks(2)/(V_impot/((d*r/k-Es_convet(1)*Es_impot(1)*V_impot*c_a/(c_a+Ks(1))/(d+k*Es_convet(1)) - Es_convet(2)*Es_impot(2)*V_impot*c_b/(c_b+Ks(2))/(d+k*Es_convet(2))   )*(d+k*Es_convet(3))/(Es_convet(3)*Es_impot(3)))-1);
metabfuns.analytical_FBGCincell_dim3={FB_cb_dim3;
    GC_cb_dim3;
    []
    };

