function metabfuns=Metab_ConvertNecessary(species_para)

metabfuns=[];
r=species_para.r;
metabfuns.r=r;

% load parameters
intercell_dim=species_para.intercell_dim;

Gamma=species_para.Gamma;
Ks=reshape(species_para.Ks,2,1);
p=species_para.p;
Es=reshape(species_para.Es,4,1);% strategy
metabfuns.Edim=length(Es);

vs=Es(1:2);
convertM=[0,Es(3);% from b to a
    -Es(4),0]; % from a to b


% for ODEs
Intakes_givenE=@(cs,xs,Es)(Es(1:2)+p).*cs-p*xs;
Growth_givenE=@(cs,xs,Es)Gamma/sum(Ks./xs);
dxdt_givenE=@(cs,xs,Es)Intakes_givenE(cs,xs,Es)-Growth_givenE(cs,xs,Es).*Ks  + sum(([0,Es(3);-Es(4),0])*xs)*[1;-1];

Intakes=@(cs,xs)Intakes_givenE(cs,xs,Es);
Growth=@(cs,xs)Growth_givenE(cs,xs,Es);
dxdt=@(cs,xs)dxdt_givenE(cs,xs,Es);

metabfuns.Intakes=Intakes;
metabfuns.Growth=Growth;
metabfuns.withincell_dxdt=dxdt;
metabfuns.intercell_dim=intercell_dim;

% for the analytical steady states given environment 
K_1=Ks(1);
K_2=Ks(2);
E_1=Es(1);
E_2=Es(2);
E_3=Es(3);
E_4=Es(4);

 growth_givencab=@(c_a,c_b)(K_2*c_a*p^4 - p^2*(E_1^2*E_3^2*K_2^2*c_a^2 + 2*E_1^2*E_3*E_4*K_1*K_2*c_a^2 + 2*E_1^2*E_3*Gamma*K_2^2*c_a^2 + 2*E_1^2*E_3*K_2^2*c_a^2*p + E_1^2*E_4^2*K_1^2*c_a^2 - 2*E_1^2*E_4*Gamma*K_1*K_2*c_a^2 + 2*E_1^2*E_4*K_1*K_2*c_a^2*p + E_1^2*Gamma^2*K_2^2*c_a^2 + 2*E_1^2*Gamma*K_2^2*c_a^2*p + E_1^2*K_2^2*c_a^2*p^2 + 2*E_1*E_2*E_3^2*K_2^2*c_a*c_b + 4*E_1*E_2*E_3*E_4*K_1*K_2*c_a*c_b - 2*E_1*E_2*E_3*Gamma*K_1*K_2*c_a*c_b + 2*E_1*E_2*E_3*Gamma*K_2^2*c_a*c_b + 2*E_1*E_2*E_3*K_1*K_2*c_a*c_b*p + 2*E_1*E_2*E_3*K_2^2*c_a*c_b*p + 2*E_1*E_2*E_4^2*K_1^2*c_a*c_b + 2*E_1*E_2*E_4*Gamma*K_1^2*c_a*c_b - 2*E_1*E_2*E_4*Gamma*K_1*K_2*c_a*c_b + 2*E_1*E_2*E_4*K_1^2*c_a*c_b*p + 2*E_1*E_2*E_4*K_1*K_2*c_a*c_b*p - 2*E_1*E_2*Gamma^2*K_1*K_2*c_a*c_b - 4*E_1*E_2*Gamma*K_1*K_2*c_a*c_b*p + 2*E_1*E_2*K_1*K_2*c_a*c_b*p^2 + 2*E_1*E_3^2*K_2^2*c_a^2*p + 2*E_1*E_3^2*K_2^2*c_a*c_b*p + 4*E_1*E_3*E_4*K_1*K_2*c_a^2*p + 4*E_1*E_3*E_4*K_1*K_2*c_a*c_b*p - 2*E_1*E_3*Gamma*K_1*K_2*c_a*c_b*p + 4*E_1*E_3*Gamma*K_2^2*c_a^2*p + 2*E_1*E_3*Gamma*K_2^2*c_a*c_b*p + 2*E_1*E_3*K_1*K_2*c_a*c_b*p^2 + 4*E_1*E_3*K_2^2*c_a^2*p^2 + 2*E_1*E_3*K_2^2*c_a*c_b*p^2 + 2*E_1*E_4^2*K_1^2*c_a^2*p + 2*E_1*E_4^2*K_1^2*c_a*c_b*p + 2*E_1*E_4*Gamma*K_1^2*c_a*c_b*p - 4*E_1*E_4*Gamma*K_1*K_2*c_a^2*p - 2*E_1*E_4*Gamma*K_1*K_2*c_a*c_b*p + 2*E_1*E_4*K_1^2*c_a*c_b*p^2 + 4*E_1*E_4*K_1*K_2*c_a^2*p^2 + 2*E_1*E_4*K_1*K_2*c_a*c_b*p^2 - 2*E_1*Gamma^2*K_1*K_2*c_a*c_b*p + 2*E_1*Gamma^2*K_2^2*c_a^2*p - 4*E_1*Gamma*K_1*K_2*c_a*c_b*p^2 + 4*E_1*Gamma*K_2^2*c_a^2*p^2 + 2*E_1*K_1*K_2*c_a*c_b*p^3 + 2*E_1*K_2^2*c_a^2*p^3 + E_2^2*E_3^2*K_2^2*c_b^2 + 2*E_2^2*E_3*E_4*K_1*K_2*c_b^2 - 2*E_2^2*E_3*Gamma*K_1*K_2*c_b^2 + 2*E_2^2*E_3*K_1*K_2*c_b^2*p + E_2^2*E_4^2*K_1^2*c_b^2 + 2*E_2^2*E_4*Gamma*K_1^2*c_b^2 + 2*E_2^2*E_4*K_1^2*c_b^2*p + E_2^2*Gamma^2*K_1^2*c_b^2 + 2*E_2^2*Gamma*K_1^2*c_b^2*p + E_2^2*K_1^2*c_b^2*p^2 + 2*E_2*E_3^2*K_2^2*c_a*c_b*p + 2*E_2*E_3^2*K_2^2*c_b^2*p + 4*E_2*E_3*E_4*K_1*K_2*c_a*c_b*p + 4*E_2*E_3*E_4*K_1*K_2*c_b^2*p - 2*E_2*E_3*Gamma*K_1*K_2*c_a*c_b*p - 4*E_2*E_3*Gamma*K_1*K_2*c_b^2*p + 2*E_2*E_3*Gamma*K_2^2*c_a*c_b*p + 2*E_2*E_3*K_1*K_2*c_a*c_b*p^2 + 4*E_2*E_3*K_1*K_2*c_b^2*p^2 + 2*E_2*E_3*K_2^2*c_a*c_b*p^2 + 2*E_2*E_4^2*K_1^2*c_a*c_b*p + 2*E_2*E_4^2*K_1^2*c_b^2*p + 2*E_2*E_4*Gamma*K_1^2*c_a*c_b*p + 4*E_2*E_4*Gamma*K_1^2*c_b^2*p - 2*E_2*E_4*Gamma*K_1*K_2*c_a*c_b*p + 2*E_2*E_4*K_1^2*c_a*c_b*p^2 + 4*E_2*E_4*K_1^2*c_b^2*p^2 + 2*E_2*E_4*K_1*K_2*c_a*c_b*p^2 + 2*E_2*Gamma^2*K_1^2*c_b^2*p - 2*E_2*Gamma^2*K_1*K_2*c_a*c_b*p + 4*E_2*Gamma*K_1^2*c_b^2*p^2 - 4*E_2*Gamma*K_1*K_2*c_a*c_b*p^2 + 2*E_2*K_1^2*c_b^2*p^3 + 2*E_2*K_1*K_2*c_a*c_b*p^3 + E_3^2*K_2^2*c_a^2*p^2 + 2*E_3^2*K_2^2*c_a*c_b*p^2 + E_3^2*K_2^2*c_b^2*p^2 + 2*E_3*E_4*K_1*K_2*c_a^2*p^2 + 4*E_3*E_4*K_1*K_2*c_a*c_b*p^2 + 2*E_3*E_4*K_1*K_2*c_b^2*p^2 - 2*E_3*Gamma*K_1*K_2*c_a*c_b*p^2 - 2*E_3*Gamma*K_1*K_2*c_b^2*p^2 + 2*E_3*Gamma*K_2^2*c_a^2*p^2 + 2*E_3*Gamma*K_2^2*c_a*c_b*p^2 + 2*E_3*K_1*K_2*c_a*c_b*p^3 + 2*E_3*K_1*K_2*c_b^2*p^3 + 2*E_3*K_2^2*c_a^2*p^3 + 2*E_3*K_2^2*c_a*c_b*p^3 + E_4^2*K_1^2*c_a^2*p^2 + 2*E_4^2*K_1^2*c_a*c_b*p^2 + E_4^2*K_1^2*c_b^2*p^2 + 2*E_4*Gamma*K_1^2*c_a*c_b*p^2 + 2*E_4*Gamma*K_1^2*c_b^2*p^2 - 2*E_4*Gamma*K_1*K_2*c_a^2*p^2 - 2*E_4*Gamma*K_1*K_2*c_a*c_b*p^2 + 2*E_4*K_1^2*c_a*c_b*p^3 + 2*E_4*K_1^2*c_b^2*p^3 + 2*E_4*K_1*K_2*c_a^2*p^3 + 2*E_4*K_1*K_2*c_a*c_b*p^3 + Gamma^2*K_1^2*c_b^2*p^2 - 2*Gamma^2*K_1*K_2*c_a*c_b*p^2 + Gamma^2*K_2^2*c_a^2*p^2 + 2*Gamma*K_1^2*c_b^2*p^3 - 4*Gamma*K_1*K_2*c_a*c_b*p^3 + 2*Gamma*K_2^2*c_a^2*p^3 + K_1^2*c_b^2*p^4 + 2*K_1*K_2*c_a*c_b*p^4 + K_2^2*c_a^2*p^4)^(1/2) + K_1*c_b*p^4 - E_3*p*(E_1^2*E_3^2*K_2^2*c_a^2 + 2*E_1^2*E_3*E_4*K_1*K_2*c_a^2 + 2*E_1^2*E_3*Gamma*K_2^2*c_a^2 + 2*E_1^2*E_3*K_2^2*c_a^2*p + E_1^2*E_4^2*K_1^2*c_a^2 - 2*E_1^2*E_4*Gamma*K_1*K_2*c_a^2 + 2*E_1^2*E_4*K_1*K_2*c_a^2*p + E_1^2*Gamma^2*K_2^2*c_a^2 + 2*E_1^2*Gamma*K_2^2*c_a^2*p + E_1^2*K_2^2*c_a^2*p^2 + 2*E_1*E_2*E_3^2*K_2^2*c_a*c_b + 4*E_1*E_2*E_3*E_4*K_1*K_2*c_a*c_b - 2*E_1*E_2*E_3*Gamma*K_1*K_2*c_a*c_b + 2*E_1*E_2*E_3*Gamma*K_2^2*c_a*c_b + 2*E_1*E_2*E_3*K_1*K_2*c_a*c_b*p + 2*E_1*E_2*E_3*K_2^2*c_a*c_b*p + 2*E_1*E_2*E_4^2*K_1^2*c_a*c_b + 2*E_1*E_2*E_4*Gamma*K_1^2*c_a*c_b - 2*E_1*E_2*E_4*Gamma*K_1*K_2*c_a*c_b + 2*E_1*E_2*E_4*K_1^2*c_a*c_b*p + 2*E_1*E_2*E_4*K_1*K_2*c_a*c_b*p - 2*E_1*E_2*Gamma^2*K_1*K_2*c_a*c_b - 4*E_1*E_2*Gamma*K_1*K_2*c_a*c_b*p + 2*E_1*E_2*K_1*K_2*c_a*c_b*p^2 + 2*E_1*E_3^2*K_2^2*c_a^2*p + 2*E_1*E_3^2*K_2^2*c_a*c_b*p + 4*E_1*E_3*E_4*K_1*K_2*c_a^2*p + 4*E_1*E_3*E_4*K_1*K_2*c_a*c_b*p - 2*E_1*E_3*Gamma*K_1*K_2*c_a*c_b*p + 4*E_1*E_3*Gamma*K_2^2*c_a^2*p + 2*E_1*E_3*Gamma*K_2^2*c_a*c_b*p + 2*E_1*E_3*K_1*K_2*c_a*c_b*p^2 + 4*E_1*E_3*K_2^2*c_a^2*p^2 + 2*E_1*E_3*K_2^2*c_a*c_b*p^2 + 2*E_1*E_4^2*K_1^2*c_a^2*p + 2*E_1*E_4^2*K_1^2*c_a*c_b*p + 2*E_1*E_4*Gamma*K_1^2*c_a*c_b*p - 4*E_1*E_4*Gamma*K_1*K_2*c_a^2*p - 2*E_1*E_4*Gamma*K_1*K_2*c_a*c_b*p + 2*E_1*E_4*K_1^2*c_a*c_b*p^2 + 4*E_1*E_4*K_1*K_2*c_a^2*p^2 + 2*E_1*E_4*K_1*K_2*c_a*c_b*p^2 - 2*E_1*Gamma^2*K_1*K_2*c_a*c_b*p + 2*E_1*Gamma^2*K_2^2*c_a^2*p - 4*E_1*Gamma*K_1*K_2*c_a*c_b*p^2 + 4*E_1*Gamma*K_2^2*c_a^2*p^2 + 2*E_1*K_1*K_2*c_a*c_b*p^3 + 2*E_1*K_2^2*c_a^2*p^3 + E_2^2*E_3^2*K_2^2*c_b^2 + 2*E_2^2*E_3*E_4*K_1*K_2*c_b^2 - 2*E_2^2*E_3*Gamma*K_1*K_2*c_b^2 + 2*E_2^2*E_3*K_1*K_2*c_b^2*p + E_2^2*E_4^2*K_1^2*c_b^2 + 2*E_2^2*E_4*Gamma*K_1^2*c_b^2 + 2*E_2^2*E_4*K_1^2*c_b^2*p + E_2^2*Gamma^2*K_1^2*c_b^2 + 2*E_2^2*Gamma*K_1^2*c_b^2*p + E_2^2*K_1^2*c_b^2*p^2 + 2*E_2*E_3^2*K_2^2*c_a*c_b*p + 2*E_2*E_3^2*K_2^2*c_b^2*p + 4*E_2*E_3*E_4*K_1*K_2*c_a*c_b*p + 4*E_2*E_3*E_4*K_1*K_2*c_b^2*p - 2*E_2*E_3*Gamma*K_1*K_2*c_a*c_b*p - 4*E_2*E_3*Gamma*K_1*K_2*c_b^2*p + 2*E_2*E_3*Gamma*K_2^2*c_a*c_b*p + 2*E_2*E_3*K_1*K_2*c_a*c_b*p^2 + 4*E_2*E_3*K_1*K_2*c_b^2*p^2 + 2*E_2*E_3*K_2^2*c_a*c_b*p^2 + 2*E_2*E_4^2*K_1^2*c_a*c_b*p + 2*E_2*E_4^2*K_1^2*c_b^2*p + 2*E_2*E_4*Gamma*K_1^2*c_a*c_b*p + 4*E_2*E_4*Gamma*K_1^2*c_b^2*p - 2*E_2*E_4*Gamma*K_1*K_2*c_a*c_b*p + 2*E_2*E_4*K_1^2*c_a*c_b*p^2 + 4*E_2*E_4*K_1^2*c_b^2*p^2 + 2*E_2*E_4*K_1*K_2*c_a*c_b*p^2 + 2*E_2*Gamma^2*K_1^2*c_b^2*p - 2*E_2*Gamma^2*K_1*K_2*c_a*c_b*p + 4*E_2*Gamma*K_1^2*c_b^2*p^2 - 4*E_2*Gamma*K_1*K_2*c_a*c_b*p^2 + 2*E_2*K_1^2*c_b^2*p^3 + 2*E_2*K_1*K_2*c_a*c_b*p^3 + E_3^2*K_2^2*c_a^2*p^2 + 2*E_3^2*K_2^2*c_a*c_b*p^2 + E_3^2*K_2^2*c_b^2*p^2 + 2*E_3*E_4*K_1*K_2*c_a^2*p^2 + 4*E_3*E_4*K_1*K_2*c_a*c_b*p^2 + 2*E_3*E_4*K_1*K_2*c_b^2*p^2 - 2*E_3*Gamma*K_1*K_2*c_a*c_b*p^2 - 2*E_3*Gamma*K_1*K_2*c_b^2*p^2 + 2*E_3*Gamma*K_2^2*c_a^2*p^2 + 2*E_3*Gamma*K_2^2*c_a*c_b*p^2 + 2*E_3*K_1*K_2*c_a*c_b*p^3 + 2*E_3*K_1*K_2*c_b^2*p^3 + 2*E_3*K_2^2*c_a^2*p^3 + 2*E_3*K_2^2*c_a*c_b*p^3 + E_4^2*K_1^2*c_a^2*p^2 + 2*E_4^2*K_1^2*c_a*c_b*p^2 + E_4^2*K_1^2*c_b^2*p^2 + 2*E_4*Gamma*K_1^2*c_a*c_b*p^2 + 2*E_4*Gamma*K_1^2*c_b^2*p^2 - 2*E_4*Gamma*K_1*K_2*c_a^2*p^2 - 2*E_4*Gamma*K_1*K_2*c_a*c_b*p^2 + 2*E_4*K_1^2*c_a*c_b*p^3 + 2*E_4*K_1^2*c_b^2*p^3 + 2*E_4*K_1*K_2*c_a^2*p^3 + 2*E_4*K_1*K_2*c_a*c_b*p^3 + Gamma^2*K_1^2*c_b^2*p^2 - 2*Gamma^2*K_1*K_2*c_a*c_b*p^2 + Gamma^2*K_2^2*c_a^2*p^2 + 2*Gamma*K_1^2*c_b^2*p^3 - 4*Gamma*K_1*K_2*c_a*c_b*p^3 + 2*Gamma*K_2^2*c_a^2*p^3 + K_1^2*c_b^2*p^4 + 2*K_1*K_2*c_a*c_b*p^4 + K_2^2*c_a^2*p^4)^(1/2) - E_4*p*(E_1^2*E_3^2*K_2^2*c_a^2 + 2*E_1^2*E_3*E_4*K_1*K_2*c_a^2 + 2*E_1^2*E_3*Gamma*K_2^2*c_a^2 + 2*E_1^2*E_3*K_2^2*c_a^2*p + E_1^2*E_4^2*K_1^2*c_a^2 - 2*E_1^2*E_4*Gamma*K_1*K_2*c_a^2 + 2*E_1^2*E_4*K_1*K_2*c_a^2*p + E_1^2*Gamma^2*K_2^2*c_a^2 + 2*E_1^2*Gamma*K_2^2*c_a^2*p + E_1^2*K_2^2*c_a^2*p^2 + 2*E_1*E_2*E_3^2*K_2^2*c_a*c_b + 4*E_1*E_2*E_3*E_4*K_1*K_2*c_a*c_b - 2*E_1*E_2*E_3*Gamma*K_1*K_2*c_a*c_b + 2*E_1*E_2*E_3*Gamma*K_2^2*c_a*c_b + 2*E_1*E_2*E_3*K_1*K_2*c_a*c_b*p + 2*E_1*E_2*E_3*K_2^2*c_a*c_b*p + 2*E_1*E_2*E_4^2*K_1^2*c_a*c_b + 2*E_1*E_2*E_4*Gamma*K_1^2*c_a*c_b - 2*E_1*E_2*E_4*Gamma*K_1*K_2*c_a*c_b + 2*E_1*E_2*E_4*K_1^2*c_a*c_b*p + 2*E_1*E_2*E_4*K_1*K_2*c_a*c_b*p - 2*E_1*E_2*Gamma^2*K_1*K_2*c_a*c_b - 4*E_1*E_2*Gamma*K_1*K_2*c_a*c_b*p + 2*E_1*E_2*K_1*K_2*c_a*c_b*p^2 + 2*E_1*E_3^2*K_2^2*c_a^2*p + 2*E_1*E_3^2*K_2^2*c_a*c_b*p + 4*E_1*E_3*E_4*K_1*K_2*c_a^2*p + 4*E_1*E_3*E_4*K_1*K_2*c_a*c_b*p - 2*E_1*E_3*Gamma*K_1*K_2*c_a*c_b*p + 4*E_1*E_3*Gamma*K_2^2*c_a^2*p + 2*E_1*E_3*Gamma*K_2^2*c_a*c_b*p + 2*E_1*E_3*K_1*K_2*c_a*c_b*p^2 + 4*E_1*E_3*K_2^2*c_a^2*p^2 + 2*E_1*E_3*K_2^2*c_a*c_b*p^2 + 2*E_1*E_4^2*K_1^2*c_a^2*p + 2*E_1*E_4^2*K_1^2*c_a*c_b*p + 2*E_1*E_4*Gamma*K_1^2*c_a*c_b*p - 4*E_1*E_4*Gamma*K_1*K_2*c_a^2*p - 2*E_1*E_4*Gamma*K_1*K_2*c_a*c_b*p + 2*E_1*E_4*K_1^2*c_a*c_b*p^2 + 4*E_1*E_4*K_1*K_2*c_a^2*p^2 + 2*E_1*E_4*K_1*K_2*c_a*c_b*p^2 - 2*E_1*Gamma^2*K_1*K_2*c_a*c_b*p + 2*E_1*Gamma^2*K_2^2*c_a^2*p - 4*E_1*Gamma*K_1*K_2*c_a*c_b*p^2 + 4*E_1*Gamma*K_2^2*c_a^2*p^2 + 2*E_1*K_1*K_2*c_a*c_b*p^3 + 2*E_1*K_2^2*c_a^2*p^3 + E_2^2*E_3^2*K_2^2*c_b^2 + 2*E_2^2*E_3*E_4*K_1*K_2*c_b^2 - 2*E_2^2*E_3*Gamma*K_1*K_2*c_b^2 + 2*E_2^2*E_3*K_1*K_2*c_b^2*p + E_2^2*E_4^2*K_1^2*c_b^2 + 2*E_2^2*E_4*Gamma*K_1^2*c_b^2 + 2*E_2^2*E_4*K_1^2*c_b^2*p + E_2^2*Gamma^2*K_1^2*c_b^2 + 2*E_2^2*Gamma*K_1^2*c_b^2*p + E_2^2*K_1^2*c_b^2*p^2 + 2*E_2*E_3^2*K_2^2*c_a*c_b*p + 2*E_2*E_3^2*K_2^2*c_b^2*p + 4*E_2*E_3*E_4*K_1*K_2*c_a*c_b*p + 4*E_2*E_3*E_4*K_1*K_2*c_b^2*p - 2*E_2*E_3*Gamma*K_1*K_2*c_a*c_b*p - 4*E_2*E_3*Gamma*K_1*K_2*c_b^2*p + 2*E_2*E_3*Gamma*K_2^2*c_a*c_b*p + 2*E_2*E_3*K_1*K_2*c_a*c_b*p^2 + 4*E_2*E_3*K_1*K_2*c_b^2*p^2 + 2*E_2*E_3*K_2^2*c_a*c_b*p^2 + 2*E_2*E_4^2*K_1^2*c_a*c_b*p + 2*E_2*E_4^2*K_1^2*c_b^2*p + 2*E_2*E_4*Gamma*K_1^2*c_a*c_b*p + 4*E_2*E_4*Gamma*K_1^2*c_b^2*p - 2*E_2*E_4*Gamma*K_1*K_2*c_a*c_b*p + 2*E_2*E_4*K_1^2*c_a*c_b*p^2 + 4*E_2*E_4*K_1^2*c_b^2*p^2 + 2*E_2*E_4*K_1*K_2*c_a*c_b*p^2 + 2*E_2*Gamma^2*K_1^2*c_b^2*p - 2*E_2*Gamma^2*K_1*K_2*c_a*c_b*p + 4*E_2*Gamma*K_1^2*c_b^2*p^2 - 4*E_2*Gamma*K_1*K_2*c_a*c_b*p^2 + 2*E_2*K_1^2*c_b^2*p^3 + 2*E_2*K_1*K_2*c_a*c_b*p^3 + E_3^2*K_2^2*c_a^2*p^2 + 2*E_3^2*K_2^2*c_a*c_b*p^2 + E_3^2*K_2^2*c_b^2*p^2 + 2*E_3*E_4*K_1*K_2*c_a^2*p^2 + 4*E_3*E_4*K_1*K_2*c_a*c_b*p^2 + 2*E_3*E_4*K_1*K_2*c_b^2*p^2 - 2*E_3*Gamma*K_1*K_2*c_a*c_b*p^2 - 2*E_3*Gamma*K_1*K_2*c_b^2*p^2 + 2*E_3*Gamma*K_2^2*c_a^2*p^2 + 2*E_3*Gamma*K_2^2*c_a*c_b*p^2 + 2*E_3*K_1*K_2*c_a*c_b*p^3 + 2*E_3*K_1*K_2*c_b^2*p^3 + 2*E_3*K_2^2*c_a^2*p^3 + 2*E_3*K_2^2*c_a*c_b*p^3 + E_4^2*K_1^2*c_a^2*p^2 + 2*E_4^2*K_1^2*c_a*c_b*p^2 + E_4^2*K_1^2*c_b^2*p^2 + 2*E_4*Gamma*K_1^2*c_a*c_b*p^2 + 2*E_4*Gamma*K_1^2*c_b^2*p^2 - 2*E_4*Gamma*K_1*K_2*c_a^2*p^2 - 2*E_4*Gamma*K_1*K_2*c_a*c_b*p^2 + 2*E_4*K_1^2*c_a*c_b*p^3 + 2*E_4*K_1^2*c_b^2*p^3 + 2*E_4*K_1*K_2*c_a^2*p^3 + 2*E_4*K_1*K_2*c_a*c_b*p^3 + Gamma^2*K_1^2*c_b^2*p^2 - 2*Gamma^2*K_1*K_2*c_a*c_b*p^2 + Gamma^2*K_2^2*c_a^2*p^2 + 2*Gamma*K_1^2*c_b^2*p^3 - 4*Gamma*K_1*K_2*c_a*c_b*p^3 + 2*Gamma*K_2^2*c_a^2*p^3 + K_1^2*c_b^2*p^4 + 2*K_1*K_2*c_a*c_b*p^4 + K_2^2*c_a^2*p^4)^(1/2) + E_1*K_2*c_a*p^3 + 2*E_3*K_2*c_a*p^3 + E_4*K_1*c_a*p^3 + E_4*K_2*c_a*p^3 + E_2*K_1*c_b*p^3 + E_3*K_1*c_b*p^3 + E_3*K_2*c_b*p^3 + 2*E_4*K_1*c_b*p^3 + Gamma*K_2*c_a*p^3 + Gamma*K_1*c_b*p^3 + E_3^2*K_2*c_a*p^2 + E_4^2*K_1*c_a*p^2 + E_3^2*K_2*c_b*p^2 + E_4^2*K_1*c_b*p^2 + 2*E_1*E_3*K_2*c_a*p^2 + E_1*E_4*K_1*c_a*p^2 + E_1*E_3^2*K_2*c_a*p + E_1*E_4^2*K_1*c_a*p + E_1*E_4*K_2*c_a*p^2 + E_3*E_4*K_1*c_a*p^2 + E_3*E_4*K_2*c_a*p^2 + E_2*E_3*K_1*c_b*p^2 + E_2*E_3*K_2*c_b*p^2 + 2*E_2*E_4*K_1*c_b*p^2 + E_2*E_3^2*K_2*c_b*p + E_2*E_4^2*K_1*c_b*p + E_3*E_4*K_1*c_b*p^2 + E_3*E_4*K_2*c_b*p^2 + E_1*Gamma*K_2*c_a*p^2 + E_3*Gamma*K_2*c_a*p^2 + 2*E_4*Gamma*K_1*c_a*p^2 + E_4*Gamma*K_2*c_a*p^2 + E_2*Gamma*K_1*c_b*p^2 + E_3*Gamma*K_1*c_b*p^2 + 2*E_3*Gamma*K_2*c_b*p^2 + E_4*Gamma*K_1*c_b*p^2 + 2*E_1*E_3*E_4*Gamma*K_1*c_a + 2*E_1*E_3*E_4*Gamma*K_2*c_a + 2*E_2*E_3*E_4*Gamma*K_1*c_b + 2*E_2*E_3*E_4*Gamma*K_2*c_b + E_1*E_3*E_4*K_1*c_a*p + E_1*E_3*E_4*K_2*c_a*p + E_2*E_3*E_4*K_1*c_b*p + E_2*E_3*E_4*K_2*c_b*p + E_1*E_3*Gamma*K_2*c_a*p + 2*E_1*E_4*Gamma*K_1*c_a*p + E_1*E_4*Gamma*K_2*c_a*p + 2*E_3*E_4*Gamma*K_1*c_a*p + 2*E_3*E_4*Gamma*K_2*c_a*p + E_2*E_3*Gamma*K_1*c_b*p + 2*E_2*E_3*Gamma*K_2*c_b*p + E_2*E_4*Gamma*K_1*c_b*p + 2*E_3*E_4*Gamma*K_1*c_b*p + 2*E_3*E_4*Gamma*K_2*c_b*p)/(2*(p^2*(E_3*K_2^2 + E_4*K_1^2 + 3*E_3*K_1*K_2 + 3*E_4*K_1*K_2) + Gamma*(K_1*p + E_3*K_1 + E_3*K_2)*(K_2*p + E_4*K_1 + E_4*K_2) + 2*K_1*K_2*p^3 + p*(E_3*K_2 + E_4*K_1)*(E_3 + E_4)*(K_1 + K_2)));
 xs_givencabg=@(c_a,c_b,g)[(c_a*p^2+ E_1*E_3*c_a+ E_2*E_3*c_b- E_3*K_1*g - E_3*K_2*g + E_1*c_a*p + E_3*c_a*p + E_3*c_b*p - K_1*g*p)/(p*(E_3+ E_4+ p) );
     (c_b*p^2+ E_1*E_4*c_a+ E_2*E_4*c_b- E_4*K_1*g - E_4*K_2*g + E_4*c_a*p + E_2*c_b*p + E_4*c_b*p - K_2*g*p)/(p*(E_3+ E_4+ p) )];
 
 analytical_steadygivenenv=@(env)[growth_givencab(env(1),env(2));
 xs_givencabg(env(1),env(2),growth_givencab(env(1),env(2)))];
 metabfuns.analytical_steadygivenenv=analytical_steadygivenenv;

 % if E is unknown
 
%Growth_givenE=@(cs,xs,Es)G(cs(1),cs(2),xs,Es);
%dxdt_givenE=@(cs,xs,Es)dP_dt(cs(1),cs(2),xs,Es);


% FB depends on d
 if E_3>10^-4
     GC_W_dim2=@(c_a,g)(p^2*(E_1^2*Gamma^2*c_a^2 + 2*E_1*E_3*Gamma*K_2*c_a*g - 2*E_1*E_4*Gamma*K_1*c_a*g - 2*E_1*Gamma^2*K_1*c_a*g + 2*E_1*Gamma^2*c_a^2*p - 2*E_1*Gamma*K_1*c_a*g*p + E_3^2*K_2^2*g^2 + 2*E_3*E_4*K_1*K_2*g^2 - 2*E_3*Gamma*K_1*K_2*g^2 + 2*E_3*Gamma*K_2*c_a*g*p + 2*E_3*K_1*K_2*g^2*p + E_4^2*K_1^2*g^2 + 2*E_4*Gamma*K_1^2*g^2 - 2*E_4*Gamma*K_1*c_a*g*p + 2*E_4*K_1^2*g^2*p + Gamma^2*K_1^2*g^2 - 2*Gamma^2*K_1*c_a*g*p + Gamma^2*c_a^2*p^2 + 2*Gamma*K_1^2*g^2*p - 2*Gamma*K_1*c_a*g*p^2 + K_1^2*g^2*p^2)^(1/2) + E_3*p*(E_1^2*Gamma^2*c_a^2 + 2*E_1*E_3*Gamma*K_2*c_a*g - 2*E_1*E_4*Gamma*K_1*c_a*g - 2*E_1*Gamma^2*K_1*c_a*g + 2*E_1*Gamma^2*c_a^2*p - 2*E_1*Gamma*K_1*c_a*g*p + E_3^2*K_2^2*g^2 + 2*E_3*E_4*K_1*K_2*g^2 - 2*E_3*Gamma*K_1*K_2*g^2 + 2*E_3*Gamma*K_2*c_a*g*p + 2*E_3*K_1*K_2*g^2*p + E_4^2*K_1^2*g^2 + 2*E_4*Gamma*K_1^2*g^2 - 2*E_4*Gamma*K_1*c_a*g*p + 2*E_4*K_1^2*g^2*p + Gamma^2*K_1^2*g^2 - 2*Gamma^2*K_1*c_a*g*p + Gamma^2*c_a^2*p^2 + 2*Gamma*K_1^2*g^2*p - 2*Gamma*K_1*c_a*g*p^2 + K_1^2*g^2*p^2)^(1/2) + E_4*p*(E_1^2*Gamma^2*c_a^2 + 2*E_1*E_3*Gamma*K_2*c_a*g - 2*E_1*E_4*Gamma*K_1*c_a*g - 2*E_1*Gamma^2*K_1*c_a*g + 2*E_1*Gamma^2*c_a^2*p - 2*E_1*Gamma*K_1*c_a*g*p + E_3^2*K_2^2*g^2 + 2*E_3*E_4*K_1*K_2*g^2 - 2*E_3*Gamma*K_1*K_2*g^2 + 2*E_3*Gamma*K_2*c_a*g*p + 2*E_3*K_1*K_2*g^2*p + E_4^2*K_1^2*g^2 + 2*E_4*Gamma*K_1^2*g^2 - 2*E_4*Gamma*K_1*c_a*g*p + 2*E_4*K_1^2*g^2*p + Gamma^2*K_1^2*g^2 - 2*Gamma^2*K_1*c_a*g*p + Gamma^2*c_a^2*p^2 + 2*Gamma*K_1^2*g^2*p - 2*Gamma*K_1*c_a*g*p^2 + K_1^2*g^2*p^2)^(1/2) - Gamma*c_a*p^3 + K_1*g*p^3 - E_1*Gamma*c_a*p^2 - E_3*Gamma*c_a*p^2 - E_4*Gamma*c_a*p^2 + E_3*K_1*g*p^2 + E_3*K_2*g*p^2 + 2*E_4*K_1*g*p^2 + E_3^2*K_2*g*p + E_4^2*K_1*g*p + Gamma*K_1*g*p^2 - E_1*E_3*Gamma*c_a*p - E_1*E_4*Gamma*c_a*p - 2*E_3*E_4*Gamma*c_a*p + E_3*E_4*K_1*g*p + E_3*E_4*K_2*g*p + E_3*Gamma*K_1*g*p + 2*E_3*Gamma*K_2*g*p + E_4*Gamma*K_1*g*p - 2*E_1*E_3*E_4*Gamma*c_a + 2*E_3*E_4*Gamma*K_1*g + 2*E_3*E_4*Gamma*K_2*g)/(2*E_3*Gamma*(E_2 + p)*(E_4 + p));
 
 else
     GC_W_dim2=@(c_a,g)-(K_2/(K_1/(p*(E_1*c_a - K_1*g + c_a*p)) - Gamma/(g*p*(E_4 + p))) + E_1*E_4*c_a - E_4*K_1*g - E_4*K_2*g + E_4*c_a*p - K_2*g*p)/(E_2*p + E_4*p + p^2 + E_2*E_4);
 end
 
 incellvfun=[];
 incellvfun=@(c_a,c_b,g)[(c_a*p^2 + E_1*E_3*c_a + E_2*E_3*c_b - E_3*K_1*g - E_3*K_2*g + E_1*c_a*p + E_3*c_a*p + E_3*c_b*p - K_1*g*p)/(E_3*p + E_4*p + p^2);
    (c_b*p^2 + E_1*E_4*c_a + E_2*E_4*c_b - E_4*K_1*g - E_4*K_2*g + E_4*c_a*p + E_2*c_b*p + E_4*c_b*p - K_2*g*p)/(E_3*p + E_4*p + p^2)
 ];

metabfuns.analytical_FBGCincell_dim2={[],GC_W_dim2,incellvfun};
  
species_id=0;
if E_2<10^-4 && E_3<10^-4
    species_id=1;
elseif E_1<10^-4 && E_4<10^-4
    species_id=2;
elseif E_3<10^-4 && E_4<10^-4
    species_id=3;
end
 

% general, no short cut
metabfuns.Growth_givenE=Growth_givenE;
metabfuns.dxdt_givenE=dxdt_givenE;
    
if species_id==1 % a->b converter
    
   
elseif species_id==2 % b->a converter
    
elseif species_id==3 % importer
    
    maxg_givenenv=@(c_a,c_b)(Gamma*c_a*c_b*(2*p + 1))/(Gamma*K_2*c_a + Gamma*K_1*c_b + K_2*c_a*p + K_1*c_b*p + 2*K_1^(1/2)*K_2^(1/2)*c_a^(1/2)*c_b^(1/2)*p);
    optE1=@(c_a,c_b)(Gamma*K_1*c_b + K_1*c_b*p - K_2*c_a*p^2 + K_1*c_b*p^2 - Gamma*K_2*c_a*p + Gamma*K_1*c_b*p + K_1^(1/2)*K_2^(1/2)*c_a^(1/2)*c_b^(1/2)*p)/(Gamma*K_2*c_a + Gamma*K_1*c_b + K_2*c_a*p + K_1*c_b*p + 2*K_1^(1/2)*K_2^(1/2)*c_a^(1/2)*c_b^(1/2)*p);

    analytical_opt_givenenv=@(env)[maxg_givenenv(env(1),env(2));
        optE1(env(1),env(2));
        1-optE1(env(1),env(2));
        0;
        0;];
    
    metabfuns.analytical_opt_givenenv=analytical_opt_givenenv;
end
