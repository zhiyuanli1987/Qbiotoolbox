function metabfuns=Metab_crossfeeding(species_para)

zerothresh=10^-3;

metabfuns=[];
% load parameters

r=species_para.r;
metabfuns.r=r;

V_1=species_para.V_1;
V_2=species_para.V_2;
k=species_para.k;
V_4=species_para.V_4;

Es=species_para.Es;
K_1=species_para.Ks(1);
K_2=species_para.Ks(2);
K_3=species_para.Ks(3);
K_4=species_para.Ks(4);
K_5=species_para.Ks(5);
K_6=species_para.Ks(6);
K_7=species_para.Ks(7);



% one within cell varible: P
% define functions
J_1=@(S,W,P,Es)max(0,Es(1)*V_1*(S-P/K_3 )/(K_1+S+P/K_5 )); % first reaction 
J_2=@(S,W,P,Es)Es(2)*V_2*P/(P+K_2 ); % second reaction
J_3=@(S,W,P,Es)Es(3)*k*(P-W); % diffuse out
J_4=@(S,W,P,Es)V_4*Es(4)*(W-P/K_6)/(W+K_4+P/K_7 ); % active intake

Intake_S=@(S,W,P,Es)J_1(S,W,P,Es);
Intake_W=@(S,W,P,Es)J_4(S,W,P,Es)-J_3(S,W,P,Es);
G=@(S,W,P,Es)J_1(S,W,P,Es)+J_2(S,W,P,Es);
dP_dt=@(S,W,P,Es)J_1(S,W,P,Es)-J_2(S,W,P,Es)+J_4(S,W,P,Es)-J_3(S,W,P,Es);%-G(S,W,P,Es)*P;

% there is only one within cell dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define: Intake rates as a vector
Intakes=@(cs,xs)[Intake_S(cs(1),cs(2),xs,Es);
    Intake_W(cs(1),cs(2),xs,Es);];

% Define: Growth  rate:
Growth=@(cs,xs)G(cs(1),cs(2),xs,Es);
dxdt=@(cs,xs)dP_dt(cs(1),cs(2),xs,Es);

Growth_givenE=@(cs,xs,Es)G(cs(1),cs(2),xs,Es);
dxdt_givenE=@(cs,xs,Es)dP_dt(cs(1),cs(2),xs,Es);

if Es(2)<zerothresh && Es(4)<zerothresh
    speciesid=1; % G1
elseif Es(3)<zerothresh && Es(4)<zerothresh
    speciesid=2; % G
elseif Es(1)<zerothresh && Es(3)<zerothresh
    speciesid=3; % G
else
    speciesid=0; % not a particular species 
end

%%%%%%%%%%%%%% FLUX balance and growth contours %%%%%%%%%%%%%%%%%%%
if speciesid>0
    E_1=Es(1);
    E_2=Es(2);
    E_3=Es(3);
    E_4=Es(4);


  
    switch speciesid
        case 1 % S1
            % W as a function of S
            steadyP_givenSW=@(S,W,d)-(Es(1)*K_4*K_5*V_1 - Es(3)*K_3*W^2*k - (Es(1)^2*K_4^2*K_5^2*V_1^2 + 2*Es(1)^2*K_4*K_5^2*V_1^2*W + Es(1)^2*K_5^2*V_1^2*W^2 + 2*Es(1)*Es(3)*K_1*K_3*K_4^2*K_5^2*V_1*k + 4*Es(1)*Es(3)*K_1*K_3*K_4*K_5^2*V_1*W*k + 2*Es(1)*Es(3)*K_1*K_3*K_5^2*V_1*W^2*k + 4*Es(1)*Es(3)*K_3^2*K_4^2*K_5*S*V_1*k + 8*Es(1)*Es(3)*K_3^2*K_4*K_5*S*V_1*W*k + 4*Es(1)*Es(3)*K_3^2*K_5*S*V_1*W^2*k + 2*Es(1)*Es(3)*K_3*K_4^2*K_5^2*S*V_1*k - 2*Es(1)*Es(3)*K_3*K_4^2*K_5*V_1*W*k + 4*Es(1)*Es(3)*K_3*K_4*K_5^2*S*V_1*W*k - 4*Es(1)*Es(3)*K_3*K_4*K_5*V_1*W^2*k + 2*Es(1)*Es(3)*K_3*K_5^2*S*V_1*W^2*k - 2*Es(1)*Es(3)*K_3*K_5*V_1*W^3*k - 2*Es(1)*Es(4)*K_3*K_4*K_5*V_1*V_4*W - 2*Es(1)*Es(4)*K_3*K_5*V_1*V_4*W^2 + Es(3)^2*K_1^2*K_3^2*K_4^2*K_5^2*k^2 + 2*Es(3)^2*K_1^2*K_3^2*K_4*K_5^2*W*k^2 + Es(3)^2*K_1^2*K_3^2*K_5^2*W^2*k^2 + 2*Es(3)^2*K_1*K_3^2*K_4^2*K_5^2*S*k^2 + 2*Es(3)^2*K_1*K_3^2*K_4^2*K_5*W*k^2 + 4*Es(3)^2*K_1*K_3^2*K_4*K_5^2*S*W*k^2 + 4*Es(3)^2*K_1*K_3^2*K_4*K_5*W^2*k^2 + 2*Es(3)^2*K_1*K_3^2*K_5^2*S*W^2*k^2 + 2*Es(3)^2*K_1*K_3^2*K_5*W^3*k^2 + Es(3)^2*K_3^2*K_4^2*K_5^2*S^2*k^2 + 2*Es(3)^2*K_3^2*K_4^2*K_5*S*W*k^2 + Es(3)^2*K_3^2*K_4^2*W^2*k^2 + 2*Es(3)^2*K_3^2*K_4*K_5^2*S^2*W*k^2 + 4*Es(3)^2*K_3^2*K_4*K_5*S*W^2*k^2 + 2*Es(3)^2*K_3^2*K_4*W^3*k^2 + Es(3)^2*K_3^2*K_5^2*S^2*W^2*k^2 + 2*Es(3)^2*K_3^2*K_5*S*W^3*k^2 + Es(3)^2*K_3^2*W^4*k^2 + 2*Es(3)*Es(4)*K_1*K_3^2*K_4*K_5*V_4*W*k + 2*Es(3)*Es(4)*K_1*K_3^2*K_5*V_4*W^2*k + 2*Es(3)*Es(4)*K_3^2*K_4*K_5*S*V_4*W*k + 2*Es(3)*Es(4)*K_3^2*K_4*V_4*W^2*k + 2*Es(3)*Es(4)*K_3^2*K_5*S*V_4*W^2*k + 2*Es(3)*Es(4)*K_3^2*V_4*W^3*k + Es(4)^2*K_3^2*V_4^2*W^2)^(1/2) + Es(1)*K_5*V_1*W - Es(4)*K_3*V_4*W - Es(3)*K_3*K_4*W*k + Es(3)*K_1*K_3*K_4*K_5*k + Es(3)*K_3*K_4*K_5*S*k + Es(3)*K_1*K_3*K_5*W*k + Es(3)*K_3*K_5*S*W*k)/(2*Es(3)*K_3*k*(K_4 + W));

            GC_W_dim2=@(S,g) -(K_3*g^2 + E_1*K_5*V_1*g + E_3*K_1*K_3*K_5*g*k + E_3*K_3*K_5*S*g*k - E_1*E_3*K_3*K_5*S*V_1*k)/(E_3*k*(K_3*g + E_1*K_5*V_1));
            FB_W_dim2=@(S,supplys) supplys(1) + supplys(2) - S;

        case 2 % G
            GC_W_dim2=@(S,g) 100*(S<((g*(K_2*K_3*g + 2*E_1*K_2*K_5*V_1 - K_1*K_3*K_5*g + 2*E_2*K_1*K_3*K_5*V_2))/(K_3*K_5*(g - 2*E_1*V_1)*(g - 2*E_2*V_2)))) +        (-0.0001)*(S>((g*(K_2*K_3*g + 2*E_1*K_2*K_5*V_1 - K_1*K_3*K_5*g + 2*E_2*K_1*K_3*K_5*V_2))/(K_3*K_5*(g - 2*E_1*V_1)*(g - 2*E_2*V_2)))) ;
            FB_W_dim2=@(S,supplys)supplys(2);
            steadyP_givenSW=@(S,W,d)-(Es(1)*K_2*K_5*V_1 - K_5^(1/2)*(K_5*Es(1)^2*K_2^2*V_1^2 + 2*K_5*Es(1)^2*K_2*K_3*S*V_1^2 + K_5*Es(1)^2*K_3^2*S^2*V_1^2 + 2*K_5*Es(1)*Es(2)*K_1*K_2*K_3*V_1*V_2 - 2*K_5*Es(1)*Es(2)*K_1*K_3^2*S*V_1*V_2 + 4*Es(1)*Es(2)*K_2*K_3^2*S*V_1*V_2 + 2*K_5*Es(1)*Es(2)*K_2*K_3*S*V_1*V_2 - 2*K_5*Es(1)*Es(2)*K_3^2*S^2*V_1*V_2 + K_5*Es(2)^2*K_1^2*K_3^2*V_2^2 + 2*K_5*Es(2)^2*K_1*K_3^2*S*V_2^2 + K_5*Es(2)^2*K_3^2*S^2*V_2^2)^(1/2) + Es(2)*K_1*K_3*K_5*V_2 - Es(1)*K_3*K_5*S*V_1 + Es(2)*K_3*K_5*S*V_2)/(2*(Es(1)*K_5*V_1 + Es(2)*K_3*V_2));

        case 3 % S2
            GC_W_dim2=@(S,g) (g*(K_2*K_6*g + E_4*K_2*K_7*V_4 - K_4*K_6*K_7*g + E_2*K_4*K_6*K_7*V_2))/(K_6*K_7*(g - E_2*V_2)*(g - E_4*V_4));
            FB_W_dim2=@(S,supplys)100*(S<supplys(1));
            steadyP_givenSW=@(S,W,d)-(E_4*K_2*K_7*V_4 - K_7^(1/2)*(K_7*E_2^2*K_4^2*K_6^2*V_2^2 + 2*K_7*E_2^2*K_4*K_6^2*V_2^2*W + K_7*E_2^2*K_6^2*V_2^2*W^2 + 2*K_7*E_2*E_4*K_2*K_4*K_6*V_2*V_4 + 4*E_2*E_4*K_2*K_6^2*V_2*V_4*W + 2*K_7*E_2*E_4*K_2*K_6*V_2*V_4*W - 2*K_7*E_2*E_4*K_4*K_6^2*V_2*V_4*W - 2*K_7*E_2*E_4*K_6^2*V_2*V_4*W^2 + K_7*E_4^2*K_2^2*V_4^2 + 2*K_7*E_4^2*K_2*K_6*V_4^2*W + K_7*E_4^2*K_6^2*V_4^2*W^2)^(1/2) + E_2*K_4*K_6*K_7*V_2 + E_2*K_6*K_7*V_2*W - E_4*K_6*K_7*V_4*W)/(2*(E_2*K_6*V_2 + E_4*K_7*V_4)); 

    end
    metabfuns.analytical_FBGCincell_dim2={FB_W_dim2,GC_W_dim2,steadyP_givenSW};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metabfuns.Intakes=Intakes;
metabfuns.Growth=Growth;
metabfuns.withincell_dxdt=dxdt;

metabfuns.Growth_givenE=Growth_givenE;
metabfuns.dxdt_givenE=dxdt_givenE;



 
metabfuns.intercell_dim=species_para.intercell_dim;
metabfuns.Edim=length(species_para.Es);

%% steady state given current enzyme allocation
% distinguish the species for S1,S2 and P

if speciesid>0 % there is analytical solution only when species belongs to one of them
        
    analytical_steadygivenenv=@(env)[G(env(1),env(2),steadyP_givenSW(env(1),env(2),nan),Es);
        steadyP_givenSW(env(1),env(2),nan)];
    
    metabfuns.analytical_steadygivenenv=analytical_steadygivenenv;
end

%% optimal enzyme allocation strategy
clear Es
%%%%%%%%%%%%%% given S and W
if speciesid>0 % there is analytical solution only when species belongs to one of them
    switch speciesid
        case 1 % S1
            optE1_givenSW=@(S,W) (K_3^2*W^2*k^2 + K_1^2*K_3^2*K_5^2*k^2 + K_3^2*K_5^2*S^2*k^2 + K_3*K_5*V_1*W*k + 2*K_1*K_3^2*K_5^2*S*k^2 - K_1*K_3*K_5^2*V_1*k - K_3*K_5^2*S*V_1*k - 2*K_3^2*K_5*S*V_1*k + 2*K_1*K_3^2*K_5*W*k^2 + 2*K_3^2*K_5*S*W*k^2 + K_3^(1/2)*K_5^(3/2)*V_1^(3/2)*k^(1/2)*(K_1*K_5 + K_3*S + K_5*S)^(1/2) - K_1*K_3^(3/2)*K_5^(3/2)*V_1^(1/2)*k^(3/2)*(K_1*K_5 + K_3*S + K_5*S)^(1/2) - K_3^(3/2)*K_5^(3/2)*S*V_1^(1/2)*k^(3/2)*(K_1*K_5 + K_3*S + K_5*S)^(1/2) - K_3^(3/2)*K_5^(1/2)*V_1^(1/2)*W*k^(3/2)*(K_1*K_5 + K_3*S + K_5*S)^(1/2))/(K_1^2*K_3^2*K_5^2*k^2 + 2*K_1*K_3^2*K_5^2*S*k^2 + 2*K_1*K_3^2*K_5*W*k^2 - 2*K_1*K_3*K_5^2*V_1*k + K_3^2*K_5^2*S^2*k^2 - 4*K_3^2*K_5*S*V_1*k + 2*K_3^2*K_5*S*W*k^2 + K_3^2*W^2*k^2 - 2*K_3*K_5^2*S*V_1*k + 2*K_3*K_5*V_1*W*k + K_5^2*V_1^2);
            maxg_givenSW=@(S,W) -(K_5*V_1*k*(W - K_3*S)*(K_5*V_1 + K_3*W*k + K_1*K_3*K_5*k + K_3*K_5*S*k - 2*K_3^(1/2)*K_5^(1/2)*V_1^(1/2)*k^(1/2)*(K_1*K_5 + K_3*S + K_5*S)^(1/2)))/(K_1^2*K_3^2*K_5^2*k^2 + 2*K_1*K_3^2*K_5^2*S*k^2 + 2*K_1*K_3^2*K_5*W*k^2 - 2*K_1*K_3*K_5^2*V_1*k + K_3^2*K_5^2*S^2*k^2 - 4*K_3^2*K_5*S*V_1*k + 2*K_3^2*K_5*S*W*k^2 + K_3^2*W^2*k^2 - 2*K_3*K_5^2*S*V_1*k + 2*K_3*K_5*V_1*W*k + K_5^2*V_1^2);
            optEs=@(S,W)[optE1_givenSW(S,W),0,1-optE1_givenSW(S,W),0]';
        case 2 % G
            optE1_givenSW=@(S,W) -(2*K_1^3*K_3^3*K_5^2*V_1*V_2^4 - 2*K_2^(1/2)*V_1^(1/2)*V_2^(1/2)*(K_2^2*K_5^2*V_1^3*V_2 + K_3^3*K_5*S^2*V_1*V_2^3 + K_1^2*K_3^2*K_5^2*V_1*V_2^3 + K_3^2*K_5^2*S^2*V_1*V_2^3 + K_3^3*K_5*S^2*V_1^2*V_2^2 + K_3^2*K_5^2*S^2*V_1^2*V_2^2 + K_1*K_3^3*K_5*S*V_1*V_2^3 + K_2*K_3*K_5^2*S*V_1^3*V_2 - 2*K_1*K_2*K_3*K_5^2*V_1^2*V_2^2 + 2*K_1*K_3^2*K_5^2*S*V_1*V_2^3 - 2*K_2*K_3*K_5^2*S*V_1^2*V_2^2 - 3*K_2*K_3^2*K_5*S*V_1^2*V_2^2 + K_1*K_3^2*K_5^2*S*V_1^2*V_2^2 + 2*K_2^(1/2)*K_3^(3/2)*K_5^(3/2)*S*V_1^(5/2)*V_2^(3/2)*(K_1*K_5 + K_3*S + K_5*S)^(1/2) - 2*K_2^(1/2)*K_3^(5/2)*K_5^(1/2)*S*V_1^(3/2)*V_2^(5/2)*(K_1*K_5 + K_3*S + K_5*S)^(1/2))^(1/2)*(K_1^3*K_3^3*K_5^2*V_1*V_2^3 - 2*K_1^2*K_2*K_3^2*K_5^2*V_1^2*V_2^2 + K_2*K_3^2*K_5^2*S^2*V_1^2*V_2^2 + K_1^2*K_3^3*K_5^2*S*V_1^2*V_2^2 + K_1*K_2^2*K_3*K_5^2*V_1^3*V_2 + K_2^2*K_3*K_5^2*S*V_1^3*V_2 + K_2^2*K_3^2*K_5*S*V_1^3*V_2 + K_2*K_3^3*K_5*S^2*V_1^3*V_2 + K_2*K_3^2*K_5^2*S^2*V_1^3*V_2 + K_2*K_3^3*K_5*S^2*V_1^2*V_2^2 + K_1^2*K_3^3*K_5^2*S*V_1*V_2^3 - 2*K_1*K_2*K_3^2*K_5^2*S*V_1^2*V_2^2 + K_1*K_2*K_3^2*K_5^2*S*V_1^3*V_2 - 3*K_1*K_2*K_3^3*K_5*S*V_1^2*V_2^2 - 2*K_2^(3/2)*K_3^(3/2)*K_5^(3/2)*S*V_1^(5/2)*V_2^(3/2)*(K_1*K_5 + K_3*S + K_5*S)^(1/2) - 2*K_2^(3/2)*K_3^(5/2)*K_5^(1/2)*S*V_1^(5/2)*V_2^(3/2)*(K_1*K_5 + K_3*S + K_5*S)^(1/2) + 2*K_1*K_2^(1/2)*K_3^(5/2)*K_5^(3/2)*S*V_1^(3/2)*V_2^(5/2)*(K_1*K_5 + K_3*S + K_5*S)^(1/2) + 2*K_1*K_2^(1/2)*K_3^(5/2)*K_5^(3/2)*S*V_1^(5/2)*V_2^(3/2)*(K_1*K_5 + K_3*S + K_5*S)^(1/2))^(1/2) - 4*K_1^2*K_2*K_3^2*K_5^2*V_1^2*V_2^3 + 2*K_1*K_3^3*K_5^2*S^2*V_1^2*V_2^3 + 2*K_1^2*K_3^3*K_5^2*S*V_1^2*V_2^3 + 2*K_1*K_2^2*K_3*K_5^2*V_1^3*V_2^2 + 2*K_2^2*K_3^2*K_5*S*V_1^3*V_2^2 + 2*K_1*K_3^3*K_5^2*S^2*V_1*V_2^4 + 2*K_2*K_3^3*K_5*S^2*V_1^2*V_2^3 + 2*K_2*K_3^3*K_5*S^2*V_1^3*V_2^2 + 4*K_1^2*K_3^3*K_5^2*S*V_1*V_2^4 - 4*K_1*K_2*K_3^2*K_5^2*S*V_1^2*V_2^3 + 2*K_1*K_2*K_3^2*K_5^2*S*V_1^3*V_2^2 - 6*K_1*K_2*K_3^3*K_5*S*V_1^2*V_2^3 - 4*K_2^(3/2)*K_3^(5/2)*K_5^(1/2)*S*V_1^(5/2)*V_2^(5/2)*(K_1*K_5 + K_3*S + K_5*S)^(1/2) + 4*K_1*K_2^(1/2)*K_3^(5/2)*K_5^(3/2)*S*V_1^(5/2)*V_2^(5/2)*(K_1*K_5 + K_3*S + K_5*S)^(1/2))/(2*K_5*V_1*V_2*(K_2*V_1 - K_1*K_3*V_2)*(K_5*K_1^2*K_3^2*V_2^2 - 2*K_5*K_1*K_2*K_3*V_1*V_2 + 2*K_5*K_1*K_3^2*S*V_1*V_2 + 2*K_5*K_1*K_3^2*S*V_2^2 + K_5*K_2^2*V_1^2 - 4*K_2*K_3^2*S*V_1*V_2 + 2*K_5*K_2*K_3*S*V_1^2 - 2*K_5*K_2*K_3*S*V_1*V_2 + K_5*K_3^2*S^2*V_1^2 + 2*K_5*K_3^2*S^2*V_1*V_2 + K_5*K_3^2*S^2*V_2^2));
            maxg_givenSW=@(S,W) (2*S*(K_2*K_3*K_5*V_1^2*V_2 + K_1*K_3^2*K_5*V_1*V_2^2 + K_3^2*K_5*S*V_1*V_2^2 + K_3^2*K_5*S*V_1^2*V_2 - 2*K_2^(1/2)*K_3^(3/2)*K_5^(1/2)*V_1^(3/2)*V_2^(3/2)*(K_1*K_5 + K_3*S + K_5*S)^(1/2)))/(K_5*K_1^2*K_3^2*V_2^2 - 2*K_5*K_1*K_2*K_3*V_1*V_2 + 2*K_5*K_1*K_3^2*S*V_1*V_2 + 2*K_5*K_1*K_3^2*S*V_2^2 + K_5*K_2^2*V_1^2 - 4*K_2*K_3^2*S*V_1*V_2 + 2*K_5*K_2*K_3*S*V_1^2 - 2*K_5*K_2*K_3*S*V_1*V_2 + K_5*K_3^2*S^2*V_1^2 + 2*K_5*K_3^2*S^2*V_1*V_2 + K_5*K_3^2*S^2*V_2^2);
            optEs=@(S,W)[optE1_givenSW(S,W),1-optE1_givenSW(S,W),0,0]';
        case 3 % S2  
            
            maxg_givenSW=@(S,W)(W*(K_2*K_6*K_7*V_2*V_4^2 + K_4*K_6^2*K_7*V_2^2*V_4 + K_6^2*K_7*V_2*V_4^2*W + K_6^2*K_7*V_2^2*V_4*W - 2*K_2^(1/2)*K_6^(3/2)*K_7^(1/2)*V_2^(3/2)*V_4^(3/2)*(K_6*W + K_7*W + K_4*K_7)^(1/2)))/(K_7*K_2^2*V_4^2 - 2*K_7*K_2*K_4*K_6*V_2*V_4 - 4*K_2*K_6^2*V_2*V_4*W - 2*K_7*K_2*K_6*V_2*V_4*W + 2*K_7*K_2*K_6*V_4^2*W + K_7*K_4^2*K_6^2*V_2^2 + 2*K_7*K_4*K_6^2*V_2^2*W + 2*K_7*K_4*K_6^2*V_2*V_4*W + K_7*K_6^2*V_2^2*W^2 + 2*K_7*K_6^2*V_2*V_4*W^2 + K_7*K_6^2*V_4^2*W^2);
            optE4_givenSW=@(S,W)-(K_4^3*K_6^3*K_7^2*V_2^4*V_4 - K_2^(1/2)*V_2^(1/2)*V_4^(1/2)*(K_2^2*K_7^2*V_2*V_4^3 + K_6^3*K_7*V_2^3*V_4*W^2 + K_4^2*K_6^2*K_7^2*V_2^3*V_4 + K_6^2*K_7^2*V_2^3*V_4*W^2 + K_6^3*K_7*V_2^2*V_4^2*W^2 + K_6^2*K_7^2*V_2^2*V_4^2*W^2 + K_2*K_6*K_7^2*V_2*V_4^3*W + K_4*K_6^3*K_7*V_2^3*V_4*W - 2*K_2*K_4*K_6*K_7^2*V_2^2*V_4^2 - 2*K_2*K_6*K_7^2*V_2^2*V_4^2*W - 3*K_2*K_6^2*K_7*V_2^2*V_4^2*W + 2*K_4*K_6^2*K_7^2*V_2^3*V_4*W + K_4*K_6^2*K_7^2*V_2^2*V_4^2*W + 2*K_2^(1/2)*K_6^(3/2)*K_7^(3/2)*V_2^(3/2)*V_4^(5/2)*W*(K_6*W + K_7*W + K_4*K_7)^(1/2) - 2*K_2^(1/2)*K_6^(5/2)*K_7^(1/2)*V_2^(5/2)*V_4^(3/2)*W*(K_6*W + K_7*W + K_4*K_7)^(1/2))^(1/2)*(K_4^3*K_6^3*K_7^2*V_2^3*V_4 - 2*K_2*K_4^2*K_6^2*K_7^2*V_2^2*V_4^2 + K_2*K_6^2*K_7^2*V_2^2*V_4^2*W^2 + K_4^2*K_6^3*K_7^2*V_2^2*V_4^2*W + K_2^2*K_4*K_6*K_7^2*V_2*V_4^3 + K_2^2*K_6*K_7^2*V_2*V_4^3*W + K_2^2*K_6^2*K_7*V_2*V_4^3*W + K_2*K_6^3*K_7*V_2*V_4^3*W^2 + K_2*K_6^2*K_7^2*V_2*V_4^3*W^2 + K_2*K_6^3*K_7*V_2^2*V_4^2*W^2 + K_4^2*K_6^3*K_7^2*V_2^3*V_4*W - 2*K_2*K_4*K_6^2*K_7^2*V_2^2*V_4^2*W + K_2*K_4*K_6^2*K_7^2*V_2*V_4^3*W - 3*K_2*K_4*K_6^3*K_7*V_2^2*V_4^2*W - 2*K_2^(3/2)*K_6^(3/2)*K_7^(3/2)*V_2^(3/2)*V_4^(5/2)*W*(K_6*W + K_7*W + K_4*K_7)^(1/2) - 2*K_2^(3/2)*K_6^(5/2)*K_7^(1/2)*V_2^(3/2)*V_4^(5/2)*W*(K_6*W + K_7*W + K_4*K_7)^(1/2) + 2*K_2^(1/2)*K_4*K_6^(5/2)*K_7^(3/2)*V_2^(3/2)*V_4^(5/2)*W*(K_6*W + K_7*W + K_4*K_7)^(1/2) + 2*K_2^(1/2)*K_4*K_6^(5/2)*K_7^(3/2)*V_2^(5/2)*V_4^(3/2)*W*(K_6*W + K_7*W + K_4*K_7)^(1/2))^(1/2) - 2*K_2*K_4^2*K_6^2*K_7^2*V_2^3*V_4^2 + K_4*K_6^3*K_7^2*V_2^3*V_4^2*W^2 + K_4^2*K_6^3*K_7^2*V_2^3*V_4^2*W + K_2^2*K_4*K_6*K_7^2*V_2^2*V_4^3 + K_2^2*K_6^2*K_7*V_2^2*V_4^3*W + K_2*K_6^3*K_7*V_2^2*V_4^3*W^2 + K_2*K_6^3*K_7*V_2^3*V_4^2*W^2 + K_4*K_6^3*K_7^2*V_2^4*V_4*W^2 + 2*K_4^2*K_6^3*K_7^2*V_2^4*V_4*W + K_2*K_4*K_6^2*K_7^2*V_2^2*V_4^3*W - 2*K_2*K_4*K_6^2*K_7^2*V_2^3*V_4^2*W - 3*K_2*K_4*K_6^3*K_7*V_2^3*V_4^2*W - 2*K_2^(3/2)*K_6^(5/2)*K_7^(1/2)*V_2^(5/2)*V_4^(5/2)*W*(K_6*W + K_7*W + K_4*K_7)^(1/2) + 2*K_2^(1/2)*K_4*K_6^(5/2)*K_7^(3/2)*V_2^(5/2)*V_4^(5/2)*W*(K_6*W + K_7*W + K_4*K_7)^(1/2))/(K_7*V_2*V_4*(K_2*V_4 - K_4*K_6*V_2)*(K_7*K_2^2*V_4^2 - 2*K_7*K_2*K_4*K_6*V_2*V_4 - 4*K_2*K_6^2*V_2*V_4*W - 2*K_7*K_2*K_6*V_2*V_4*W + 2*K_7*K_2*K_6*V_4^2*W + K_7*K_4^2*K_6^2*V_2^2 + 2*K_7*K_4*K_6^2*V_2^2*W + 2*K_7*K_4*K_6^2*V_2*V_4*W + K_7*K_6^2*V_2^2*W^2 + 2*K_7*K_6^2*V_2*V_4*W^2 + K_7*K_6^2*V_4^2*W^2));
            optEs=@(S,W)[0,1-optE4_givenSW(S,W),0,optE4_givenSW(S,W)]';
    end
    analytical_opt_givenenv=@(env)[maxg_givenSW(env(1),env(2));
        optEs(env(1),env(2))];
    metabfuns.analytical_opt_givenenv=analytical_opt_givenenv;
end


% %%%%%%%%%%%%%% given D
if speciesid>0 % there is analytical solution only when species belongs to one of them
    switch speciesid
    case 1 % S1
        optE1_S1_givenDW=@(D,W)(K_3*V_1*D - V_1^(1/2)*(K_5*V_1 + K_3*D)^(1/2)*(K_3*V_1*D + K_5*V_1*D + K_3*V_1*W*k + K_5*V_1*W*k - K_3*W*D*k - K_5*W*D*k + K_1*K_3*K_5*V_1*k - K_1*K_3*K_5*D*k)^(1/2) + K_3*V_1*W*k + K_5*V_1*W*k + K_1*K_3*K_5*V_1*k)/(V_1*(K_3*W*k - K_5*V_1 + K_5*W*k + K_1*K_3*K_5*k));
        minS_S1_givenDW= @(D,W)-((W - D/(k*((K_3*V_1*D - V_1^(1/2)*(K_5*V_1 + K_3*D)^(1/2)*(K_3*V_1*D + K_5*V_1*D + K_3*V_1*W*k + K_5*V_1*W*k - K_3*W*D*k - K_5*W*D*k + K_1*K_3*K_5*V_1*k - K_1*K_3*K_5*D*k)^(1/2) + K_3*V_1*W*k + K_5*V_1*W*k + K_1*K_3*K_5*V_1*k)/(V_1*(K_3*W*k - K_5*V_1 + K_5*W*k + K_1*K_3*K_5*k)) - 1)))/K_3 + (D*(K_1 + (W - D/(k*((K_3*V_1*D - V_1^(1/2)*(K_5*V_1 + K_3*D)^(1/2)*(K_3*V_1*D + K_5*V_1*D + K_3*V_1*W*k + K_5*V_1*W*k - K_3*W*D*k - K_5*W*D*k + K_1*K_3*K_5*V_1*k - K_1*K_3*K_5*D*k)^(1/2) + K_3*V_1*W*k + K_5*V_1*W*k + K_1*K_3*K_5*V_1*k)/(V_1*(K_3*W*k - K_5*V_1 + K_5*W*k + K_1*K_3*K_5*k)) - 1)))/K_5)*(K_3*W*k - K_5*V_1 + K_5*W*k + K_1*K_3*K_5*k))/(K_3*V_1*D - V_1^(1/2)*(K_5*V_1 + K_3*D)^(1/2)*(K_3*V_1*D + K_5*V_1*D + K_3*V_1*W*k + K_5*V_1*W*k - K_3*W*D*k - K_5*W*D*k + K_1*K_3*K_5*V_1*k - K_1*K_3*K_5*D*k)^(1/2) + K_3*V_1*W*k + K_5*V_1*W*k + K_1*K_3*K_5*V_1*k))/((D*(K_3*W*k - K_5*V_1 + K_5*W*k + K_1*K_3*K_5*k))/(K_3*V_1*D - V_1^(1/2)*(K_5*V_1 + K_3*D)^(1/2)*(K_3*V_1*D + K_5*V_1*D + K_3*V_1*W*k + K_5*V_1*W*k - K_3*W*D*k - K_5*W*D*k + K_1*K_3*K_5*V_1*k - K_1*K_3*K_5*D*k)^(1/2) + K_3*V_1*W*k + K_5*V_1*W*k + K_1*K_3*K_5*V_1*k) - 1);
        
        optEs=@(D,W)[optE1_S1_givenDW(D,W),0,1-optE1_S1_givenDW(D,W),0 ]';
        minSW=@(D,W)minS_S1_givenDW(D,W);
    case 2 % G
        optE1_G_givenD=@(D)-(2*K_1*K_3*K_5*V_1*V_2^2 - K_2^(1/2)*V_1^(1/2)*V_2^(1/2)*(D*K_3*V_2 - D*K_5*V_1 + 2*K_5*V_1*V_2)^(1/2)*(D*K_2*K_3*V_1 + D*K_2*K_5*V_1 - D*K_1*K_3*K_5*V_1 - D*K_1*K_3*K_5*V_2 + 2*K_1*K_3*K_5*V_1*V_2)^(1/2) + D*K_2*K_3*V_1*V_2 - D*K_1*K_3*K_5*V_1*V_2)/(2*K_5*V_1*V_2*(K_2*V_1 - K_1*K_3*V_2));
        minS_G_givenD=@(D)-(D*K_2^(1/2)*(K_2*V_1 - K_1*K_3*V_2)*(D*K_3*V_2 - D*K_5*V_1 + 2*K_5*V_1*V_2)^(1/2)*(D*K_2*K_3*V_1 + D*K_2*K_5*V_1 - D*K_1*K_3*K_5*V_1 - D*K_1*K_3*K_5*V_2 + 2*K_1*K_3*K_5*V_1*V_2)^(1/2))/(K_3*V_2^(1/2)*(D + (2*K_1*K_3*K_5*V_1*V_2^2 - K_2^(1/2)*V_1^(1/2)*V_2^(1/2)*(D*K_3*V_2 - D*K_5*V_1 + 2*K_5*V_1*V_2)^(1/2)*(D*K_2*K_3*V_1 + D*K_2*K_5*V_1 - D*K_1*K_3*K_5*V_1 - D*K_1*K_3*K_5*V_2 + 2*K_1*K_3*K_5*V_1*V_2)^(1/2) + D*K_2*K_3*V_1*V_2 - D*K_1*K_3*K_5*V_1*V_2)/(K_5*V_2*(K_2*V_1 - K_1*K_3*V_2)))*(2*K_2*K_5*V_1^(3/2)*V_2 - K_2^(1/2)*V_2^(1/2)*(D*K_3*V_2 - D*K_5*V_1 + 2*K_5*V_1*V_2)^(1/2)*(D*K_2*K_3*V_1 + D*K_2*K_5*V_1 - D*K_1*K_3*K_5*V_1 - D*K_1*K_3*K_5*V_2 + 2*K_1*K_3*K_5*V_1*V_2)^(1/2) - D*K_2*K_5*V_1^(3/2) + D*K_2*K_3*V_1^(1/2)*V_2));
        optEs=@(D,W)[optE1_G_givenD(D),1-optE1_G_givenD(D),0,0]';
        minSW=@(D)minS_G_givenD(D);
    case 3 % S2  
        optE4_givenD=@(g)-(K_2*K_6*V_2*V_4*g + K_4*K_6*K_7*V_2^2*V_4 - K_2^(1/2)*V_2^(1/2)*V_4^(1/2)*(K_7*V_2*V_4 + K_6*V_2*g - K_7*V_4*g)^(1/2)*(K_2*K_6*V_4*g + K_2*K_7*V_4*g + K_4*K_6*K_7*V_2*V_4 - K_4*K_6*K_7*V_2*g - K_4*K_6*K_7*V_4*g)^(1/2) - K_4*K_6*K_7*V_2*V_4*g)/(K_7*V_2*V_4*(K_2*V_4 - K_4*K_6*V_2));
 
        minW_S2_givenD=@(g)-(K_2^(1/2)*g*(K_2*V_4 - K_4*K_6*V_2)*(K_7*V_2*V_4 + K_6*V_2*g - K_7*V_4*g)^(1/2)*(K_2*K_6*V_4*g + K_2*K_7*V_4*g + K_4*K_6*K_7*V_2*V_4 - K_4*K_6*K_7*V_2*g - K_4*K_6*K_7*V_4*g)^(1/2))/(K_6*V_2^(1/2)*(g + (K_2*K_6*V_2*V_4*g + K_4*K_6*K_7*V_2^2*V_4 - K_2^(1/2)*V_2^(1/2)*V_4^(1/2)*(K_7*V_2*V_4 + K_6*V_2*g - K_7*V_4*g)^(1/2)*(K_2*K_6*V_4*g + K_2*K_7*V_4*g + K_4*K_6*K_7*V_2*V_4 - K_4*K_6*K_7*V_2*g - K_4*K_6*K_7*V_4*g)^(1/2) - K_4*K_6*K_7*V_2*V_4*g)/(K_7*V_2*(K_2*V_4 - K_4*K_6*V_2)))*(K_2*K_7*V_2*V_4^(3/2) - K_2*K_7*V_4^(3/2)*g - K_2^(1/2)*V_2^(1/2)*(K_7*V_2*V_4 + K_6*V_2*g - K_7*V_4*g)^(1/2)*(K_2*K_6*V_4*g + K_2*K_7*V_4*g + K_4*K_6*K_7*V_2*V_4 - K_4*K_6*K_7*V_2*g - K_4*K_6*K_7*V_4*g)^(1/2) + K_2*K_6*V_2*V_4^(1/2)*g));

        optEs=@(D,W)[0,1-optE4_givenD(D),0,optE4_givenD(D)];
        minSW=@(D)minW_S2_givenD(D);
    end
    
    
    metabfuns.analytical_optE_givenDW = optEs;
    metabfuns.analytical_minSW_givenDW= minSW;
    
end
 
metabfuns.optS1_E=@(Ss,g)(K_3^2*K_5*V_1*g - (K_3*K_5*V_1*(K_5*V_1 + K_3*g)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k))^(1/2) + K_3*K_5*V_1*g + K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k + K_1*K_3^2*K_5^2*V_1*k)/(K_1*K_3*K_5^2*V_1*k - K_3*K_5^2*V_1^2 - K_5^2*V_1^2 + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k + K_1*K_3^2*K_5^2*V_1*k);
metabfuns.optS1_S=@(Ss,g)(K_5*V_1*((K_1*K_3*g*k*(K_5^2*V_1^2 + K_3*K_5^2*V_1^2 + K_3^2*K_5*V_1*g + K_3*K_5*V_1*g - K_3^(1/2)*K_5^(1/2)*V_1^(1/2)*(K_5*V_1 + K_3*g)^(1/2)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k)^(1/2)))/(V_1*(K_3^2*Ss*k - K_3*K_5*V_1 - K_5*V_1 + K_1*K_3^2*K_5*k + K_1*K_3*K_5*k + K_3*K_5*Ss*k)) - (g*(K_3^2*K_5*V_1*g + K_3*K_5*V_1*g + K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k + K_1*K_3^2*K_5^2*V_1*k - K_3^(1/2)*K_5^(1/2)*V_1^(1/2)*(K_5*V_1 + K_3*g)^(1/2)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k)^(1/2)))/(K_3^2*Ss*k - K_3*K_5*V_1 - K_5*V_1 + K_1*K_3^2*K_5*k + K_1*K_3*K_5*k + K_3*K_5*Ss*k) - K_3*g^2 + (Ss*k*(K_5^2*V_1^2 + K_3*K_5^2*V_1^2 + K_3^2*K_5*V_1*g + K_3*K_5*V_1*g - K_3^(1/2)*K_5^(1/2)*V_1^(1/2)*(K_5*V_1 + K_3*g)^(1/2)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k)^(1/2))*(K_3^2*K_5*V_1*g + K_3*K_5*V_1*g + K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k + K_1*K_3^2*K_5^2*V_1*k - K_3^(1/2)*K_5^(1/2)*V_1^(1/2)*(K_5*V_1 + K_3*g)^(1/2)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k)^(1/2)))/(K_5*V_1*(K_3^2*Ss*k - K_3*K_5*V_1 - K_5*V_1 + K_1*K_3^2*K_5*k + K_1*K_3*K_5*k + K_3*K_5*Ss*k)^2) + (K_3*Ss*g*k*(K_5^2*V_1^2 + K_3*K_5^2*V_1^2 + K_3^2*K_5*V_1*g + K_3*K_5*V_1*g - K_3^(1/2)*K_5^(1/2)*V_1^(1/2)*(K_5*V_1 + K_3*g)^(1/2)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k)^(1/2)))/(K_5*V_1*(K_3^2*Ss*k - K_3*K_5*V_1 - K_5*V_1 + K_1*K_3^2*K_5*k + K_1*K_3*K_5*k + K_3*K_5*Ss*k)))*(K_3^2*Ss*k - K_3*K_5*V_1 - K_5*V_1 + K_1*K_3^2*K_5*k + K_1*K_3*K_5*k + K_3*K_5*Ss*k)^2)/(k*(K_5^2*V_1^2 + K_3*K_5^2*V_1^2 + K_3^2*K_5*V_1*g + K_3*K_5*V_1*g - K_3^(1/2)*K_5^(1/2)*V_1^(1/2)*(K_5*V_1 + K_3*g)^(1/2)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k)^(1/2))*(K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_3^3*K_5*V_1*g + K_3^3*Ss*g*k + K_3^2*K_5^2*V_1*g + K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k + K_3^3*K_5*Ss*V_1*k + K_1*K_3^2*K_5*g*k + K_1*K_3^3*K_5*g*k + K_3^2*K_5*Ss*g*k - K_3^3*K_5*Ss*g*k + 2*K_1*K_3^2*K_5^2*V_1*k + K_1*K_3^3*K_5^2*V_1*k + K_3^2*K_5^2*Ss*V_1*k - K_1*K_3^2*K_5^2*g*k - K_1*K_3^3*K_5^2*g*k - K_3^2*K_5^2*Ss*g*k - K_3^(1/2)*K_5^(1/2)*V_1^(1/2)*(K_5*V_1 + K_3*g)^(1/2)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k)^(1/2) - K_3^(3/2)*K_5^(1/2)*V_1^(1/2)*(K_5*V_1 + K_3*g)^(1/2)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k)^(1/2)));
metabfuns.optS1_W=@(Ss,g)(K_5*V_1*(K_3*g^2 + (g*(K_3^2*K_5*V_1*g + K_3*K_5*V_1*g + K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k + K_1*K_3^2*K_5^2*V_1*k - K_3^(1/2)*K_5^(1/2)*V_1^(1/2)*(K_5*V_1 + K_3*g)^(1/2)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k)^(1/2)))/(K_3^2*Ss*k - K_3*K_5*V_1 - K_5*V_1 + K_1*K_3^2*K_5*k + K_1*K_3*K_5*k + K_3*K_5*Ss*k) - (K_1*K_3*g*k*(K_5^2*V_1^2 + K_3*K_5^2*V_1^2 + K_3^2*K_5*V_1*g + K_3*K_5*V_1*g - K_3^(1/2)*K_5^(1/2)*V_1^(1/2)*(K_5*V_1 + K_3*g)^(1/2)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k)^(1/2)))/(V_1*(K_3^2*Ss*k - K_3*K_5*V_1 - K_5*V_1 + K_1*K_3^2*K_5*k + K_1*K_3*K_5*k + K_3*K_5*Ss*k)) - (K_3*Ss*g*k*(K_5^2*V_1^2 + K_3*K_5^2*V_1^2 + K_3^2*K_5*V_1*g + K_3*K_5*V_1*g - K_3^(1/2)*K_5^(1/2)*V_1^(1/2)*(K_5*V_1 + K_3*g)^(1/2)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k)^(1/2)))/(V_1*(K_3^2*Ss*k - K_3*K_5*V_1 - K_5*V_1 + K_1*K_3^2*K_5*k + K_1*K_3*K_5*k + K_3*K_5*Ss*k)) + (K_3*Ss*k*(K_5^2*V_1^2 + K_3*K_5^2*V_1^2 + K_3^2*K_5*V_1*g + K_3*K_5*V_1*g - K_3^(1/2)*K_5^(1/2)*V_1^(1/2)*(K_5*V_1 + K_3*g)^(1/2)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k)^(1/2))*(K_3^2*K_5*V_1*g + K_3*K_5*V_1*g + K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k + K_1*K_3^2*K_5^2*V_1*k - K_3^(1/2)*K_5^(1/2)*V_1^(1/2)*(K_5*V_1 + K_3*g)^(1/2)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k)^(1/2)))/(K_5*V_1*(K_3^2*Ss*k - K_3*K_5*V_1 - K_5*V_1 + K_1*K_3^2*K_5*k + K_1*K_3*K_5*k + K_3*K_5*Ss*k)^2))*(K_3^2*Ss*k - K_3*K_5*V_1 - K_5*V_1 + K_1*K_3^2*K_5*k + K_1*K_3*K_5*k + K_3*K_5*Ss*k)^2)/(k*(K_5^2*V_1^2 + K_3*K_5^2*V_1^2 + K_3^2*K_5*V_1*g + K_3*K_5*V_1*g - K_3^(1/2)*K_5^(1/2)*V_1^(1/2)*(K_5*V_1 + K_3*g)^(1/2)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k)^(1/2))*(K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_3^3*K_5*V_1*g + K_3^3*Ss*g*k + K_3^2*K_5^2*V_1*g + K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k + K_3^3*K_5*Ss*V_1*k + K_1*K_3^2*K_5*g*k + K_1*K_3^3*K_5*g*k + K_3^2*K_5*Ss*g*k - K_3^3*K_5*Ss*g*k + 2*K_1*K_3^2*K_5^2*V_1*k + K_1*K_3^3*K_5^2*V_1*k + K_3^2*K_5^2*Ss*V_1*k - K_1*K_3^2*K_5^2*g*k - K_1*K_3^3*K_5^2*g*k - K_3^2*K_5^2*Ss*g*k - K_3^(1/2)*K_5^(1/2)*V_1^(1/2)*(K_5*V_1 + K_3*g)^(1/2)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k)^(1/2) - K_3^(3/2)*K_5^(1/2)*V_1^(1/2)*(K_5*V_1 + K_3*g)^(1/2)*(K_5^2*V_1*g + K_3*K_5^2*V_1*g + K_3^2*K_5*V_1*g + K_1*K_5^2*V_1*k + K_5^2*Ss*V_1*k + K_3^2*Ss*g*k + K_3*K_5*V_1*g + K_3*K_5*Ss*V_1*k + K_1*K_3*K_5*g*k + K_3*K_5*Ss*g*k + 2*K_1*K_3*K_5^2*V_1*k + K_3*K_5^2*Ss*V_1*k + K_3^2*K_5*Ss*V_1*k - K_1*K_3*K_5^2*g*k + K_1*K_3^2*K_5*g*k - K_3*K_5^2*Ss*g*k - K_3^2*K_5*Ss*g*k + K_1*K_3^2*K_5^2*V_1*k - K_1*K_3^2*K_5^2*g*k)^(1/2)));
