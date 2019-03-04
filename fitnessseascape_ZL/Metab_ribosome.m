function metabfuns=Metab_ribosome(species_para)

metabfuns=[];

m=species_para.m; % protein concentration
% import
K_c=species_para.K_c;
K_p=species_para.K_p;
V_c=species_para.V_c;
V_p=species_para.V_p;

% growth
Gamma=species_para.Gamma;
K_a=species_para.K_a;

% ribosome synthesis
synRform=species_para.synRform;% 1-3
k=species_para.k;
K_rp=species_para.K_rp;
K_rr=species_para.K_rr;
w_r=species_para.w_r;
w_p=species_para.w_p;

% form of protein allocation 
fform=species_para.fform;
fs=species_para.Es;
f_c=fs(1);
f_pv=fs(2);
f_rv=fs(3);

% form of degradation
dform=species_para.dform;




% all varibles:[a,rp,rr,r]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% define functions, there are several ways
g=@(a,rp,rr,r)Gamma*r*a/(a+K_a);

jc=@(c)V_c*c/(c+K_c);
jp=@(p)V_p*p/(p+K_p);

if synRform==1
    synR=@(a,rp,rr,r)k*rp*rr;
elseif synRform==2
    synR=@(a,rp,rr,r)k*rp/(rp+K_rp)*rr/(rr+K_rr);
else % 3
    synR=@(a,rp,rr,r)k*rp/(rp+K_rp)*rr;
end
    
if fform==1 % static value
    f_p=@(p,a,rp,rr,r)f_pv;
    f_r=@(p,a,rp,rr,r)f_rv;
elseif fform==2 % regulate by the value of g and influx
    
    f_r=@(p,a,rp,rr,r)(1-f_c )*(jp(p)/w_r )/(g(a,rp,rr,r)/w_p +jp(p)/w_r );
    f_p=@(p,a,rp,rr,r)(1-f_c-f_r(p,a,rp,rr,r));

end

if dform==1 % all dilute
    da=@(a,rp,rr,r)g(a,rp,rr,r);
    drp=@(a,rp,rr,r)g(a,rp,rr,r);
    drr=@(a,rp,rr,r)g(a,rp,rr,r);
    dr=@(a,rp,rr,r)g(a,rp,rr,r);
elseif dform==2
    % no degradation for some of them
    da=@(a,rp,rr,r)g(a,rp,rr,r);
    drp=@(a,rp,rr,r)0;
    drr=@(a,rp,rr,r)0;
    dr=@(a,rp,rr,r)g(a,rp,rr,r);
elseif dform==3
    da=@(a,rp,rr,r)0.01;
    drp=@(a,rp,rr,r)0.01;
    drr=@(a,rp,rr,r)0.01;
    dr=@(a,rp,rr,r)g(a,rp,rr,r);
end

Intake_c=@(c, a,rp,rr,r)f_c*jc(c);
Intake_p=@(p, a,rp,rr,r)f_p(p,a,rp,rr,r)*jp(p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dadt=  @(c,p, a,rp,rr,r)m*Intake_c(c, a,rp,rr,r)-m*g(a,rp,rr,r)-da(a,rp,rr,r)*a;
drpdt= @(c,p, a,rp,rr,r)m*f_r(p,a,rp,rr,r)*g(a,rp,rr,r)-m*w_p*synR(a,rp,rr,r)-drp(a,rp,rr,r)*rp;
drrdt= @(c,p, a,rp,rr,r)m*Intake_p(p, a,rp,rr,r)-m*w_r*synR(a,rp,rr,r)-drr(a,rp,rr,r)*rr;
drdt=  @(c,p, a,rp,rr,r)m*synR(a,rp,rr,r)-dr(a,rp,rr,r)*r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define: Intake rates as a vector
Intakes=@(cs,xs)[Intake_c(cs(1), xs(1),xs(2),xs(3),xs(4));
    Intake_p(cs(2), xs(1),xs(2),xs(3),xs(4))];

% Define: Growth  rate:
Growth=@(cs,xs)g(xs(1),xs(2),xs(3),xs(4));
dxdt=@(cs,xs)[dadt(cs(1),cs(2),xs(1),xs(2),xs(3),xs(4));
    drpdt(cs(1),cs(2),xs(1),xs(2),xs(3),xs(4));
    drrdt(cs(1),cs(2),xs(1),xs(2),xs(3),xs(4));
    drdt(cs(1),cs(2),xs(1),xs(2),xs(3),xs(4))];

%%%%%%%%%%%%%% OPTIMIZATION %%%%%%%%%%%%%%%%%%%
analytical_opt_Es_g_givenenv=[];
analytical_opt_Es_env_givenD=[];
threeF_opt_Es_g_givenenv=[];
threeF_opt_Es_env_givenD=[];
if dform==2 && fform==2 % artificial balance of ribosome synthesis, and some degradation for a
    % use the three-way optimization method
    % given environment, what's the optimal allocation of f
    F1=@(f_c,g,a,j_c,j_p)Gamma*(1-f_c )*(j_p/w_r *1/w_p )/(g/w_p +j_p/w_r )*m*a/(a+K_a )-g;
    F2=@(f_c,g,a,j_c,j_p)m*f_c*j_c-m*g-g*a;
    F3=@(f_c,g,a,j_c,j_p)(a*(K_a+ a))/(K_a*(f_c- 1) )+(j_c*m)/g; % F1_a/F1_f-F2_a/F2_f
    
    threeF_opt_Es_g_givenenv.threewayf=@(x,env)[F1(x(1),x(2),x(3),jc(env(1)),jp(env(2)));
        F2(x(1),x(2),x(3),jc(env(1)),jp(env(2)));
        F3(x(1),x(2),x(3),jc(env(1)),jp(env(2)))];
   threeF_opt_Es_g_givenenv.Edim=1; % which E (or f) is optimized
   threeF_opt_Es_g_givenenv.xdim=1;% which within cell varible appear
   
   
   threeF_opt_Es_env_givenD.threewayf=@(x,D,env1)[F1(x(1),D,x(2),jc(env1),x(3));
        F2(x(1),D,x(2),jc(env1),x(3));
        F3(x(1),D,x(2),jc(env1),x(3))];
   threeF_opt_Es_env_givenD.Edim=1; % which E (or f) is optimized
   threeF_opt_Es_env_givenD.xdim=1;% which within cell varible appear
   
   Es_given_fcgjp=@(f_c,g,jp)[f_c;
       (1-f_c )*(g/w_p) /(g/w_p +jp/w_r); % f_p
       (1-f_c )*(jp/w_r)/(g/w_p +jp/w_r);% f_r
    ];
   
   Es_givensolution=@(solution,env)Es_given_fcgjp(solution(1),solution(2),jp(env(2)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metabfuns.Intakes=Intakes;
metabfuns.Growth=Growth;
metabfuns.withincell_dxdt=dxdt;

metabfuns.analytical_FBGCincell_dim2={[]
    [];
    [];
    };
 
metabfuns.intercell_dim=species_para.intercell_dim;
metabfuns.Edim=length(species_para.Es);

metabfuns.analytical_opt_Es_g_givenenv=analytical_opt_Es_g_givenenv;
metabfuns.analytical_opt_Es_env_givenD=analytical_opt_Es_env_givenD;

metabfuns.threeF_opt_Es_g_givenenv=threeF_opt_Es_g_givenenv;
metabfuns.threeF_opt_Es_env_givenD=threeF_opt_Es_env_givenD;
metabfuns.threeF_opt_Es_g_givenenv.Es_givensolution=Es_givensolution;
