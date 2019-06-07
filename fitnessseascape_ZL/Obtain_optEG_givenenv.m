function [optEs,optg,fval]=Obtain_optEG_givenenv(env,metabfun,inis,nonzeroEdim)
% given environment, what's the optimal allocation of resouces and
% correspoinding growth
options = optimset('Display','off');

 E_dim=metabfun.Edim;

 if nargin > 3
     ;% get nonzeroEdim from inputs
 else
      nonzeroEdim=1:E_dim;
 end

if isfield(metabfun,'analytical_opt_givenenv')
    analytical_opt_givenenv=metabfun.analytical_opt_givenenv;
    gEs=analytical_opt_givenenv(env);
    optg=gEs(1);
    optEs=gEs(2:end);
    fval=0;
    
elseif isfield(metabfun,'threeF_opt_Es_g_givenenv')
    
    threeF_opt_Es_g_givenenv=metabfun.threeF_opt_Es_g_givenenv;
    threewayf=threeF_opt_Es_g_givenenv.threewayf;
    
    Es_givensolution=threeF_opt_Es_g_givenenv.Es_givensolution;
    tozerofun=@(x)threewayf(x,env);
    
    if nargin < 3 % no iniv
        iniv=rand(1,3);
    else
        iniEs=inis.iniEs;
        inig= inis.inig;
        inixs=inis.inixs;
        
        Edim=threeF_opt_Es_g_givenenv.Edim; % which E (or f) is optimized
        xdim=threeF_opt_Es_g_givenenv.xdim;% which within cell varible appear
        iniv=[iniEs(Edim),inig,inixs(xdim)];
    end

    % f_c,g,a,j_c,j_p
    % x is always in the format of [Ei,g,intercell j]
    lb=zeros(3,1);
    ub=inf*ones(3,1);
    ub(1)=1;
    
    [solution,fval]=lsqnonlin(tozerofun,iniv,lb,ub,options);
    optg=solution(2);
    optEs=Es_givensolution(solution,env);
    
    
    
    
else % no analytical solution
    % first is the dimension for intercell varibles
    intercell_dim=metabfun.intercell_dim;

    % nonzeroEdim
    realEnum=length(nonzeroEdim);
    convm=zeros(4,realEnum);
    for i=1:realEnum
        convm(nonzeroEdim(i),i)=1;
    end
    
    realEtoEs=@(realE)convm*realE;
    
    tominfun=@(x)-1*metabfun.Growth_givenE(env,x(1:intercell_dim),realEtoEs(x((intercell_dim+1):(realEnum+intercell_dim))));%(env,xs,Es)
    dxdt_givenE=@(x)metabfun.dxdt_givenE(env,x(1:intercell_dim),realEtoEs(x((intercell_dim+1):(realEnum+intercell_dim))));
    
    if nargin < 3 % no iniv
        iniv=rand(1,intercell_dim+realEnum);
    else
        iniEs=reshape(inis.iniEs(nonzeroEdim),1,realEnum);
        inixs=reshape(inis.inixs,1,intercell_dim);
        iniv=[inixs,iniEs];
    end
    
    
    lb=zeros(1,intercell_dim+realEnum);
    ub=inf*ones(1,intercell_dim+realEnum);
    ub((intercell_dim+1):(realEnum+intercell_dim))=1; % for enzyme
    A=[];
    b=[];
    Aeq=[zeros(1,intercell_dim),ones(1,realEnum)];
    beq=1;
   
    nonlcon = @(x)deal(0,dxdt_givenE(x));

    iniv=iniv';
    [xfinal,fval] = fmincon(tominfun,iniv,A,b,Aeq,beq,lb,ub,nonlcon,options);
   
    
    optrealEs=xfinal((intercell_dim+1):end);
    optEs=realEtoEs(optrealEs);
    optg=-1*fval;
    fval=max(abs(dxdt_givenE(xfinal)));
    
    
end

end

