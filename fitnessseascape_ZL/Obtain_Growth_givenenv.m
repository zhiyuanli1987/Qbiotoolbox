function [growthrate,withincellsteady,fval]=Obtain_Growth_givenenv(env,metabfun,odepara)
% given enviornment, what's the value for the growth rate of one species
growthrate=nan;
withincellsteady=[];




env=reshape(env,length(env),1);
Growthfunn=metabfun.Growth;% @(cs,xs)
% Define: Within cell varibles:
withincell_dxdt=@(xs)metabfun.withincell_dxdt(env,xs); %@(cs,xs)[];

if metabfun.intercell_dim==0 % no within cell varibles
    growthrate=Growthfunn(env,nan);
    withincellsteady=nan;
    fval=0;
elseif isfield(metabfun,'analytical_steadygivenenv')
    Growth_givenenv=metabfun.analytical_steadygivenenv;
    gandx=Growth_givenenv(env);
    
    growthrate=gandx(1);
    withincellsteady=gandx(2:end);
    fval=withincell_dxdt(withincellsteady);
    
else % first need to get the within cell varible to steady state, then compute the growth rate
    % run ODE for approaching the steady state, then solve for steady state
    
    intercell_dim=metabfun.intercell_dim;
    
    if nargin > 2
        iniv=odepara.iniv;
        timespan_runode=odepara.timespan_runode;
        runodeornot=odepara.runodeornot;
        plotresult=odepara.plotresult;
    else % use default  value
        iniv=rand(intercell_dim,1);
        timespan_runode=[0,500];
        runodeornot=1;
        plotresult=0;
    end
    
    iniv=reshape(iniv,length(iniv),1);
    
    if runodeornot==1 
      
        intercell_ode=@(t,x)withincell_dxdt(x);
       
        [timelist,trajecties]=ode23(intercell_ode,timespan_runode,iniv);
        iniv=trajecties(end,:)';
        if plotresult
            plot(timelist,trajecties);
            set(gca,'fontsize',22)
            xlabel('Time')
            ylabel('Withincellvarible');
            legend(strsplit(num2str(1:intercell_dim)));
        end
    end
    
    options = optimset('Display','off');
    lb=zeros(intercell_dim,1);
    ub=inf*ones(intercell_dim,1);
    
    
    [withincellsteady,fval]=lsqnonlin(withincell_dxdt,iniv,lb,ub,options);
    growthrate=metabfun.Growth(env,withincellsteady);
    
end
