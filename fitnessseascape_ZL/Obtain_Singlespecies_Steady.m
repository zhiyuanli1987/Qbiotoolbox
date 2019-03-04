function [steadyv,fval]=Obtain_Singlespecies_Steady(chemostat_para,ones_metabfun,iniv)
% in a chemostat, what's the steady state environment constructed by one
% species

c_supplys=chemostat_para.c_supplys;
env_dim=length(c_supplys);
intercell_dim=ones_metabfun.intercell_dim;
varible_num=env_dim+1+intercell_dim;



dxdt=@(x)ODE_onechemostat(x,chemostat_para,{ones_metabfun});

% if there is no initial value, run the ODE for a while for steady states
if nargin <3
    odetouse=@(t,x)dxdt(x);
    timespan_runode=[0,300];
    randini=rand(varible_num,1);
    
   
    [timelist,trajecties]=ode23(odetouse,timespan_runode,randini);
    iniv=trajecties(end,:)';
end


options = optimset('Display','off');
lb=zeros(varible_num,1);
ub=inf*ones(varible_num,1);
ub(1:env_dim)=max(c_supplys);



[steadyv,fval]=lsqnonlin(dxdt,iniv,lb,ub,options);



