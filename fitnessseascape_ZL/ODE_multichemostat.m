function dxdt=ODE_multichemostat(x,chemostat_para,species_metabfuns,leakrate,chemostat_network)

% chemostat_network: i->j means transportation of materials from i to j
species_num=length(species_metabfuns);
intercell_dim=species_metabfuns{1}.intercell_dim;

c_supplys=chemostat_para.c_supplys;
env_dim=length(c_supplys);

chemostat_num=size(chemostat_network,1);
for i=1:chemostat_num % cannot leak to self
    chemostat_network(i,i)=0;
end
eachchemo_variblenum=env_dim+species_num+species_num*intercell_dim;
 

exchangedim=1:(env_dim+species_num);
% varibles are arranged in this way: all varibles in chemostat 1.... all varibles in chemostat 2...
x_varible_chemo=reshape(x,eachchemo_variblenum,chemostat_num);
 
%%%%%%%%%%%%%%%%
% % the changing rate for each chemostate if they were seperated
changerate_varible_chemo=0*x_varible_chemo;
for c=1:chemostat_num
     changerate_varible_chemo(:,c)=ODE_onechemostat(x_varible_chemo(:,c),chemostat_para,species_metabfuns);
end

% exchange between connected chemostat
[fromcs,tocs]=find(chemostat_network==1);
for i=1:length(fromcs)
    cid_i=fromcs(i);
    cid_j=tocs(i);
    
    % flow of nutirents and biomass from i to j
    leakageflux_i2j=x_varible_chemo(exchangedim,cid_i)*leakrate;
    changerate_varible_chemo(exchangedim,cid_i)=changerate_varible_chemo(exchangedim,cid_i)-leakageflux_i2j;
    changerate_varible_chemo(exchangedim,cid_j)=changerate_varible_chemo(exchangedim,cid_j)+leakageflux_i2j; 
end
 
dxdt=reshape(changerate_varible_chemo,eachchemo_variblenum*chemostat_num,1);