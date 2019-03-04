function dxdt=ODE_onechemostat(x,chemostat_para,species_metabfuns)
% get the changing rate for multiple species in the same environment
% x and dxdt should be the format of: env1,env2....M1,M2....withincellspecies....

species_num=length(species_metabfuns);
intercell_dim=species_metabfuns{1}.intercell_dim;

c_supplys=chemostat_para.c_supplys;
d=chemostat_para.d;
env_dim=length(c_supplys);


% environment is the same for all species

d_cs_dt=d*(c_supplys-x(1:env_dim)); % environments
d_m_dt=zeros(species_num,1);
withincell_dxdt=zeros(intercell_dim,species_num);

for s=1:species_num
    % for this species, the envionrment, m and within cell varibles
    metabfuns=species_metabfuns{s};
    
    cs=x(1:env_dim);
    m=x((env_dim)+s);
    
    xs=x((env_dim+species_num)+(s-1)*intercell_dim+(1:intercell_dim));
    

    d_cs_dt=d_cs_dt-m*metabfuns.Intakes(cs,xs); % intake by species
    
    
    d_m_dt(s)=m*(metabfuns.Growth(cs,xs)-d);% growth rate of this species
    
    if intercell_dim>0
        withincell_dxdt(:,s)=metabfuns.withincell_dxdt(cs,xs);
    end
end

dxdt=[d_cs_dt;d_m_dt;reshape(withincell_dxdt,intercell_dim*species_num,1)];
