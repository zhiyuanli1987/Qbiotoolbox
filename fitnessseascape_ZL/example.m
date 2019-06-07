clear
clc

% this .m file provide example and explaination for using this set of
% programs

addpath(genpath('External_codes'))

%%
% explaination for programs:
% figuregenerating.m generates figures for the manuscript

% Metab_xx.m are functions define the import rate, growth rate and changing rate of intracelluar metabolites for a type of metabolic models, its inputs are
% parameters(including metabolic strategy) for this type of metabolic model, and output "metabfuns" that are specific functions for metabolic model given this parameter

% Nutrientspace_onespecies.m: given the metabolic model with certain
% parameter and metabolic strategy, what's the growth rate in the nutrient
% space

% Obtain_Growth_givenenv.m: given the metabolic model with certain
% parameter and metabolic strategy, and given the environment, what's the
% growth rate of the cell

% Nutrientspace_optspecies.m: given the metabolic model with certain
% parameter, what's the optimal metabolic strategy that gives the maximal
% growth rate in the nutrient space

% Obtain_optEG_givenenv.m: given the metabolic model with certain
% parameter, and given the environment, what's the optimal metabolic strategy that gives the maximal
% growth rate

% Analysis_OptM.m: given the results obtained from
% Nutrientspace_optspecies.m, obtain the maximal growth contour and the
% maximal flux balance curve

% Obtain_Singlespecies_Steady.m: given chemostat parameter, what's the steady state environment constructed by one
% species

% Obtain_supplyvector.m: given the steady state environment and the species
% parameter, what's the chemostat supplies that lead to this steady state
% environment

% ODE_onechemostat.m:obtain the ODEs for species dynamics in one chemostat

% ODE_multichemostat.m: obtain the ODEs for species dynamics in multiple
% linked chemostat with same chemostat parameter. Additional to the
% parameters for species and chemostat, one need to supply the system with
% a leakage rate and a matrix "chemostat_network" to describe the
% connectivities between chemostats

% Visualize_chemostatdynamic.m: Visualize the population dynamics, and visualize simulation dynamics in the nutrient space, 

% Visualize_chemostatsteady.m: Visualize the growth rate (optional) and GC for one
% species in the chemostat 



%% example of ultilizing this toolbox for investigating competition in a chemostat
% example actions to be performed
a1_generatenutrientspace=1;
a2_obtainGCFB=1;
a3_obtainsteadyandSL=1;
a4_competition=1;
a5_optimalstrategy=0;
%%%%%%%%%%%%%%%%%%%%
% step 1: choose a metabolic model
MetabModel=@Metab_MinImport;

% step 2: define a parameter set
minformul_species_para=[];
minformul_species_para.intercell_dim=0;
minformul_species_para.Ks=[0.5,0.5]';
minformul_species_para.Gamma=10;
minformul_species_para.r=100;

minformul_species_para.Es=[0.3;0.7]; % metabolic strategy. It is Alpha in the manuscript.

metabfun_min=MetabModel(minformul_species_para); % functions for one species generated by this set of parameter

% define one set of chemostat parameters
chemostat_para=[];
chemostat_para.c_supplys=[1,0.7]';
chemostat_para.d=1;


    ranges=[0,1;
        0,1]; % ranges in the nutrient space
% possible action 1: generate the nutirent space
if a1_generatenutrientspace

    pnumbydim=[20,20];% number of points to sample in each nutrient space
    showprogress=0;
    
    nspace_drawpara=[];
    nspace_drawpara.colortable=colormap('summer');
    nspace_drawpara.drawcontour=0;
    
    %figure;
    hold on
    [nutreint_boarders_m,nutreint_values_m,growthrates_m,consumptions_m,withincellv_m,ch]=Nutrientspace_onespecies(metabfun_min,ranges,pnumbydim,0,nspace_drawpara);
    title('Growth rates in nutrient space')
    
end

% possible action 2: obatain growth contour and flux balance curve 
% line
if a2_obtainGCFB
    
    %figure;
    hold on;
    % Visualize flux balance curve and growth contour
    steady_drawpara=[];
    steady_drawpara.pnum=100;% number of points for GC and FB
    steady_drawpara.growthcolor=[1,0,0];
    steady_drawpara.fluxcolor=[0,0,1];
    steady_drawpara.linewidth=1.5;
    Visualize_chemostatsteady(chemostat_para,metabfun_min,steady_drawpara);
end


% possible action 3: obtain steady state and supply line
if a3_obtainsteadyandSL
    
    % single species steady state 
    steadyv=Obtain_Singlespecies_Steady(chemostat_para,metabfun_min);
    env=steadyv(1:2);
    
    % supply line for this steady state
    [FBv_givenIa,onelargesupplypoint]=Obtain_supplyvector(metabfun_min,env);
    
%    figure;
    hold on;   
    scatter(env(1),env(2),100,[0,0,0],'filled');
    plot([env(1),onelargesupplypoint(1)],[env(2),onelargesupplypoint(2)],'-- k','linewidth',1);        
    axis(reshape(ranges',1,4))
    set(gca,'fontsize',24)
    xlabel('c_a');
    ylabel('c_b');
end

% possible action 4: For multi-species competition, draw GCs, get fitness landscape, run simulation
if a4_competition
    
    % define two species by their metabolic strategies
    Ea_byspecies=[0.3,0.6]; % The alpha_a for two species

    % assign each species a defined color
    species_num=length(Ea_byspecies);    
    species_colors=distinguishable_colors(2, {'w','k'});
 
    steady_drawpara=[];
    steady_drawpara.showfluxcurve=0;
    
    metabfuns_consortia=[]; % record the corresponding equations for each species in the consortia
    steadyenv_byspecies=[];
    figure; hold on;
    for s=1:species_num
        
        Es=[Ea_byspecies(s);1-Ea_byspecies(s)];
        this_species=minformul_species_para;
        this_species.Es=Es;
        
        metabfuns_consortia{s}=MetabModel(this_species);
        
        % GC and steady states for each species
        steady_drawpara.growthcolor=species_colors(s,:);
        Visualize_chemostatsteady(chemostat_para,metabfuns_consortia{s},steady_drawpara);
        [steadyv,fval]=Obtain_Singlespecies_Steady(chemostat_para,metabfuns_consortia{s});
        scatter(steadyv(1),steadyv(2),100,steady_drawpara.growthcolor,'filled');
   
        % record the species-specific steady state
        steadyenv_byspecies(s,:)=steadyv(1:2);
    end
    
    % simulate the population dynamics
    odepara=[];
    odepara.timespan_runode=[0,100];
    
    dynamic_drawpara=[];
    dynamic_drawpara.linewidths=[0.5,1];
    dynamic_drawpara.speciescolor=species_colors;
        
    Visualize_chemostatdynamic(chemostat_para,metabfuns_consortia,dynamic_drawpara,odepara);

    % obtain fitness landscape for the two steady environment
    pnum=100;
    strategy_list=linspace(0,1,pnum);

    figure;hold on;
    for e=1:size(steadyenv_byspecies,1)
        env=steadyenv_byspecies(e,:)';
        fitness_list=nan*strategy_list;
        for i=1:pnum
            tempspecies=minformul_species_para;
            tempspecies.Es=[strategy_list(i),1-strategy_list(i)]';
            tempmetabfun=MetabModel(tempspecies);
            fitness_list(i)=Obtain_Growth_givenenv(env,tempmetabfun);
        end
        % fitness landscape under this enviroment
        
        plot(strategy_list,fitness_list,'color',species_colors(e,:),'linewidth',1);
        if ~(e>2)
            scatter(Ea_byspecies(e),chemostat_para.d,150,species_colors(e,:),'d','filled')
        end
    end
    title('Fitness landscape')
    set(gca,'fontsize',22);
    xlabel('Strategies (Alpha_a)');
    ylabel('Fitness')
    
end

% possible action 5: obtain non-invasible (optimal) strategies/growth rates in the nutrient
% space given a supply condition
if a5_optimalstrategy
    
    pnumbydim=[30,30];
    
    [nutreint_boarders_m,nutreint_values_m,optgrowthrates_m,opt_Es,consumptions_m,withincellv_m]=Nutrientspace_optspecies(minformul_species_para,MetabModel,ranges,pnumbydim);
    resultsfromoptspecies={nutreint_boarders_m,nutreint_values_m,optgrowthrates_m,opt_Es,consumptions_m,withincellv_m};
    
    opt_drawpara=[];
    opt_drawpara.colornutrientspace=1;
    opt_drawpara.grwoth_colormap=colormap('bone');
    opt_drawpara.plotGC=3;
    opt_drawpara.plotGCcolor=colormap('jet');
    opt_drawpara.minmaxv_colors=[0.2,0.8];
    
    figure;hold on;
    [alongGCs,alongFBs]=Analysis_OptM(metabfun_min,resultsfromoptspecies,chemostat_para.d,[],opt_drawpara);
    shading interp
   
    % choose some species on the supply curve
    lw=2;
    steady_drawpara=[];
    steady_drawpara.pnum=1000;
    steady_drawpara.linewidth=lw;
    steady_drawpara.showfluxcurve=0;
    
    alongGC=alongGCs{1};
    teststrategies=[0.3,0.5,0.7];
    for i=1:length(teststrategies)
        [minv,minloc]=min(abs(alongGC.Es_optGC(1,:)-teststrategies(i)));
        env=alongGC.envs(:,minloc);
        Es=alongGC.Es_optGC(:,minloc);
        
        species_one=minformul_species_para;
        species_one.Es=Es;
        
        metabfun=MetabModel(species_one);
        
        colorid=alongGC.colorids(minloc);
        
        steady_drawpara.growthcolor=opt_drawpara.plotGCcolor(colorid,:);
        Visualize_chemostatsteady(chemostat_para,metabfun,steady_drawpara);
        
        % supply vector
        lw=1;
        [FBv_givenIa,onelargesupplypoint]=Obtain_supplyvector(metabfun,env);
        plot([env(1),onelargesupplypoint(1)],[env(2),onelargesupplypoint(2)],'-- k','linewidth',lw);
    end
    
    axis([0 0.5 0 0.5])

    
end