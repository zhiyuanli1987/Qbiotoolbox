% this program is for generating figures for the manuscript
clear
clc

addpath(genpath('External_codes'))


defaultfs=22;


%%%%%% tasks to perform
F1DE_demonstration=0; % Figure 1 D-E: demonstrate growth contour and flux balance curve the nutrient space

F2BC_invasion_rule=0; % Figure 2 B-C: invasion rule
F2DE_mutualinvasion=0;% Figure 2 D-E; unlimited co-existence

F3_oscillation=0; %     Figure 3: oscillation

F4_multistable=0; %     Figure 4: multistability, chain of invasion and non-invasible strategy
FS2_spatial=0;  %       Figure S2: spantanious spatial gradient

F5_optinvasion=0;%      Figure 5: two mutually invading species constructed an non-invasible envionrment
F6_newdim=1;    %       Figure 6: species creating a new nutrient dimension

%%
%%%%%% define species parameters for different metabolic models
% necessary nutirents
min_species_para=[];
min_species_para.intercell_dim=0;
min_species_para.Ks=[0.7;1.3];
min_species_para.Es=[0.3;0.7];
min_species_para.Gamma=10;
min_species_para.r=100;

% substitutable nutirents
sum_species_para=[];
sum_species_para.intercell_dim=0;
sum_species_para.Ks=[1.2;0.8];
sum_species_para.Vs=[3;3];
sum_species_para.Es=[0.4;0.6];
sum_species_para.r=100;

%
species_iandc_subs=[];
species_iandc_subs.Ks=[0.5,0.5,0.5]';
species_iandc_subs.k=nan;
species_iandc_subs.Es=nan;
species_iandc_subs.intercell_dim=length(species_iandc_subs.Ks);
species_iandc_subs.r=100;
species_iandc_subs.V_import=1000;


%% Figure 1 D-E: demonstrate growth contour and flux balance curve the nutrient space
if F1DE_demonstration
    
    % basic definitions:
    MetabModel=@Metab_MinImport;
    
    % define chemostat parameter
    chemostat_para=[];
    chemostat_para.c_supplys=[1,0.7]';
    chemostat_para.d=1;
    
    csupply_a_list=[0.6,0.8,1];% the possible supply concentration in a-dimension
    dlist=[0.5,1,1.6]; % the possible growth rates to test
    
    % parameters for picture drawing
    nspace_drawpara=[];nspace_drawpara.colortable=colormap('summer');
    ranges=[0,1;
        0,1];
    pnumbydim=[50,50];
    pnum=200;
    
    GClw=1;
    GCcolor=[1,0.8,0;
        0.8,0,0;
        0.3,0,0];
    ms=200;
    
    supplyl_lw=0.8;
    supplyl_color=[0,0,0];
    
    FBlw=1;
    FB_colors=[       0.75, 0, 0.75;
        0, 0.75, 1
        0, 0, 1
        ];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generating plot
    metabfun=MetabModel(min_species_para); % choose a metabolic model
    
    % obtain cell growth in the nutreint space;
    figure;hold on
    [nutreint_boarders_m,nutreint_values_m,growthrates_m,consumptions_m,withincellv_m,ch]=Nutrientspace_onespecies(metabfun,ranges,pnumbydim,0,nspace_drawpara);
    set(ch,'ytick',[0.5,1,1.5])
    
    % the steady state by one supply condition
    steadyv=Obtain_Singlespecies_Steady(chemostat_para,metabfun);
    env=steadyv(1:2);
    
    % supply line for this steady state
    [FBv_givenIa,onelargesupplypoint]=Obtain_supplyvector(metabfun,env);
    
    %%%%%%%% obtain growth contour, flux balance curve, supply lines
    FB_cb_dim2=metabfun.analytical_FBGCincell_dim2{1};
    GC_cb_dim2=metabfun.analytical_FBGCincell_dim2{2};
    alist=linspace(ranges(1,1),ranges(1,2),pnum);
    
    % Growth contour
    GC_blist=nan*alist;
    for j=1:pnum
        GC_blist(j)=GC_cb_dim2(alist(j),chemostat_para.d);
    end
    goodlist=find(GC_blist>ranges(2,1) & GC_blist<ranges(2,2));
    
    plot(alist,GC_blist,'linewidth',GClw,'color',GCcolor(2,:));
    scatter(env(1),env(2),ms,GCcolor(2,:),'filled');
    plot([env(1),onelargesupplypoint(1)],[env(2),onelargesupplypoint(2)],':','linewidth',supplyl_lw,'color',supplyl_color);
    % supply line
    
    % FB for different supply conditions
    for i=1:length(csupply_a_list)
        c_supply= FBv_givenIa(csupply_a_list(i));
        
        % flux balance curve based on this supply condition
        alist=linspace(ranges(1,1),c_supply(1),pnum);
        FB_blist=nan*alist;
        for j=1:pnum
            FB_blist(j)=FB_cb_dim2(alist(j),c_supply);
        end
        goodlist=find(FB_blist>ranges(2,1) & FB_blist<c_supply(2));
        plot(alist(goodlist),FB_blist(goodlist),'linewidth',FBlw,'color',FB_colors(i,:));
        scatter(c_supply(1),c_supply(2),ms,FB_colors(i,:));
    end
    axis(reshape(ranges',1,4))
    set(gca,'fontsize',defaultfs,'xtick',[0:0.5:1],'ytick',[0,0.5,1])
    xlabel('');
    ylabel('');
    title('');
    axis square
    box on;
    
    % how change of growth rate may cause change of nutirent limitation
    figure; hold on;
    c_supply= chemostat_para.c_supplys;
    
    % GC for different d
    for d=1:length(dlist)
        GC_blist=nan*alist;
        for j=1:pnum
            GC_blist(j)=GC_cb_dim2(alist(j),dlist(d));
        end
        goodlist=find(GC_blist>ranges(2,1) & GC_blist<ranges(2,2));
        plot(alist,GC_blist,'linewidth',GClw,'color',GCcolor(d,:));
        
        chemostat_para.d=dlist(d);
        steadyv=Obtain_Singlespecies_Steady(chemostat_para,metabfun);
        scatter(steadyv(1),steadyv(2),ms,GCcolor(d,:),'filled');
    end
    
    % flux balance curve based on this supply condition
    alist=linspace(ranges(1,1),c_supply(1),pnum);
    FB_blist=nan*alist;
    for j=1:pnum
        FB_blist(j)=FB_cb_dim2(alist(j),c_supply);
    end
    goodlist=find(FB_blist>ranges(2,1) & FB_blist<c_supply(2));
    plot(alist(goodlist),FB_blist(goodlist),'linewidth',FBlw,'color',FB_colors(3,:));
    scatter(c_supply(1),c_supply(2),ms,FB_colors(3,:));
    
    axis(reshape(ranges',1,4))
    set(gca,'fontsize',defaultfs,'xtick',[0:0.5:1],'ytick',[0,0.5,1])
    xlabel('');
    ylabel('');
    title('');
    axis square
    box on;
    
end

%%
if F2BC_invasion_rule
    
    MetabModel=@Metab_substitutable;
    species_colors=distinguishable_colors(2, {'w','k'});
    
    % two species, the second is invader
    E1_lists=[0.2,0.6];
    supplylist=[0.5,1
        1,0.5;
        1,1];
    iniinvader_m_percent=0.2;% fraction of the invader in the initial condition
    
    % drawing parameters
    ranges=[0,1;0,1];
    pnumbydim=2*[100,100];
    
    steady_drawpara=[];
    steady_drawpara.pnum=100;
    steady_drawpara.ranges=ranges;
    steady_drawpara.showfluxcurve=0;
    steady_drawpara.linewidth=1;
    
    ms=200;
    
    metabfuns_twospecies=[];
    for s=1:length(E1_lists)
        Es=[E1_lists(s);1-E1_lists(s)];
        this_species=sum_species_para;
        this_species.Es=Es;
        
        metabfuns_twospecies{s}=MetabModel(this_species);
    end
    
    chemostat_para=[];
    chemostat_para.d=1;
    for i=1:size(supplylist)
        figure;
        chemostat_para.c_supplys=supplylist(i,:)';
        % background of the invader
        [nutreint_boarders_m,nutreint_values_m,growthrates_m,consumptions_m,withincellv_m,ch]=Nutrientspace_onespecies(metabfuns_twospecies{2},ranges,pnumbydim,0);
        tocolorm=growthrates_m;
        tocolorm(tocolorm<chemostat_para.d)=nan;
        
        colormap('summer')
        h=pcolor(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),tocolorm);
        h.LineStyle='none';
        ch=colorbar;
        caxis([min(growthrates_m(:)),max(growthrates_m(:))]);
        hold on;
        
        % growth contour for each species
        inisteady=[];
        for s=1:length(E1_lists)
            metabfun=metabfuns_twospecies{s};
            steady_drawpara.growthcolor=species_colors(s,:);
            
            Visualize_chemostatsteady(chemostat_para,metabfun,steady_drawpara);
            [steadyv,fval]=Obtain_Singlespecies_Steady(chemostat_para,metabfun);
            scatter(steadyv(1),steadyv(2),ms,steady_drawpara.growthcolor,'filled');
            
            if s==1
                inisteady=[steadyv;iniinvader_m_percent*steadyv(end)]';
            end
        end
        % draw the "invasion zone" of the second species
        axis square
        scatter(supplylist(i,1),supplylist(i,2),ms,[0,0,0])
        set(gca,'fontsize',defaultfs,'xtick',[0,0.5,1],'ytick',[0.5,1])
        set(ch,'ytick',[0,1])
        xlabel('')
        ylabel('')
        
        % run the invasion simulation
        odepara=[];
        odepara.timespan_runode=[0,100];
        odepara.inistates=inisteady;
        dynamic_drawpara=[];
        dynamic_drawpara.arrow_num=1;
        dynamic_drawpara.linewidths=[0.5,1];
        dynamic_drawpara.arrow_percent=1/30;
        dynamic_drawpara.speciescolor=species_colors;
        Visualize_chemostatdynamic(chemostat_para,metabfuns_twospecies,dynamic_drawpara,odepara);
        set(gca,'fontsize',40,'xtick',odepara.timespan_runode,'ytick',0:100:300)
        box on;
    end
    
end

%%
if F2DE_mutualinvasion
    
    MetabModel=@Metab_substitutable;
    
    % parameter for species
    E1_lists=[0.2,0.6,0.1,0.8,0.3,0.5,0.7,0.9];
    ranges=[0,1;0,1];
    ms=200;
    finalcolor=[1,0,1];
    
    chemostat_para=[];
    chemostat_para.d=1;
    chemostat_para.c_supplys=[1;1];
    
    
    species_num=length(E1_lists);
    species_colors=distinguishable_colors(species_num, {'w','k'});
    
    steady_drawpara=[];
    steady_drawpara.pnum=100;
    steady_drawpara.ranges=ranges;
    steady_drawpara.showfluxcurve=0;
    steady_drawpara.linewidth=1;
    
    % fitness landscape formed by each species, and the final stage
    for showid=1:2
        figure;hold on;
        
        if showid==1
            speciesid=setdiff(1:species_num,4);
        else
            speciesid=1:2;
        end
        
        % for each species present
        metabfun_consortia=[];
        steadyenv_byspecies=[];
        h=[];
        for s=1:length(speciesid)
            sid=speciesid(s);
            
            Es=[E1_lists(sid);1-E1_lists(sid)];
            this_species=sum_species_para;
            this_species.Es=Es;
            
            metabfun=MetabModel(this_species);
            metabfun_consortia{s}=metabfun;
            
            %GC
            steady_drawpara.growthcolor=species_colors(sid,:);
            [alongGC,h(s)]=Visualize_chemostatsteady(chemostat_para,metabfun,steady_drawpara);
            
            % steady state
            [steadyv,fval]=Obtain_Singlespecies_Steady(chemostat_para,metabfun);
            scatter(steadyv(1),steadyv(2),ms,steady_drawpara.growthcolor,'filled');
            
            steadyenv_byspecies(sid,:)=steadyv(1:2);
        end
        scatter(chemostat_para.c_supplys(1),chemostat_para.c_supplys(2),ms,[0,0,0])
        axis square;
        set(gca,'fontsize',defaultfs,'xtick',[0:0.5:1],'ytick',[0,0.5,1]);
        xlabel(''); ylabel('')
        %legend(h,strsplit(num2str(E1_lists(speciesid))))
        box on
        
        odepara=[];
        odepara.timespan_runode=[0,100];
        odepara.inistates=[0.2,0.6,200/length(speciesid)*ones(1,length(speciesid))];
        dynamic_drawpara=[];
        dynamic_drawpara.linewidths=[0.5,1];
        dynamic_drawpara.speciescolor=species_colors(speciesid,:);
        dynamic_drawpara.arrow_sefrac=[0.3,0.8];
        if showid==2
            scatter(finalstate(1),finalstate(2),ms,finalcolor,'filled')
        end
        
        finalstate=Visualize_chemostatdynamic(chemostat_para,metabfun_consortia,dynamic_drawpara,odepara);
        set(gca,'fontsize',40,'xtick',odepara.timespan_runode,'ytick',0:50:200);
        box on
        
        % the fitness landscape
    end % for showid=1:2
    
    
    % all recorded envs
    pnum=100;
    ms=200;
    lw=1;
    
    allenvs=[steadyenv_byspecies;finalstate(1:2)];
    colorstouse=[species_colors(speciesid,:);
        finalcolor];
    strategy_list=linspace(0,0.8,pnum);
    
    figure;hold on;
    for e=1:size(allenvs,1)
        env=allenvs(e,:)';
        fitness_list=nan*strategy_list;
        for i=1:pnum
            tempspecies=sum_species_para;
            tempspecies.Es=[strategy_list(i),1-strategy_list(i)]';
            tempmetabfun=MetabModel(tempspecies);
            fitness_list(i)=Obtain_Growth_givenenv(env,tempmetabfun);
        end
        % fitness landscape under this enviroment
        
        plot(strategy_list,fitness_list,'color',colorstouse(e,:),'linewidth',lw);
        if ~(e>2)
            scatter(E1_lists(e),chemostat_para.d,ms*1.5,colorstouse(e,:),'d','filled')
        end
    end
    xlim([min(strategy_list),max(strategy_list)])
    ylim([0.9,1.15])
    axis square
    set(gca,'fontsize',40,'ytick',[0.9,1,1.1],'xtick',0.2:0.4:1);
    %xlabel('E_a')
    % ylabel('g')
    box on;
end

%%

if F3_oscillation
    Es=[];
    MetabModel=@Metab_assimilation;
    
    species_num=3;
    kvalue=11;% k=11: stable, k=10: unstable
    
    
    species_colors=[1,0,0;
        0,0.7,0;
        0,0,1];
    
    % chemsotat parameters
    chemostat_para=[];
    chemostat_para.c_supplys=[1,1,1]';
    chemostat_para.d=1;
    
    
    % for k=10;E1s=[0.15,0.2; 0.1, 0.25; 0.26,0.04];
    E1s=[0.15,0.2; % for nutrient a
        0.1, 0.25;  % for nutrient b
        0.26,0.04];
    
    ranges=[0.1,0.6;
        0.1,0.6;
        0.1,0.6;];
    
    
    
    if kvalue==11
        simulationtime=650;
        inistate=[
            0.3490
            0.3765
            0.2483
            0.6984
            0.1228
            0.2183
            19.2551
            11.4387
            59.8347
            73.7865
            20.0625
            8.8176
            10.9874
            78.0682
            15.5982]';
        arrowpercent=20;
    else
        simulationtime=20000;
        inistate=[        0.4923
            0.2223
            0.3557
            0.0460
            0.0003
            0.9537
            24.6703
            8.7548
            76.2771
            93.3764
            15.4867
            11.9423
            14.1813
            57.1879
            20.7934]';
        arrowpercent=20;
    end
    
    arrow_sefrac=[0.1,0.9];
    
    viewangle=[139,30];
    timetick=[0:50:200];
    
    
    
    E2s=[E1s(3,:);
        E1s(1,:);
        E1s(2,:)];
    
    E3s=[E1s(2,:);
        E1s(3,:);
        E1s(1,:)];
    
    E1s=E1s/sum(E1s(:));
    E2s=E2s/sum(E2s(:));
    E3s=E3s/sum(E3s(:));
    Es(1,:)=E1s(:)';
    Es(2,:)=E2s(:)';
    Es(3,:)=E3s(:)';
    
    
    
    
    figure;hold on;
    
    h=[];
    metabfun_consortia=[];
    species_metabfuns=[];
    envs_bys=[];
    steadys=[];
    for s=1:size(Es)
        drawpara=[];
        drawpara.pnum=100;
        drawpara.ranges=ranges;
        drawpara.growthcolor=species_colors(s,:);
        drawpara.fluxcolor=drawpara.growthcolor;
        drawpara.linewidth=1;
        drawpara.showfluxcurve=1;
        drawpara.surfacealpha=0.3;
        
        species=species_iandc_subs;
        species.Es=Es(s,:)';
        species.k=kvalue;
        metabfun=MetabModel(species);
        
        Visualize_chemostatsteady(chemostat_para,metabfun,drawpara);
        [steady,fval]=Obtain_Singlespecies_Steady(chemostat_para,metabfun);
        
        envs_bys(s,:)=steady(1:3);
        metabfun_consortia{s}=metabfun;
        
        h(s)=scatter3(steady(1),steady(2),steady(3),50,species_colors(s,:),'filled');
        steadys(s,1:length(steady))=steady;
    end
    hold on;
    
    set(gca,'fontsize',defaultfs,'xtick',[0.03,0.05],'ytick',[0.03,0.05],'ztick',[0.03,0.05]);
    axis(reshape(ranges',1,6))
    view(viewangle);
    xlabel('{\it c}_a')
    ylabel('{\it c}_b')
    zlabel('{\it c}_c')
    
    %print('-clipboard','-dbitmap','-r600')
    
    % dynamics
    drawpara=[];
    drawpara.arrow_sefrac=arrow_sefrac;
    drawpara.arrow_percent=arrowpercent;
    drawpara.speciescolor=species_colors;
    drawpara.arrow_num=3;
    odepara=[];
    odepara.inistates=inistate;
    odepara.timespan_runode=[1,simulationtime];
    [finalstate,traject]=Visualize_chemostatdynamic(chemostat_para,metabfun_consortia,drawpara,odepara);
    
    
    % fitness by species;
    fitness_s_ss=[];
    for s=1:species_num
        env=envs_bys(s,:);
        for ss=1:species_num% species
            [growthrate,withincellsteady,fval]=Obtain_Growth_givenenv(env,metabfun_consortia{ss});
            fitness_s_ss(s,ss)=growthrate;
        end
    end
    figure;
    Easypcolor(fitness_s_ss,'','');
    axis square
    set(gca,'fontsize',defaultfs);
    caxis([min(fitness_s_ss(:)),max(fitness_s_ss(:))])
    ch=colorbar;
    set(gca,'fontsize',defaultfs);
    
    
    figure;hold on;
    if kvalue==11
        subplot(2,1,1); hold on;
    end
    env_dim=3;
    for s=1:species_num
        plot(traject(:,1),traject(:,env_dim+s+1),'linewidth',1,'color',species_colors(s,:))
    end
    set(gca,'fontsize',defaultfs)
    
    ylabel('Biomass');
    xlim(odepara.timespan_runode)
    box on;
    
    if kvalue==11
        
        % the ever-changing fitness landscape
        fitness_bys_alongtraject=[];
        for t=1:size(traject,1)
            env=traject(t,2:4)';
            if kvalue==1||mod(t,10)==0
                fitnessthree=[];
                for ss=1:species_num
                    fitnessthree(ss)=Obtain_Growth_givenenv(env,metabfun_consortia{ss});
                end
                fitness_bys_alongtraject(:,end+1)=fitnessthree;
                
            end
        end
        subplot(2,1,2)
        Easypcolor(fitness_bys_alongtraject,'','');
        set(gca,'fontsize',defaultfs);
        caxis([min(fitness_bys_alongtraject(:)),max(fitness_bys_alongtraject(:))])
        xlabel('Time')
        %ch=colorbar;set(ch,'ytick',[0.95,1,1.05,1.1])
        % xlim(odepara.timespan_runode)
        
    end
    
    
    
end

%%
if F4_multistable
    minformul_species_para=[];
    minformul_species_para.intercell_dim=0;
    minformul_species_para.Ks=[0.5,0.5]';
    minformul_species_para.Gamma=10;
    minformul_species_para.Es=[nan,nan]';
    minformul_species_para.r=100;
    
    MetabModel=@Metab_MinImport;
    
    metabfun_general=MetabModel(minformul_species_para);
    
    startstrategies=[0.35,0.65];
    progressive_eachround=5;
    ms=200;
    lw=1;
    ranges=[0,0.3;
        0,0.3];
    % parameter for chemostats
    chemostat_para=[];
    chemostat_para.c_supplys=[1,1]';
    chemostat_para.d=1;
    
    % parameter for species
    species_num=10;
    species_colors=distinguishable_colors(species_num, {'w','k'});
    
    %%%%%%%%%%
    % generate the bistability for the initial species
    E1_lists=startstrategies;
    
    steady_drawpara=[];
    steady_drawpara.pnum=1000;
    steady_drawpara.ranges=ranges;
    steady_drawpara.linewidth=lw;
    steady_drawpara.showfluxcurve=0;
    
    
    metabfun_consortia=[];
    steadyenv_byspecies=[];
    figure;
    h=[];
    
    hold on;
    for s=1:length(E1_lists)
        E1=E1_lists(s);
        species_one=minformul_species_para;
        species_one.Es=[E1;1-E1];
        
        metabfun=MetabModel(species_one);
        metabfun_consortia{s}=metabfun;
        
        
        steady_drawpara.growthcolor=species_colors(s,:);
        [gc,h(s)]=Visualize_chemostatsteady(chemostat_para,metabfun,steady_drawpara);
        [steady,fval]=Obtain_Singlespecies_Steady(chemostat_para,metabfun);
        scatter(steady(1),steady(2),ms,steady_drawpara.growthcolor,'filled');
        steadyenv_byspecies(s,:)=steady(1:2);
        
    end
    
    dynamic_drawpara=[];
    dynamic_drawpara.showtimecourse=0;
    dynamic_drawpara.axislim=ranges;
    dynamic_drawpara.linewidths=[0.5,1];
    dynamic_drawpara.arrow_num=1;
    dynamic_drawpara.se_frac=[0.2,0.2];
    
    odepara=[];
    odepara.timespan_runode=[0,100];
    odepara.inistates=[0.25,0.05,500,100
        0.05,0.25,100,500];
    Visualize_chemostatdynamic(chemostat_para,metabfun_consortia,dynamic_drawpara,odepara);
    axis square
    set(gca,'fontsize',defaultfs,'xtick',[0:0.1:0.5],'ytick',[0:0.1:0.5])
    box on;
    axis(reshape(ranges',1,4))
    legend(h,strsplit(num2str(E1_lists)))
    xlabel('')
    ylabel('')
    
    %%%% the fitness landscape in the environment formed by the bistability  species
    figure;hold on;
    strategy_list=linspace(0.3,0.7,1000);
    for se=1:size(steadyenv_byspecies,1)
        env=steadyenv_byspecies(se,:)';
        fitness_list=nan*strategy_list;
        for i=1:length(strategy_list)
            tempspecies=minformul_species_para;
            tempspecies.Es=[strategy_list(i),1-strategy_list(i)]';
            tempmetabfun=MetabModel(tempspecies);
            fitness_list(i)=Obtain_Growth_givenenv(env,tempmetabfun);
        end
        plot(strategy_list,fitness_list,'linewidth',lw,'color',species_colors(se,:));
        scatter(E1_lists(se),chemostat_para.d,3*ms,species_colors(se,:),'d','filled');
    end
    axis([min(strategy_list) max(strategy_list) 0.7,max(fitness_list)*1.05])
    set(gca,'fontsize',40,'xtick',[min(strategy_list) max(strategy_list)],'ytick',1)
    box on
    
    % the chain of invasion
    figure; hold on;
    pnum=1001;
    strategy_list=linspace(startstrategies(1),startstrategies(2),pnum);
    for s=1:2
        current_Ea=startstrategies(s);
        for r=1:(progressive_eachround+1)
            
            drawspecies=1;
            currentspecies=minformul_species_para;
            if r>progressive_eachround
                if s==2
                    current_Ea=0.5;
                    colorcode=[0,0,0];
                else
                    drawspecies=0;
                end
            else
                
                colorid=r*2+(s-2);
                colorcode=species_colors(colorid,:);
            end
            
            if drawspecies
                currentspecies.Es=[current_Ea,1-current_Ea]';
                currentmetabfun=MetabModel(currentspecies);
                % the environment constructed by the current species
                [steady,fval]=Obtain_Singlespecies_Steady(chemostat_para,currentmetabfun);
                currentenv=steady(1:2);
                
                % fitness landscape under current envirnoment
                fitness_list=nan*strategy_list;
                for i=1:pnum
                    tempspecies=minformul_species_para;
                    tempspecies.Es=[strategy_list(i),1-strategy_list(i)]';
                    tempmetabfun=MetabModel(tempspecies);
                    fitness_list(i)=Obtain_Growth_givenenv(currentenv,tempmetabfun);
                end
                
                
                
                plot(strategy_list,fitness_list,'color',colorcode,'linewidth',lw); hold on;
                scatter(current_Ea,chemostat_para.d,1.5*ms,colorcode,'d','filled');
                
                % choose the next strategy
                [maxv,maxloc]=max(fitness_list);
                current_Ea=strategy_list(maxloc);
            end
        end % for r
    end % for s
    axis([min(strategy_list),max(strategy_list),0.95,1.08])
    set(gca,'fontsize',defaultfs,'xtick',0:0.1:1,'ytick',0.95:0.05:2);
    box on;
    
    %%%%%% optimization
    
    pnumbydim=3*[20,20];
    
    
    [nutreint_boarders_m,nutreint_values_m,optgrowthrates_m,opt_Es,consumptions_m,withincellv_m]=Nutrientspace_optspecies(minformul_species_para,MetabModel,ranges,pnumbydim);
    resultsfromoptspecies={nutreint_boarders_m,nutreint_values_m,optgrowthrates_m,opt_Es,consumptions_m,withincellv_m};
    
    opt_drawpara=[];
    opt_drawpara.colornutrientspace=1;
    opt_drawpara.grwoth_colormap=colormap('bone');
    opt_drawpara.plotGC=3;
    opt_drawpara.plotGCcolor=colormap('jet');
    opt_drawpara.minmaxv_colors=[0.2,0.8];
    
    figure;hold on;
    [alongGCs,alongFBs]=Analysis_OptM(metabfun_general,resultsfromoptspecies,chemostat_para.d,[],opt_drawpara);
    shading interp
    % choose some species on the supply curve
    
    lw=2;
    steady_drawpara=[];
    steady_drawpara.pnum=1000;
    steady_drawpara.linewidth=lw;
    steady_drawpara.showfluxcurve=0;
    
    alongGC=alongGCs{1};
    teststrategies=[startstrategies,0.5];
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
    
    
    set(gca,'fontsize',defaultfs,'xtick',0:0.1:1,'ytick',0.1:0.1:1);
    axis(reshape(ranges',1,4)+[0,-0.01,0,-0.01])
    axis square;
    xlabel('')
    ylabel('')
    
    
end




%%
if FS2_spatial
    minformul_species_para=[];
    minformul_species_para.intercell_dim=0;
    minformul_species_para.Ks=[0.5,0.5]';
    minformul_species_para.Gamma=10;
    minformul_species_para.Es=[nan,nan]';
    minformul_species_para.r=100;
    
    MetabModel=@Metab_MinImport;
    
    metabfun_general=MetabModel(minformul_species_para);
    
    startstrategies=[0.35,0.65];
    
    
    leakrate=1; % 10, 0.5
    chemostatnum=20;
    chemostat_network=zeros(chemostatnum);
    % a "connect to nearest neighbor" pattern
    for c=1:(chemostatnum-1)
        chemostat_network(c,c+1)=1;
        chemostat_network(c+1,c)=1;
    end
    
    % parameter for chemostats
    chemostat_para=[];
    chemostat_para.c_supplys=[1,1]';
    chemostat_para.d=1;
    
    
    
    Eslegend=[];
    E1list=linspace(startstrategies(1),startstrategies(2),2);
    for s=1:length(E1list)
        E1=E1list(s);
        species_one=minformul_species_para;
        species_one.Es=[E1;1-E1];
        species{s}=species_one;
        Eslegend{s}=num2str(E1,'%0.2f');
    end
    
    species_num=length(species);
    species_colors=distinguishable_colors(species_num, {'w','k'});
    
    
    species_metabfuns=[];
    
    drawpara=[];
    drawpara.pnum=1000;
    drawpara.ranges(1,:)=[0.01,1];
    drawpara.ranges(2,:)=[0.01,1];
    drawpara.growthcolor=species_colors;%localpoints(k,:);
    
    drawpara.ms=200;
    drawpara.linewidth=1;
    drawpara.showfluxcurve=0;
    
    
    figure;
    hold on;
    h=[];
    steadyenv_byspecies=zeros(length(species),2);
    fitness_species_inenvs=zeros(length(species));
    for s=1:length(species)
        metabfun=MetabModel(species{s});
        species_metabfuns{s}=metabfun;
        drawpara.growthcolor=species_colors(s,:);%localpoints(k,:);
        drawpara.fluxcolor=drawpara.growthcolor;
        drawpara.surfacealpha=0.5;
        Visualize_chemostatsteady(chemostat_para,metabfun,drawpara);
        
        [steady,fval]=Obtain_Singlespecies_Steady(chemostat_para,metabfun);
        h(s)=scatter(steady(1),steady(2),50,species_colors(s,:),'filled');
        
        steadyenv_byspecies(s,:)=steady(1:2);
        % the steady-state fitness of each species under the envionrment
        % constructed by diferent species
    end
    
    
    dxdt_noleak=@(t,x)ODE_multichemostat(x,chemostat_para,species_metabfuns,0,chemostat_network);
    dxdt_withleak=@(t,x)ODE_multichemostat(x,chemostat_para,species_metabfuns,leakrate,chemostat_network);
    
    
    inistate_byc=zeros(length(chemostat_para.c_supplys)+length(species_metabfuns),chemostatnum);
    %for c=1:chemostatnum
    %    inistate_byc(length(chemostat_para.c_supplys)+1,c)=0.1;
    %end
    inistate_byc(length(chemostat_para.c_supplys)+1,1)=0.1;
    inistate_byc(length(chemostat_para.c_supplys)+2,chemostatnum)=0.2;
    ini=reshape(inistate_byc,length(inistate_byc(:)),1);
    
    % run to steady state with zero leakage, then start to connect
    [timelist1,trajecties1]=ode23(dxdt_noleak,[0,100],ini);
    [timelist2,trajecties2]=ode23(dxdt_withleak,timelist1(end)+[0,5000],trajecties1(end,:)');
    timelist=[timelist1;timelist2];
    trajecties=[trajecties1;trajecties2];
    
    env_byc=[];
    spe_byc=[];
    figure;
    hold on;
    for c=1:chemostatnum
        envdim=(c-1)*(size(inistate_byc,1))+(1:length(chemostat_para.c_supplys));
        speciesdim=(c-1)*(size(inistate_byc,1))+length(chemostat_para.c_supplys)+(1:length(species_metabfuns));
        env_byc(c,:)=trajecties(end,envdim);
        spe_byc(c,:)=trajecties(end,speciesdim);
        
        plot(timelist,trajecties(:,envdim(1)));
    end
    
    % get a value on information
    threshhold=0.01;
    env_byint=round(env_byc/threshhold);
    [distinctstates,rep,ids]=unique(env_byint,'rows');
    tabfreq=tabulate(ids);
    prob=tabfreq(:,3)/100;
    mulinfor=-1*sum(prob.*log2(prob)); % in bit
    
    
    figure;hold on; % the species
    lw=1;
    for s=1:species_num
        plot(spe_byc(:,s),'linewidth',lw,'color',species_colors(s,:));
    end
    plot(sum(spe_byc,2),'--','linewidth',lw,'color',[0,0,0])
    set(gca,'fontsize',22);
    ylabel('Biomass Density')
    xlabel('Chemostat ID');
    xlim([1,chemostatnum])
    
    
    figure;hold on; % the concentration
    env_colors=[0,0.8,0;
        0,0.8,0.8];
    env_num=2;
    for e=1:env_num
        plot(env_byc(:,e),'linewidth',lw,'color',env_colors(e,:));
    end
    set(gca,'fontsize',22);
    ylabel('Concentration')
    xlabel('Chemostat ID');
    xlim([1,chemostatnum])
    title(mulinfor)
    
    
    % for each enviornment, what's the fitness values
    strategy_list=linspace(0.3,0.7,50);
    
    fitness_cid_strategy=zeros(chemostatnum,length(strategy_list));
    
    for c=1:chemostatnum
        env=env_byc(c,:)';
        for i=1:length(strategy_list)
            tempspecies=minformul_species_para;
            tempspecies.Es=[strategy_list(i),1-strategy_list(i)]';
            tempmetabfun=MetabModel(tempspecies);
            fitness_cid_strategy(c,i)=Obtain_Growth_givenenv(env,tempmetabfun);
        end
    end
    figure;
    set(gca,'fontsize',defaultfs)
    [X,Y]=meshgrid(1:chemostatnum,strategy_list);
    h=surface(X,Y,fitness_cid_strategy');
    h.LineStyle='none';
    %shading interp
    xlim([1,chemostatnum])
    ch=colorbar;
    xlabel('Chemostat ID')
    ylabel('Strategies')
    ylabel(ch,'Fitness');
    
end

%%
if F5_optinvasion
    ms=200;
    contour_ms=100;
    contour_lw=1;
    supplyv_lw=1;
    fitness_lw=1;
    ranges=[0,1;
        0,1];
    patchcolor=0.8*[1,1,1];
    
    species_colors=[0,0.8,0;0,0,1;0.9,0,0];
    env_dim=2;
    
    MetabModel=@Metab_ConvertNecessary;
    
    chemostat_para=[];
    chemostat_para.d=0.2;
    chemostat_para.c_supplys=[0.5;1];
    
    % test t model
    Gamma=1;
    K_1=1;
    K_2=1;
    p=0.2;
    E_1=0.7;
    E_2=0;
    E_3=0;
    E_4=1-E_1-E_2-E_3;
    
    general_species.Gamma=Gamma;
    general_species.Ks=[K_1;
        K_2;
        ];
    general_species.p=p;
    general_species.Es=[E_1;
        E_2;
        E_3;
        E_4;];
    general_species.intercell_dim=2;
    general_species.r=100;
    general_metabfun=MetabModel(general_species);
    
    
    %% different
    % test
    convertera2b_species_para=general_species;
    convertera2b_species_para.Es=[0.5,0,0,0.5]';
    
    converterb2a_species_para=general_species;
    converterb2a_species_para.Es=[0,0.5,0.5,0]';
    
    importer_species_para=general_species;
    importer_species_para.Es=[0.5,0.5,0,0]';
    
    Ca2b_id=1;
    Cb2a_id=2;
    Impo_id=3;
    species_consortia={convertera2b_species_para,converterb2a_species_para,importer_species_para};
    importer_otherE_id=[1,4;
        2,3;
        1,2;];
    
    % for optimization
    
    pnumbydim=10*[10,10];
    
    if exist('optgrowth_m','var')==0
        optresults_path=['./output/optinvasion',strjoin(strsplit(num2str(pnumbydim)),'-'),'.mat'];
        if exist(optresults_path,'file')
            load(optresults_path);
        else
            showprogress=1;
            [nutreint_boarders_m,nutreint_values_m,optgrowth_m,optEs_m,consumptions_m,withincellv_m]=Nutrientspace_optspecies(general_species,MetabModel,ranges,pnumbydim,showprogress);
            
            save(optresults_path,'nutreint_boarders_m','nutreint_values_m','optgrowth_m','optEs_m','consumptions_m','withincellv_m');
        end
    end
    
    % determine the species type for each optimal solution
    optsid=nan*optgrowth_m;
    colormatrix=nan*zeros(pnumbydim(1),pnumbydim(1),3);
    
    for  i=1:pnumbydim(1)
        for j=1:pnumbydim(2)
            optE_thisenv=optEs_m(i,j,:);
            [maxv,maxloc]=max(sum(optE_thisenv(importer_otherE_id),2));
            if maxv>0.9
                sid=maxloc;
                colormatrix(i,j,:)=species_colors(sid,:);
            else
                sid=0;
            end
            optsid(i,j)=sid;
        end
    end
    
    contourlevels=0.1:0.1:0.4;
    
    figure;
    hold on;
    ph=surface(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),0*optsid,colormatrix);
    ph.LineStyle='none';
    
    bonecolortable=colormap('bone');
    colorids=round(linspace(1,size(bonecolortable,1),length(contourlevels)));
    for i=1:length(contourlevels)
        contour(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optgrowth_m,contourlevels(i)*[1,1],'color',bonecolortable(colorids(i),:),'linewidth',contour_lw)
    end
    axis square
    set(gca,'fontsize',defaultfs,'xtick',0:0.5:1,'ytick',0:0.5:1);
    ch=colorbar();
    caxis([min(contourlevels),max(contourlevels)])
    %caxis([min(optgrowth_m(:)),max(optgrowth_m(:))]);
    set(ch,'ytick',0:0.2:1)
    box on;
    
    
    % strategy in the nutrient space
    figure;
    for e=1:4
        subplot(2,2,e); hold on;
        
        colormap('default')
        ph=pcolor(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optEs_m(:,:,e));
        ph.LineStyle='none';
        ch=colorbar;
        caxis([0 1])
        set(ch,'ytick',[0,1])
        
        for i=1:length(contourlevels)
            contour(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optgrowth_m,contourlevels(i)*[1,1],'color',bonecolortable(colorids(i),:),'linewidth',contour_lw)
        end
        
        axis square
        set(gca,'fontsize',18,'xtick',[0,1],'ytick',[0,1]);
        box on
    end
    
    % extract the environment at the conjuction point
    % obtain the optimal species of two types, in that environment, draw their
    % cgrowth contour and supply line
    figure;
    optdrawpara=[];
    optdrawpara.colornutrientspace=0;
    optdrawpara.plotGC=0;
    
    resultsfromoptspecies={nutreint_boarders_m,nutreint_values_m,optgrowth_m,optEs_m,consumptions_m,withincellv_m};
    alongGCs=Analysis_OptM(general_metabfun,resultsfromoptspecies,chemostat_para.d,[],optdrawpara);
    alongGC=alongGCs{1};
    
    % find there the strategy switch type
    sid_alongGC=[];
    for i=1:size(alongGC.envs,2)
        sid_alongGC(i)=optsid(alongGC.D_ijs(1,i),alongGC.D_ijs(2,i));
    end
    loc_beforeswitch=find(sid_alongGC(1:(end-1))~=sid_alongGC(2:end));
    
    % locate the environment for switching to happen
    for sw=1
        loc_before=loc_beforeswitch(sw);
        loc_after=loc_beforeswitch(sw)+1;
        env_before=reshape(nutreint_values_m(alongGC.D_ijs(1,loc_before),alongGC.D_ijs(2,loc_before),:),env_dim,1);
        env_after=reshape(nutreint_values_m(alongGC.D_ijs(1,loc_after),alongGC.D_ijs(2,loc_after),:),env_dim,1);
        
        env_middle=mean([env_before,env_after],2);
        
        % optimal species for this environment in involving type
        optconsortia_paras=[];
        stypelist=sid_alongGC([loc_before,loc_after]);
        for s=stypelist
            this_speciespara=species_consortia{s};
            this_metabfun=MetabModel(this_speciespara);
            nonzeroEdim=importer_otherE_id(s,:);
            
            inis=[];
            inis.iniEs=this_speciespara.Es;
            inis.inig= chemostat_para.d;
            inis.inixs=[0.1,0.1]';
            
            [optE_thiss,optg_thiss,fval]=Obtain_optEG_givenenv(env_middle,this_metabfun,inis,nonzeroEdim);
            
            optspecies=this_speciespara;
            optspecies.Es=optE_thiss;
            optconsortia_paras{end+1}=optspecies;
        end
        
        figure;
        hold on;
        steady_drawpara=[];
        steady_drawpara.showfluxcurve=0;
        
        % the optimal growth contour
        % for each species, how it maps to colors
        minmaxEimporter_bys=zeros(3,2);
        for s=1:3
            importerdim=importer_otherE_id(s,1);
            Esforthis=alongGC.Es_optGC(importerdim,sid_alongGC==s);
            minmaxEimporter_bys(s,:)=[min(Esforthis),max(Esforthis)];
        end
        
        GC_ijs=[];
        GC_envs=[];
        GC_colors=[];
        for k=1:size(alongGC.envs,2)
            if mod(k,4)==1 &&(k==1 || sum(abs(alongGC.D_ijs(:,k)-GC_ijs(:,end)))>1)
                sid=sid_alongGC(k);
                GC_ijs(:,end+1)=alongGC.D_ijs(:,k);
                GC_envs(:,end+1)=alongGC.envs(:,k);
                
                importerdim=importer_otherE_id(sid,1);
                localE=alongGC.Es_optGC(importerdim,k);
                GC_colors(end+1,:)=species_colors(sid,:)+[1,1,1]*(minmaxEimporter_bys(sid,2)-localE)/(0.1+minmaxEimporter_bys(sid,2)-minmaxEimporter_bys(sid,1));
            end
        end
        
        scatter(GC_envs(1,:)',GC_envs(2,:)',contour_ms,GC_colors,'s')
        
        envlists=[];env_byspecies=[];
        optconsortia_metabfuns=[];
        coexistregion_xy=env_middle;
        envcolors=[];
        for s=1:length(optconsortia_paras)
            sid=stypelist(s);
            this_speciespara=optconsortia_paras{s};
            this_metabfun=MetabModel(this_speciespara);
            optconsortia_metabfuns{end+1}=this_metabfun;
            
            steady_drawpara.growthcolor=species_colors(sid,:);
            steady_drawpara.ranges=ranges;
            Visualize_chemostatsteady(chemostat_para,this_metabfun,steady_drawpara);
            
            steadyv=Obtain_Singlespecies_Steady(chemostat_para,this_metabfun);
            scatter(steadyv(1),steadyv(2),ms,steady_drawpara.growthcolor,'filled');
            envcolors(end+1,:)=steady_drawpara.growthcolor;
            % supply vector
            [FBv_givenIa,onelargesupplypoint]=Obtain_supplyvector(this_metabfun,env_middle);
            plot([env_middle(1),onelargesupplypoint(1)],[env_middle(2),onelargesupplypoint(2)],'k :','linewidth',supplyv_lw) ;
            
            envlists(:,end+1)=steadyv(1:2);
            env_byspecies{end+1}=sid;
            
            coexistregion_xy(:,end+1)=onelargesupplypoint;
        end
        
        % patch the vector
        ph=patch(coexistregion_xy(1,:),coexistregion_xy(2,:),patchcolor);
        ph.LineStyle='none';
        ph.FaceAlpha=0.5;
        
        scatter(chemostat_para.c_supplys(1),chemostat_para.c_supplys(2),ms,'k')
        axis square
        box on;
        set(gca,'fontsize',defaultfs,'xtick',0:0.5:1,'ytick',0:0.5:1)
        xlabel('')
        ylabel('')
        
        
        drawpara=[];
        drawpara.axislim=ranges;
        drawpara.speciescolor=species_colors(stypelist,:);
        drawpara.arrow_num=1;
        drawpara.showtimecourse=0;
        odepara=[];
        odepara.inistates=[0.75,0.6,0.05,0.05,5*ones(1,4)];
        odepara.timespan_runode=[0,500];
        [finalstate,traject]=Visualize_chemostatdynamic(chemostat_para,optconsortia_metabfuns,drawpara,odepara);
        scatter(finalstate(1),finalstate(2),ms,sum(species_colors(stypelist,:)),'filled')
        
        % time course for simulating lots of species on the optimal curve
        morespecie_num=10;
        otherstrategies=(unique(alongGC.Es_optGC','rows'))';
        manyspecies_metabfun=optconsortia_metabfuns;
        toextract_ids=round(linspace(1,size(otherstrategies,2),morespecie_num));
        for i=1:morespecie_num
            Es=otherstrategies(:,toextract_ids(i));
            
            if sum(abs(Es-optconsortia_paras{1}.Es))>0.01 & sum(abs(Es-optconsortia_paras{2}.Es))>0.01
                onespecies=general_species;
                onespecies.Es=Es;
                manyspecies_metabfun{end+1}=MetabModel(onespecies);
            end
        end
        length(manyspecies_metabfun)
        
        showcompetition=1;
        if showcompetition==1
            timespan_runode= [0,100000];
            ini=[0,0,0.01*ones(1,length(manyspecies_metabfun)),5*ones(1,2*length(manyspecies_metabfun))]';
            odetouse=@(t,x)ODE_onechemostat(x,chemostat_para,manyspecies_metabfun);
            [timelist,trajecties]=ode23(odetouse,timespan_runode,ini);
            
            colortable=distinguishable_colors(length(manyspecies_metabfun)+3, {'k'});
            colortable(1:3,:)=[];
            figure;hold on;
            for s=1:length(manyspecies_metabfun)
                if s<3
                    colorcode=species_colors(stypelist(s),:);
                else
                    colorcode=colortable(s,:);
                end
                plot(timelist,trajecties(:,env_dim+s),'linewidth',1,'color',colorcode)
            end
            xlim([10,max(timelist)])
            set(gca,'xscale','log')
            set(gca,'fontsize',40,'xtick',[10,1000,max(timespan_runode)])
            box on;
        end
        
        
        
        
        
        
        
        
        % fitnes landscape at this final environment
        envlists(:,end+1)=[finalstate(1:env_dim)'];
        env_byspecies{end+1}=stypelist;
        
        envcolors(end+1,:)=sum(species_colors(stypelist,:));
        figure;
        ms=150;
        for e=1:size(envlists,2)
            
            
            subplot(3,2,e*2-1); hold on;
            
            env=envlists(:,e);
            strategy_list=linspace(0.1,0.99,500);
            maxf=0;
            for s=1:3
                fitness_list=nan*strategy_list;
                for i=1:length(strategy_list)
                    Es=zeros(4,1);
                    Es(importer_otherE_id(s,1))=strategy_list(i);
                    Es(importer_otherE_id(s,2))=1-strategy_list(i);
                    
                    tempspecies=this_speciespara;
                    tempspecies.Es=Es;
                    tempmetabfun=MetabModel(tempspecies);
                    [gr,withincellsteady]=Obtain_Growth_givenenv(env,tempmetabfun);
                    
                    if max(withincellsteady)>0
                        fitness_list(i)=gr;
                    end
                    
                    maxf=max(maxf,max(fitness_list));
                end
                plot(strategy_list,fitness_list,'linewidth',fitness_lw,'color',species_colors(s,:));
                if ismember(s,env_byspecies{e})
                    loc=find(stypelist==s);
                    Es=optconsortia_paras{loc}.Es;
                    scatter(Es(importer_otherE_id(s,1)),chemostat_para.d,ms,species_colors(s,:),'d','filled')
                end
                
                xlim([min(strategy_list),max(strategy_list)])
                ylim([min(fitness_list),maxf*1.05])
                %axis square
                box on
                set(gca,'fontsize',15,'xtick',[0.1,0.9],'ytick',0:0.2:1)
                ax = gca;
                ax.XColor = envcolors(e,:);
                ax.YColor = envcolors(e,:);
            end
        end % for e=1:size(envlists,2)
        
    end % for sw=1
    
end






%%
if F6_newdim
    paren = @(x, varargin) x(varargin{:});
    
    
    env_dim=2;
    max_speciesnum=3;
    species_colors=[0,0,1;
        0,0.8,0;
        0.9,0,0];
    
    
    MetabModel=@Metab_crossfeeding;
    contour_lw=1;
    
    % get the optimal species and corresponding growth rates
    ranges=[0,1.5;%1.5
        0,0.3];%0.25
    pnumbydim=[100,100]*2;
    drawpara=[];
    showprogress=1;
    
    % parameter for chemostats
    chemostat_para=[];
    chemostat_para.c_supplys=[1.8,0]'; % supply 1.5
    chemostat_para.d=0.6; %0.3, 0.6
    
    % parameter for species
    newdim_species_para=[];
    newdim_species_para.intercell_dim=1;
    newdim_species_para.r=100;
    
    
    newdim_species_para.Es= [0.2,0.1,0.6,0.1];%[0.2,0.3,0.4,0.1];
    
    newdim_species_para.V_1=5;
    newdim_species_para.V_2=1;
    newdim_species_para.k=8;
    newdim_species_para.V_4=10;
    
    newdim_species_para.Ks=[0.5, 0.5, 0.5, 0.1, 0.5, 15,10];
    newdim_species_metabfun=MetabModel(newdim_species_para);
    
    %%% define three species
    % G
    species_special=newdim_species_para;
    species_special.Es=[0.5,0.5,0,0];
    metabfun=MetabModel(species_special);
    optEs=metabfun.analytical_optE_givenDW;% (D,W);
    optG_Es=optEs(chemostat_para.d,0);
    species_optG=species_special;
    species_optG.Es=optG_Es;
    
    % S1
    species_special.Es=[0.5,0,0.5,0];
    metabfun=MetabModel(species_special);
    optEs=metabfun.analytical_optE_givenDW;% (D,W);
    optS1_Es=optEs(chemostat_para.d,0);% note that it is optimized for 0
    species_goodS1=species_special;
    species_goodS1.Es=optS1_Es;
    
    %S2
    species_special.Es=[0,0.5,0,0.5];
    metabfun=MetabModel(species_special);
    optEs=metabfun.analytical_optE_givenDW;% (D,W);
    minSW=metabfun.analytical_minSW_givenDW;
    optS2_Es=optEs(chemostat_para.d,0);
    
    species_optS2=species_special;
    species_optS2.Es=optS2_Es;
    
    species_consortia={species_goodS1,species_optS2,species_optG};
    S1id=1;
    S2id=2;
    Gid=3;
    
    %%%%%%%%%%%%%%%%%%%%
    importer_otherE_id=[1,3;
        4,2;
        1,2;];
    
    if ~exist('optEs_all','var')
        
        optresults_path=['./output/newdim',strjoin(strsplit(num2str(pnumbydim)),'-'),'.mat'];
        if exist(optresults_path,'file')
            load(optresults_path);
        else
            
            optgrowth_bys=zeros(pnumbydim(1),pnumbydim(2),3);
            optEs_bys=zeros(pnumbydim(1),pnumbydim(2),4,3);
            consumptions_bys=zeros(pnumbydim(1),pnumbydim(2),2,3);
            withincellv_bys=zeros(pnumbydim(1),pnumbydim(2),3);
            for s=1:3
                speciespara=species_consortia{s};
                [nutreint_boarders_m,nutreint_values_m,optgrowth_bys(:,:,s),optEs_bys(:,:,:,s),consumptions_bys(:,:,:,s),withincellv_bys(:,:,s)]=Nutrientspace_optspecies(speciespara,MetabModel,ranges,pnumbydim,showprogress,drawpara);
            end
            
            save(optresults_path,'nutreint_boarders_m','nutreint_values_m','optgrowth_bys','optEs_bys','consumptions_bys','withincellv_bys');
        end
        
        
        % compare and get the optimal species on every point
        optEs_all=nan*zeros(pnumbydim(1),pnumbydim(2),4);
        optgrowth_all=nan*zeros(pnumbydim(1),pnumbydim(2));
        consumptions_all=nan*zeros(pnumbydim(1),pnumbydim(2),2);
        optsid=nan*zeros(pnumbydim(1),pnumbydim(2));
        for i=1:pnumbydim(1)
            for j=1:pnumbydim(2)
                
                [maxv,sid]=max(optgrowth_bys(i,j,:));
                if maxv>0
                    optsid(i,j)=sid;
                    optgrowth_all(i,j)=maxv;
                    optEs_all(i,j,:)=optEs_bys(i,j,:,sid);
                    consumptions_all(i,j,:)=consumptions_bys(i,j,:,sid);
                end
            end
        end
    end
    
    contourlevels=[0.2:0.2:0.8];
    % the growth rates
    
    colormatrix=nan*zeros(pnumbydim(1),pnumbydim(1),3);
    for i=1:pnumbydim(1)
        for j=1:pnumbydim(2)
            if ~isnan(optsid(i,j))
                
                % get the importer E value
                
                sid=optsid(i,j);
                speciescolor=species_colors(sid,:);
                importerEvalue=optEs_all(i,j,importer_otherE_id(sid,1));
                colormatrix(i,j,:)= importerEvalue*speciescolor+[1,1,1]*(1-importerEvalue);
                
                colormatrix(i,j,:)=species_colors(optsid(i,j),:);
            end
        end
    end
    
    bonecolortable=colormap('bone');
    colorids=round(linspace(1,size(bonecolortable,1),length(contourlevels)));
    % species id in  utrient space
    
    figure;hold on;
    ph=surface(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),0*optsid,colormatrix);
    ph.LineStyle='none';
    hold on;
    %xlabel('S')
    %ylabel('W')
    for i=1:length(contourlevels)
        contour(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optgrowth_all,contourlevels(i)*[1,1],'color',bonecolortable(colorids(i),:),'linewidth',contour_lw)
    end
    colormap(bonecolortable)
    ch=colorbar;
    caxis([min(contourlevels),max(contourlevels)])
    set(ch,'ytick',contourlevels)
    axis square
    axis(reshape(ranges',1,4))
    set(gca,'fontsize',defaultfs,'xtick',[0:0.5:1.5],'ytick',0.1:0.1:0.3)
    
    
    
    % the strategy
    figure;
    for e=1:4
        subplot(2,2,e); hold on;
        colormap('default')
        ch=colorbar;
        ph=pcolor(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optEs_all(:,:,e));
        ph.LineStyle='none';
        caxis([0 1])
        
        for i=1:length(contourlevels)
            contour(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optgrowth_all,contourlevels(i)*[1,1],'color',bonecolortable(colorids(i),:),'linewidth',contour_lw)
        end
        
        set(ch,'ytick',[0,1])
        set(gca,'fontsize',15,'xtick',[0.5,1.5],'ytick',[0,0.2]);
        axis square
    end
    
    
    species_num=3;
    optmetabfun_consortia=[];
    optspecies_consortia=[];
    
    
    Dlist=[0.4,0.6];
    supplybyd=[];
    supplybyd(:,1)=[1;
        0];
    supplybyd(:,2)=[1.8;
        0];
    
    ini_mass=1;
    add_mass=5;
    
    
    contour_ms=100;
    steady_ms=200;
    fitness_ms=150;
    supplyv_lw=1;
    contour_lw=1;
    fitness_lw=1;
    patchcolor=0.5*ones(1,3);
    
    for i=1:length(Dlist)
        
        % define controllable under this condition.
        D=Dlist(i);
        chemostat_para.d=D;
        chemostat_para.c_supplys=supplybyd(:,i);
        axisrange=[ranges(1,1),max(ranges(1,2),supplybyd(1,i)); ranges(2,:)];
        
        if i==1
            axisrange=[0 1.05;0,0.2];
        else
             axisrange=[0 1.2;0,0.2];
        end
        
        
        % the maximal growth contour
        optdrawpara=[];
        optdrawpara.colornutrientspace=0;
        optdrawpara.plotGC=0;
        
        resultsfromoptspecies={nutreint_boarders_m,nutreint_values_m,optgrowth_all,optEs_all,consumptions_all,withincellv_bys};
        alongGCs=Analysis_OptM(newdim_species_metabfun,resultsfromoptspecies,chemostat_para.d,[],optdrawpara);
        alongGC=alongGCs{1};
        
        sid_alongGC=[];
        for k=1:size(alongGC.envs,2)
            sid_alongGC(k)=optsid(alongGC.D_ijs(1,k),alongGC.D_ijs(2,k));
        end
        loc_beforeswitch=find(sid_alongGC(1:(end-1))~=sid_alongGC(2:end));
        
        % draw the maximal growth contour
        minmaxEimporter_bys=zeros(3,2); % for each species, how it maps to colors
        for s=1:3
            importerdim=importer_otherE_id(s,1);
            thisslist=find(sid_alongGC==s);
            if ~isempty(thisslist)
                Esforthis=alongGC.Es_optGC(importerdim,thisslist);
                minmaxEimporter_bys(s,:)=[min(Esforthis),max(Esforthis)];
            end
        end
        GC_ijs=[];
        GC_envs=[];
        GC_colors=[];
        for k=1:size(alongGC.envs,2)
            if mod(k,8)==1 &&(k==1 || sum(abs(alongGC.D_ijs(:,k)-GC_ijs(:,end)))>1)
                sid=sid_alongGC(k);
                GC_ijs(:,end+1)=alongGC.D_ijs(:,k);
                GC_envs(:,end+1)=alongGC.envs(:,k);
                
                importerdim=importer_otherE_id(sid,1);
                localE=alongGC.Es_optGC(importerdim,k);
                GC_colors(end+1,:)=species_colors(sid,:)+[1,1,1]*(minmaxEimporter_bys(sid,2)-localE)/(0.01+minmaxEimporter_bys(sid,2)-minmaxEimporter_bys(sid,1));
            end
        end
        
        figure;hold on;
        box on;
        axis square
        scatter(GC_envs(1,:)',GC_envs(2,:)',contour_ms,GC_colors,'s')
        
        steady_drawpara=[];
        steady_drawpara.ranges=[ranges(1,1),max(ranges(1,2),supplybyd(1,i));
            ranges(2,:)];
        steady_drawpara.showfluxcurve=0;
        steady_drawpara.linewidth=contour_lw;
        steady_drawpara.pnum=1000;
        
        
        for sw=1:length(loc_beforeswitch)
            loc_before=loc_beforeswitch(sw);
            loc_after=loc_beforeswitch(sw)+1;
            env_before=reshape(nutreint_values_m(alongGC.D_ijs(1,loc_before),alongGC.D_ijs(2,loc_before),:),env_dim,1);
            env_after=reshape(nutreint_values_m(alongGC.D_ijs(1,loc_after),alongGC.D_ijs(2,loc_after),:),env_dim,1);
            env_middle=mean([env_before,env_after],2);
            
            % optimal species at this point
            optconsortia_paras=[];
            stypelist=sid_alongGC([loc_before,loc_after]);
            
            % determine if this is the point to show dynamics
            showsimulation=0;
            if ismember(1,stypelist)
                showsimulation=1;
            end
            
            coexistregion_xy=env_middle;
            for s=stypelist
                this_speciespara=species_consortia{s};
                this_metabfun=MetabModel(this_speciespara);
                nonzeroEdim=importer_otherE_id(s,:);
                
                inis=[];
                inis.iniEs=this_speciespara.Es;
                inis.inig= chemostat_para.d;
                inis.inixs=[0.1,0.1]';
                
                [optE_thiss,optg_thiss,fval]=Obtain_optEG_givenenv(env_middle,this_metabfun,inis,nonzeroEdim);
                
                optspecies=this_speciespara;
                optspecies.Es=optE_thiss;
                optconsortia_paras{end+1}=optspecies;
                opt_metabfun=MetabModel(optspecies);
                
                % the supply vector at this point
                [FBv_givenIa,onelargesupplypoint]=Obtain_supplyvector(opt_metabfun,env_middle);
                plot([env_middle(1),onelargesupplypoint(1)],[env_middle(2),onelargesupplypoint(2)],'k :','linewidth',supplyv_lw) ;
                
                coexistregion_xy(:,end+1)=onelargesupplypoint;
                
                optmetabfun_consortia{s}=opt_metabfun;
                optspecies_consortia{s}=optspecies;
                % show species GC if necessary
                if showsimulation
                    
                    
                    steady_drawpara.growthcolor=species_colors(s,:);
                    Visualize_chemostatsteady(chemostat_para,opt_metabfun,steady_drawpara);
                    
                    steadyv=Obtain_Singlespecies_Steady(chemostat_para,opt_metabfun);
                    if 0%steadyv(3)>10^-2% non zero biomass
                        scatter(steadyv(1),steadyv(2),steady_ms,steady_drawpara.growthcolor,'filled')
                    end
                    axisrange(:,2)=max(axisrange(:,2),steadyv(1:2)*1.1);
                end
                
            end
            ph=patch(coexistregion_xy(1,:),coexistregion_xy(2,:),patchcolor);
            ph.LineStyle='none';
            ph.FaceAlpha=0.5;
        end % for sw=1:length(loc_beforeswitch)
        
        scatter(chemostat_para.c_supplys(1), chemostat_para.c_supplys(2)+(axisrange(2,1)),steady_ms,'k')
        axis(reshape(axisrange',1,4));
        set(gca,'fontsize',defaultfs)
        set(gca,'xtick',0:0.5:1.5,'ytick',0.1:0.1:0.4);
        xlabel('')
        ylabel('')
        
        % simulation
        if showsimulation
            % species to start
            if i==1
                firstround_species=S1id;
                secondround_species=[S1id,S2id];
                timerun_12=[100,200];
            else
                firstround_species=[S1id,S2id];
                secondround_species=[S1id,Gid];
                timerun_12=[100,200];
            end
            
            odetouse=@(t,x)ODE_onechemostat(x,chemostat_para,optmetabfun_consortia);
            
            timespan1=[0,timerun_12(1)];
            timespan2=timerun_12(1)+[0,timerun_12(2)];
            inienv=[0.1,0.1];
            inibiomass=zeros(1,length(optmetabfun_consortia));
            inibiomass(firstround_species)=ini_mass;
            ini_incell=0*ones(1,length(optmetabfun_consortia));
            if i==1
                ini1=[inienv,inibiomass,ini_incell]';
            else
                ini1=[    1.4364
                    0.2501
                    0.3636
                    0.1134
                    0
                    0.4205
                    1.5167
                    0.3948];
            end
            
            
            % run for the first round
            [timelist1,trajecties1]=ode23(odetouse,timespan1,ini1);
            % change the initial and run for the second round
            ini2=trajecties1(end,:);
            addloc=env_dim+(setdiff(secondround_species,firstround_species));
            ini2(addloc)=ini2(addloc)+add_mass;
            %[timelist2,trajecties2]=ode23(odetouse,timespan2,ini2);
            %timelist=[timelist1;timelist2];
            %trajecties=[trajecties1;trajecties2];
            
            dynamic_drawpara=[];
            dynamic_drawpara.showtimecourse=0;
            dynamic_drawpara.arrow_num=1;
            dynamic_drawpara.arrow_angle=40;
            if i==1
                dynamic_drawpara.arrow_sefrac=[0.25,0.26];
            else
                dynamic_drawpara.arrow_sefrac=[0.7,0.71];
            end
            dynamic_drawpara.axislim=axisrange;
            dynamic_drawpara.arrow_percent=1;
            dynamic_drawpara.linewidths=[0.5,0.5];
            
            odepara=[];
            odepara.timespan_runode=timespan2;
            odepara.inistates=ini2;
            [finalstate,traject]=Visualize_chemostatdynamic(chemostat_para,optmetabfun_consortia,dynamic_drawpara,odepara);
            
            timelist=[timelist1;traject(:,1)];
            trajecties=[trajecties1;traject(:,2:end)];
            
            % two environments
            
            % before/after
            envcolors=[];
            envcolors(1,:)=sum(species_colors(firstround_species,:),1);
            envcolors(2,:)=sum(species_colors(secondround_species,:),1);
            
            
            env1=trajecties1(end,1:2);
            env2=finalstate(end,1:2);
            scatter(env1(1),env1(2),steady_ms,envcolors(1,:),'filled');
            scatter(env2(1),env2(2),steady_ms,envcolors(2,:),'filled');
            
            
            %%% plot trajectories
            
            
            figure; hold on;
            box on;
            for sid=1:length(optmetabfun_consortia)
                plot(timelist,trajecties(:,env_dim+sid),'linewidth',1,'color',species_colors(sid,:))
            end
            set(gca,'fontsize',40,'xtick',[0,timerun_12(1),sum(timerun_12)],'ytick',20:20:500);
            ylim([0.001,1.05*max(max(trajecties(:,env_dim+(1:species_num))))])
            
            
            
            
        end
        
        % fitness landscape
        figure;
        allenvs=[env1',env2'];
        strategy_list=linspace(0.1,0.7,100);
        for e=1:2
            subplot(2,2,2*e-1);hold on;
            box on
            
            env=allenvs(:,e);
            % fitness landscape under this enviroment
            if e==1
                toshowspecies=firstround_species;
            else
                toshowspecies=secondround_species;
            end
            
            for s=1:3
                fitness_list=nan*strategy_list;
                importerdim=importer_otherE_id(s,1);
                otherdim=importer_otherE_id(s,2);
                for j=1:length(strategy_list)
                    Es=zeros(4,1);
                    Es(importerdim)=strategy_list(j);
                    Es(otherdim)=1-strategy_list(j);
                    tempspecies=optspecies_consortia{1};
                    tempspecies.Es=Es;
                    tempmetabfun=MetabModel(tempspecies);
                    fitness_list(j)=Obtain_Growth_givenenv(env,tempmetabfun);
                    
                end
                plot(strategy_list,fitness_list,'color',species_colors(s,:),'linewidth',fitness_lw);
                if ismember(s,toshowspecies)
                    scatter(optspecies_consortia{s}.Es(importerdim),chemostat_para.d,fitness_ms,species_colors(s,:),'d','filled');
                    
                end
            end
            
            set(gca,'fontsize',20,'xtick',[min(strategy_list),max(strategy_list)],'ytick',0.2:0.2:2)
            xlim([min(strategy_list),max(strategy_list)])
            ax = gca;
            ax.XColor = envcolors(e,:);
            ax.YColor = envcolors(e,:);
        end
        
    end % different D
end