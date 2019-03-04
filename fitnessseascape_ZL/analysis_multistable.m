clear
clc
% for analyzing bistability, multistability

% choices to make
modelid=1; % 1: Min, 2: 1/(sum(1/Import))
progressive_eachround=5;
max_speciesnum=progressive_eachround*2;% maximal number of species possiblely to show
species_colors=distinguishable_colors(max_speciesnum, {'w','k'});

startstrategies=[0.35,0.65];
fs=24;


% tasks 
twospecies_bistability=0;
opt_species=1;
progressiveevo=0;

opt_show2=1;


%%
if modelid==1
    MetabModel=@Metab_MinImport;
    
    % parameter for chemostats
    chemostat_para=[];
    chemostat_para.c_supplys=[1,1]';
    chemostat_para.d=1;
    
    % parameter for species
    species_general=[];
    species_general.intercell_dim=0;
    species_general.Ks=[0.5,0.5]';
    species_general.Gamma=10;
    species_general.Es=[nan,nan]';
elseif modelid==2
    MetabModel=@Metab_RwImport;
    
    % parameter for chemostats
    chemostat_para=[];
    chemostat_para.c_supplys=[1,1]';
    chemostat_para.d=1;
    
    % parameter for species
    species_general=[];
    species_general.intercell_dim=0;
    species_general.Ks=[0.5,0.5]';
    species_general.Gamma=25;
    species_general.Es=[nan,nan]';
end

metabfun_general=MetabModel(species_general);

%%
% first, the two strategies form bistability
if twospecies_bistability
    colortouse=species_colors(1:2,:);
    
    strategy_legend=[];
    fa_list=startstrategies; % 0.3,0.7, 9; 0.1,0.9, 9; 0.1,0.9, 17;
    for s=1:length(fa_list)
        E1=fa_list(s);
        species_one=species_general;
        species_one.Es=[E1;1-E1];
        species{s}=species_one;
        strategy_legend{s}=['f_a=',num2str(E1,'%0.2f')];
    end
    species_num=length(species);

    drawpara_steady=[];
    drawpara_steady.pnum=1000;
    drawpara_steady.ranges=[0.01,1;
        0.01,1]*0.3;
    drawpara_steady.linewidth=2;
    drawpara_steady.showfluxcurve=0;

    species_metabfuns=[];
    steadyenv_byspecies=zeros(species_num,2);

    figure;
    hold on;
    h=[];
    % GC
    for s=1:species_num
        metabfun=MetabModel(species{s});
        species_metabfuns{s}=metabfun;
        drawpara_steady.growthcolor=colortouse(s,:);%localpoints(k,:);
        Visualize_chemostatsteady(chemostat_para,metabfun,drawpara_steady);

        [steady,fval]=Obtain_Singlespecies_Steady(chemostat_para,metabfun);
        h(s)=scatter(steady(1),steady(2),150,colortouse(s,:),'filled');
        steadyenv_byspecies(s,:)=steady(1:2);
    end
    
    % trajectory
    drawpara_dynamic=[];
    drawpara_dynamic.linewidths=[1,1];
    drawpara_dynamic.speciescolor=colortouse;
    drawpara_dynamic.showtimecourse=0;
    drawpara_dynamic.arrow_num=1;
    drawpara_dynamic.arrow_sefrac=[0.1,0.3];
    drawpara_dynamic.arrow_percent=1/40; % percent of the arrow for each axis;
    drawpara_dynamic.axislim=drawpara_steady.ranges;
    drawpara_dynamic.timespan_runode=[0,100];
    drawpara_dynamic.arrow_angle=20;

    odepara=[];
    odepara.timespan_runode=[0,100];
    odepara.inistates=[0.15,0.05,8,1
        0.05,0.1,10,15];
    Visualize_chemostatdynamic(chemostat_para,species_metabfuns,drawpara_dynamic,odepara);
    axis square  
    %title('Two species bistability');
    
    %legend(h,strategy_legend)
    set(gca,'fontsize',24,'xtick',[0.1 0.3],'ytick',[0.1,0.3])
    box on;
    
    %%%% the fitness landscape in the environment formed by one species
    figure;hold on;
    h=[];
    flist=linspace(0.3,0.7,1000);
    for se=1:size(steadyenv_byspecies,1)
        env=steadyenv_byspecies(se,:)';
        glist=nan*flist;
        for i=1:length(flist)
             tempspecies=species_general;
             tempspecies.Es=[flist(i),1-flist(i)]';
             tempmetabfun=MetabModel(tempspecies);
             glist(i)=Obtain_Growth_givenenv(env,tempmetabfun);
        end
        plot(flist,glist,'linewidth',1.5,'color',colortouse(se,:));
        h(se)=scatter(species{se}.Es(1),chemostat_para.d,200,colortouse(se,:),'filled');
    end
    axis([min(flist) max(flist) 0.6 1.1])
    set(gca,'fontsize',24,'xtick',[0 0.3],'ytick',[0,0.3])
    %xlabel('f_a')
    %ylabel('g')
    %title('fitness landscape')
    %legend(h,strategy_legend)
    axis square
    box on
    
    
end

%% optimal species and other properties for this system
if opt_species
    colors=[     1     0     0;
             0,0,0
     0     0     1
];
    ranges=[0.01,1;
        0.01,1]*0.3;
    pnumbydim=[20,20];
    clistpnum=101;
    fs=24;
    
    [nutreint_boarders_m,nutreint_values_m,optgrowthrates_m,opt_Es,consumptions_m,withincellv_m]=Nutrientspace_optspecies(species_general,MetabModel,ranges,pnumbydim);
    resultsfromoptspecies={nutreint_boarders_m,nutreint_values_m,optgrowthrates_m,opt_Es,consumptions_m,withincellv_m};
    
    % optimal growth rate
    figure;hold on;
    colormap(bone)
    ph=pcolor(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optgrowthrates_m);
    shading interp
    ph.LineStyle='none';
    ch=colorbar;
    set(gca,'fontsize',fs);
    
    % optimal growth contours and the vecors on it
    % get optimal growth contour and the optimal strategies along it 
    clist=linspace(ranges(1,1),ranges(1,2),clistpnum);
    Es_optGC=nan*zeros(2,clistpnum);

    unknowndim=2;    
    if isfield(metabfun_general,'opt_GC')
        opt_GC=metabfun_general.opt_GC;
        unknowndim=opt_GC.unkowndim;
        GCfun=@(x)opt_GC.GC(x,chemostat_para.d);
        optGC_list=nan*clist;
        for i=1:length(clist)
            otherv=GCfun(clist(i));
            if otherv>0
                optGC_list(i)=otherv;
                % obtian the enzyme strategy
                Es_optGC(:,i)=opt_GC.Es(clist(i),chemostat_para.d);
            end
        end
    end
    
    if unknowndim==2
        x=clist;
        y=optGC_list;
    else
        y=clist;
        x=optGC_list;
    end
    
    if 0
        freezeColors
        colormap('default')
        caxis([0,1])
        z = zeros(size(clist));
        col = Es_optGC(1,:);
        surface([x;x],[y;y],[z;z],[col;col],'facecol','no', 'edgecol','interp', 'linew',4);
    end
    
      % for some points along the contour, generate the enzyme strategy &
    % species, and their corresponding flux balance curve restriction
    if 0
        drawpara_steady=[];
        drawpara_steady.pnum=1000;
        drawpara_steady.ranges=[0.01,1;
            0.01,1]*0.3;
        drawpara_steady.linewidth=2.5;
        drawpara_steady.showfluxcurve=0;


        selectaroundenv=[0.1,0.1248,0.2];
        strategys=[];
        for i=[];%1:3%1:length(selectaroundenv)
            [minv,minloc]=min(abs(x-selectaroundenv(i)));
            env=[x(minloc);y(minloc)];
            strategys(:,end+1)=Es_optGC(:,minloc);

            species_on_optGC=species_general;
            species_on_optGC.Es=strategys(:,i);
            metabfun=MetabModel(species_on_optGC);
            drawpara_steady.growthcolor=colors(i,:);%localpoints(k,:);
            if i==3
                h=Visualize_chemostatsteady(chemostat_para,metabfun,drawpara_steady);
            

            [FBv_givenIa,onelargesupplypoint,fval]=Obtain_supplyvector(metabfun,env);

            plot([env(1),onelargesupplypoint(1)],[env(2),onelargesupplypoint(2)],'-- k','linewidth',1)
        
            end
        end
    end
    axis([ranges(1,1),ranges(1,2)-0.025,ranges(2,1),ranges(2,2)-0.025])
    set(gca,'xtick',[0.1,0.2],'ytick',[0.1,0.2]);
    set(ch,'ytick',[0,1])
    axis square
    xlabel('')
    ylabel('')
end


% given enviornment and species, what's the supply that stablize this
% enviornment






% 
% fitness_species_inenvs=zeros(species_num);
% for s=1:length(species)
%     for ss=1:length(species)
%         fitness_species_inenvs(s,ss)=Obtain_Growth_givenenv(steadyenv_byspecies(s,:),species_metabfuns{ss});
%     end
% end
% 
% flist=linspace(0.25,0.75,101);
% figure;hold on;
% h=[];
% invadenum=zeros(species_num,1);
% for s=1:species_num
%     invadenum(s)=sum(fitness_species_inenvs(s,:)>fitness_species_inenvs(s,s));
%     s
%     glist=0*flist;
%     for i=1:length(flist)
%         tempspecies=species_general;
%         tempspecies.Es=[flist(i),1-flist(i)]';
%         tempmetabfun=MetabModel(tempspecies);
%         glist(i)=Obtain_Growth_givenenv(steadyenv_byspecies(s,:),tempmetabfun);
%     end
%     plot(flist,glist,'color',species_colors(s,:));
%     %plot(fa_list,fitness_species_inenvs(s,:),'color',species_colors(s,:),'linewidth',1);
%     h(s)=scatter(fa_list(s),fitness_species_inenvs(s,s),100,species_colors(s,:),'filled');
% end
% legend(h,strategy_legend);
% set(gca,'fontsize',22);
% 
% xlim([0.25 0.75])
% how many species can invade it
% %%

if progressiveevo
pnum=501;

flist=linspace(startstrategies(1),startstrategies(2),pnum);

for s=1:2
    currentfa=startstrategies(s);
    for r=1:progressive_eachround
        currentfa
        currentspecies=species_general;
        currentspecies.Es=[currentfa,1-currentfa]';
        currentmetabfun=MetabModel(currentspecies);
        % the environment constructed by the current species
        [steady,fval]=Obtain_Singlespecies_Steady(chemostat_para,currentmetabfun);
        currentenv=steady(1:2);
        
        % fitness landscape under current envirnoment
        
        glist=0*flist;
        for i=1:pnum
             tempspecies=species_general;
                tempspecies.Es=[flist(i),1-flist(i)]';
                tempmetabfun=MetabModel(tempspecies);
            glist(i)=Obtain_Growth_givenenv(currentenv,tempmetabfun);
        end
        
        colorid=r*2+(s-2);
        
        plot(flist,glist,'color',species_colors(colorid,:),'linewidth',2); hold on;
        scatter(currentfa,chemostat_para.d,100,species_colors(colorid,:),'d','filled');
        
        % choose the next strategy
        [maxv,maxloc]=max(glist);
        currentfa=flist(maxloc);
        
    end
end
% and the optimal one
if 1
currentfa=0.5;
optcolor=[0,0,0];
currentspecies=species_general;
        currentspecies.Es=[currentfa,1-currentfa]';
        currentmetabfun=MetabModel(currentspecies);
        % the environment constructed by the current species
        [steady,fval]=Obtain_Singlespecies_Steady(chemostat_para,currentmetabfun);
        currentenv=steady(1:2);
        
        % fitness landscape under current envirnoment
        
        glist=0*flist;
        for i=1:pnum
             tempspecies=species_general;
                tempspecies.Es=[flist(i),1-flist(i)]';
                tempmetabfun=MetabModel(tempspecies);
            glist(i)=Obtain_Growth_givenenv(currentenv,tempmetabfun);
        end
        
        plot(flist,glist,'color',optcolor,'linewidth',2); hold on;
        scatter(currentfa,chemostat_para.d,100,optcolor,'d','filled'); 

      
end
set(gca,'fontsize',fs,'xtick',[0.35,0.5,0.65],'ytick',[1]);
        %xlabel('Strategies')
        %ylabel('Growth rate')
        axis([0.35 0.65 0.9 1.08])
end
  

if opt_show2==2
    
end

% question for next: space factor