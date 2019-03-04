%clear
clc
paren = @(x, varargin) x(varargin{:});

%%%
getopot=0;



%%%
env_dim=2;
max_speciesnum=3;
species_colors=[0,0,1;0,0.8,0;1,0,0];

MetabModel=@Metab_crossfeeding;

% parameter for chemostats
chemostat_para=[];
chemostat_para.c_supplys=[1.8,0]'; % supply 1.5
chemostat_para.d=0.6; %0.3, 0.6

% parameter for species
species_general=[];
species_general.intercell_dim=1;


species_general.Es= [0.2,0.1,0.6,0.1];%[0.2,0.3,0.4,0.1];

species_general.V_1=5;
species_general.V_2=1;
species_general.k=8;
species_general.V_4=10;

species_general.Ks=[0.5, 0.5, 0.5, 0.1, 0.5, 15,10];


%%
% three species
fprintf('G\n');
%G
species_special=species_general;
species_special.Es=[0.5,0.5,0,0];
metabfun=MetabModel(species_special);
optEs=metabfun.analytical_optE_givenDW;% (D,W);
minSW=metabfun.analytical_minSW_givenDW;
optG_Es=optEs(chemostat_para.d,0);
minSW(chemostat_para.d)

species_optG=species_special;
species_optG.Es=optG_Es;

%


% generate some representititve species
% S1
fprintf('S1\n');
species_special.Es=[0.5,0,0.5,0];
metabfun=MetabModel(species_special);
optEs=metabfun.analytical_optE_givenDW;% (D,W);
minSW=metabfun.analytical_minSW_givenDW;
optS1_Es=optEs(chemostat_para.d,0);% note that it is optimized for 0
minSW(chemostat_para.d,0)

species_goodS1=species_special;
species_goodS1.Es=optS1_Es;

% S2
fprintf('S2\n');
species_special.Es=[0,0.5,0,0.5];
metabfun=MetabModel(species_special);
optEs=metabfun.analytical_optE_givenDW;% (D,W);
minSW=metabfun.analytical_minSW_givenDW;
optS2_Es=optEs(chemostat_para.d,0);
minSW(chemostat_para.d)

species_optS2=species_special;
species_optS2.Es=optS2_Es;
metabfun=MetabModel(species_optS2);


%%
%gcfun(1,chemostat_para.d)
%
importer_otherdims=[1,] % 
% get the optimal species and corresponding growth rates
ranges=[0,1.5;%1.5
    0,0.3];%0.25
pnumbydim=[100,100]*2;
drawpara=[];
showprogress=1;

species_consortia={species_goodS1,species_optS2,species_optG};
S1id=1;
S2id=2;
Gid=3;

importer_otherE_id=[1,3;
    4,2;
    1,2;]
if ~exist('optEs_all','var')
    optgrowth_bys=zeros(pnumbydim(1),pnumbydim(2),3);
    optEs_bys=zeros(pnumbydim(1),pnumbydim(2),4,3);
    consumptions_bys=zeros(pnumbydim(1),pnumbydim(2),2,3);
    withincellv_bys=zeros(pnumbydim(1),pnumbydim(2),3);
    for s=1:3
        speciespara=species_consortia{s};
        [nutreint_boarders_m,nutreint_values_m,optgrowth_bys(:,:,s),optEs_bys(:,:,:,s),consumptions_bys(:,:,:,s),withincellv_bys(:,:,s)]=Nutrientspace_optspecies(speciespara,MetabModel,ranges,pnumbydim,showprogress,drawpara);
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


% the growth rates

    figure;
    hold on;
    colormap('default');
    colormatrix=nan*zeros(pnumbydim(1),pnumbydim(1),3);
    for i=1:pnumbydim(1)
        for j=1:pnumbydim(2)
            if ~isnan(optsid(i,j))
            colormatrix(i,j,:)=species_colors(optsid(i,j),:);
            end
        end
    end
    ph=surface(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optsid,colormatrix);
    ph.LineStyle='none';
    xlabel('S')
    ylabel('W')
    contour3(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optgrowth_all+100,'k','linewidth',1)

    
    axis square
    set(gca,'fontsize',24,'xtick',[0.5,1.5],'ytick',[0,0.2]);
    freezeColors
    colortable=flipud(colormap(hot));
    colormap(colortable)
    contour(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optgrowth_all,'linewidth',2,'ShowText','on');
    


figure;
    hold on;
    colormap('bone');
    ph=pcolor(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optgrowth_all);

% the strategy
figure;
for i=1:4
    subplot(2,2,i); hold on;
    colormap('default')
    ch=colorbar;
    ph=pcolor(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optEs_all(:,:,i));
    ph.LineStyle='none';
    
    caxis([0 1])
    set(ch,'ytick',[0,1])
    set(gca,'fontsize',15,'xtick',[0.5,1.5],'ytick',[0,0.2]);
    contour(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optgrowth_all,'r','linewidth',1);
    axis square
end


% the very acturate "optimal contour" based on the above knowledge
Dlist=[0.4,0.6];
Wlistnum=[15,30];

% not yet fully optimized
ini_mass=0.1;
add_mass=0.1;
ini_ms=[ini_mass,0,0];
ini_previous=[ones(1,env_dim), ini_ms,0.1*ones(1,species_general.intercell_dim*length(ini_ms))]';
for i=1:length(Dlist)
    Wlist=linspace(ranges(2,1),ranges(2,2),Wlistnum(i));

    D=Dlist(i);
    chemostat_para.d=D;
    
    %%%%%%%%% estblish the optimal species in every time
    metabfun_consortia=[];
    for s=1:max_speciesnum
        tempspecies=species_consortia{s};
        tempmetabfun=MetabModel(tempspecies);
        optEs=tempmetabfun.analytical_optE_givenDW(chemostat_para.d,0);% (D,W);

        % the optimal one
        opt_thiss=tempspecies;
        opt_thiss.Es=optEs;
        optmetabfun=MetabModel(opt_thiss);

                
        species_consortia{s}=opt_thiss;
        metabfun_consortia{s}=optmetabfun;
    end %

    
    
    figure;
    hold on;
    lw=2;
    
    S2_minWvalue=metabfun_consortia{S2id}.analytical_minSW_givenDW(D);
    G_minSvalue=metabfun_consortia{Gid}.analytical_minSW_givenDW(D);
    S1GC_SgivenW=@(W)metabfun_consortia{S1id}.analytical_minSW_givenDW(D,W);
    
    % estimate the intersection point
    
    S1GC_Slist=0*Wlist;
    for j=1:length(Wlist)
        S1GC_Slist(j)=S1GC_SgivenW(Wlist(j));
    end
    nearsolutionloc=find(S1GC_Slist(1:(end-1))<G_minSvalue & S1GC_Slist(2:end)>G_minSvalue);
    
    tosolvefun=@(w)S1GC_SgivenW(w)-G_minSvalue;
    S1_G_intersect_W=fzero(tosolvefun,Wlist(nearsolutionloc:(nearsolutionloc+1)));
    S1_G_intersect_S=G_minSvalue;
    
    S1_S2_intersect_W=S2_minWvalue;
    S1_S2_intersect_S=S1GC_SgivenW(S2_minWvalue);
    
    % whichever the smallest W point will show
    [minWend,GorS2]=min([S1_G_intersect_W,S1_S2_intersect_W]);
    boarder_Wlist=linspace(ranges(2,1),minWend,length(Wlist));
    boarder_S1GC_Slist=0*boarder_Wlist;
    optE1_list=0*boarder_Wlist;
    % where to show the supply vector
    supplyshowpoint=3;
    ptoshow=round(linspace(1,Wlistnum(i),supplyshowpoint+1));ptoshow(1)=[];
    vector_endx=[];
    for j=1:length(boarder_Wlist)
        boarder_S1GC_Slist(j)=S1GC_SgivenW(boarder_Wlist(j));
        % optimal strategy in this environment
        optS1_Es_thisW=metabfun_consortia{S1id}.analytical_optE_givenDW(chemostat_para.d,boarder_Wlist(j));
        optE1_list(j)=optS1_Es_thisW(1);
        
        if ismember(j,ptoshow)
            env=[boarder_S1GC_Slist(j),boarder_Wlist(j)]';
            optspecies_thisenv=species_general;
            optspecies_thisenv.Es=optS1_Es_thisW;
            tempmetabfun=MetabModel(optspecies_thisenv);
            FBv_givenIa=Obtain_supplyvector(tempmetabfun,env);
            Ib_givenIa=@(x)paren(FBv_givenIa(x),2);
            vector_endx(end+1)=fzero(Ib_givenIa,[0,chemostat_para.c_supplys(1)*1.1]);
        end
    end
    S1copoint=[boarder_S1GC_Slist(end),minWend];
    
    % plot contours that would show
    % G
    if GorS2==1 % G coexist with S1
        optpair=[S1id,Gid];
        plot([G_minSvalue,G_minSvalue],[S1_G_intersect_W,S2_minWvalue],'color',species_colors(Gid,:), 'linewidth',lw);
    else
        optpair=[S1id,S2id];
    end
    % S2
    plot([0,min(G_minSvalue,S1_S2_intersect_S)],[S2_minWvalue,S2_minWvalue],'color',species_colors(S2id,:), 'linewidth',lw);
    % S1
    ms=50;
    [colortable,colorids]= Colormap_expression(optE1_list-min(optE1_list),[0,0,0;0.8,0.8,1;0,0,1],0);
    scatter(boarder_S1GC_Slist,boarder_Wlist,ms,colortable(colorids,:),'s','filled')
    %plot(boarder_S1GC_Slist,boarder_Wlist,'color',species_colors(S1id,:), 'linewidth',lw)
    
    % plot the supply vector if needed
   for p=[];%1:length(ptoshow)
       j=ptoshow(p);
       env=[boarder_S1GC_Slist(j),boarder_Wlist(j)]';
       plot([env(1),vector_endx(p)],[env(2),0],'-- k','linewidth',0.8);
   end
    
    %  extract the optimal S1 strategy at the intersection point
    optS1_Es=metabfun_consortia{S1id}.analytical_optE_givenDW(chemostat_para.d,minWend);% (D,W)
    species_optS1=species_special;
    species_optS1.Es=optS1_Es;
    metabfun_optS1=MetabModel(species_optS1);
    
    species_consortia{S1id}=species_optS1;
    metabfun_consortia{S1id}=metabfun_optS1;

    shownutrientspace=1;
    
  
    

    
    % growth contour for each individual species, and the steady state they
    % created
    
    % dynamical simulations
    % for S1-S2 optimal coexistence: first run S1, then add S2
    % for S1-G optimal coexistence first run S1 and S2, then add G
    
    species_num=3;
    if GorS2==2 % S1 with S2
        add_sid=S2id;
        runtime1=50;
        runtime2=100;
    else
        add_sid=Gid;
        ini_previous=[ 1.4359    0.2500    0.3641    0.1141   0    0.4203    1.5151    0.3947];
        runtime1=50;
        runtime2=200;
    end
    
    
    odetouse=@(t,x)ODE_onechemostat(x,chemostat_para,metabfun_consortia);
    
    timespan1=[0,runtime1];
    timespan2=runtime1+[0,runtime2];
    ini1=ini_previous;
    % run for the first round
    [timelist1,trajecties1]=ode23(odetouse,timespan1,ini1);
    % change the initial and run for the second round
    ini2=trajecties1(end,:);
    addloc=env_dim+add_sid;
    ini2(addloc)=ini2(addloc)+add_mass;
    [timelist2,trajecties2]=ode23(odetouse,timespan2,ini2);
    timelist=[timelist1;timelist2];
    trajecties=[trajecties1;trajecties2];
    

      if shownutrientspace
            hold on;
            drawpara=[];
            drawpara.ranges=[1.1,1.6;0, 0.3];
            drawpara.showfluxcurve=0;
            drawpara.pnum=1000;
           % growth contour and self-steady of the optimal pair right at the conjuction point
           
           slist=1:3;
           for s=slist
                
               ms=150;
               sid=s;%optpair(s);
               drawpara.growthcolor=species_colors(sid,:);
               drawpara.linewidth=2;
               Visualize_chemostatsteady(chemostat_para,metabfun_consortia{sid},drawpara);
               hold on;
               if s~=2
                   steadyv=Obtain_Singlespecies_Steady(chemostat_para,metabfun_consortia{sid});
                   scatter(steadyv(1),steadyv(2),ms,species_colors(sid,:),'filled')
               end
           end
            
            axis square
            %axis([drawpara.ranges(1,1),drawpara.ranges(1,2),drawpara.ranges(2,1),drawpara.ranges(2,2)])
            %set(gca,'fontsize',24,'xtick',[1.1,1.6],'ytick',[0,0.3]);
      
            
            drawpara.showtimecourse=0;
            drawpara.axislim=drawpara.ranges;
            drawpara.arrow_num=1;
            odepara=[];
            odepara.inistates=ini2;
            odepara.timespan_runode=timespan2;
            Visualize_chemostatdynamic(chemostat_para,metabfun_consortia,drawpara,odepara)
            xlabel('')
            ylabel('')
      end
          axis square;
    axis([ranges(1,1),ranges(1,2),ranges(2,1),ranges(2,2)])
    set(gca,'fontsize',24,'xtick',[0.5,1.5],'ytick',[0,0.2]);
    
    
    
    figure;
    subplot(2,1,1);hold on;
    for sid=1:species_num
        plot(timelist,trajecties(:,env_dim+sid),'linewidth',2,'color',species_colors(sid,:))
    end
    title(D)
    set(gca,'fontsize',22,'xtick',[0,runtime1,runtime1+runtime2],'ytick',[0.5]);
    ylim([0.001,1.05*max(max(trajecties(:,env_dim+(1:species_num))))])
    
    ini_previous=trajecties(end,:);
    
    
    %%%%% fitness landscape in some enviornment 
    
    for envid=1:2
        if envid==1
            givenenv=trajecties1(end,1:2);
        else
            givenenv=S1copoint;
        end
        
        origs=[];
        figure;hold on;
        title(envid)
        %givenenv=trajecties2(end,1:2);
        

        importE_list=linspace(0.01,0.99,300);
        for s=1:species_num
            fitness_list=nan*importE_list;
            for j=1:length(importE_list)
                importerE=importE_list(j);
                Es=zeros(1,4);

                importerdim=importer_otherE_id(s,1);
                otherdim=importer_otherE_id(s,2);

                Es(importerdim)=importerE;
                Es(otherdim)=1-importerE;

                tempspecies=species_general;
                tempspecies.Es=Es;
                tempmetabfun=MetabModel(tempspecies);
                [growthrate,withincellsteady,fval]=Obtain_Growth_givenenv(givenenv,tempmetabfun);
                if growthrate>0 && abs(fval)<0.01
                    fitness_list(j)=growthrate;
                end 
            end
            plot(importE_list,fitness_list,'color',species_colors(s,:),'linewidth',2);

            actualg=Obtain_Growth_givenenv(givenenv,metabfun_consortia{s});
            ms=150;
            scatter(species_consortia{s}.Es(importerdim),actualg,ms,species_colors(s,:),'d','filled');
            %[maxv,maxloc]=max(fitness_list);
            %scatter(importE_list(maxloc),maxv,ms,species_colors(s,:),'o');

            origs(end+1)=actualg;
        end
        
        xlimvalue=[0.1,0.7];
        ylim([min(origs)*0.9,max(origs)*1.05]);
        xlim(xlimvalue)
        set(gca,'fontsize',24);
        set(gca,'xtick',0.5,'ytick',chemostat_para.d)
    end
    
end



%


% figure;hold on;
% ph=pcolor(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optgrowth_all);
% ph.LineStyle='none';
% title('opt G')
% contour(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optgrowth_all,'k');

% figure;
% % plot the consumption vectors along the optimal growth contours
% resultsfromoptspecies=[];
%
% resultsfromoptspecies{1}=nutreint_boarders_m;
% resultsfromoptspecies{2}=nutreint_values_m;
% resultsfromoptspecies{3}=optgrowth_all;
% resultsfromoptspecies{4}=optEs_all;
% resultsfromoptspecies{5}=consumptions_all;
% resultsfromoptspecies{6}=[];
%
% dlist=[0.3,0.7];
% [alongGCs,alongFBs]=Analysis_OptM(metabfun,resultsfromoptspecies,dlist,[]);
% for i=1:length(alongGCs)
%     alongGC=alongGCs{i};
%     envs=alongGC.envs;
%     Es_optGC=alongGC.Es_optGC;
%
%
%     % find some points
%     for j=1:1:(size(envs,2)-1)
%         env=envs(:,j);
%         Es=Es_optGC(:,j);
%
%         optspecies_thisenv=species_general;
%         optspecies_thisenv.Es=Es;
%
%        optmetabfun=MetabModel(optspecies_thisenv);
%        [FBv_givenIa,onelargesupplypoint]=Obtain_supplyvector(optmetabfun,env);% [FBv_givenIa,onelargesupplypoint,fval]=Obtain_supplyvector(metabfun,env)
%         plot([env(1),onelargesupplypoint(1)],[env(2),onelargesupplypoint(2)],'g');hold on;
%
%
%     end
% end
% axis([ranges(1,1),ranges(1,2),ranges(2,1),ranges(2,2)])
% %
%
% %%
% % run simulation of three optimized species at the given points:
% species_consortia={species_goodS1,species_optS2,species_optG};
% species_num=length(species_consortia);
% env_dim=2;
%
% timespan_runode=[0,500];
% ini=[0,0, 0.01*ones(1,species_num),0.1*rand(1,species_general.intercell_dim*species_num)]';
%
%
% species_metabfuns=[];
% for s=1:species_num
%     species_metabfuns{s}=MetabModel(species_consortia{s});
% end
%
% odetouse=@(t,x)ODE_onechemostat(x,chemostat_para,species_metabfuns);
% [timelist1,trajecties1]=ode23(odetouse,timespan_runode,ini);
% %figure;
% %plot(timelist1,trajecties1,'linewidth',2);
% %legend({'S','W','S1','S2','P_s1','P_s2'})
% trajecties1(end,1:2)
% nextini=trajecties1(end,:);
% nextini(4)=0.1;% for S2
% [timelist2,trajecties2]=ode23(odetouse,timespan_runode+max(timespan_runode),nextini);
% timelist=[timelist1;timelist2];
% trajecties=[trajecties1;trajecties2];
%
% plot(trajecties1(end,1),trajecties1(end,2),'r o')
% plot(trajecties2(end,1),trajecties2(end,2),'m o')
%
%
% figure;hold on;
% for s=1:species_num
%     plot(timelist,trajecties(:,env_dim+s),'linewidth',1,'color',species_colors(s,:))
% end
% set(gca,'fontsize',20);
% xlabel('Time');
% ylabel('Density');
% legend({'S1','S2','G'})




