%clear
clc
species_colors=[0,0,1;0,0.8,0;1,0,0];

envdim=3;


MetabModel=@Metab_importandconvert_substituable;

if envdim==2
    
    chemostat_para=[];
    chemostat_para.c_supplys=[1,1]';
    chemostat_para.d=1;
    
    species_general=[];
    species_general.V_import=120;
    species_general.Ks=[0.5,0.5]';
    species_general.k=1;
    species_general.Es=[0.2,0.3,   0.2,0.3]; % first 3: importer, second three: converter   
    species_general.intercell_dim=length(species_general.Ks);
    
    metabfun=MetabModel(species_general);
    
    
    if 1     % for mutual invasion
    E1s=[0.3,0.4; % for nutrient a 
         0.2,0.1;  % for nutrient b
         ];
        
    E2s=[0.2,0.1;
        0.3,0.4];
    end
    
    if 0 % for bistability
        E1s=[0.23,0.4; % for nutrient a 
         0.27,0.1;  % for nutrient b
         ];
        
    E2s=[0.27,0.1;
        0.23,0.4];
    end
    

    
    E3s=[0.25,0.25;
        0.25,0.25];
    
    E1s=E1s/sum(E1s(:));
    E2s=E2s/sum(E2s(:));
    E3s=E3s/sum(E3s(:));
    Es(1,:)=E1s(:)';
    Es(2,:)=E2s(:)';
    Es(3,:)=E3s(:)';
    
    figure;hold on;
    
    h=[];
    for s=1:size(Es)
        drawpara=[];
        drawpara.pnum=100;
        drawpara.ranges=[0,1;0,1];
        drawpara.growthcolor=species_colors(s,:);
        drawpara.fluxcolor=drawpara.growthcolor;
        drawpara.linewidth=1;
        drawpara.showfluxcurve=1;
        drawpara.surfacealpha=1;
        
        species=species_general;
        species.Es=Es(s,:)';
        metabfun=MetabModel(species);
        
        
        Visualize_chemostatsteady(chemostat_para,metabfun,drawpara);
        [steady,fval]=Obtain_Singlespecies_Steady(chemostat_para,metabfun);
        h(s)=scatter(steady(1),steady(2),50,species_colors(s,:),'filled');
    end
    legend(h,strsplit(num2str(1:size(Es))));
    

elseif envdim==3
%%    
chemostat_para=[];     
chemostat_para.c_supplys=[1,1,1]';     
chemostat_para.d=1;          
species_general=[];     
species_general.V_import=50;     
species_general.Ks=[0.5,0.5,0.5]';     
species_general.k=1;     
species_general.Es=nan;    
species_general.intercell_dim=length(species_general.Ks);


E1s=[0.15,0.4; % for nutrient a           
    0.1, 0.05;  % for nutrient b 
    0.25,0.05];



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
    species_metabfuns=[];
    envs_bys=[];
    steadys=[];
    for s=1:size(Es)
        drawpara=[];
        drawpara.pnum=100;
        drawpara.ranges=[0.1,0.4;
            0.1,0.4;
            0.1,0.4;];
        drawpara.growthcolor=species_colors(s,:);
        drawpara.fluxcolor=drawpara.growthcolor;
        drawpara.linewidth=1;
        drawpara.showfluxcurve=1;
        drawpara.surfacealpha=0.3;
        
        species=species_general;
        species.Es=Es(s,:)';
        metabfun=MetabModel(species);
        
        Visualize_chemostatsteady(chemostat_para,metabfun,drawpara);
        [steady,fval]=Obtain_Singlespecies_Steady(chemostat_para,metabfun);
        
        envs_bys(s,:)=steady(1:3);
        species_metabfuns{s}=metabfun;
        
        h(s)=scatter3(steady(1),steady(2),steady(3),50,species_colors(s,:),'filled');
        steadys(s,1:length(steady))=steady;
    end
    hold on;
    
    set(gca,'fontsize',20,'xtick',[],'ytick',[],'ztick',[]);
    axis([0.1 0.4 0.1 0.4 0.1 0.4]) 
    view([142.8,18]);
    
 
    iniv=[ 0.2954
    0.1873
    0.2860
    0.0151
    0.0640
    0.1945
    1.9226
    1.2568
    4.1956
    4.6514
    1.5319
    1.8487
    1.7211
    3.1982
    1.9292];


    odetouse=@(t,x)ODE_onechemostat(x,chemostat_para,species_metabfuns);
    [timelist,trajecties]=ode23(odetouse,[1,143.65],iniv);
   
      v = VideoWriter('nutirentspace.mp4','MPEG-4');
      open(v);
     for t=1:length(timelist)
         if 1%mod(t,3)==1
            plot3(trajecties(1:t,1),trajecties(1:t,2),trajecties(1:t,3),'k','linewidth',3);
            %scatter3(trajecties(t,1),trajecties(t,2),trajecties(t,3),10,'k','filled')
            title(['T=',num2str(t)])
           frame = getframe;
            writeVideo(v,frame);
            
         end
     end
     
    close(v);
    
     
     figure;
    v = VideoWriter('fitnesslandscape.mp4','MPEG-4');
    open(v);
    clear M
    count=0;
    for t=1:length(timelist)
        
        if 1%mod(t,1)==1
            env=trajecties(t,1:3);
            glist=[];
            clf
            for ss=1:size(Es)% species
                [growthrate,withincellsteady,fval]=Obtain_Growth_givenenv(env,species_metabfuns{ss});
                glist(ss)=growthrate;
                bh=bar(ss,growthrate);hold on;
                bh.FaceColor=species_colors(ss,:);
                bh.LineStyle='none';
            end

            xlabel('Species')
            ylabel('Fitness')
            title(['T=',num2str(t)])
            axis([0,4,0.92,1.1]);
            set(gca,'fontsize',22,'xtick',[], 'ytick',[1]);
             count=count+1;
            frame = getframe;
            writeVideo(v,frame);
        end
    end
    close(v);
    
    
    if 0
        drawpara=[];
        odepara=[];
        %odepara.inistates=finalstate
        odepara.timespan_runode=[1,1000];
        %odepara.inistates=finalstate;
        finalstate=Visualize_chemostatdynamic(chemostat_para,species_metabfuns,drawpara,odepara)
        legend(strsplit(num2str(1:size(Es))))

        % the fitness of each species under the envionment they created
    
    end
    
       fitness_s_ss=[];
    for s=1:size(Es)
        env=envs_bys(s,:);
        for ss=1:size(Es)% species
            [growthrate,withincellsteady,fval]=Obtain_Growth_givenenv(env,species_metabfuns{ss});
            fitness_s_ss(s,ss)=growthrate;
        end
    end
    
%    figure;
%    Easypcolor(fitness_s_ss,'','');
%    caxis([min(fitness_s_ss(:)),max(fitness_s_ss(:))])
%    colorbar

    
%     % two chemostats
%     leakrate=0.00001;
%     chemostatnum=3;
%     chemostat_network=ones(chemostatnum);
%             % a "connect to nearest neighbor" pattern
%             for c=1:(chemostatnum-1)
%                 chemostat_network(c,c+1)=1;
%                 chemostat_network(c+1,c)=1;
%             end
%             dxdt_noleak=@(t,x)ODE_multichemostat(x,chemostat_para,species_metabfuns,0,chemostat_network);
%             dxdt_withleak=@(t,x)ODE_multichemostat(x,chemostat_para,species_metabfuns,leakrate,chemostat_network);
%             
%             
%             ini=[envs_bys(1,:),1,0,0,ones(1,9), envs_bys(2,:),0,1,0,ones(1,9), envs_bys(3,:),0,0,1,ones(1,9)]';
%             %ini=[testini_1;testini_2];
%             
%             % run to steady state with zero leakage, then start to connect
%             [timelist1,trajecties1]=ode23(dxdt_noleak,[0,1000],ini);
%             [timelist2,trajecties2]=ode23(dxdt_withleak,timelist1(end)+[0,500000],trajecties1(end,:)');
%             timelist=[timelist1;timelist2];
%             trajecties=[trajecties1;trajecties2];
%             
%              env_byc=[];
%              spe_byc=[];
%             figure;
%             hold on;
%             for c=1:chemostatnum
%                 envdim=(c-1)*15+(1:length(chemostat_para.c_supplys));
%                 speciesdim=(c-1)*15+length(chemostat_para.c_supplys)+(1:length(species_metabfuns));
%                 env_byc(c,:)=trajecties(end,envdim);
%                 spe_byc(c,:)=trajecties(end,speciesdim);
%                 
%                 subplot(3,1,c); hold on;
%                 for s=1:3
%                     plot(timelist,trajecties(:,speciesdim(s)),'color',species_colors(s,:));
%                 end
%                 ylim([0,0.1])
%             end
%             
            
            %figure;
            %Easypcolor(spe_byc,'','')
%figure;
%plot(trajecties)
end






