function [alongGCs,alongFBs]=Analysis_OptM(metabfun,resultsfromoptspecies,Dlist,supplies,drawpara,clist)
% analyze the optimal growth matrix and strategies for growth contour (different D), flux
% balance curve (different supply);

nutreint_boarders_m=resultsfromoptspecies{1};
nutreint_values_m=resultsfromoptspecies{2};
optgrowthrates_m=resultsfromoptspecies{3};
opt_Es=resultsfromoptspecies{4};
consumptions_m=resultsfromoptspecies{5};
withincellv_m=resultsfromoptspecies{6};

envdim=size(nutreint_boarders_m,3);


ranges=[min(nutreint_boarders_m(1,:)),max(nutreint_boarders_m(1,:));
    min(nutreint_boarders_m(:,1)),max(nutreint_boarders_m(:,1))];

if ~(nargin>5)
    clist=linspace(min(nutreint_boarders_m(:)),max(nutreint_boarders_m(:)),100);
end
clistpnum=length(clist);



% default drawing parameters
colornutrientspace=1;
grwoth_colormap=gray;
plotGC=1; % 1: black, 2: vector; 3: colorful; 0: not to plot
plotGClw=1;
plotFB=1;
plotFBlw=1;
minmaxv_colors=[];
if nargin > 4
    if isfield(drawpara,'colornutrientspace')
        colornutrientspace=drawpara.colornutrientspace;
    end
    if isfield(drawpara,'grwoth_colormap')
        grwoth_colormap=drawpara.grwoth_colormap;
    end
    if isfield(drawpara,'plotGC')
        plotGC=drawpara.plotGC;
    end
    if isfield(drawpara,'plotGCcolor')
        plotGCcolor=drawpara.plotGCcolor; % will be a table if plotGC=3
    end
    if isfield(drawpara,'plotGClw')
        plotGClw=drawpara.plotGClw;
    end
    
    if isfield(drawpara,'plotFB')
        plotFB=drawpara.plotFB;
    end
    if isfield(drawpara,'plotFBlw')
        plotFBlw=drawpara.plotFBlw;
    end
    if isfield(drawpara,'minmaxv_colors')
        minmaxv_colors=drawpara.minmaxv_colors;
    end
end


if colornutrientspace
    colormap(grwoth_colormap); hold on;
    ph=pcolor(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optgrowthrates_m);
    ph.LineStyle='none';
    
    set(gca,'fontsize',20);
    xlabel('c_a');
    ylabel('c_b');
    
    
end

% obtain growth contours, and the environment with strategies, steady
% states linked to it
alongGCs=[];

for i=1:length(Dlist)
    D=Dlist(i);
    
    alongGC=[];
    alongGC.d=D;
    % if there is an analytical solution
    
    D_ijs=round(contourc(optgrowthrates_m,[D,D]));
        D_ijs(:,1)=[];
    if isfield(metabfun,'opt_GC')
        Es_optGC=nan*zeros(2,clistpnum);
        
        opt_GC=metabfun.opt_GC;
        unknowndim=opt_GC.unkowndim;
        GCfun=@(x)opt_GC.GC(x,D);
        optGC_list=nan*clist;
        for j=1:length(clist)
            otherv=GCfun(clist(j));
            if otherv>0
                optGC_list(j)=otherv;
                % obtian the enzyme strategy
                Es_optGC(:,j)=opt_GC.Es(clist(j),D);
            end
        end
        if unknowndim==2
            envs=[clist;
                optGC_list];
        else
            envs=[optGC_list;
                clist;];
        end
        todel=(isnan(optGC_list) );
        envs(:,todel)=[];
        Es_optGC(:,todel)=[];
        
    else % find it by contour
        %contour(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optgrowthrates_m,[D,D],'g')
        
        envs=[];
        Es_optGC=[];
        for j=1:size(D_ijs,2)
            envs(:,j)=nutreint_values_m(D_ijs(2,j),D_ijs(1,j),:);
            Es_optGC(:,j)=opt_Es(D_ijs(2,j),D_ijs(1,j),:);
        end
    end
    % remove the nan li
    alongGC.envs=envs;
    alongGC.Es_optGC=Es_optGC;
    alongGC.D_ijs=flipud(D_ijs);
    
    if plotGC>0
        colorcode=[];
        if plotGC==1
            colorcode=[1,0,0];
        elseif plotGC==2
            colorcode=plotGCcolor(i,:);
        end
        
        if ~isempty(colorcode)
            
            plot(envs(1,:),envs(2,:),'color',colorcode,'linewidth',plotGClw);

        elseif plotGC==3
            % plotGCcolor is a table for mapping colors
            
            x=envs(1,:);
            y=envs(2,:);
            z = zeros(size(x));
            
            colorids=[];
            colorm=[];
            if isempty(minmaxv_colors)
                maxv=max(Es_optGC(1,:));
                minv=min(Es_optGC(1,:));
            else
                maxv=minmaxv_colors(2);
                minv=minmaxv_colors(1);
            end
            for j=1:length(x)
                v=Es_optGC(1,j);
                colorid=min(size(plotGCcolor,1),max(1,round(size(plotGCcolor,1)*(v-minv)/(maxv-minv))));
                colorm(1,j,:)=plotGCcolor(colorid,:);
                colorm(2,j,:)=plotGCcolor(colorid,:);
                colorids(j)=colorid;
            end
            alongGC.colorids=colorids;
            surface([x;x],[y;y],[z;z],colorm,'facecol','no', 'edgecol','interp', 'linew',1);
            
            freezeColors
            colormap(plotGCcolor);
            ch=colorbar;
            set(ch,'ytick',[min(optgrowthrates_m(:)),max(optgrowthrates_m(:))]);
            set(ch,'yticklabel',{num2str(minv,2),num2str(maxv,2)});
        end
        %axis([ranges(1,1),ranges(1,2),ranges(2,1),ranges(2,2)])
    end
    
    alongGCs{i}=alongGC;
    % if one need to plot this
end % i=1:length(Dlist)

alongFBs=[];
% get the flux balance points
for k=1:size(supplies,2)
    
    onesupply=supplies(:,k);
    FB=[];
    FB.supply=onesupply;
    
    
    densitybydim=zeros(size(optgrowthrates_m,1),size(optgrowthrates_m,2),envdim);
    for i=1:size(optgrowthrates_m,1)
        for j=1:size(optgrowthrates_m,2)
            env=reshape(nutreint_values_m(i,j,:),envdim,1);
            densitybydim(i,j,:)=(onesupply-env)./reshape(consumptions_m(i,j,:),envdim,1);
        end
    end
    
    mdiff=densitybydim(:,:,1)-densitybydim(:,:,2);
    
    
    
    
    % get contours for this one
    nutreint_values_m(:,1,1);
    zerijs=round(contourc(mdiff,[0,0]));
    zerijs(:,1)=[];
    [value,order]=sort(zerijs(1,:));
    zerijs=zerijs(:,order);
    
    
    envs=[];
    Es_alongFB=[];
    for i=1:size(zerijs,2)
        envs(:,i)=nutreint_values_m(zerijs(1,i),zerijs(2,i),:);
        Es_alongFB(:,i)=opt_Es(zerijs(1,i),zerijs(2,i),:);
    end
    
    
    
    FB.ijs=zerijs;
    FB.envs=envs;
    FB.Es_alongFB=Es_alongFB;
    
    alongFBs{k}=FB;
    
    if plotFB
        plot(envs(1,:),envs(2,:),'color',[0.8,0,0],'linewidth',plotFBlw);
    end
    
end