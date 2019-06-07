function [alongGC,h]=Visualize_chemostatsteady(chemostat_para,metabfun,drawpara)

c_supplys=chemostat_para.c_supplys;
d=chemostat_para.d;
env_dim=length(c_supplys);

pnum=100;
ranges=[0.01*zeros(env_dim,1),c_supplys];
growthcolor=[0,0,0];%localpoints(k,:);
fluxcolor=[0,0,0];
lw=1;
showfluxcurve=1;
surfacealpha=0.5;


% parameters for drawing
if nargin > 2
    if isfield(drawpara,'pnum')
        pnum=drawpara.pnum;
    end
    if isfield(drawpara,'ranges')
        ranges=drawpara.ranges;
    end
    if isfield(drawpara,'growthcolor')
        growthcolor=drawpara.growthcolor;
    end
    if  isfield(drawpara,'fluxcolor')
        fluxcolor=drawpara.fluxcolor;
    end
    if isfield(drawpara,'linewidth')
        lw=drawpara.linewidth;
    end
    if isfield(drawpara,'showfluxcurve')
        showfluxcurve=drawpara.showfluxcurve;
    end
    if isfield(drawpara,'surfacealpha')
        surfacealpha=drawpara.surfacealpha;
        
    end
    
else % if there is no inputs for drawing
    
end


if env_dim==2
    
    alist=linspace(ranges(1,1),ranges(1,2),pnum);
    GCFB_b=nan*zeros(2,pnum);
    
    analytical_FBGCmin=metabfun.analytical_FBGCincell_dim2;
    if ~isempty(analytical_FBGCmin) % there are analytical solution
        
        
        
        incellvfun=analytical_FBGCmin{3};
        
        for l=1:2
            tousefun=[];
            if l==1
                if ~isempty(analytical_FBGCmin{2})
                    tousefun=@(c_a)analytical_FBGCmin{2}(c_a,d);
                end
            else
                if ~isempty(analytical_FBGCmin{1})
                    tousefun=@(c_a)analytical_FBGCmin{1}(c_a,c_supplys);
                end
            end
            
            if ~isempty(tousefun)
                for p=1:pnum

                    c_a=alist(p);
                    c_b=tousefun(c_a);
                    
                    xs=1;
                    
                    if ~isempty(incellvfun)
                        xs=incellvfun(c_a,c_b,d);
                    end

                    if isempty(incellvfun) || (~(c_b<-0.01) & all(~(xs<0) & ~isinf(xs) & ~isnan(xs)))
                        GCFB_b(l,p)=c_b;
                    end
                end
            end
        end
    end
    
    goodlist=find(~(GCFB_b(1,:)<-0.01) & ~isinf(GCFB_b(1,:)) & ~isnan(GCFB_b(1,:)));
    % plot the GC and FB
    alongGC.envs=[alist(goodlist);GCFB_b(1,goodlist)];
    h=plot(alist(goodlist),GCFB_b(1,goodlist),'color',growthcolor,'linewidth',lw);hold on;
    
    if showfluxcurve
        h=plot(alist,GCFB_b(2,:),'color',fluxcolor,'linewidth',lw);hold on;
    end
    

    
    axis([ranges(1,:),ranges(2,:)])
    set(gca,'fontsize',20)
    xlabel('c_a');
    ylabel('c_b');
    
    
elseif env_dim==3
    
    alist=linspace(ranges(1,1),ranges(1,2),pnum);inter_a=alist(2)-alist(1);
    blist=linspace(ranges(2,1),ranges(2,2),pnum);inter_b=blist(2)-blist(1);
    
    analytical_FBGCmin=metabfun.analytical_FBGCincell_dim3;
    
    if ~isempty(analytical_FBGCmin{1})
        FB_cbcc_dim3 = @(c_a)analytical_FBGCmin{1}(c_a,c_supplys);% @(c_a,c_supplys)
    else
        FB_cbcc_dim3=[];
    end
    
    if ~isempty(analytical_FBGCmin{2})
        GC_cc_dim3 = @(c_a,c_b)analytical_FBGCmin{2}(c_a,c_b,d);% @(c_a,g)
    else
        GC_cc_dim3=[];
    end
    
    incellvfun=analytical_FBGCmin{3};
    
    % for the growth contour
    aboarder_m=zeros(pnum+1);
    bboarder_m=zeros(pnum+1);
    c_GC_m=nan*zeros(pnum+1);
    b_FB_list=zeros(1,pnum);
    c_FB_list=zeros(1,pnum);
    for i=1:(pnum+1)
        for j=1:(pnum+1)
            
            if ~(i>pnum)
                c_a=alist(i);
                aborder=c_a-1/2*inter_a;
            else
                aborder=c_a+1/2*inter_a;
            end
            
            if ~(j>pnum)
                c_b=blist(j);
                bborder=c_b-1/2*inter_b;
            else
                bborder=c_b+1/2*inter_b;
            end
            
            aboarder_m(i,j)=aborder;
            bboarder_m(i,j)=bborder;
            
            if ~(i>pnum) && ~(j>pnum)
                % compute c_c value for this c_a and c_b
                
                if ~isempty(GC_cc_dim3)
                    c_GC_m(i,j)=GC_cc_dim3(c_a,c_b);
                end
            end
        end
    end
    
    c_GC_m(c_GC_m<0)=nan;
    
    h=surf(aboarder_m,bboarder_m,c_GC_m);hold on;
    h.LineStyle='none';
    
    if all(growthcolor==[0,0,0])
        ; % do nothing
    else
        h.FaceColor=growthcolor;
    end
    h.FaceAlpha=surfacealpha;
    
    if showfluxcurve
        for i=1:pnum
            if ~isempty(FB_cbcc_dim3)
                bc=FB_cbcc_dim3(alist(i));
                b_FB_list(i)=bc(1);
                c_FB_list(i)=bc(2);
            end
        end
        badlist=b_FB_list<0 | c_FB_list<0;
        b_FB_list(badlist)=nan;
        c_FB_list(badlist)=nan;
        plot3(alist,b_FB_list,c_FB_list,'linewidth',lw,'color',fluxcolor);
    end
    axis([ranges(1,:),ranges(2,:),ranges(3,:)])
    
else
    msg = 'env dim should be 2 or 3';
    error(msg);
end

