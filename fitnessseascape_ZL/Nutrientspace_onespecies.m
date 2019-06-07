function [nutreint_boarders_m,nutreint_values_m,growthrates_m,consumptions_m,withincellv_m,ch]=Nutrientspace_onespecies(metabfun,ranges,pnumbydim,showprogress,drawpara,maxrepeat)
ch=nan;
if  nargin <4
    showprogress=0;
end

if nargin<5
    drawpara=[];
end

if nargin<6
    maxrepeat=100;
end



Intakes=metabfun.Intakes;
pnumbydim_cell=num2cell(pnumbydim);
% the growth rate and consumption at every point of the nutrient space
env_dim=size(ranges,1);
intercell_dim=metabfun.intercell_dim;

growthrates_m=zeros(pnumbydim);
nutreint_boarders_m=nan*zeros(pnumbydim_cell{:},env_dim);
nutreint_values_m=nan*zeros(pnumbydim_cell{:},env_dim);
consumptions_m=nan*zeros(pnumbydim_cell{:},env_dim);
withincellv_m=nan*zeros(pnumbydim_cell{:},intercell_dim);

envboarderlist=[];
for e=1:env_dim
    envboarderlist(e,1:pnumbydim(e))=linspace(ranges(e,1),ranges(e,2),pnumbydim(e));
end

% for each environment
totalrecombnum= prod(pnumbydim);
for i=1:totalrecombnum
    
    % which environment it is now
    subs=[];
    if env_dim==2
        [subs(1),subs(2)]=ind2sub(pnumbydim,i);
    elseif env_dim==3
        [subs(1),subs(2),subs(3)]=ind2sub(pnumbydim,i);
    else
        msg='nutirent dimension should be 2 or 3';
        perror(msg)
    end
    
    cellsub=num2cell(subs);
    
    env=zeros(env_dim,1);
    
    for e=1:env_dim
        % assign to boarders
        nutreint_boarders_m(cellsub{:},e)=envboarderlist(e,subs(e));
        
        % extract the environment if it is not the last one
    end
    
    if all(subs<pnumbydim)
        for e=1:env_dim
            env(e)=mean(envboarderlist(e,(subs(e):subs(e)+1)));
            nutreint_values_m(cellsub{:},e)= env(e);
        end
        
        
        % compute the steady state of the system
        % collect some possible initial values from surrounding dots
        possible_inis=[];
        if intercell_dim~=0
            for e=1:env_dim
                nearbysub=subs;
                if subs(e)>1
                    nearbysub(e)=nearbysub(e)-1;
                    cellx=num2cell(nearbysub);
                    iniv=[];
                    for d=1:intercell_dim
                        iniv(d)=withincellv_m(cellx{:},d);
                    end
                    possible_inis(end+1,:)=iniv;
                end
            end
        end
        
        
        runtime=1;
        foundsolution=0;
        
        odepara=[];
        odepara.timespan_runode=[0,100];
        odepara.plotresult=0;
        while runtime<maxrepeat && foundsolution==0
            if ~(runtime>size(possible_inis,1))
                iniv=possible_inis(runtime,:);
                odepara.runodeornot=0;
            else
                iniv=rand(1,intercell_dim);
                odepara.runodeornot=1;
            end
            odepara.iniv=iniv;
            
            
            [growthrate,withincellsteady,fval]=Obtain_Growth_givenenv(env,metabfun,odepara);

            
            if ~(growthrate<0) && all(~(withincellsteady<-0.001)) && abs(fval)<0.001
                foundsolution=1;
                growthrates_m(i)=growthrate;
                withincellv_m(cellsub{:},:)=withincellsteady; % nan*zeros(pnumbydim,intercell_dim);
                
                % compute the consumption rate
                consumptions_m(cellsub{:},:)=Intakes(env,withincellsteady);
            end
            
            runtime=runtime+1;
            % run till it find the right solution, or break    
        end
            
        if showprogress
            if mod(i,10)==0
                fprintf('Constructing nutreint space: %f\n',100*i/totalrecombnum);
            end
        end
    end
end

if ~isempty(drawpara)
    
    if isfield(drawpara,'colortable')
        colormap(drawpara.colortable)
    end
    
    h=pcolor(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),growthrates_m);
    h.LineStyle='none';
    ch=colorbar;
    xlabel('c_a')
    ylabel('c_b');
    title('Growth rate'); 
    set(gca,'fontsize',22);
    
    if isfield(drawpara,'drawcontour') && drawpara.drawcontour==1
        
        
        contour(nutreint_values_m(:,:,1),nutreint_values_m(:,:,2),growthrates_m,'k','ShowText','on')
    end
    
    
end


