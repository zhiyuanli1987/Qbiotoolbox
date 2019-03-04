function [nutreint_boarders_m,nutreint_values_m,optgrowthrates_m,opt_Es,consumptions_m,withincellv_m]=Nutrientspace_optspecies(speciespara,MetabModel,ranges,pnumbydim,showprogress,drawpara,runpara)

metabfun=MetabModel(speciespara);
if  nargin <5
    showprogress=0;
end

if nargin<6
    drawpara=[];
end

maxrepeat=100;
additionalrun=3;
if nargin>6
    if isfield(runpara,'maxrepeat')
        maxrepeat=runpara.maxrepeat;
    end
    if isfield(runpara,'additionalrun')
        additionalrun=runpara.additionalrun;
    end
end


pnumbydim_cell=num2cell(pnumbydim);
env_dim=size(ranges,1);
intercell_dim=metabfun.intercell_dim;
f_dim=metabfun.Edim;

nutreint_boarders_m=nan*zeros(pnumbydim_cell{:},env_dim);
nutreint_values_m=nan*zeros(pnumbydim_cell{:},env_dim);
optgrowthrates_m=nan*zeros(pnumbydim_cell{:});
opt_Es=nan*zeros(pnumbydim_cell{:},f_dim);
consumptions_m=nan*zeros(pnumbydim_cell{:},env_dim);
withincellv_m=nan*zeros(pnumbydim_cell{:},intercell_dim);

envboarderlist=[];
for e=1:env_dim
    envboarderlist(e,1:pnumbydim(e))=linspace(ranges(e,1),ranges(e,2),pnumbydim(e));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    
    % assign to boarders
    for e=1:env_dim
        nutreint_boarders_m(cellsub{:},e)=envboarderlist(e,subs(e));
    end
    
    if all(subs<pnumbydim)
        for e=1:env_dim
            env(e)=mean(envboarderlist(e,(subs(e):subs(e)+1)));
            nutreint_values_m(cellsub{:},e)= env(e);
        end
        
        % compute the steady state of the system
        % collect some possible initial values from surrounding dots
        possible_inis=[];
        for e=1:env_dim
            nearbysub=subs;
            if subs(e)>1
                nearbysub(e)=nearbysub(e)-1;
                cellx=num2cell(nearbysub);
                inis=[];
                inis.iniEs=opt_Es(cellx{:},:);
                inis.inig=optgrowthrates_m(cellx{:});
                inis.inixs=withincellv_m(cellx{:},:);
                possible_inis{end+1}=inis;
            end
        end
        
        
        % obtain optimal strategy under this environment
        runtime=1;
        foundsolution=0;
        while runtime<maxrepeat && foundsolution==0
            
            optEs_candidate=[];
            optg_candidate=[];
            fval_candidate=[];
            for rr=1:additionalrun
                if ~(runtime>length(possible_inis))
                    % the initial value obtined from surroundings
                    [optEs_candidate(:,rr),optg_candidate(rr),fval_candidate(rr)]=Obtain_optEG_givenenv(env,metabfun,possible_inis{runtime});
                else
                    [optEs_candidate(:,rr),optg_candidate(rr),fval_candidate(rr)]=Obtain_optEG_givenenv(env,metabfun);
                end
                runtime=runtime+1;
            end
            % choose the optimal one
            [optg,optloc]=max(optg_candidate);
            fval=fval_candidate(optloc);
            optEs=optEs_candidate(:,optloc);
            
            
            
        
            % is a valid solution, add it
            if ~(optg<0) && all(~(optEs<-0.001)) && all(~(optEs>1)) && abs(fval)<0.001
                foundsolution=1;
                optgrowthrates_m(i)=optg;
                opt_Es(cellsub{:},:)=optEs;
                
                % within cell steady
                optspecies=speciespara;
                optspecies.Es=optEs;
                
                optmetabfun=MetabModel(optspecies);

                if intercell_dim>0
                    [g2,withincellsteady,fval2]=Obtain_Growth_givenenv(env,optmetabfun);
                    withincellv_m(cellsub{:},:)=withincellsteady; % nan*zeros(pnumbydim,intercell_dim);
                else
                    withincellsteady=nan;
                end
                % compute the consumption rate
                Intakes=optmetabfun.Intakes;
                consumptions_m(cellsub{:},:)=Intakes(env,withincellsteady);
            end
            
            
            % run till it find the right solution, or break    
        end % while runtime<maxrepeat && foundsolution==0
            
        
        if showprogress
            if mod(i,5)==0
                fprintf('Constructing nutreint space: %f\n',100*i/totalrecombnum);
            end
        end
    end
end

if ~isempty(drawpara)
    figure;
    h=pcolor(nutreint_boarders_m(:,:,1),nutreint_boarders_m(:,:,2),optgrowthrates_m);
    h.LineStyle='none';
    colorbar;
    xlabel('c_a')
    ylabel('c_b');
    title('Growth rate'); 
    hold on;
    contour(nutreint_values_m(:,:,1),nutreint_values_m(:,:,2),optgrowthrates_m,'k')
    set(gca,'fontsize',22);
    
end


