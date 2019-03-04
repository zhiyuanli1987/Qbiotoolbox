function pointsid=Visualize_supplylines(alongGC,metabfun,drawpara,additionalenv)
% draw the supply vector aline the line, and extract the species, if not
% yet

% default
selectpnum=3;
linestr='--';
linecolor=[0,0,0];
linewidth=1;
if nargin > 2
    
    if isfield(drawpara,'selectpnum')
        selectpnum=drawpara.selectpnum;
    end
    
    if isfield(drawpara,'linestr')
        linestr=drawpara.linestr;
    end
    
    if isfield(drawpara,'linewidth')
        linewidth=drawpara.linewidth;
    end
    
    if isfield(drawpara,'linecolor')
        linecolor=drawpara.linecolor;
    end
    
end

additionalenvnum=0;
if nargin>3
    additionalenvnum=size(additionalenvnum,2);
    
end

envs=alongGC.envs;

% pick points on the GC

% pick points on the GC
distalongenv=zeros(size(envs,2),1);
for i=2:(size(envs,2))
    distalongenv(i)=distalongenv(i-1)+norm(envs(:,i)-envs(:,i-1));
end
distb=linspace(min(distalongenv),max(distalongenv),(selectpnum+2));
distb(1)=[];
distb(end)=[];

pointsid=[];
for i=1:length(distb)
    loc=find(distalongenv(1:(end-1))<distb(i) & distalongenv(2:end)>distb(i));
    if ~isempty(loc)
        pointsid(end+1)=loc;
    end
end
pointsid=unique(pointsid);


for i=1:(length(pointsid)+additionalenvnum)
    if i>length(pointsid)
        
        env=additionalenv(:,i-length(pointsid));
    else
        env=envs(:,pointsid(i));
    end
    [FBv_givenIa,onelargesupplypoint,fval]=Obtain_supplyvector(metabfun,env);
    
    if size(linecolor,1)==1
        colortouse=linecolor;
    else
        colortouse=linecolor(i,:);
    end
    
    plot([env(1),onelargesupplypoint(1)],[env(2),onelargesupplypoint(2)],linestr,'color',linecolor,'linewidth',linewidth);
end
