function [finalstate,traject]=Visualize_chemostatdynamic(chemostat_para,species_metabfuns,drawpara,odepara)
% run the chemostat ODE, obtain :
% 1. the trajectories in nutirent space
% 2. time-course

c_supplys=chemostat_para.c_supplys;
env_dim=length(c_supplys);
species_num=length(species_metabfuns);
intercell_dim=species_metabfuns{1}.intercell_dim;

% parameters in running ODEs
timespan_runode=[0,100];
inistates=[rand(1,env_dim),0.1*rand(1,species_num),0.1*rand(1,intercell_dim*species_num)];
if nargin > 3
    if isfield(odepara,'timespan_runode')
        timespan_runode=odepara.timespan_runode;
    end
    if isfield(odepara,'inistates')
        inistates=odepara.inistates;
    end    
end


% parameters for drawing
lw1=1;
lw2=1;

showtimecourse=1;
arrow_num=2;
arrow_sefrac=[0.1,0.5];
arrow_percent=1/50;
axislim=[zeros(env_dim,1),c_supplys];
arrow_angle=20;
if nargin>2
    if isfield(drawpara,'linewidths')
        lw1=drawpara.linewidths(1);
        lw2=drawpara.linewidths(2);
    end
    if isfield(drawpara,'speciescolor')
        speciescolor=drawpara.speciescolor;
    else
        speciescolor=colormap('lines');
    end
    if isfield(drawpara,'showtimecourse')
        showtimecourse=drawpara.showtimecourse;
    end
    if isfield(drawpara,'arrow_num')
        arrow_num=drawpara.arrow_num;
    end
    if isfield(drawpara,'arrow_sefrac')
        arrow_sefrac=drawpara.arrow_sefrac;
    end
    if isfield(drawpara,'arrow_percent')
        arrow_percent=drawpara.arrow_percent; % percent of the arrow for each axis;
    end
    if isfield(drawpara,'axislim')
        axislim=drawpara.axislim;
    end
    if isfield(drawpara,'arrow_angle')
        arrow_angle=drawpara.arrow_angle;
    end
end


% define ODEs
odetouse=@(t,x)ODE_onechemostat(x,chemostat_para,species_metabfuns);

% run odes
for i=1:size(inistates,1)
    
    ini=inistates(i,:);
    [timelist,trajecties]=ode23(odetouse,timespan_runode,ini);
    

    % generate arrow plot on the nutirent dim
    
    % the length traveled in the nutrient dimension
    envtraject=trajecties(:,1:env_dim);
    lengthfromstart=zeros(length(timelist),1);
    for k=2:length(timelist)
        lengthfromstart(k)=lengthfromstart(k-1)+sum(abs(envtraject(k,:)-envtraject(k-1,:)));
    end
    lengthfromstart(1)=[];
    
   
    
    % find the distance to cut
    lenlist=linspace(arrow_sefrac(1),arrow_sefrac(2),arrow_num)*lengthfromstart(end);
    
    
    arrowloc=zeros(arrow_num,1);
    for a=1:arrow_num
        arrowloc(a)=max(find(lengthfromstart<lenlist(a)));
    end

    
    % construct a triangle for arrow
    if env_dim==2
        plot(envtraject(:,1),envtraject(:,2),'k','linewidth',lw1); 
        for a=1:length(arrowloc)
            startp=envtraject(arrowloc(a),:)';
%             vectordir=(envtraject(arrowloc(a)+1,:)-envtraject(arrowloc(a),:));
%             dirtriangle=repmat(startp,1,3)+ZL_Triangle((axislim(:,2)-axislim(:,1))*arrow_percent,vectordir,arrow_angle);
%             
%             ZL_Triangle((axislim(:,2)-axislim(:,1))*arrow_percent,vectordir,arrow_angle)
%             
%             
%             
%             patch(dirtriangle(1,:),dirtriangle(2,:),[0,0,0])
%             
            
    timeloc=min(find(timelist>arrow_percent+timelist(arrowloc(a))));
            endp=envtraject(timeloc,:)';
            
            arrow(startp,endp)
            % use arrow percent as the interval time
           

            %mArrow3([startp;0],[endp;0]);
            
        end
    elseif env_dim==3
        
        plot3(envtraject(:,1),envtraject(:,2),envtraject(:,3),'k','linewidth',lw1); 
        for a=1:length(arrowloc)
            startp=envtraject(arrowloc(a),:)';
            
            % use arrow percent as the interval time
            timeloc=min(find(timelist>arrow_percent+timelist(arrowloc(a))));
            endp=envtraject(timeloc,:)';
            
            %vectordir=(envtraject(arrowloc(a)+1,:)-envtraject(arrowloc(a),:))./(timelist(arrowloc(a)+1)-timelist(arrowloc(a)));
            %endp=startp+arrow_percent*50*vectordir';

            mArrow3(startp,endp);
        end
    end
    
   
    
    
    
end

finalstate=trajecties(end,:);

traject=[timelist,trajecties];
% generate time-course of populations for the last run
if showtimecourse
    figure;hold on;
    for s=1:species_num
        plot(timelist,trajecties(:,env_dim+s),'linewidth',lw2,'color',speciescolor(s,:))
    end
    set(gca,'fontsize',20);
    xlabel('Time');
    ylabel('Biomass');
end
