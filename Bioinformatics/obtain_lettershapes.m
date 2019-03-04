%clear
clc
% this program generate the structure "letterclass", only need to run it for once

%%%%%%%%%%%%%% customizable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% assign colors for nucleotides and amino acids, can be changed for different visualization
ACGT={'A','C','G','T'};
ACGTcolor={[0,0.8,0], % color for nucleotides, following the sequence in ACGT
    [0,0,1],
    [244, 220, 66]/255,
    [1,0,0]};
nuc_colormaping=containers.Map(ACGT,ACGTcolor);

Astrs={'A',
    'R'
    'N'
    'D'
    'C'
    'Q'
    'E'
    'G'
    'H'
    'I'
    'L'
    'K'
    'M'
    'F'
    'P'
    'S'
    'T'
    'W'
    'Y'
    'V'
    'B'
    'Z'
    'X'};
Acolors={[0,0,0]; % colors for amino acids, following the sequence in Astrs
    [0,0,1];
    [20, 255, 255]/255;
    [1,0,0];
    [1,1,0];
    [20, 255, 255]/255;
    [1,0,0];
    [1,1,1]*0.5;
    [137, 165, 237]/255;
    [0,1,0];
    [0,1,0];
    [0,0,1];
    [1,1,0];
    [5,86,178]/255;
    [219, 138, 144]/255;
    [249, 179, 0]/255;
    [249, 179, 0]/255;
    [234, 82, 219]/255;
    [0.0196    0.3373    0.6980];
    [0,1,0];
    [211, 182, 124]/255;
    [211, 182, 124]/255;
    [211, 182, 124]/255;
    };
aa_colormaping=containers.Map(Astrs,Acolors);
%%%%%%%%%%%%%%%%%%%%%%%%%/customizable%%%%%%%%%%%%%%%%%%%%%



%%%%%%%% parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% some parameters for image extraction, may be changed, but be careful%%%%%
plottest=1;inputcolor=[1,1,0];
dist_thresh=0.05;
resize_percent=0.3; % how much resolution to shrink.maximal is 1. smaller value leads to more rough edges but also smaller image size
%%%%%%%% /parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%% DON'T CHANGE THINGS BLOW!!!%%%%%%%%%%%%%%%%%
% use command "close all" for closing all ploting windows
letterclass=[];
letterclass.alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
letterclass.map=containers.Map(cellstr(letterclass.alphabet'),1:length(letterclass.alphabet));
for l=1:length(letterclass.alphabet)
    letter=letterclass.alphabet(l);
    
    fs=500;
    text(0,0, letter, 'HorizontalAlignment','center','fontsize',fs);
    axis([-1 1 -1 1])
    axis off
    
    
    F = getframe(gcf);
    [threeDim, Map] = frame2im(F);
    I_ori=mat2gray(threeDim);I_ori=I_ori(:,:,1);
    I=imresize(I_ori,resize_percent);
    
    

    edgeM = edge(I,'Canny');
    [edge_js,edge_is]=find(edgeM>0.99);
    
    movev=min(edge_is)*size(I,1)/size(I,2);
    edge_js=-1*edge_js;
    edge_js=edge_js-min(edge_js);
    edge_is=(edge_is-min(edge_is))*size(I,1)/size(I,2);% this is monotonically increasing
    
    maxv=max([edge_is;edge_js]);
    edge_is=edge_is/maxv;movev=movev/maxv;
    edge_js=edge_js/maxv;
    
    datanum=length(edge_js);
    dot_counted=zeros(1,datanum);
    curvesdots=[];
    curves_dotnum=[];
    while sum(dot_counted)<datanum
        
        notyetcounted=find(dot_counted==0);
        stack=[];
        nextneighbor=notyetcounted(1);
        while ~isempty(nextneighbor)
            % find the one that is within a distance and closet
            stack=[stack,nextneighbor];
            dot_counted(nextneighbor)=1;
            notyetcounted=find(dot_counted==0);
            
            currentdot=nextneighbor;
            nextneighbor=[];
            
            dist_to_current=max([abs(edge_is(notyetcounted)-edge_is(currentdot)),abs(edge_js(notyetcounted)-edge_js(currentdot))],[],2);
            candidatesloc=find(dist_to_current<dist_thresh);
            if ~isempty(candidatesloc)
                [minv,minloc]=min(dist_to_current(candidatesloc));
                nextneighbor=notyetcounted(candidatesloc(minloc));
            end
        end
        curvesdots{end+1}=stack;
        curves_dotnum(end+1)=length(stack);
        
        
    end
    
    [value,order]=sort(curves_dotnum,'descend');
    curvesdots=curvesdots(order);
    letters.curves_xy=[];
    letters.curvenum=0;
    letters.str=letter;
    for i=1:length(curvesdots)
        if length(curvesdots{i})>100
            letters.curvenum=letters.curvenum+1;
            letters.curves_xy{letters.curvenum}=[edge_is(curvesdots{i}),edge_js(curvesdots{i})];
        end
    end
    
    clf
    letterclass.curves_xy{l}=letters.curves_xy;
    letterclass.curvenum(l)=letters.curvenum;
    letterclass.strs{l}=letters.str;
    letterclass.centerx(l)=(max(edge_is)-min(edge_is))/2;
    
end

if plottest
    for j=1:length(letterclass.strs)
        figure;
        
        % plot the results
        for i=1:letterclass.curvenum(j)
            xs=letterclass.curves_xy{j}{i}(:,1)+(0.5-letterclass.centerx(j));
            ys=letterclass.curves_xy{j}{i}(:,2);
            if i==1
                tofillcolor=inputcolor;
            else
                tofillcolor=[1,1,1];
            end
            fh=fill(xs,ys,tofillcolor); hold on;
        end
        title(letters.curvenum)
        
        %plot(edge_is,edge_js,'.');
        axis([0,1,0,1])
    end
end

weight_thresh=0;
SeqLogo_ZL(Astrs,'aa',aa_colormaping,letterclass,weight_thresh);


