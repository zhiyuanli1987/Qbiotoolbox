function [colortable,colorids]= Colormap_expression(data,threecolor_negzeroposi,dataonlog)
% input:values
% linear mapping
% return: colormap for all positive value -> one, near zero-> two, negative to be three,
% and colorid for each vlaues
data(isinf(data))=nan;
colortablesize=64;
maxpos=0;
minneg=0;

if dataonlog==0
    threshv=0;
else % dataonlog==1
    threshv=1;
end
posivalues=data(data>threshv);if ~isempty(posivalues);maxpos= max(posivalues);end
negvalues=data(data<threshv);if ~isempty(negvalues);minneg= min(negvalues);end




% how many dot out of the 64 to assign to the zero or negative?
if maxpos==minneg
    colortable=repmat(threecolor_negzeroposi(2,:),colortablesize,1);
    colorids=ones(size(data));
else
    pos_num=round((colortablesize-1)*maxpos/(maxpos-minneg));
    neg_num=(colortablesize-1)-pos_num;
    

    
    % make colors from neg to zero to positive
    negdiff=threecolor_negzeroposi(2,:)-threecolor_negzeroposi(1,:);
    increlist=linspace(0,1,neg_num+1);

    
    negcolors=repmat(threecolor_negzeroposi(1,:),neg_num,1)+repmat(increlist(1:(end-1))',1,3).*repmat(negdiff,neg_num,1);
    
    posidiff=threecolor_negzeroposi(3,:)-threecolor_negzeroposi(2,:);
    increlist=linspace(0,1,pos_num+1);
    poscolors=repmat(threecolor_negzeroposi(2,:),pos_num,1)+repmat(increlist(2:end)',1,3).*repmat(posidiff,pos_num,1);
    
    colortable=[negcolors;threecolor_negzeroposi(2,:);poscolors];
    
    % get the ids for each data
    colorids=ones(size(data))*(neg_num+1);
    colorids(isnan(data))=nan;
    ispos=(data>threshv);
    isneg=(data<threshv);
    colorids(ispos)=neg_num+1+min(pos_num,round(pos_num*(data(ispos)-threshv)/(maxpos-threshv)));
    colorids(isneg)=max(1,min(neg_num,round(neg_num*(data(isneg)-minneg)/(threshv-minneg))));

end



