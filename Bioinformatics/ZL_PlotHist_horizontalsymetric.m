function ZL_PlotHist_horizontalsymetric(data,boarderlist,colorcode,xcenter,hmax)

hold on;
% get boarders for bins
%boarderlist(1)=min(min(data),boarderlist(1));
%boarderlist(end)=max(max(data),boarderlist(end));

% get count
realboarder=boarderlist;
realboarder(1)=-inf;
realboarder(end)=inf;

count=zeros(1,length(boarderlist)-1);
for i=1:(length(boarderlist)-1)
    count(i)=sum(~(data<realboarder(i)) & ~(data>realboarder(i+1)));
end
count=count/sum(count);

% plot
for i=1:(length(boarderlist)-1)
    if count(i)>0
        xposis=xcenter+[-1,1]*1/2*count(i)*hmax;
        yposis=[boarderlist(i),boarderlist(i+1)];

        h=patch([xposis(1),xposis(2),xposis(2),xposis(1)],[yposis(1),yposis(1),yposis(2),yposis(2)],colorcode);
        h.LineStyle='none';
    end
end




