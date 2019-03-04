function ch=ZL_Densityscatter(x,y,neigborcrite,colortable,dotsize,maxneighbornum,showdistdistri,squaredistthre)

abs_thresh=neigborcrite(1);
percent_thresh=neigborcrite(2);


nonnanlist=(~isnan(x) & ~isnan(y) & ~isinf(x) & ~isinf(y));
datanum=sum(nonnanlist);
x=reshape(x(nonnanlist),datanum,1);
y=reshape(y(nonnanlist),datanum,1);


% distance matrix between these points
if isempty(squaredistthre)
    datap=[x,y];

    distline=pdist(datap);
    % threshold for the first percentdist
    if isnan(abs_thresh)
        sortedvalue=sort(distline);
        abs_thresh=sortedvalue(round(length(distline)*percent_thresh));
    end


    distm=squareform(distline);
    neighbornum=sum(1.0*(distm<abs_thresh));
    
else % only based on x-y distance
    neighbornum=0*x;
    xlist=(min(x)*0.99):(squaredistthre):(max(x)*1.01);
    ylist=(min(y)*0.99):(squaredistthre):(max(x)*1.01);
    
    for i=2:length(xlist)
        [i,length(xlist)]
        for j=2:length(ylist)
            datahere=(~(x<xlist(i-1)) & ~(x>xlist(i)) & ~(y<ylist(j-1)) & ~(y>ylist(j)) );
            neighbornum(datahere)=sum(datahere);
        end
    end
end

max(neighbornum)
    


    % map the neighbor number to color
    minv=0;
    if maxneighbornum>0 % not adjustable

        maxv=maxneighbornum;

    else
        maxv=max(neighbornum);

    end


 colorids=min(max(round(size(colortable,1)*(neighbornum-minv)/(maxv-minv)),1),size(colortable,1));

    


scatter(x,y,dotsize,colortable(colorids,:),'filled')
colormap(colortable)
ch=colorbar;
set(ch,'YTick',[0,1],'YTickLabel',{num2str(minv),num2str(maxv)})
set(gca,'fontsize',20)


if showdistdistri
    
    
    figure;
    subplot(2,1,1)
    hist(distline,100)
    subplot(2,1,2)
    hist(neighbornum,10)
end

