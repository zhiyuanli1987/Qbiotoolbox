function hh=Easypcolor(data,xticklabel,yticklabel)
    li=length(data(:,1));
    lj=length(data(1,:));
    data(li+1,:)=0;
    data(:,lj+1)=0;
    hh=pcolor(data);
    %colorbar
    
    if ~isempty(xticklabel)
        set(gca,'XTick',1.5:1:lj+1.5);
    else
        set(gca,'XTick',[]);
    end
    set(gca,'XTickLabel',xticklabel);
%     
    if ~isempty(yticklabel)
        set(gca,'YTick',1.5:1:li+1.5);
    else
        set(gca,'YTick',[]);
    end
    set(gca,'YTickLabel',yticklabel);
    
    set(hh,'LineStyle','none');