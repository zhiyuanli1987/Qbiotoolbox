function densitymatrix=ZL_Densityplot(x,y,pnum,colortable,showlog)

voidvalue=nan;

nonnanlist=find(~isnan(x) & ~isnan(y) & ~isinf(x) & ~isinf(y));
x=x(nonnanlist);
y=y(nonnanlist);



densitymatrix=zeros(pnum,pnum);

xids=min(pnum,max(1,ceil(pnum*(x-min(x))/(max(x)-min(x)))));
yids=min(pnum,max(1,ceil(pnum*(y-min(y))/(max(y)-min(y)))));
mids=pnum*(xids-1)+yids;
freq=tabulate(mids);



if showlog
    densitymatrix(freq(:,1))=log10(freq(:,2));
else
    densitymatrix(freq(:,1))=(freq(:,2));
end



densitymatrix(isinf(densitymatrix))=voidvalue;
densitymatrix(densitymatrix==0)=voidvalue;

xlist=[linspace(min(x),max(x),pnum),max(x)*1.01];
ylist=[linspace(min(y),max(y),pnum),max(y)*1.01];
densitymatrix(end+1,:)=0;
densitymatrix(:,end+1)=0;
colormap(colortable);
hh=pcolor(xlist,ylist,densitymatrix);
set(hh,'LineStyle','none');
% 
%colortable=colormap;densityids=max(1,min(64,round(64*(densitymatrix-min(densitymatrix(:)))/(max(densitymatrix(:))-min(densitymatrix(:))))));colorcodes=colortable(densityids(mids),:);n=min(10000,length(x));scatter(x(1:n),y(1:n),[],colorcodes(1:n,:),'.')
% 
