function threepoints=ZL_Triangle(elip_ab,vectordir,halfangle)
% return the three points for a triangle in a ellipsis,with center in zero
% halfangle: is by degree, 90 as maximal
drawmethod=2;

drawcase=0;
%%%%
halfangle=deg2rad(halfangle);
threepoints=zeros(2,3);

% a*cos(t),b*sin(t)
a=elip_ab(1);
b=elip_ab(2);

angle=atan(vectordir(2)/vectordir(1));

% solve for the vertex
v_x=((-1)^(vectordir(1)<0))*sqrt(1/(1/a^2+(vectordir(2)/vectordir(1)/b)^2));
v_y=((-1)^(vectordir(2)<0))*sqrt(1/(1/b^2+(vectordir(1)/vectordir(2)/a)^2));

threepoints(:,1)=[v_x;v_y];

%the two other intersections
if drawmethod==1
    
    for i=1:2
        
        if i==1
            newangle=angle+halfangle;
        else
            newangle=angle-halfangle;
        end
        
        
        k=sin(newangle)/cos(newangle);
        testy1=((v_y + k*(k^2 + 1 - k^2*v_x^2 + 2*k*v_x*v_y - v_y^2)^(1/2) - k*v_x))/(k^2 + 1);
        testx1=v_x-(v_y-testy1)/k;
        
        testy2=-((k*(k^2 + 1 - k^2*v_x^2 + 2*k*v_x*v_y - v_y^2)^(1/2) - v_y + k*v_x))/(k^2 + 1);
        testx2=v_x-(v_y-testy2)/k;
        
        newtheta(1)=atan(testy1/testx1);if testx1<0; newtheta(1)=newtheta(1)+pi; end
        newtheta(2)=atan(testy2/testx2);if testx2<0; newtheta(2)=newtheta(2)+pi; end
        
        %
        new_x(1)=a*cos(newtheta(1));
        new_y(1)=b*sin(newtheta(1));
        
        new_x(2)=a*cos(newtheta(2));
        new_y(2)=b*sin(newtheta(2));
        
        
        [maxv,maxloc]=max(abs(new_x-v_x)+abs(new_y-v_y));
        
        threepoints(:,i+1)=[new_x(maxloc);new_y(maxloc)];
        
    end
    
    
else
    
    angle=atan(vectordir(2)/vectordir(1));
    % the two other intersections
    for i=1:2
        
        if i==1
            newangle=angle+halfangle;
        else
            newangle=angle-halfangle;
        end
        
        
        k=sin(newangle)/cos(newangle);
        
        new_y(1)=(b*(b*v_y + a*k*(a^2*k^2 + b^2 - k^2*v_x^2 + 2*k*v_x*v_y - v_y^2)^(1/2) - b*k*v_x))/(a^2*k^2 + b^2);
        new_x(1)=v_x-(v_y-new_y(1))/k;
        
        new_y(2)=-(b*(a*k*(a^2*k^2 + b^2 - k^2*v_x^2 + 2*k*v_x*v_y - v_y^2)^(1/2) - b*v_y + b*k*v_x))/(a^2*k^2 + b^2);
        new_x(2)=v_x-(v_y-new_y(2))/k;
        
        
        [maxv,maxloc]=max(abs(new_x-v_x)+abs(new_y-v_y));
        
        threepoints(:,i+1)=[new_x(maxloc);new_y(maxloc)];
    end
end
%
%
if drawcase

    figure;
    anglelist=linspace(0,2*pi,50);
    xlist=a*cos(anglelist);
    ylist=b*sin(anglelist);
    plot(xlist,ylist);hold on;

    patch(threepoints(1,:),threepoints(2,:),[0,0,0])



end

