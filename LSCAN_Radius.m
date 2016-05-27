function R=LSCAN_Radius(a,Area)
%return the radius of curvature of a portion of circle, knowing its area
%a is the radius of the section cutting the circle
%(=rc the radius of the contractile ring)
%a is a minimal value for the radius of curvature
m=length(a);
R=zeros(m,1);

for i=1:m
    if(isnan(Area(i))) R(i)=NaN;
    else
        if(Area(i)<=pi/2*a(i)^2)
            %circle is smaller than an hemisphere
            R(i)=fzero(@(x)SmallCircleArea(a(i),x)-Area(i),[a(i),300*a(i)]);

        else
            %circle is larger than an hemisphere
            R(i)=fzero(@(x)CircleArea(a(i),x)-Area(i),[a(i),300*a(i)]);
        end
    end
end

end

%area of the smaller segment of a circle, a is the radius of the secant
function Area=SmallCircleArea(a,R)
Area=R.^2.*(asin(a./R)-a./R.*sqrt(1-a.^2./R.^2));
end

%area of the larger segment of a circle, a is the radius of the secant
function Area=CircleArea(a,R)
Area=R.^2.*(pi-asin(a./R)+a./R.*sqrt(1-a.^2./R.^2));
end