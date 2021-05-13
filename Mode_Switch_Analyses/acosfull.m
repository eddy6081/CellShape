function theta = acosfull(x,y)
    r=sqrt(x.^2 + y.^2);
    z = x./r;
    theta = acos(z);
    for el=1:length(y)
        if y(el)<0
            theta(el) = -theta(el);
        end
    end
end