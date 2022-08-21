function [c, s] = orthogonalization(x, y)
    if (norm(x) == norm(y))
        c = 1/sqrt(2);
        s = 1/sqrt(2);
    else
        t = 2*x'*y/(norm(y)^2-norm(x)^2);
        c = sqrt((1+1/sqrt(1+t^2))/2);
        s = sqrt((1-1/sqrt(1+t^2))/2);
    end
end