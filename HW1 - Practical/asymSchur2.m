function [c1, s1, c2, s2] = asymSchur2(A, p, q)
    a = A([p, q], [p, q]);
    if (a(1, 2) == a(2, 1))
        c = 1; 
        s = 0;
    else
        t = (a(2, 1)-a(1, 2))/(a(1, 1)+a(2, 2));
        c = 1/(sqrt(1+t^2));
        s = t*c;
    end
    temp = [c, s; -s, c]*a;
    [c2, s2] = symSchur2(temp, 1, 2);
    c1 = c * c2 + s * s2;
    s1 = c * s2 - s * c2;
end