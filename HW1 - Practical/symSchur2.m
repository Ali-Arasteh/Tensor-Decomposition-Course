function [c, s] = symSchur2(A, p, q)
    a = A([p, q], [p, q]);
    if (a(1, 2) == 0)
        c = 1;
        s = 0;
    else
        tau = (a(2, 2)-a(1, 1))/(2*a(1, 2));
        if (tau >= 0)
            t_min = 1/(tau+sqrt(1+tau^2));
        else
            t_min = 1/(tau-sqrt(1+tau^2));
        end
        c = 1/(sqrt(1+t_min^2));
        s = t_min*c;
    end
end