function s = off(A)
    s = norm(A, 'fro')^2;
    for i=1:min(size(A))
        s = s-A(i, i)^2;
    end
    s = sqrt(s);
end