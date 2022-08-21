function [jU2, jS2, jV2] = Jacobi_svd_2sided(A)
    tol = 0.0001;
    delta = tol*norm(A, 'fro');
    [m, n] = size(A);
    jU2 = eye(m);
    jS2 = A;
    jV2 = eye(n);
    while (off(jS2) > delta)
        for p=1:min(m, n)-1
            for q=p+1:min(m, n)
                [c1, s1, c2, s2] = asymSchur2(jS2, p, q);
                J1 = eye(m);
                J1([p, q], [p, q]) = [c1, s1; -s1, c1];
                J2 = eye(n);
                J2([p, q], [p, q]) = [c2, s2; -s2, c2];
                jU2 = jU2*J1;
                jS2 = J1'*jS2*J2;
                jV2 = jV2*J2;
            end
        end
        if m < n
            for p=1:m
                for q=m+1:n
                    if jS2(p, p) == 0
                        c2 = 0;
                        s2 = 1;
                    else
                        t = -jS2(p, q)/jS2(p, p);
                        c2 = 1/sqrt(1+t^2);
                        s2 = t*c2;
                    end
                    J2 = eye(n);
                    J2([p, q], [p, q]) = [c2, s2; -s2, c2];
                    jS2 = jS2*J2;
                    jV2 = jV2*J2;
                end
            end
        elseif m > n
            for p=n+1:m
                for q=1:n
                    if jS2(q, q) == 0
                        c1 = 0;
                        s1 = 1;
                    else
                        t = -jS2(p, q)/jS2(q, q);
                        c1 = 1/sqrt(1+t^2);
                        s1 = t*c1;
                    end
                    J1 = eye(m);
                    J1([q, p], [q, p]) = [c1, s1; -s1, c1];
                    jU2 = jU2*J1;
                    jS2 = J1'*jS2;
                end
            end   
        end
    end
    for i=1:min(m, n)
        if jS2(i, i) < 0
            jS2(i, i) = -jS2(i, i);
            jU2(:, i) = -jU2(:, i);
        end
    end
    jS2 = diag(jS2);
    for i=1:min(m, n)-1
        for j=1:min(m, n)-i
            if jS2(j) < jS2(j+1)
                temp = jU2(:, j);
                jU2(:, j) = jU2(:, j+1);
                jU2(:, j+1) = temp;
                temp = jS2(j);
                jS2(j) = jS2(j+1);
                jS2(j+1) = temp;
                temp = jV2(:, j);
                jV2(:, j) = jV2(:, j+1);
                jV2(:, j+1) = temp;
            end
        end
    end
    temp = jS2;
    jS2 = zeros(m, n);
    for i=1:min(m, n)
        jS2(i, i) = temp(i);
    end
end