function [jV, jD] = Jacobi_eig(A)
    tol = 0.0001;
    delta = tol*norm(A, 'fro');
    n = size(A, 1);
    jV = eye(n);
    jD = A;
    while (off(jD) > delta)
        for p=1:n-1
            for q=p+1:n
                [c, s] = symSchur2(jD, p, q);
                J = eye(n);
                J([p, q], [p, q]) = [c, s; -s, c];
                jV = jV*J;
                jD = J'*jD*J;
            end
        end
    end
    jD = diag(jD);
    for i=1:n-1
        for j=1:n-i
            if jD(j) > jD(j+1)
                temp = jD(j);
                jD(j) = jD(j+1);
                jD(j+1) = temp;
                temp = jV(:, j);
                jV(:, j) = jV(:, j+1);
                jV(:, j+1) = temp;
            end
        end
    end
    jD = diag(jD);
end