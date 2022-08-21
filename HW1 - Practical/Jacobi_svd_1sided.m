function [jU1, jS1, jV1] = Jacobi_svd_1sided(A)
    [m, n] = size(A);
    tol = 0.0001;
    if m <= n    
        delta = tol * norm(A*A', 'fro');
        D = A';
        jV1 = eye(m);
        while (off(D'*D) > delta)
            for p=1:m-1
                for q=p+1:m
                    [c, s] = orthogonalization(D(:, p), D(:, q));
                    J = eye(m);
                    J([p, q], [p, q]) = [c, s; -s, c];
                    D = D*J;
                    jV1 = jV1*J;
                end
            end
        end
        jU1 = zeros(n, m);
        jS1 = zeros(n, m);
        for i = 1:m
            jU1(:, i) = D(:, i)/norm(D(:, i));
            jS1(i, i) = norm(D(:, i));
        end
        temp = jU1;
        jU1 = jV1;
        jS1 = jS1';
        jV1 = temp;
    else
        delta = tol * norm(A'*A, 'fro');
        D = A;
        jV1 = eye(n);
        while (off(D'*D) > delta)
            for p=1:n-1
                for q=p+1:n
                    [c, s] = orthogonalization(D(:, p), D(:, q));
                    J = eye(n);
                    J([p, q], [p, q]) = [c, s; -s, c];
                    D = D*J;
                    jV1 = jV1*J;
                end
            end
        end
        jU1 = zeros(m, n);
        jS1 = zeros(m, n);
        for i = 1:n
            jU1(:, i) = D(:, i)/norm(D(:, i));
            jS1(i, i) = norm(D(:, i));
        end
    end
        jS1 = diag(jS1);
    for i=1:min(m, n)-1
        for j=1:min(m, n)-i
            if jS1(j) < jS1(j+1)
                temp = jU1(:, j);
                jU1(:, j) = jU1(:, j+1);
                jU1(:, j+1) = temp;
                temp = jS1(j);
                jS1(j) = jS1(j+1);
                jS1(j+1) = temp;
                temp = jV1(:, j);
                jV1(:, j) = jV1(:, j+1);
                jV1(:, j+1) = temp;
            end
        end
    end
    temp = jS1;
    jS1 = zeros(min(m, n), min(m, n));
    for i=1:min(m, n)
        jS1(i, i) = temp(i);
    end
end