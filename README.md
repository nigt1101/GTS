A=[10.1000 0.3000 0.1000 0.3000 0.5000 0.3000 0.4000 0.2000 0.3000 0.5000;		
0.4000 20.1000 0.5000 0.5000 0.4000 0.1000 0.2000 0.2000 0.2000 0.4000;		
0.2000 0.2000 30.1000 0.1000 0.2000 0.5000 0.2000 0.5000 0.4000 0.3000;		
0.3000 0.5000 0.4000 40.2000 0.3000 0.5000 0.3000 0.1000 0.1000 0.3000;		
0.1000 0.1000 0.5000 0.1000 50.3000 0.3000 0.1000 0.1000 0.4000 0.3000;		
0.4000 0.5000 0.5000 0.1000 0.1000 60.3000 0.1000 0.1000 0.1000 0.2000;		
0.2000 0.3000 0.1000 0.5000 0.2000 0.2000 70.5000 0.4000 0.2000 0.3000;		
0.4000 0.5000 0.2000 0.3000 0.1000 0.5000 0.5000 80.4000 0.4000 0.3000;		
0.4000 0.1000 0.2000 0.3000 0.1000 0.2000 0.3000 0.4000 90.4000 0.5000;		
0.4000 0.3000 0.5000 0.1000 0.2000 0.1000 0.1000 0.3000 0.1000 100.4000];

% Jacobi
    n = size(A, 1);
    T = diag(diag(A));
    Tinv = diag(diag(A).^(-1));
    
    B = eye(n) - Tinv*A;
    q = norm(B, 1);
    ma = max(abs(diag(A)));
    mi = min(abs(diag(A)));
    landa = ma/mi;
    
    X0 = zeros(n, n);
    X = B*X0 + Tinv;
    while norm(X - X0, 1) > e*(1-q)*(1/q)*(1/landa)
        X0 = X;
        X = B*X0 + Tinv;
    end

% Newton neu cheo troi cot

    n = size(A, 1);
    %T = diag(diag(A));
    Tinv = diag(diag(A).^(-1));
    E = eye(n);
    Xdau = Tinv;
    X0 = Tinv;
    B = eye(n) - Tinv*A;
    q = norm(E - A*X0);
    e = 0.001;
    k = 1
    X = X0*(2*E - A*X0);
    while norm(Xdau, 1)*(q^(2^k))/(1-q) > e
        X0 = X;
        X = X0*(2*E - A*X0);
        k=k+1;
    end

% Gauss - Seidal
n = size(A, 1);
   
    T = diag(diag(A));
    U = triu(A) - T;
    L = tril(A) - T; 
    Tinv = diag(diag(A).^(-1));
    E = eye(n);

    X0 = Tinv;
    
    s = 0;
    for j = 1 : n
        tmp = 0
        for i = j+1 : n
            tmp = tmp + abs(A(i, j));
        end
        if s < tmp
            s = tmp
        end
    end

    l = 0;
    for j = 1 : n 
        tmp1 = 0;
        tmp2 = 0;
        for i = 1 : j 
            tmp1 = tmp1 + abs(A(i,j))
        end
        for i = j+1 : n 
            tmp2 = tmp2 + abs(A(i,j))
        end
        if l < tmp1/(1-tmp2)
            l = tmp1/(1-tmp2)
        end
    end
    X = -Tinv*L*X0 - Tinv*U*X + Tinv;
    e = 0.00001
    while norm(X - X0, 1)*s/((1-s)*(1-l)) > e
        X0 = X;
        X = -Tinv*L*X0 - Tinv*U*X + Tinv;
    end

    

