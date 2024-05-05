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
% Giai dung ma tran nghich dao bang Cholesky
    A=[1 1 -3 2;
        1 -2 0 -1;
        0 1 1 3;
        2 -3 2 0]
    B= eye(4)
    function isSymmetric = isMatrixSymmetric(matrix)
        % Kiểm tra số hàng và số cột của ma trận
        [rows, cols] = size(matrix);
    
    % Kiểm tra nếu ma trận không phải là ma trận vuông
    if rows ~= cols
        disp('Ma trận không phải là ma trận vuông');
        isSymmetric = false;
        return;
    end
    
    % Kiểm tra ma trận đối xứng
    isSymmetric = isequal(matrix, matrix');
    end
    
    function X = Cholesky(A,B)
        n = size(A, 1)
        Q = zeros(n)
        Q(1,1) = sqrt(A(1,1))
        Q(1,2:n) = A(1,2:n)/Q(1,1)
        for i = 2:n
            for j = i:n
                tmp = 0
                tmp2 = 0
                for k = 1:i-1
                    tmp = tmp + Q(k,i)^2
                    tmp2 = tmp2 + Q(k,i)*Q(k,j)
                end
            Q(i,i) = sqrt(A(i,i) - tmp)
            Q(i,j) = (A(i,j) - tmp2)/Q(i,i)
            end
        end
        QT = transpose(Q)
        Y = zeros(n,1)
        Y(1,1) = B(1,1)/Q(1,1)
        for i = 2:n
            tmp=0
            for k = 1:i-1
                tmp = tmp + Q(k,i)*Y(k,1)
            end
            Y(i,1) = (B(i,1)-tmp)/Q(i,i)
        end
        X = zeros(n,1)
        X(n,1) = Y(n,1)/Q(n,n)
        for i = n-1:-1:1
            tmp1 = 0
            for k = i+1:n
                tmp1 = tmp1 + Q(i,k)*X(k,1)
            end
            X(i,1) = (Y(i,1) - tmp1)/Q(i,i)
        end
    end
    
    if isMatrixSymmetric(A)
        X1 = Cholesky(A,B(1:4,1))
        X2 = Cholesky(A,B(1:4,2))
        X3 = XCholesky(A,B(1:4,3))
        X4 = Cholesky(A,B(1:4,4))
        X = [X1 X2 X3 X4]
    else
        M = transpose(A)*A 
        X1 = Cholesky(M,transpose(A)*B(1:4,1));
        X2 =  Cholesky(M,transpose(A)*B(1:4,2));
        X3 =  Cholesky(M,transpose(A)*B(1:4,3));
        X4 =  Cholesky(M,transpose(A)*B(1:4,4));
        X = [X1 X2 X3 X4]
    end

% Tim nghiem dung bang Cholesky
    A=[1 1 -3 2;
    1 -2 0 -1;
    0 1 1 3;
    2 -3 2 0]
    B=[6;-6;16;6]
    function isSymmetric = isMatrixSymmetric(matrix)
        % Kiểm tra số hàng và số cột của ma trận
        [rows, cols] = size(matrix);
        
    % Kiểm tra nếu ma trận không phải là ma trận vuông
    if rows ~= cols
        disp('Ma trận không phải là ma trận vuông');
        isSymmetric = false;
        return;
    end
    
    % Kiểm tra ma trận đối xứng
    isSymmetric = isequal(matrix, matrix');
    end
    
    function Cholesky(A,B)
        n = size(A, 1)
        Q = zeros(n)
        Q(1,1) = sqrt(A(1,1))
        Q(1,2:n) = A(1,2:n)/Q(1,1)
        for i = 2:n
            for j = i:n
                tmp = 0
                tmp2 = 0
                for k = 1:i-1
                    tmp = tmp + Q(k,i)^2
                    tmp2 = tmp2 + Q(k,i)*Q(k,j)
                end
            Q(i,i) = sqrt(A(i,i) - tmp)
            Q(i,j) = (A(i,j) - tmp2)/Q(i,i)
            end
        end
        QT = transpose(Q)
        Y = zeros(n,1)
        Y(1,1) = B(1,1)/Q(1,1)
        for i = 2:n
            tmp=0
            for k = 1:i-1
                tmp = tmp + Q(k,i)*Y(k,1)
            end
            Y(i,1) = (B(i,1)-tmp)/Q(i,i)
        end
        X = zeros(n,1)
        X(n,1) = Y(n,1)/Q(n,n)
        for i = n-1:-1:1
            tmp1 = 0
            for k = i+1:n
                tmp1 = tmp1 + Q(i,k)*X(k,1)
            end
            X(i,1) = (Y(i,1) - tmp1)/Q(i,i)
        end
    end
    
    if isMatrixSymmetric(A)
        Cholesky(A)
    else
        M = transpose(A)*A
        N = transpose(A)*B 
        Cholesky(M,N);
    end

% giai dung nghiem cua hpt bang gauss
    function B = GaussThuan(A, b)
        n = size(A, 1);
        B=[A, b];
        for k = 1: n-1
            if B(k, k) == 0
                disp('loi chia cho 0')
                return
            else
                for i = k+1:n
                    B(i,:) = B(i, :) - B(k, :)*B(i,k)/B(k, k);
                end
            end
        end
    end
    function nghiem = Gaussnghich(B)
        %B la ma tran tam giac tren
        n = size(B, 1);
        nghiem = zeros(n, 1);
         for k =n:-1:1
                if k == n
                    nghiem(n) = B(n,n+1)/B(n,n);
                else
                    S = B(k,n+1)-sum(B(k,1:n)*nghiem(1:n,1))
                    nghiem(k,1) = S/B(k,k);
                end   
         end
     end
    
    C = GaussThuan(A, B)
    X = zeros(size(A,1),1)
    X = Gaussnghich(C)
