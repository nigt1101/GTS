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

% phan bu dai so tim dung ma tran nghich dao

    A=[1 5;
        -2 4]
    
    n = size(A,1)
    C = zeros(n)
    for i = 1:n
        for j = 1:n
            Aij = [A(1:i-1,1:j-1) A(1:i-1,j+1:n);
                    A(i+1:n,1:j-1) A(i+1:n,j+1:n)]
            
            C(i,j) = power(-1,i+j)*det(Aij)
        end
    end
    Anghich = transpose(C)/det(A)

% vien quanh tim dung ma tran nghich dao

    A = [50 107 36;
        25 54 20;
        31 66 21]
    n = size(A,1)
    function Ainv = VienQuanh(A)
        a = 1/A(1,1)
        n = size(A,1);
        for i = 2:n
                a = VienQuanh1(A(1:i,1:i),a);
        end
        Ainv = a;
    end
    function S = VienQuanh1(A, a)
        %A la ma tran can tim nghich dao
        %a la nghich dao cua ma tran (n-1 x n-1)
        n =size(A,1);
        %A1 = A(1:n-1,1:n-1);
        U = A(n,1:n-1);
        V = A(1:n-1,n);
        ann = A(n,n);
        beta12 = a*V*(U*a*V - ann)^(-1);
        beta11 = (eye(n-1) - beta12*U)*a;
        beta22 = (ann -U*a*V)^(-1);
        beta21 = -beta22*U*a;
        S=[beta11, beta12; beta21, beta22];
    end
    X = zeros(n)
    X = VienQuanh(A)

% giai dung hpt bang gaussjordan

    A=[1 1 -3 2 6;
    1 -2 0 -1 -6;
    0 1 1 3 16;
    2 -3 2 0 6] %ma tran bo sung
    n = size(A,1)
    function x = GJordan(A)
    
    [m,n] = size(A);
    
    i = 1;
    j = 1;

    while i <= m && j <= n
       [p, k] = max(abs(A(i:m,j)));
       k = k+i-1;
       % Swap i-th and k-th rows.
       A([i k],j:n) = A([k i],j:n);
       % Divide the pivot row by the pivot element.
       A(i,j:n) = A(i,j:n)./A(i,j);
       % Subtract multiples of the pivot row from all the other rows.
       for k = [1:i-1 i+1:m]
           A(k,j:n) = A(k,j:n) - A(k,j).*A(i,j:n);
       end
       i = i + 1;
       j = j + 1;
       %disp(A);
    end
    x = A(1:m,n);
    end
    
    X = zeros(n,1)
    X = GJordan(A)

% thuật toán tổng quát sử dụng phương pháp lũy thừa và xuống thang

    function check = KTsongsong(u,v,tol)
    check = (norm( u./norm(u,2)-v./norm(v,2), 2) <= tol) || (norm( u./norm(u,2)+v./norm(v,2), 2) <= tol);
    end
    
    function [lambda, v] = powerMethod(A,tol,M)
    if nargin == 1
        tol = 1e-10;
        M = 200;
    end
    if nargin == 2
        M = 200;
    end
    
    %tol là giá trị hiệu chỉnh
    %M số lần lặp nhiều nhất cho phép
    
    n = size(A,1);
    x = ones(n,1);
    m = 1;
    check = false;
    lambda = [];
    v =zeros(n,0);
    % truong hop 1 gtr troi
    y1 = x;
    while check == false && (m <= M)
        y1 = A*y1;
        y2 = A*y1;
        m= m+1;
        check = KTsongsong(y1,y2,tol);
    end
    if check == true
        disp('y1 song song voi y2')
        lambda = [lambda; mean(y2(y2~=0)./y1(y1~=0))];
        v1 = y1./norm(y1,2);
        v = [v,v1];
        disp(lambda)
        disp(v)
        return
    end
    %truong hop 2 nghiem doi nhau
    if check == false
        m = 1;
        y1 = x;
        while check ==false && (m <= M)
            y1 = A*y1;
            y2 = A*y1;
            y3 = A*y2;
            m= m+1;
            check = KTsongsong(y1,y3,tol);
        end
    if check == true
        disp('y1 song song voi y3')
        lambda1 = sqrt(mean(y3(y3~=0)./y1(y1~=0)));
        lambda = [lambda;lambda1;-lambda1];
        v1 = y2+lambda1*y1;
        v1 = v1/norm(v1,2);
        v2 = y2-lambda1*y1;
        v2 = v2./norm(v2,2);
        v=[v,v1,v2];
        disp(lambda)
        disp(v)
        return
    end
    end
    %truong hop 2 nghiem phuc lien hop
    if check == false
        disp('Khong co vecto nao song song')
         % Tính y1, y2 và y3 bằng cách lũy thừa ma trận và nhân với vector x
        y1 = ((A^(2*m+2)) * x);
        y2 = ((A^(2*m)) * x);
        y3 = ((A^(2*m+1)) * x);
    
    % Tìm các chỉ số của các phần tử khác 0 trong y1
    index = find(y1 ~= 0);
    
    % Kiểm tra xem index có đủ phần tử không
    if length(index) < 2
        error('Không tìm thấy đủ các phần tử khác 0 trong y1.');
    end
    
    % Lấy các chỉ số đầu tiên và thứ hai
    r = index(1);
    s = index(2);

    end
        syms z
        p = det([1, y2(r), y2(s);...
              z, y3(r), y3(s); z^2 y1(r) y1(s)]);
        lambda = double(solve(p,z));
        v1 = y3-lambda(2)*y1;
        v1 = v1/norm(v1,2);
        v2 = y3-lambda(1)*y1;
        v2 = v2/norm(v2,2);
        v= [v1,v2];
        disp(lambda)
        disp(v)
    end
    
    
    function [A1] = XuongThang(A,v)
    %tra ra ma tran A1 da xuong thang
    [~,i] = max(abs(v));
    v = v/v(i);
    n = size(A);
    theta = eye(n);
    theta(:,i) = theta(:,i) - v;
    A1 = theta*A;
    end
    
    function [eigvalue, v] = ttTongQuat(A)
        eigvalue = [];
        v = [];
        n = size(A, 1);
    
    
    i = 1;
    while i <= n
        disp('Lan lap thu :')
        disp(i)
        [eigv, x] = powerMethod(A);
        
        cnt = numel(x);
        if cnt == 2
            i = i + 2;
            A = XuongThang(A, x(:, 1));
            A = XuongThang(A, x(:, 2));
            
            disp('Ma tran sau khi xuong thang')
            disp(A)
        else
            i = i + 1;
            A = XuongThang(A, x(:, 1));
            
            disp('Ma tran sau khi xuong thang')
            disp(A)
        end
        v = [v x];
        eigvalue = [eigvalue; eigv];
    end
    end
    
    % Định nghĩa ma trận A
    A = [[2 3 2],
        [4 3 5],
        [3 2 9]];
    
    % Gọi hàm powerMethod để tính trị riêng và vector riêng
    [eigvalue, v] = ttTongQuat(A);
    disp(eigvalue)
    disp(v)
    
    % Tính giá trị riêng và vector riêng của ma trận A
    [V, D] = eig(A);
    
    % D là ma trận đường chéo chứa các giá trị riêng
    disp('Ma trận giá trị riêng:');
    disp(D);
    
    % V là ma trận chứa các vector riêng ứng với các giá trị riêng tương ứng
    disp('Ma trận vector riêng:');
    disp(V);


    
