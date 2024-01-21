function K = myRandomPlace(A,B,P)
    n = size(B,2);
    tn = 0;
    t0 = 1:n;
    t0 = t0(randperm(length(t0)));
    T = zeros(n,n);
    for i = 1:n
        while T(t0(i),i) == 0
            T(t0(i),i) = abs(rand(1))+0.1;
        end
    end
    %T
    while n ~= tn
        T2 = rand(n,n);
        [Q,R] = qr(T2);
        T2 = Q*R;
        tn = rank(T2);
    end

    for i = 2:n
        for j = 1:i-1
            T2(i,j) = 0;%转换为上三角矩阵
        end
    end
    T = T*T2;
    B2 = B*T;
    K2 = myPlace(A,B2,P);
    K = T*K2;
end