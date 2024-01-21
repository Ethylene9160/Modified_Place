function K = myAdvancedRandomPlace(A,B,P)
    n = size(B,2);
    T = zeros(n,n);
    while rank(T) ~= n
        for i = 1:n
            for j = 1:n
                T(i,j) = getRandomNumber(0.6);
            end
        end
    end
    
    
    
    B2 = B*T;
    K2 = myPlace(A,B2,P);
    K = T*K2;
end

function rd = getRandomNumber(threhold)
    rd = rand(1);
    if rd >= threhold
        rd = rd/threhold-1+0.05;
    else
        rd = 0;
    end
end