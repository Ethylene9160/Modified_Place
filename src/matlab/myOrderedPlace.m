
function K = myOrderedPlace(A,B,P,order)
    n = size(B,2);
    
    
    if size(B,2) ~= size(order,2)
        error("order haven't been given!");
    end
    T0 = zeros(n,n);
    for i = 1:n
        T0(order(i),i) = 1;
    end
    
    B = B*T0;
    K = myPlace(A,B,P);    
    K = T0*K;
end