% Just a beta version!
function K = myPoleMIMOPlacementRecursion(A,B,P)
    if size(B,2) == 1
        K = myPoleSISOPlacement(A,B,P)
    else
        disp(B)
        if size(B,2) == 1
            Bi = B
        else
            Bi = B(:,1:1)
        end
        [T, A_c, B_c, A_uc, B_uc] = controllabilityDecomposition(A,Bi)
        r = size(A_c, 1)
        if r ~= 0
            K_up = myPoleSISOPlacement(A_c, B_c, P(1:r))
        else
            K_up=[]
        end
        A = A_uc
        %B = B_uc(:,2:end)
        B = B(r+1:end,2:end);
        P = P(r+1:end);
        r;
        s = size(A,1);
        if size(A, 1) == 0
            str = "over!"
            %k_up_size= size(K_up,2)
            %K = [K_up;zeros(1,size(K_up,2))]
            K_up = [K_up;zeros(1,size(K_up,2))]
            K = K_up*inv(T)
            
        else
            
            str = "continue!"
            K_down = myPoleMIMOPlacementRecursion(A,B,P)
            
            K = [K_up, zeros(size(K_up,1), size(K_down,2));
                zeros(size(K_down,1), size(K_up, 2)), K_down]
            
            K=K*inv(T)
        end
    end
end