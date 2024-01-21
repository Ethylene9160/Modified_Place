function K = myBestRisePlace(A,B,P,fIp)
    iSize = size(B,2);
    if iSize < fIp
        error("cannot get access to the colun you refer!");
    end
    P_abs = abs(P)
    index = 0;
    pmax = 0;
    for i = 1:iSize
        if pmax < P_abs(i)
            pmax = P_abs(i);
            index = i;
        end
    end
    pt = P(index);
    P(index) = P(fIp);
    P(fIp) = pt;
    K = myPlace(A,B,P);
end