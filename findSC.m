% This file is found in the netalign package but there's a small bug which has
% been corrected here.
function [SC] = findSC(S, li, lj)
    n = max(li);
    m = max(lj);
    ne = size(S, 1);
    %flag = zeros(ne, 1);  % Fix by Eric.
    flag = zeros(max(max(li), max(lj)), 1);
    [rp ci] = sparse_to_csr(S);
    sind = zeros(size(ci));
    ri = zeros(size(ci));
    SC = zeros(size(ci, 1), 4);
    tot = 0;
    for i = 1:ne
        for j = rp(i):rp(i+1) - 1
            ri(j) = i;
            k = li(ci(j));
            if (flag(k) == 0)
                tot = tot + 1;
                flag(k) = tot;
            end
            sind(j) = flag(k);
        end
        for j = rp(i):rp(i+1) - 1
            k = li(ci(j));
            flag(k) = 0;
        end
    end
    SC(:,1) = sind;
    tot = 0;
    for i = 1:ne
        for j = rp(i):rp(i+1) - 1
            k = lj(ci(j));
            if (flag(k) == 0)
                tot = tot + 1;
                flag(k) = tot;
            end
            sind(j) = flag(k);
        end
        for j = rp(i):rp(i+1) - 1
            k = lj(ci(j));
            flag(k) = 0;
        end
    end
    SC(:,2) = sind;
    St = sparse(ri, ci, SC(:,1));
    [i j k] = find(St);
    SC(:,3) = k;
    St = sparse(ri, ci, SC(:,2));
    [i j k] = find(St);
    SC(:,4) = k;
end
