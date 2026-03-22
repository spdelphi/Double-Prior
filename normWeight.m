function [w]  = normWeight(w)
    scale = sqrt(sum(w.^2));
    idx = scale ~= 0;
    w(:,idx) = w(:,idx)./scale(idx);
end