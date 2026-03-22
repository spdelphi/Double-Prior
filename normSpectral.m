function [normLib,max_le,min_le] = normSpectral(lib,stretch)
    if nargin<2 || isempty(stretch)
        stretch = 0;
    else
        stretch = 1;
    end
    
    [ bands,edNum]  = size( lib );
    para_length = zeros(1,edNum);
    for i = 1:edNum
        para_length(i) = norm(lib(:,i)); 
    end
    normLib = lib./ para_length;
    if stretch == 0
        normLib = normLib;
    else 
        normLib = normLib * max(para_length);
    end
    max_le = max(para_length);
    min_le = min(para_length);
end