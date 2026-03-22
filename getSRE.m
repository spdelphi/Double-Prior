function sre = getSRE(x,abundance)
    % sre = 10*log10(norm(x,'fro')/norm(abundance-x,'fro'));
    sre = 10*log10( sum(x(:).^2) / sum( (x(:) - abundance(:)).^2 ) );
end