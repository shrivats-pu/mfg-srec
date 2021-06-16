function yt = ComputeYt(x, ytp1, mu, sigma)

    % internal points
    alpha = (ytp1(2:end,:)-ytp1(1:end-1,:))./(x(2:end,:)-x(1:end-1,:));
    beta = ytp1(1:end-1,:) - alpha .* x(1:end-1,:);

    f = (beta + alpha .* mu(1:end-1,:) ) .* ( normcdf( (x(2:end,:)-mu(2:end,:) )/sigma) ...
                                            - normcdf( (x(1:end-1,:)-mu(1:end-1,:) )/sigma) ) ...
            - sigma*alpha.*( normpdf( (x(2:end,:)-mu(2:end,:) )/sigma) ...
                           - normpdf( (x(1:end-1,:)-mu(1:end-1,:) )/sigma) );
        
    yt = sum(f,1);
    
    % below x(1)
    alpha = (ytp1(2,:)-ytp1(1,:))./(x(2,:)-x(1,:));
    beta = ytp1(2,:) - alpha.*x(2,:);
    yt = yt + (beta + alpha .*mu(1,:) ) .*  normcdf( (x(1,:)-mu(1,:) )/sigma)  ...
        - sigma*alpha .* normpdf( (x(1,:)-mu(1,:) )/sigma);
    
    % above x(n)
    alpha = (ytp1(end,:)-ytp1(end-1,:))./(x(end,:)-x(end-1,:));
    beta = ytp1(end-1,:) - alpha .* x(end-1,:);
    yt = yt + (beta + alpha .*mu(end,:) ) .*  ( 1-normcdf( (x(end,:)-mu(end,:) )/sigma)) ...
        + sigma*alpha .* normpdf( (x(end,:)-mu(end,:) )/sigma);

    yt = yt';
end