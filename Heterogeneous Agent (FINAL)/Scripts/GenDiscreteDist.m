function p = GenDiscreteDist(x, mu, sigma)

    p = zeros(size(x));
    
    lowx_mid_pt = 0.5*(x(1)+x(2));
    p(1) = normcdf( lowx_mid_pt, mu, sigma);

    mid_pt = 0.5*(x(2:end) + x(1:end-1));
    
    p(2:end-1) = normcdf( mid_pt(2:end), mu, sigma) ...
        -  normcdf( mid_pt(1:end-1), mu, sigma);
    p(end) = 1-normcdf( mid_pt(end), mu, sigma);

end

