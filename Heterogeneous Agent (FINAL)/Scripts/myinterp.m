function y = myinterp(xi, yi, x)

    y =zeros(size(x));
    
    dx = xi(2)-xi(1);
    
    xmin = min(xi);
    xmax = max(xi);

    idx = (x>xmin) & (x<xmax);
    n = floor( ( x(idx) - xmin)/dx) + 1;
    
    y(idx) = (yi(n+1)-yi(n))./(xi(n+1)-xi(n)) .* ( x(idx) - xi(n)) + yi(n);
    
    idx = (x<=xmin);
    y(idx) = (yi(1)-yi(2))./(xi(1)-xi(2)) .*(x(idx)-xi(2)) + yi(2);
    
    idx = (x>=xmax);
    y(idx) = (yi(end)-yi(end-1))./(xi(end)-xi(end-1)) .*(x(idx)-xi(end-1)) + yi(end-1);
    

end

