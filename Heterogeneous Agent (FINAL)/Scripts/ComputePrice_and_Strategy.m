function [S, g, Gamma] =  ComputePrice_and_Strategy(x, p, Y, zeta, gamma, h, P, t, K, pi_k)
% computes price and strategy, as well as numerous explanatory and
% descriptive plots
S = NaN(size(Y,2), 1);
g = NaN(length(x), size(Y,2), K);
Gamma = NaN(length(x), size(Y,2), K);



for n = 1 : size(Y,2)
    S(n) = sum(P ./ gamma .* pi_k .* reshape(sum( p(:,n, :).* Y(:,n, :), 1), [1,K])) / sum(pi_k ./ gamma);
    for k = 1 : K
        g(:,n, k) = h(k) + P/zeta(k) * Y(:,n,k);
        Gamma(:,n, k) = 1/gamma(k) *( P * Y(:,n,k) - S(n) );
    end
end



%save2pdf("summary_of_controls.pdf", f1, 600)
%save2pdf("equilibrium_S.pdf", f2, 600)

%disp(S);