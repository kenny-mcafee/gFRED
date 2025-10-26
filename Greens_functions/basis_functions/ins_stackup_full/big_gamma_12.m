function out = big_gamma_12(N,L_s,k_RTV,k_s)
    out = (L_s.^(2*(1:N)'-2)).*((2*(1:N)'-2)./((1:N)')).*(1-k_s/k_RTV);

end