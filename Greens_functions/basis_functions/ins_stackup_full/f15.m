function f_out = f15(x,N,L_s,L_RTV,k_s,k_RTV)
    f_out = f1(x,N,L_s,L_RTV) + big_gamma_12(N,L_s,k_RTV,k_s).*(1-(1+L_RTV/L_s+x/L_s).^((1:N)'));
end