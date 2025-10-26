function f_out = f2(x,N,L_s,L_RTV,k_0,k_s,k_RTV)
    f_out = f15(x,N,L_s,L_RTV,k_s,k_RTV) + big_gamma_23(N,L_RTV,L_s,k_0,k_s,k_RTV).*(1-(1+x/(L_s+L_RTV)).^((1:N)'));
     
end