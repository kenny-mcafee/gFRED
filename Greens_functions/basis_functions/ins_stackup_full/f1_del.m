function f_del_out = f1_del(x,N,L_s,L_RTV)
    f_del_out = (2*(1:N)'-2).*(x+L_s+L_RTV).^(2*(1:N)'-3);
 
end