function [coeff_vect,integrand_mat,weights_vect] = gauss_quadrature(tau_upper,tau_lower,N,legendre_general)
    
    points = legendre_general{N}(1,:);
    weights_vect = legendre_general{N}(2,:);
   

    [tau_diff_mat,points_mat] = meshgrid(tau_upper-tau_lower,points);
    [tau_sum_mat,~] = meshgrid(tau_upper+tau_lower,points);
    integrand_mat = (tau_diff_mat)/2.*points_mat + (tau_sum_mat)/2.*ones(size(points_mat));
    coeff_vect = (tau_upper-tau_lower)/2;
end
