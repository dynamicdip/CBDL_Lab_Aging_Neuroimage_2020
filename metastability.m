% calculate metastabilty of n modes.
% input is Nmodes * Mphase
function m = metastability(phase_matrix)
phase_matrix = exp(1i*phase_matrix);
order_parameter = mean(phase_matrix);
m = std(abs(order_parameter));






