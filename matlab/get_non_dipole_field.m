function [X_nd, Y_nd, Z_nd] = get_non_dipole_field(year, lat, lon, alt)
if nargin < 4
    alt = 0;
end

[g, h] = get_coeffs_igrf13_pure(year);

n_max = size(g, 1) - 1;

[X_total, Y_total, Z_total] = calculate_field(g, h, lat, lon, alt, n_max);

g_dipole = zeros(size(g));
h_dipole = zeros(size(h));
g_dipole(2, 1:2) = g(2, 1:2);
h_dipole(2, 1:2) = h(2, 1:2);

[X_dipole, Y_dipole, Z_dipole] = calculate_field(g_dipole, h_dipole, lat, lon, alt, 1);

X_nd = X_total - X_dipole;
Y_nd = Y_total - Y_dipole;
Z_nd = Z_total - Z_dipole;
end

