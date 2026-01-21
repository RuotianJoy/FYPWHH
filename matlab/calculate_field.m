function [X, Y, Z] = calculate_field(g, h, lat, lon, alt, n_max)
a = 6371.2;
r = alt + a;
colat = 90.0 - lat;

phi = deg2rad(lon);
theta = deg2rad(colat);

[P, dP] = schmidt_legendre_functions(n_max, theta);

X = 0;
Y = 0;
Z = 0;

cos_m_phi = zeros(1, n_max + 1);
sin_m_phi = zeros(1, n_max + 1);
for m = 0:n_max
    cos_m_phi(m + 1) = cos(m * phi);
    sin_m_phi(m + 1) = sin(m * phi);
end

sin_theta = sin(theta);

for n = 1:n_max
    ratio = (a / r)^(n + 2);
    for m = 0:n
        g_nm = g(n + 1, m + 1);
        h_nm = h(n + 1, m + 1);

        Z = Z - (n + 1) * ratio * (g_nm * cos_m_phi(m + 1) + h_nm * sin_m_phi(m + 1)) * P(n + 1, m + 1);

        X = X + ratio * (g_nm * cos_m_phi(m + 1) + h_nm * sin_m_phi(m + 1)) * dP(n + 1, m + 1);

        if m > 0
            if abs(sin_theta) < 1e-10
            else
                Y = Y + ratio * m * (g_nm * sin_m_phi(m + 1) - h_nm * cos_m_phi(m + 1)) * P(n + 1, m + 1) / sin_theta;
            end
        end
    end
end
end

