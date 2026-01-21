function [p, dp] = schmidt_legendre_functions(n_max, theta)
costh = cos(theta);
sinth = sin(theta);

p = zeros(n_max + 2, n_max + 2);
dp = zeros(n_max + 2, n_max + 2);

p(1, 1) = 1.0;
dp(1, 1) = 0.0;

p(2, 1) = costh;
p(2, 2) = sinth;
dp(2, 1) = -sinth;
dp(2, 2) = costh;

for n = 2:n_max
    for m = 0:n
        ni = n + 1;
        mi = m + 1;
        if m > n - 2
            if m == n
                p(ni, ni) = sinth * sqrt((2 * n - 1) / (2 * n)) * p(ni - 1, ni - 1);
                dp(ni, ni) = sqrt((2 * n - 1) / (2 * n)) * (costh * p(ni - 1, ni - 1) + sinth * dp(ni - 1, ni - 1));
            else
                p(ni, ni - 1) = sqrt(2 * n - 1) * costh * p(ni - 1, ni - 1);
                dp(ni, ni - 1) = sqrt(2 * n - 1) * (costh * dp(ni - 1, ni - 1) - sinth * p(ni - 1, ni - 1));
            end
        else
            a_nm = sqrt((2 * n - 1)^2 / (n^2 - m^2));
            b_nm = sqrt(((n - 1)^2 - m^2) / (n^2 - m^2));
            p(ni, mi) = a_nm * costh * p(ni - 1, mi) - b_nm * p(ni - 2, mi);
            dp(ni, mi) = a_nm * (costh * dp(ni - 1, mi) - sinth * p(ni - 1, mi)) - b_nm * dp(ni - 2, mi);
        end
    end
end
end

