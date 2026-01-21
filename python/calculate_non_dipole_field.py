import numpy as np
from igrf13 import get_coeffs
import math
import argparse


def schmidt_legendre_functions(n_max, theta):
    costh = math.cos(theta)
    sinth = math.sin(theta)

    p = np.zeros((n_max + 2, n_max + 2))
    dp = np.zeros((n_max + 2, n_max + 2))

    p[0, 0] = 1.0
    dp[0, 0] = 0.0

    p[1, 0] = costh
    p[1, 1] = sinth
    dp[1, 0] = -sinth
    dp[1, 1] = costh

    for n in range(2, n_max + 1):
        for m in range(n + 1):
            if m > n - 2:
                if m == n:
                    p[n, n] = sinth * math.sqrt((2 * n - 1) / (2 * n)) * p[n - 1, n - 1]
                    dp[n, n] = math.sqrt((2 * n - 1) / (2 * n)) * (costh * p[n - 1, n - 1] + sinth * dp[n - 1, n - 1])
                else:
                    p[n, n - 1] = math.sqrt(2 * n - 1) * costh * p[n - 1, n - 1]
                    dp[n, n - 1] = math.sqrt(2 * n - 1) * (costh * dp[n - 1, n - 1] - sinth * p[n - 1, n - 1])
            else:
                a_nm = math.sqrt((2 * n - 1) ** 2 / (n ** 2 - m ** 2))
                b_nm = math.sqrt(((n - 1) ** 2 - m ** 2) / (n ** 2 - m ** 2))
                p[n, m] = a_nm * costh * p[n - 1, m] - b_nm * p[n - 2, m]
                dp[n, m] = a_nm * (costh * dp[n - 1, m] - sinth * p[n - 1, m]) - b_nm * dp[n - 2, m]

    return p, dp


def calculate_field(g, h, lat, lon, alt, n_max):
    a = 6371.2
    r = alt + a
    colat = 90.0 - lat

    phi = math.radians(lon)
    theta = math.radians(colat)

    P, dP = schmidt_legendre_functions(n_max, theta)

    X, Y, Z = 0, 0, 0

    cos_m_phi = [math.cos(m * phi) for m in range(n_max + 1)]
    sin_m_phi = [math.sin(m * phi) for m in range(n_max + 1)]

    sin_theta = math.sin(theta)

    for n in range(1, n_max + 1):
        ratio = (a / r) ** (n + 2)
        for m in range(n + 1):
            g_nm = g.get((n, m), 0)
            h_nm = h.get((n, m), 0)

            Z -= (n + 1) * ratio * (g_nm * cos_m_phi[m] + h_nm * sin_m_phi[m]) * P[n, m]

            X += ratio * (g_nm * cos_m_phi[m] + h_nm * sin_m_phi[m]) * dP[n, m]

            if m > 0:
                if abs(sin_theta) < 1e-10:
                    pass
                else:
                    Y += ratio * m * (g_nm * sin_m_phi[m] - h_nm * cos_m_phi[m]) * P[n, m] / sin_theta

    return X, Y, Z


def get_non_dipole_field(year, lat, lon, alt=0):
    g, h = get_coeffs(year)
    n_max = max([k[0] for k in g.keys()])

    X_total, Y_total, Z_total = calculate_field(g, h, lat, lon, alt, n_max)

    g_dipole = {k: v for k, v in g.items() if k[0] == 1}
    h_dipole = {k: v for k, v in h.items() if k[0] == 1}
    X_dipole, Y_dipole, Z_dipole = calculate_field(g_dipole, h_dipole, lat, lon, alt, 1)

    X_nd = X_total - X_dipole
    Y_nd = Y_total - Y_dipole
    Z_nd = Z_total - Z_dipole

    return X_nd, Y_nd, Z_nd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate the non-dipole magnetic field.")
    parser.add_argument('--start', type=int, default=1900, help='Start year')
    parser.add_argument('--end', type=int, default=2020, help='End year')
    parser.add_argument('--step', type=int, default=10, help='Year step')
    args = parser.parse_args()

    start_year = args.start
    end_year = args.end
    year_step = args.step
    lat_min, lat_max, lat_step = -90, 90, 5
    lon_min, lon_max, lon_step = -180, 180, 10
    altitude = 0

    lats = np.arange(lat_min, lat_max + lat_step, lat_step)
    lons = np.arange(lon_min, lon_max + lon_step, lon_step)
    years = np.arange(start_year, end_year + 1, year_step)

    print(f"Calculating non-dipole field from {start_year} to {end_year} with a step of {year_step}...")

    results = {}

    for year in years:
        print(f"Processing year: {year}")
        Z_nd_grid = np.zeros((len(lats), len(lons)))
        for i, lat in enumerate(lats):
            for j, lon in enumerate(lons):
                _, _, Z_nd = get_non_dipole_field(year, lat, lon, altitude)
                Z_nd_grid[i, j] = Z_nd
        results[year] = Z_nd_grid
        np.save(f'non_dipole_Z_{year}.npy', Z_nd_grid)

    print("Calculation complete. Data saved in .npy files.")

