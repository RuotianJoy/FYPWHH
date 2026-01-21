
# -*- coding: utf-8 -*-
"""
The International Geomagnetic Reference Field (IGRF) is a standard mathematical
description of the Earth's main magnetic field. It is updated every five years.
The 13th generation IGRF was released in December 2019.

This file provides the IGRF-13 coefficients.
Data is from the IAGA Division V, Working Group V-MOD:
https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html

The coefficients are provided for epochs from 1900 to 2020 (every 5 years),
and a predictive secular variation for 2020-2025.

Format:
Each row contains the Gauss coefficients g_n^m and h_n^m for a specific
degree n and order m.
Columns are: n, m, g_1900, h_1900, g_1905, h_1905, ..., g_2020, h_2020, g_sv, h_sv
Units are in nanoTeslas (nT) for coefficients and nT/year for secular variation.
"""

import csv


def load_coeffs_from_csv(filepath):
    coeffs = []
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)

        year_columns = [str(y) for y in range(1900, 2021, 5)]
        sv_column_name = '2025-2030'

        try:
            year_indices = [header.index(y) for y in year_columns]
            sv_index = header.index(sv_column_name)
        except ValueError as e:
            raise ValueError(f"CSV header is missing an expected column or has an unexpected format: {e}")

        for row in reader:
            if not row or not row[0]:
                continue

            g_or_h = row[0]
            n = int(row[1])
            m = int(row[2])

            new_row = [n, m, g_or_h]

            for idx in year_indices:
                new_row.append(float(row[idx]))

            new_row.append(float(row[sv_index]))

            coeffs.append(new_row)

    return coeffs


COEFFS = load_coeffs_from_csv('d:\\papertranslate\\PT\\Data\\IGRF14coeffsData.csv')

YEARS = [1900.0, 1905.0, 1910.0, 1915.0, 1920.0, 1925.0, 1930.0, 1935.0, 1940.0,
         1945.0, 1950.0, 1955.0, 1960.0, 1965.0, 1970.0, 1975.0, 1980.0, 1985.0,
         1990.0, 1995.0, 2000.0, 2005.0, 2010.0, 2015.0, 2020.0]


def get_coeffs(year):
    if not (1900 <= year <= 2025):
        raise ValueError("Year must be between 1900 and 2025")

    if year > 2020:
        epoch1_idx = len(YEARS) - 1
        t = year - YEARS[epoch1_idx]
    else:
        for i in range(len(YEARS) - 1):
            if YEARS[i] <= year < YEARS[i + 1]:
                epoch1_idx = i
                break
        else:
            epoch1_idx = len(YEARS) - 1

        if year == YEARS[epoch1_idx]:
            t = 0
            epoch2_idx = epoch1_idx
        else:
            epoch2_idx = epoch1_idx + 1
            t = (year - YEARS[epoch1_idx]) / (YEARS[epoch2_idx] - YEARS[epoch1_idx])

    g = {}
    h = {}

    for row in COEFFS:
        n, m, g_or_h = row[0], row[1], row[2]

        c1 = row[3 + epoch1_idx]

        if year > 2020:
            sv = row[-1]
            c = c1 + sv * t
        else:
            c2 = row[3 + epoch2_idx]
            c = c1 + t * (c2 - c1)

        if g_or_h == 'g':
            g[(n, m)] = c
        else:
            h[(n, m)] = c

    return g, h

