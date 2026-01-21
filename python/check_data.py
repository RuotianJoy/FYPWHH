import numpy as np
import matplotlib.pyplot as plt


def inspect_data(year):
    try:
        filepath = f'non_dipole_Z_{year}.npy'
        data = np.load(filepath)
        print(f"--- Inspection for year {year} ---")
        print(f"File: {filepath}")
        print(f"Data shape: {data.shape}")
        print(f"Min value: {np.min(data):.2f} nT")
        print(f"Max value: {np.max(data):.2f} nT")
        print(f"Mean value: {np.mean(data):.2f} nT")
        print(f"Standard deviation: {np.std(data):.2f} nT")

        plt.figure()
        plt.hist(data.flatten(), bins=50)
        plt.title(f'Histogram of Z component for {year}')
        plt.xlabel('Field Strength (nT)')
        plt.ylabel('Frequency')
        plt.grid(True)
        plt.savefig(f'histogram_{year}.png')
        plt.close()
        print(f"Saved histogram to histogram_{year}.png")
        print("-" * 30 + "\n")

    except FileNotFoundError:
        print(f"ERROR: Data file for year {year} not found!")


if __name__ == '__main__':
    inspect_data(1900)
    inspect_data(1960)
    inspect_data(2020)

