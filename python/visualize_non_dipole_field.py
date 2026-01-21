import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import imageio
import os
import argparse


def create_animation(years, gif_name='non_dipole_Z_animation.gif', duration=0.5):
    filenames = []
    vmin = float('inf')
    vmax = float('-inf')
    print("Scanning data for color scale...")
    for year in years:
        try:
            data = np.load(f'non_dipole_Z_{year}.npy')
            vmin = min(vmin, data.min())
            vmax = max(vmax, data.max())
        except FileNotFoundError:
            print(f"Data file for year {year} not found. Skipping.")
            continue

    print(f"Generating frames from {min(years)} to {max(years)}...")
    for year in years:
        try:
            data = np.load(f'non_dipole_Z_{year}.npy')
        except FileNotFoundError:
            continue

        lats = np.linspace(-90, 90, data.shape[0])
        lons = np.linspace(-180, 180, data.shape[1])
        lons, lats = np.meshgrid(lons, lats)

        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
        ax.set_global()
        ax.coastlines()
        ax.gridlines()

        contour = ax.contourf(
            lons,
            lats,
            data,
            60,
            transform=ccrs.PlateCarree(),
            cmap='coolwarm',
            vmin=vmin,
            vmax=vmax,
        )

        cbar = plt.colorbar(contour, ax=ax, orientation='vertical', shrink=0.7)
        cbar.set_label('Non-dipole Field Z-component (nT)')

        ax.set_title(f'Global Distribution of Non-dipole Field Z-component for {year}')

        filename = f'frame_{year}.png'
        plt.savefig(filename)
        plt.close(fig)
        filenames.append(filename)

    print(f"\nCreating GIF: {gif_name}...")
    with imageio.get_writer(gif_name, mode='I', duration=duration) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)
    print(f"Animation saved to {gif_name}")

    print("Cleaning up temporary files...")
    for filename in filenames:
        os.remove(filename)
    print("Done.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create an animation of the non-dipole magnetic field.")
    parser.add_argument('--start', type=int, default=1900, help='Start year')
    parser.add_argument('--end', type=int, default=2020, help='End year')
    parser.add_argument('--step', type=int, default=10, help='Year step')
    parser.add_argument('--duration', type=float, default=0.5, help='Duration of each frame in the GIF')
    parser.add_argument('--output', type=str, default='non_dipole_Z_animation.gif', help='Output GIF filename')
    args = parser.parse_args()

    years_to_process = range(args.start, args.end + 1, args.step)
    create_animation(years_to_process, gif_name=args.output, duration=args.duration)

