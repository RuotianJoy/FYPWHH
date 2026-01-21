this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
addpath(fullfile(this_dir, 'pure_matlab_igrf'));

start_year = 1900;
end_year = 2020;
year_step = 10;

lat_min = -90;
lat_max = 90;
lat_step = 5;
lon_min = -180;
lon_max = 180;
lon_step = 10;
altitude = 0;

lats = lat_min:lat_step:lat_max;
lons = lon_min:lon_step:lon_max;
years = start_year:year_step:end_year;

fprintf('Calculating non-dipole field from %d to %d with a step of %d...\n', start_year, end_year, year_step);

for y = 1:length(years)
    year = years(y);
    fprintf('Processing year: %d\n', year);
    Z_nd_grid = zeros(length(lats), length(lons));
    for i = 1:length(lats)
        lat = lats(i);
        for j = 1:length(lons)
            lon = lons(j);
            [~, ~, Z_nd] = get_non_dipole_field(year, lat, lon, altitude);
            Z_nd_grid(i, j) = Z_nd;
        end
    end
    filename = sprintf('non_dipole_Z_%d.mat', year);
    save(filename, 'Z_nd_grid');
end

fprintf('Calculation complete. Data saved in .mat files.\n');

