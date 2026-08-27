function [Hs_ERA5,Tp_ERA5,time_ERA5,D_ERA5] = Get_ERA5_Waves(lat,lon)

% Get some wave data from ERA5 reanalysis
% cdsapi (python function) is installed in a virtual environment. First
% activate the virtual environment (WaveFlood) in VSCode

% Set Python environment — path is built relative to this file's location
repo_root = fileparts(mfilename('fullpath'));
pyenv('Version', fullfile(repo_root, 'WaveFlood', 'bin', 'python'));

% Import Python module (make sure the script is in the same folder or Python path)
era5 = py.importlib.import_module('Get_ERA5_Waves');

% Define parameters to get from ERA5. Note: the CDS API can request all
% three variables in a single Python call, but they are pulled one at a
% time here to keep the MATLAB/Python handoff simple.
parameters = {'significant_height_of_combined_wind_waves_and_swell', 'peak_wave_period', 'mean_wave_direction'};

% Most recent 30 full calendar years of ERA5 data, ending with last year
% (the current year is excluded since it isn't complete yet). E.g. run on
% any date in 2038 -> years 2008-2037.
end_year   = year(datetime('now')) - 1;
start_year = end_year - 29;
years = arrayfun(@(y) num2str(y), start_year:end_year, 'UniformOutput', false);

months = {'01','02','03','04','05','06','07','08','09','10','11','12'};

days = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15',...
    '16','17','18','19','20','21','22','23','24','25','26','27','28'};

for i = 1:length(parameters)

    params = struct();
    params.variable = parameters{i};

    % Convert to Python lists (years/months/days are computed once, above)
    params.year = py.list(years);
    params.month = py.list(months);
    params.day = py.list(days);
    params.time = py.list({'00:00', '12:00'});
    % Request a 1-deg area around the target point. A single snapped point
    % risks landing on a land cell (ERA5 wave vars are NaN over land).
    % We download a small area and then select the nearest valid ocean point below.
    params.area = py.list({lat+0.5, lon-0.5, lat-0.5, lon+0.5});

    % Call Python function with filename
    era5.get_era5_data(params);

    % Find the latest .nc file in the directory (assuming Python script saves it)
    nc_files = dir('*.nc');  % Get list of all .nc files
    [~, newest_idx] = max([nc_files.datenum]);  % Find the latest file
    output_filename = nc_files(newest_idx).name;  % Get its name

    % Read grid coordinates from the netCDF
    nc_lat = double(ncread(output_filename, 'latitude'));   % [n_lat x 1]
    nc_lon = double(ncread(output_filename, 'longitude'));  % [n_lon x 1]

    % Read the downloaded data
    if strcmp(parameters{i}, 'significant_height_of_combined_wind_waves_and_swell')
        raw = ncread(output_filename, 'swh');   % [n_lon x n_lat x n_time]
        time_ERA5 = squeeze(ncread(output_filename, 'valid_time'));
    elseif strcmp(parameters{i}, 'peak_wave_period')
        raw = ncread(output_filename, 'pp1d');
    elseif strcmp(parameters{i}, 'mean_wave_direction')
        raw = ncread(output_filename, 'mwd');
    end

    % Find the nearest grid point that has valid (ocean) data.
    % Use the first timestep to identify which cells are ocean vs land (NaN).
    first_slice = raw(:,:,1);                          % [n_lon x n_lat]
    ocean_mask  = ~isnan(first_slice);                 % true = ocean cell

    % Build distance matrix (degrees) from every grid point to target lat/lon
    [LON_grid, LAT_grid] = meshgrid(nc_lon, nc_lat);   % both [n_lat x n_lon]
    dist = sqrt((LAT_grid - lat).^2 + (LON_grid - lon).^2);  % [n_lat x n_lon]

    % Mask out land cells and find the nearest ocean cell
    dist(~ocean_mask') = Inf;   % transpose because meshgrid is [lat x lon] but raw is [lon x lat]
    [~, min_idx] = min(dist(:));
    [lat_idx, lon_idx] = ind2sub(size(dist), min_idx);

    fprintf('  Nearest ocean point: lat=%.2f, lon=%.2f\n', nc_lat(lat_idx), nc_lon(lon_idx));

    % Extract the time series at that grid point
    ts = squeeze(raw(lon_idx, lat_idx, :));

    if strcmp(parameters{i}, 'significant_height_of_combined_wind_waves_and_swell')
        Hs_ERA5 = ts;
    elseif strcmp(parameters{i}, 'peak_wave_period')
        Tp_ERA5 = ts;
    elseif strcmp(parameters{i}, 'mean_wave_direction')
        D_ERA5 = ts;
    end

    % Delete the temporary .nc file after reading
    delete(output_filename);

end

disp('Done')

