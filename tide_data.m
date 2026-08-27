function [predi, obsv, d_t] = tide_data(begin_date, end_date, yrs, station)
% tide_data() Returns tide data from NOAA
%   See api info here:  https://tidesandcurrents.noaa.gov/api-helper/url-generator.html
%                       https://api.tidesandcurrents.noaa.gov/api/prod/responseHelp.html
%
%   Example api calls:
%       verified: https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?begin_date=20230701&end_date=20240629&station=8452660&product=hourly_height&datum=MTL&time_zone=lst_ldt&interval=h&units=metric&application=DataAPI_Sample&format=csv
%       prediction: https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?begin_date=20230701&end_date=20240629&station=8452660&product=predictions&datum=MTL&time_zone=lst_ldt&interval=h&units=metric&application=DataAPI_Sample&format=csv
% Format for begin_date and end_date is 'MMDD'
% Format for yrs is [yyyy:yyyy]'
% Format for station is '#########'

    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');

    options = weboptions;
    options.Timeout = 50; % Increase the default timeout from 5 to 20

    % Politeness delay [s] between successive API calls. NOAA CO-OPS
    % rate-limits their datagetter endpoint; back-to-back requests without a
    % small pause reliably trigger 403 Forbidden after a few dozen calls.
    api_pause_s = 0.3;

    filename = 'api_tidedata.csv';
    yrs = int2str(yrs);
    yrs = cellstr(yrs);

    % Retrieve verified data from API. NOAA CO-OPS recommends supplying a
    % descriptive application identifier so their throttling can distinguish
    % your traffic from the generic 'DataAPI_Sample' string (which is heavily
    % rate-limited).
    APP_ID = 'SIO_AirSea_TRIM';
    url_start = {'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?begin_date='};
    url_mid = {'&end_date='};
    url_mid2 = {'&station='};
    url_end = {['&product=hourly_height&datum=MHHW&time_zone=lst_ldt&interval=h&units=metric&application=' APP_ID '&format=csv']};
    % in url_end, adjust station code for each site and adjuist end_date value
    % to match the last recorded date for which a full 24 hrs of data is
    % available
    
    url_full = cell(length(yrs),1);
    V = table();
    for k = 1:length(yrs)
        url_Data = char(strcat(url_start, yrs{k}, begin_date, url_mid,yrs{k}, end_date, url_mid2, station, url_end));
        url = sprintf(url_Data);
        websave_retry(filename, url, options);
        pause(api_pause_s);
        data_table = readtable(filename);

        % Check for valid API response: must contain DateTime column
        if ~ismember('DateTime', data_table.Properties.VariableNames)
            warning('tide_data: no valid observed data for year %s — skipping.', yrs{k});
            continue
        end

        % Identify the water level column (first non-DateTime column) and
        % keep only DateTime + that column, renamed consistently.
        % This handles years where NOAA returns extra quality-flag columns.
        wl_col = data_table.Properties.VariableNames{2};
        data_table = data_table(:, {'DateTime', wl_col});
        data_table.Properties.VariableNames{2} = 'WaterLevel';

        V = [V; data_table];
        fprintf('Finished processing the recorded data for year %s\n',yrs{k})
    end



    % Retrieve predicted data from API
    url_start_p = {'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?begin_date='};
    url_mid_p = {'&end_date='};
    url_mid2_p = {'&station='};
    url_end_p = {['&product=predictions&datum=MHHW&time_zone=lst_ldt&interval=h&units=metric&application=' APP_ID '&format=csv']};
    url_full_p = cell(length(yrs),1);
    P = table();
    for j = 1:length(yrs)
        url_Data_p = char(strcat(url_start_p, yrs{j}, begin_date, url_mid_p, yrs{j}, end_date, url_mid2_p, station, url_end_p));
        url_p = sprintf(url_Data_p);
        websave_retry(filename, url_p, options);
        pause(api_pause_s);
        data_table_p = readtable(filename);

        % Check for valid API response: must contain DateTime column
        if ~ismember('DateTime', data_table_p.Properties.VariableNames)
            warning('tide_data: no valid predicted data for year %s — skipping.', yrs{j});
            continue
        end

        % Keep only DateTime + prediction column, renamed consistently.
        pred_col = data_table_p.Properties.VariableNames{2};
        data_table_p = data_table_p(:, {'DateTime', pred_col});
        data_table_p.Properties.VariableNames{2} = 'Prediction';

        P = [P; data_table_p];
        fprintf('Finished processing the predicted data for year %s\n',yrs{j})
    end
    
    mergedTable = outerjoin(V, P, 'Keys', 'DateTime', 'MergeKeys', true);
    
    obsv = mergedTable(:, 'WaterLevel');
    obsv = obsv{:,:}; % convert table array into vector
    
    predi = mergedTable(:,'Prediction');
    predi = predi{:,:}; % convert table array into vector
    
    dates = mergedTable(:,'DateTime');
    d_t = dates{:,:}; % convert table into a datetime table



    warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');

end


function websave_retry(filename, url, options)
% websave_retry  Wrapper around websave with exponential backoff for the
%                NOAA CO-OPS datagetter API. Retries on transient failures
%                (HTTP 403 rate-limit, 429 too-many-requests, 5xx server
%                errors, and MATLAB timeouts). Non-transient errors and
%                exhausted retries re-throw the original exception.

    max_attempts   = 5;
    base_delay_s   = 2;   % first retry waits 2 s; then 4 s, 8 s, 16 s, ...

    for attempt = 1:max_attempts
        try
            websave(filename, url, options);
            return   % success — leave the retry loop
        catch ME
            is_transient = contains(ME.message, '403') || ...
                           contains(ME.message, '429') || ...
                           contains(ME.message, '500') || ...
                           contains(ME.message, '502') || ...
                           contains(ME.message, '503') || ...
                           contains(ME.message, '504') || ...
                           contains(ME.message, 'timed out', 'IgnoreCase', true);
            if ~is_transient || attempt == max_attempts
                rethrow(ME);
            end
            delay = base_delay_s * 2^(attempt - 1);
            fprintf('  websave failed (attempt %d/%d, %s). Backing off %d s...\n', ...
                    attempt, max_attempts, strtrim(regexprep(ME.message, '\s+', ' ')), delay);
            pause(delay);
        end
    end
end