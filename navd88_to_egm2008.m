function [zOut, info] = navd88_to_egm2008(lon, lat, zNAVD88, varargin)
%NAVD88_TO_EGM2008 Convert NAVD88 orthometric heights to EGM2008(WGS84)
% ellipsoidal heights using NOAA's VDatum REST API. Replaces manually
% running points through the VDatum web GUI (vdatum.noaa.gov/vdatumweb).
%
%   zOut = NAVD88_TO_EGM2008(lon, lat, zNAVD88)
%   zOut = NAVD88_TO_EGM2008(lon, lat, zNAVD88, 'region', 'westcoast')
%   [zOut, info] = NAVD88_TO_EGM2008(...)
%
% REQUIRED INPUTS
%   lon, lat   : decimal degrees (West negative). Scalars or vectors of
%                the same size.
%   zNAVD88    : NAVD88 height(s), same size as lon/lat, or a single
%                scalar to apply to every point.
%
% OPTIONAL NAME-VALUE PAIRS
%   'region'   : VDatum region code, e.g. 'westcoast', 'contiguous',
%                'gcnmi', 'hi'. If you don't pass this, the function
%                prints VDatum's region list and asks you to pick one
%                interactively -- the wrong region for your point's
%                geography is a common cause of VDatum rejecting or
%                erroring on a request, so this is intentionally not
%                silently defaulted anymore.
%   'unit'     : 'm' or 'ft' (applied to both s_v_unit and t_v_unit).
%                Default 'm'.
%   's_v_geoid': source geoid model used to step NAVD88 -> ellipsoid.
%                Default 'geoid18' (current NGS hybrid geoid, matches
%                the VDatum GUI default).
%   'pause_s'  : seconds to pause between requests, to be polite to the
%                NOAA server on long lists of points. Default 0.1.
%
% OUTPUTS
%   zOut  : EGM2008(WGS84) height for each point, same units as input.
%           NaN for any point where the request failed -- a warning is
%           issued in that case, including VDatum's own error message
%           when one is available.
%   info  : struct array, one element per point, with fields
%           lon, lat, s_z (input height as echoed by VDatum),
%           t_z (output height), uncertainty (VDatum's +/- estimate, m
%           or ft to match 'unit'), region, h_frame (horizontal frame
%           the output lon/lat/height are referenced to -- always
%           WGS84_G1674, since that's what EGM2008 is defined against).
%
% EXAMPLE (single tide gauge, Imperial Beach area, prompts for region)
%   [z, info] = navd88_to_egm2008(-117.1611, 32.5763, 1.832);
%
% EXAMPLE (batch of stations, Guam region, region given up front)
%   lon = [144.6, 144.65, 144.7];
%   lat = [13.45, 13.47, 13.50];
%   zN  = [2.10, 1.95, 2.30];
%   [z, info] = navd88_to_egm2008(lon, lat, zN, 'region', 'gcnmi');
%
% Requires internet access to vdatum.noaa.gov. See API docs at
% https://vdatum.noaa.gov/docs/services.html for the full parameter set.

p = inputParser;
addParameter(p, 'region', '');
addParameter(p, 'unit', 'm');
addParameter(p, 's_v_geoid', 'geoid18');
addParameter(p, 'pause_s', 0.1);
parse(p, varargin{:});
region    = p.Results.region;
unit      = p.Results.unit;
sGeoid    = p.Results.s_v_geoid;
pauseSecs = p.Results.pause_s;

if isempty(region)
    region = promptForRegion();
end

lon = lon(:);
lat = lat(:);
if numel(zNAVD88) == 1
    zNAVD88 = repmat(zNAVD88, size(lon));
else
    zNAVD88 = zNAVD88(:);
end

if ~isequal(numel(lon), numel(lat), numel(zNAVD88))
    error('navd88_to_egm2008:sizeMismatch', ...
        'lon, lat, and zNAVD88 must be scalars or vectors of the same length.');
end

n = numel(lon);
zOut = nan(n, 1);
info = repmat(struct('lon', NaN, 'lat', NaN, 's_z', NaN, 't_z', NaN, ...
    'uncertainty', NaN, 'region', region, 'h_frame', 'WGS84_G1674'), n, 1);

baseURL = 'https://vdatum.noaa.gov/vdatumweb/api/convert';
opts = weboptions('Timeout', 20);

for i = 1:n
    try
        resp = webread(baseURL, ...
            's_x', num2str(lon(i), '%.8f'), ...
            's_y', num2str(lat(i), '%.8f'), ...
            's_z', num2str(zNAVD88(i), '%.4f'), ...
            'region', region, ...
            's_v_frame', 'NAVD88', ...
            's_v_unit', unit, ...
            's_v_geoid', sGeoid, ...
            't_h_frame', 'WGS84_G1674', ...
            't_v_frame', 'EGM2008', ...
            't_v_unit', unit, ...
            't_v_geoid', 'egm2008', ...
            opts);
    catch ME
        warning('navd88_to_egm2008:requestFailed', ...
            'Point %d (lon=%.5f, lat=%.5f) failed: %s', i, lon(i), lat(i), ME.message);
        continue
    end

    % VDatum returns HTTP 200 even when it can't perform the conversion --
    % in that case the body is an error object (errorCode/message) instead
    % of a result object, and there is no t_z field. Catch that here
    % rather than crashing on an unrecognized field name.
    if ~isfield(resp, 't_z')
        if isfield(resp, 'message')
            warning('navd88_to_egm2008:apiError', ...
                'Point %d (lon=%.5f, lat=%.5f) was rejected by VDatum: %s', ...
                i, lon(i), lat(i), resp.message);
        else
            warning('navd88_to_egm2008:apiError', ...
                'Point %d (lon=%.5f, lat=%.5f) returned an unexpected response with no t_z field.', ...
                i, lon(i), lat(i));
        end
        continue
    end

    tz = str2double(resp.t_z);
    if isnan(tz) || tz <= -999999
        warning('navd88_to_egm2008:noData', ...
            'Point %d (lon=%.5f, lat=%.5f) returned no-data (-999999). Likely outside the %s VDatum model bounds.', ...
            i, lon(i), lat(i), region);
        continue
    end

    zOut(i) = tz;
    info(i).lon = lon(i);
    info(i).lat = lat(i);
    info(i).s_z = str2double(resp.s_z);
    info(i).t_z = tz;
    info(i).uncertainty = str2double(resp.uncertainty);
    info(i).region = resp.region;
    info(i).h_frame = resp.t_h_frame;

    if pauseSecs > 0 && i < n
        pause(pauseSecs);
    end
end

end


function code = promptForRegion()
% Print VDatum's region list (mirrors the dropdown on vdatum.noaa.gov)
% and ask the user to pick one by number. Returns the API region code.

names = { ...
    'Alaska', ...
    'South East Alaska Tidal', ...
    'American Samoa', ...
    'Contiguous United States', ...
    'Chesapeake/Delaware Bay', ...
    'West Coast', ...
    'West Gulf Coast', ...
    'Guam and Commonwealth of Northern Mariana Islands', ...
    'Hawaii', ...
    'Puerto Rico and US Virgin Islands', ...
    'Saint George Island', ...
    'Saint Paul Island', ...
    'Saint Lawrence Island'};

codes = {'ak', 'seak', 'as', 'contiguous', 'chesapeak_delaware', ...
    'westcoast', 'wgom', 'gcnmi', 'hi', 'prvi', 'sgi', 'spi', 'sli'};

fprintf('\nSelect the VDatum region for your point(s):\n');
for k = 1:numel(names)
    fprintf('  %2d) %s\n', k, names{k});
end

sel = NaN;
while isnan(sel) || sel < 1 || sel > numel(names) || mod(sel, 1) ~= 0
    raw = input(sprintf('Enter a number (1-%d): ', numel(names)), 's');
    sel = str2double(raw);
end

code = codes{sel};
fprintf('Using region: %s (%s)\n\n', names{sel}, code);

end