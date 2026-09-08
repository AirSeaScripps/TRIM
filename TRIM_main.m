close all
clear all

% This script uses water level data from https://tidesandcurrents.noaa.gov/stations.html?type=Water+Levels
% Set units to "Metric," timezone to "LST/LDT," and datum to "NAVD88".
%
% Total Water Level (TWL) is estimated using a copula-based joint
% probability model:
%
%   CALC_WAVES = 1  (preferred):
%     A trivariate t-copula is fitted to concurrent (tide, surge, Hs)
%     observations.  Tp is estimated from a power-law regression on storm
%     peaks.  TWL is computed directly inside the Monte Carlo as:
%       TWL_i = tide_i + surge_i + Stockdon2006(Hs_i, Tp_i, beta)
%     Annual maxima of TWL are extracted and fitted with a GEV to give the
%     100-year return level.  This approach correctly captures storm-surge–
%     wave co-occurrence and avoids combining separately-derived return
%     levels (which assumes independence and overcounts joint probability).
%
%   CALC_WAVES = 0  (fallback):
%     A bivariate t-copula is fitted to (tide, surge) pairs.  Wave runup
%     is set to zero.  Useful for sites with no offshore wave data or
%     where the coast is sheltered from significant swell.
%
% Supported copula family: 't' (Student-t)
%   t-copula: symmetric upper and lower tail dependence, parameterised by
%             a correlation matrix (rho) and degrees of freedom (nu).
%             nu → inf recovers the Gaussian (no tail dependence) limit.
%   For CALC_WAVES = 0, 'Gumbel' is also accepted (upper-tail only).
%
% Requires: Statistics and Machine Learning Toolbox

% Initializations
infile = []; % Optional: set to the path of a previously-saved workspace
             % (produced by a prior run of this script) to skip
             % re-downloading NOAA/ERA5 data on repeat runs during
             % development. Leave empty ([]) to always download fresh data.

%% USER INPUT

infile_passthrough = infile; % Duplicate to pass through the load step (in case infile = [] in the loaded file)

USE_PRODUCT_COPULA = false;   % set true for independence control run (default = false)
rng(42);                       % fix seed for reproducible comparison

% TIDE SECTION: Use (for example) NOAA tide data (with observed and predicted water levels)
% Choose the best tide station here: https://tidesandcurrents.noaa.gov/map/index.html
USE_NOAA = 1; % yes = 1, no = 0
years = [1975:2025]'; % Enter the timespan to pull. Should be at least 30 years (otherwise in testing mode)
station = '9410230'; % Enter the station number (for example 
ddLat = ddm2dd("32° 52.0 N"); % Enter the coordinate string for latitude published for the tide gauge (eg: "48° 6.7 N")
ddLon = ddm2dd("117° 15.4 W"); % Enter the coordinate string for longitude published for the tide gauge (eg: "122° 45.6 W")

% Enter Datum Information (Be sure to select NAVD88 datum with [m] units from "Datum" section of NOAA tides and currents for this location).
MHHW = 	1.566;   % Enter MHHW value IN METERS above NAVD88
NAVD88 = 0;     % Enter NAVD88 IN METERS (zero if tide data is referenced to NAVD88). This is available if NAVD88 is not an available datum and a further adjustment must be made to MTL.

% Enter the transformation from NAVD88 to EGM2008
[NAVD88_to_EGM, info] = navd88_to_egm2008(ddLon, ddLat, NAVD88); % Possible to verify using the conversion tool: https://vdatum.noaa.gov/vdatumweb/

% OPTION 2: Use 30 years of wave data? If so, Enter variables required for wave analysis
% Note: Get_ERA5_Waves.m always requests the most recent 30 full calendar
% years, computed dynamically from the current date (see that file).
% Wave data from ERA: https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels
CALC_WAVES = 1; % Calculate wave runup using trivariate copula [1], or skip [0]

% --- Shoreline wave model selection ---
% Choose which model translates offshore (Hs, Tp) into the shoreline wave-
% driven contribution to TWL. Used only when CALC_WAVES = 1.
%
%   'Stockdon' (default) : Stockdon et al. (2006) for sandy beaches.
%                          Uses foreshore slope β. Applies dissipative and
%                          reflective branches based on the Iribarren number.
%   'Becker'             : Becker & Merrifield (2026) for dissipative fringing
%                          reefs similar in character to Ipan reef, Guam
%                          (δ > 1 under typical conditions). Uses the mean
%                          reef-flat bed elevation D0. Not valid for wave-like
%                          reefs (e.g. Majuro, Roi-Namur RMI) or non-reef
%                          coastlines.
WAVE_MODEL = 'Stockdon';
assert(any(strcmpi(WAVE_MODEL, {'Stockdon', 'Becker'})), ...
       'WAVE_MODEL must be ''Stockdon'' or ''Becker'', got ''%s''.', WAVE_MODEL);

% --- Stockdon (sandy beach) parameter ---
beta = 0.06;    % Foreshore beach slope [dimensionless], used when WAVE_MODEL = 'Stockdon'.

% --- Becker (dissipative reef) parameter ---
% Mean reef-flat BED ELEVATION relative to the tide gauge datum [m].
% Negative below datum (typical for a submerged reef). Extract by averaging
% a topobathymetric DEM over a window on the reef flat fronting the site.
% D0 must be in the same vertical datum as the tide gauge (convert via VDatum
% if the DEM is in ellipsoidal / MHW / a different frame).
%
% Example (Hagåtña Bay, Guam, 144.764°E 13.481°N, from 2020 NOAA NGS Topobathy
% Lidar, 100 m window on the reef flat): D0 = -0.63 m (std 0.19 m).
D0 = -0.73;     % [m], used when WAVE_MODEL = 'Becker'.

% --- Offshore wave grid point ---
lat = 52.25;    % Latitude of offshore wave calculation
lon = -176.50;  % Longitude of offshore wave calculation

% Percentile threshold for POT wave analysis [%].
% Only ERA5 records above this Hs level are used to fit the GPD tail model.
% Higher values → fewer but more extreme storm peaks, noisier GPD fit.
% Lower values  → more data but includes non-extreme events that bias ξ.
% 95 is a widely used default; 90–98 is a reasonable range.
Hs_thresh_pct = 95;

% Physical upper bound on Hs sampled in the Monte Carlo [m].
% When the fitted GPD has a positive shape parameter (ξ > 0), which is common
% for typhoon-driven sites, the distribution is unbounded above and rare deep-
% tail draws can produce Hs values well beyond anything physically possible
% at the site. This is particularly problematic for the Becker runup
% parameterization, whose setup coefficient scales linearly with Hs and can
% amplify unphysical Hs into unphysical TWL, breaking the downstream GEV fit.
% A physical cap keeps the marginal well-behaved without altering the copula
% fit or the GPD parameters themselves. Reasonable defaults:
%   - open-coast extratropical (e.g. Imperial Beach)  : Hs_max ~ 12 m
%   - Pacific typhoon-envelope (e.g. Guam, Saipan)    : Hs_max ~ 20 m
%   - hurricane-envelope Atlantic/Gulf                : Hs_max ~ 20 m
% Set to Inf to disable the cap.
Hs_max = 20;   % [m]


% --- Copula options ---
% Family: 't' (Student-t, default).  'Gumbel' is also an option when not
% using waves (and only a bivariate fit between tide and storm surge is
% needed.)
COPULA_FAMILY = 't';

% Number of Monte Carlo repetitions (each simulates 100 years)
repeats = 300;

% Surge decorrelation time [hours] — used ONLY for wave POT declustering.
% Storm surge is temporally autocorrelated — a single storm can sustain
% elevated water levels for many consecutive hours, so successive hourly
% observations are NOT independent.
%
% Previously, tau_hours was also used to define n_eff_per_year for the
% Monte Carlo (n_eff = floor(8760/tau_hours)).  That approach is sensitive
% to tau_hours: varying it from 12 to 72 h shifts the 100-yr TWL by ~0.6 m,
% making the result dependent on an arbitrary tuning parameter.
%
% The current approach instead derives n_eff_per_year directly from the
% number of POT-identified storm peaks (lambda_waves × 100 yr), exactly as
% the wave extreme-value analysis does.  tau_hours is retained only to set
% the minimum inter-storm separation during ERA5 declustering.
tau_hours = 24;

%% Extract data directly from API call to NOAA and ERA5
% MHHW_el is the height of the MHHW datum surface above NAVD88.
% It is always recomputed here from the user-specified MHHW and NAVD88
% values, regardless of whether data are loaded from a saved workspace or
% downloaded live.  This prevents a stale MHHW_el from a prior session
% being silently inherited when running in infile (debug) mode.
MHHW_el = MHHW - NAVD88;

if ~isempty(infile)
    load(infile, "-regexp", "^(?!(beta|repeats|tau_hours|Hs_thresh_pct|COPULA_FAMILY)$).+"); % Load everything except user defined tuning parameters defined above
    infile = infile_passthrough; % Reset the infile
    clearvars -except infile predi obsv d_t Hs_ERA5 Tp_ERA5 time_ERA5 D_ERA5 lat lon years station beta Hs_thresh_pct...
        MHHW NAVD88 NAVD88_to_EGM COPULA_FAMILY repeats tau_hours obsv_v predi_v USE_NOAA MHHW_el ind_nans mask CALC_WAVES...
        beta repeats tau_hours Hs_thresh_pct
    % Re-apply MHHW_el after clearvars, since it may have been overwritten
    % by the loaded workspace before clearvars ran.
    MHHW_el = MHHW - NAVD88;
    disp('Using previously downloaded data from infile (skipping NOAA/ERA5 download).')
else

    if USE_NOAA == 1
        disp('Using NOAA tides and water levels.')

        % Get the data
        disp('Getting tide and water level data...')
        [predi, obsv, d_t] = tide_data('0101','1231',years,station); % Return predicted and observed water level WRT MHHW

        % Clean the data
        d_t(end-24:end) = []; %Remove the last day (observations likely to be partially processed or incomplete)
        obsv(end-24:end) = [];
        predi(end-24:end) = [];

    elseif USE_NOAA == 0
        disp('Option to use NOAA tide and water level is turned off.')
    else
        disp('Please set USE_NOAA to be either 1 (yes) or 0 (no).')
    end

    if CALC_WAVES == 1
        % Call the function to get ERA5 data
        fprintf('Getting wave data for [%0.2f,%0.2f]\n',lat,lon)
        [Hs_ERA5,Tp_ERA5,time_ERA5,D_ERA5] = Get_ERA5_Waves(lat,lon);
        disp('NOTE: Wave data limited to recent 30 year span.')
    elseif CALC_WAVES == 0
        disp('Skipping wave data. Not using a runup model.')
    else
        disp('Please set calc_waves to be either 1 (use a runup model) or 0 (skip wave runup).')
    end

end

% Remove duplicate timestamps in d_t. NOAA CO-OPS returns data in local
% clock time (time_zone=lst_ldt), which produces a repeated hour at the
% autumn DST fall-back transition each year. Runs at DST-observing
% stations with long enough records (e.g. Alaska sites going back to the
% 1970s) can hit this, which breaks interp1 later since it requires
% unique sample points. Keep the first occurrence of each duplicate.
% Applied unconditionally (both the infile/cached-load path and the fresh
% NOAA-download path) so old cached .mat files that predate this fix are
% cleaned up on load too.
[d_t, keep_idx] = unique(d_t, 'stable');
if numel(keep_idx) < numel(predi)
    fprintf('Removed %d duplicate timestamp(s) from d_t (kept first occurrence).\n', ...
        numel(predi) - numel(keep_idx));
end
predi = predi(keep_idx);
obsv  = obsv(keep_idx);

%% Tide Analysis section

if ~isempty(infile)
    disp('Skipping tide data cleaning... this data has already been cleaned.')
else

    % Remove nan values (for visualization)
    obsv_v = obsv;                      % duplicate "observational values" vector (will have no nans)
    predi_v = predi;                    % duplicate "predicted values" vector (will have no nans)
    indo_nans = find(isnan(obsv_v));
    indp_nans = find(isnan(predi_v));
    ind_nans = [indo_nans;indp_nans];   % combined nans vector
    obsv_v(ind_nans) = [];              % remove missing observational values
    predi_v(ind_nans) = [];             % remove missing observational values

    % Make a mask for non-nan values (for visualization)
    mask = true(size(d_t));
    mask(ind_nans) = false;

    % Replace nan values with the mean for time series analysis
    % (this is fine as long as there aren't too many - the high pass error will be filtered out)
    nanpredi = sum(isnan(predi));
    predi(isnan(predi)) = nanmean(predi);
    fprintf('I replaced %d NaN values in the predicted record\n',nanpredi)

    nanobsv = sum(isnan(obsv));
    obsv(isnan(obsv)) = nanmean(obsv);
    fprintf('I replaced %d NaN values in the observed record\n',nanobsv)
end

TideExcess = obsv - predi; % This is the difference from predicted sea level

%% Wave analysis section
R2 = 0; % Initialize with zero waves (assume no waves unless otherwise told to)

if CALC_WAVES == 1
    % --- POT + GPD extreme wave analysis ---
    %
    % Fitting a GEV (or Rayleigh) to all hourly Hs is inappropriate: the GEV
    % is a limit theorem for block maxima, not raw observations. With only
    % 30 years of ERA5 data, annual block maxima yields just 30 points —
    % too few for reliable 3-parameter GEV estimation.
    %
    % Instead we use Peaks Over Threshold (POT) with a Generalized Pareto
    % Distribution (GPD). The Pickands-Balkema-de Haan theorem guarantees
    % GPD is the correct limiting distribution for threshold exceedances.
    %
    % Steps:
    %   1. Set threshold u at the Hs_thresh_pct percentile.
    %   2. Decluster: group consecutive exceedances within tau_hours into
    %      one storm event (same window as the surge independence assumption),
    %      keeping only the peak Hs of each cluster.
    %   3. Fit GPD to the storm-peak excesses (storm_peak_Hs - u).
    %   4. Compute storm rate lambda_waves [storms/yr] from the record length.
    %   5. Invert the Poisson-GPD return level equation for T = 100 yr:
    %        Hs100 = u + GPD^{-1}( 1 - 1/(lambda_waves * T) )
    %
    % Note: Hs100 is retained as a DIAGNOSTIC output only.  The 100-year
    % TWL return level is computed directly inside the trivariate Monte
    % Carlo (Step F), which avoids combining separately-derived return
    % levels and is the more physically correct approach.

    % --- Step 1: Threshold ---
    u_waves = prctile(Hs_ERA5, Hs_thresh_pct);
    fprintf('  Wave threshold u = %.2f m (%.0f%% percentile of all Hs)\n', u_waves, Hs_thresh_pct);

    % --- Step 2: Declustering ---
    % Derive ERA5 time step directly from the valid_time vector (Unix seconds).
    dt_ERA5_h = double(time_ERA5(2) - time_ERA5(1)) / 3600;   % hours between ERA5 records (valid_time is Unix seconds; cast to double to prevent int64 propagation)
    assert(dt_ERA5_h >= 1 && dt_ERA5_h <= 24, ...
           'Unexpected ERA5 time step (%.2f h). Check that valid_time is in Unix seconds.', dt_ERA5_h);
    min_sep_steps = ceil(tau_hours / dt_ERA5_h);          % timesteps required between independent storms

    exceedance_idx = find(Hs_ERA5 > u_waves);
    if isempty(exceedance_idx)
        error('No Hs exceedances above %.0f%% threshold (u = %.2f m). Lower Hs_thresh_pct.', ...
              Hs_thresh_pct, u_waves);
    end

    % A new cluster begins wherever the gap between consecutive exceedance
    % indices exceeds min_sep_steps (i.e., Hs dropped below u for long enough).
    gaps            = diff(exceedance_idx);
    cluster_starts  = [1; find(gaps > min_sep_steps) + 1];  % positions within exceedance_idx
    n_storms        = length(cluster_starts);

    storm_peak_Hs = zeros(n_storms, 1);
    storm_peak_Tp = zeros(n_storms, 1);

    for k = 1:n_storms
        if k < n_storms
            cluster_idx = exceedance_idx(cluster_starts(k) : cluster_starts(k+1) - 1);
        else
            cluster_idx = exceedance_idx(cluster_starts(k) : end);
        end
        [storm_peak_Hs(k), pk] = max(Hs_ERA5(cluster_idx));
        storm_peak_Tp(k)       = Tp_ERA5(cluster_idx(pk));
    end

    % --- Step 3: GPD fit to storm-peak excesses ---
    excesses = storm_peak_Hs - u_waves;
    pd_GPD   = fitdist(excesses, 'GeneralizedPareto', 'theta', 0);
    fprintf('  GPD fit:  sigma = %.3f m,  xi = %.3f\n', pd_GPD.sigma, pd_GPD.k);

    % --- Step 4: Storm rate ---
    n_years_waves  = length(Hs_ERA5) * dt_ERA5_h / (365.25 * 24);
    lambda_waves   = n_storms / n_years_waves;
    fprintf('  %d independent storm peaks over %.1f years  (lambda = %.2f storms/yr)\n', ...
            n_storms, n_years_waves, lambda_waves);

    % --- Step 5: 100-year return level (diagnostic only) ---
    % From the Poisson-GPD model, the annual exceedance probability is:
    %   P(Hs_annual_max > x) ≈ 1 - exp(-lambda * P_GPD(excess > x - u))
    % For large T (T=100), this simplifies to the quantile inversion:
    %   Hs_T = u + GPD^{-1}( 1 - 1/(lambda * T) )
    % Guard: the Poisson-GPD return level formula requires at least one
    % expected exceedance in the return period (lambda * T >= 1). If this
    % fails, the threshold is too high — lower Hs_thresh_pct so more storms
    % are included and lambda_waves increases.
    if lambda_waves * 100 < 1
        error(['Too few storms to estimate a 100-yr return level: ' ...
               'lambda_waves = %.3f storms/yr gives only %.2f expected exceedances ' ...
               'in 100 years.\nLower Hs_thresh_pct (currently %.0f%%) to include ' ...
               'more storm peaks and increase lambda_waves.'], ...
               lambda_waves, lambda_waves * 100, Hs_thresh_pct);
    end
    p_GPD = 1 - 1 / (lambda_waves * 100);
    Hs100 = u_waves + icdf(pd_GPD, p_GPD);
    fprintf('  Hs100 (POT+GPD, diagnostic) = %.2f m\n', Hs100);

    % Tp100: mean Tp of the most extreme declustered storm peaks (top 10%).
    % Tp is not formally modelled — using the conditional mean of the upper
    % tail of the observed storm-peak Tp is standard practice.
    Tp100_thresh = prctile(storm_peak_Hs, 90);
    Tp100 = mean(storm_peak_Tp(storm_peak_Hs > Tp100_thresh));
    fprintf('  Tp100 (conditional mean, diagnostic) = %.1f s\n', Tp100);

    % --- Plot 1: GPD fit to storm-peak excesses ---
    x_exc  = linspace(0, max(excesses), 200);
    y_gpd  = pdf(pd_GPD, x_exc);
    figure;
    hold on
    histogram(excesses, 30, 'Normalization', 'pdf', 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none')
    plot(x_exc, y_gpd, 'r-', 'LineWidth', 2.5)
    xline(Hs100 - u_waves, 'k--', 'LineWidth', 1.5, ...
          'Label', sprintf('H_{s100}−u = %.2f m', Hs100 - u_waves), ...
          'LabelHorizontalAlignment', 'left','fontsize',12,'fontname','times')
    xlabel('$H_s - u \quad [\mathrm{m}]$ (excess above threshold)', 'interpreter', 'latex', 'fontsize', 18)
    ylabel('Probability Density', 'fontsize', 18)
    title(sprintf('POT storm peaks — GPD fit  (u = %.2f m, %d storms, %.1f yr)', ...
          u_waves, n_storms, n_years_waves), 'fontsize', 14)
    legend('Declustered storm peaks', sprintf('GPD  (\\sigma=%.2f, \\xi=%.2f)', pd_GPD.sigma, pd_GPD.k), ...
           'Location', 'northeast', 'fontsize', 18)
    set(gca,'fontsize',18)
    box on; grid on;
    hold off

    % --- Plot 2: Hs-Tp scatter (storm peaks only, Hs100 marked) ---
    figure
    scatter(storm_peak_Hs, storm_peak_Tp, 20, 'filled', 'k')
    hold on
    scatter(Hs100, Tp100, 120, 'filled', 'r')
    grid on
    xlabel('$H_s \quad [\mathrm{m}]$', 'interpreter', 'latex', 'fontsize', 14)
    ylabel('$T_p \quad [\mathrm{s}]$', 'interpreter', 'latex', 'fontsize', 14)
    title(sprintf('Storm peaks at [%0.2f,%0.2f]', lat, lon), 'fontsize', 14, 'fontname', 'times')
    legend('Storm peak', '$H_{s100},\,T_{p100}$ (diagnostic)', 'interpreter', 'latex', 'fontsize', 12)
    hold off

elseif CALC_WAVES == 0
    disp('Waves were not calculated. R2 = 0 m.')
else
    error('calc_waves was not set properly. Please choose either 1 (calculate wave runup) or 0 (no waves are at this site).')
end


%% Data exploration and plots
% Total Water Level = Mean Tide Level + Sea Level Rise + Tide + Storm Surge and setup + Waves + Datum adjustment

% MTL = 0 (zero if data is extracted with reference to MTL)
% SLR = Low pass filtered TideExcess
% Tide = Predicted tide (from NOAA or harmonic constituents)
% Storm Surge and setup = TideExcess - SLR
% Waves = Extra component (not from this data set)

% Low pass filter the tide data to find sea level rise component
Fs = 1; % Sampling frequency (samples per hour)
Fc = 1/(24*365*12); % Cutoff frequency (cycles per hour)
% Normalize the cutoff frequency
Wn = Fc / (Fs / 2); % Normalized cutoff frequency (Nyquist frequency is Fs/2)

% Design a 4th-order Butterworth filter using Second Order Sections (to allow for a long Fc)
[z, p, k] = butter(4, Wn, 'low');
[sos, g] = zp2sos(z, p, k);
LPC = filtfilt(sos, g, TideExcess);

SS = TideExcess - LPC; % Storm surge and setup
yrTWL = LPC + predi + SS; % Reconstructed Total Water Level

SS(ind_nans) = []; % take out times where TideExcess calculation has no observation or predicted values
yrTWL(ind_nans) = [];

% Time series of raw data
figure;
hold on
plot(d_t,predi,'b')
plot(d_t(mask), yrTWL,'r')
grid on
legend('Predicted (tidal)','Observed')
ylabel('Depth [m]')

% Functional relationship between observations and predictions
figure;
scatter(obsv_v,predi_v,'filled')
grid on

% Distribution of tidal elevation
pd_tide_kde = fitdist(predi_v,'Kernel');
x = linspace(min(predi_v),max(predi_v),100);
y = pdf(pd_tide_kde,x);
figure;
hold on
plot(x,y,'LineWidth',3)
histogram(predi_v,100,'Normalization','pdf')
grid on
xlabel('$\eta_{\mathrm{tide}} \quad [\mathrm{m}]$','interpreter','latex')
ylabel('Probability Density')
hold off
p = 0.99;
valueAt99Percent = icdf(pd_tide_kde, p);

% Distribution of surge + wave setup
pd1 = fitdist(SS,'Kernel');
x1 = linspace(min(SS),max(SS),100);
y1 = pdf(pd1,x1);
figure;
hold on
plot(x1,y1,'LineWidth',3)
histogram(SS,100,'Normalization','pdf')
grid on
xlabel('$\eta_{\mathrm{local}} \quad [\mathrm{m}]$','interpreter','latex')
ylabel('Probability Density')
hold off
valueAt99Percent1 = icdf(pd1, p);

% Time series of \Delta \eta
TideExcess2 = TideExcess;
TideExcess2(ind_nans) = NaN; % Don't plot times where there were no observations
figure;
hold on
plot(d_t, TideExcess2,'b')
plot(d_t,LPC,'r','linewidth',2)
ylabel('Depth [m]')
legend('$\eta_{\mathrm{res}}$','$\eta_{\mathrm{SL}}$','interpreter','latex','fontsize',16)
box on
grid on

%% Estimate SLR component

x = 1:length(d_t)';
P = polyfit(x,LPC,1);
yfit = polyval(P,x);

SLR = yfit(end); % This is the latest estimate of sea level rise (reference to the most recent datum)

% Option to predict what SLR will be in the future
yr_predict = 20; % Set the number of years out to make a prediction (refine this possibly with IPCC data?)

dt = hours(d_t(2) - d_t(1))/(365.25 *24); % Time step in years
ntsteps = yr_predict./dt;
SLR_predict = P(1)*ntsteps + yfit(end); % slope * years + current SLR (with reference to last datum)


%% Copula-based 100-year return level
%
% When CALC_WAVES = 1 (preferred):
%   A trivariate (tide, surge, Hs) t-copula captures the joint dependence
%   of all three components.  TWL is computed directly as
%     TWL = tide + surge + Stockdon(Hs, Tp, beta)
%   inside each Monte Carlo draw.  Annual maxima are fitted with a GEV to
%   give the 100-year TWL.  value_99 therefore already includes wave runup.
%
% When CALC_WAVES = 0 (fallback):
%   A bivariate (tide, surge) copula is used.  Wave runup = 0 (R2 = 0).
%   value_99 is the 100-year tide + surge return level.
%
% In either case value_99 is the 100-year combined return level relative
% to MHHW, and the EGM2008 result adds SLR + datum offsets on top.

if USE_PRODUCT_COPULA
    fprintf(['CAUTION!! Running independent copula for test purposes...\n']);
else
    fprintf(['Running the fitted copula...\n']);
end


if CALC_WAVES == 1

    %% ── STEP A: DATA ALIGNMENT ───────────────────────────────────────────────
    % Subsample NOAA hourly (tide, surge) to the 12-hourly ERA5 timestamps.
    % 12-hour sampling covers all tidal phases over ~15-day beat cycles and
    % gives ~21 000 concurrent observations over 30 years — ample for copula
    % fitting.  The copula captures rank dependence at storm time scales
    % (days), which is well resolved at 12-hour resolution.

    % Reconstruct full-length surge aligned with d_t (before ind_nans removal)
    SS_full = TideExcess - LPC;   % [n_full × 1]

    % Convert time vectors to Unix seconds for matching
    t_noaa_unix = posixtime(d_t);          % [n_full × 1]
    t_era5_unix = double(time_ERA5);       % [n_ERA5 × 1]  (cast avoids int64 issues)

    % Nearest NOAA index for each ERA5 timestamp (O(n_ERA5 log n_noaa))
    noaa_idx_all = round(interp1(t_noaa_unix, (1:numel(d_t))', ...
                                  t_era5_unix, 'nearest', 'extrap'));
    noaa_idx_all = max(1, min(numel(d_t), noaa_idx_all));   % clamp to valid range

    % Reject matches more than 1 hour apart (guards against timezone offsets)
    time_offset  = abs(t_noaa_unix(noaa_idx_all) - t_era5_unix);
    era5_keep    = time_offset < 3600;                        % logical, length n_ERA5
    noaa_idx_conc = noaa_idx_all(era5_keep);
    era5_idx_conc = find(era5_keep);

    % Extract concurrent triples
    tide_conc  = predi(noaa_idx_conc);
    surge_conc = SS_full(noaa_idx_conc);
    Hs_conc    = double(Hs_ERA5(era5_idx_conc));
    Tp_conc    = double(Tp_ERA5(era5_idx_conc));

    % Remove rows flagged as bad in the NOAA record (ind_nans) or with NaN
    noaa_is_bad           = false(numel(d_t), 1);
    noaa_is_bad(ind_nans) = true;
    row_bad = noaa_is_bad(noaa_idx_conc) | isnan(tide_conc) | ...
              isnan(surge_conc) | isnan(Hs_conc);

    tide_conc  = tide_conc(~row_bad);
    surge_conc = surge_conc(~row_bad);
    Hs_conc    = Hs_conc(~row_bad);
    Tp_conc    = Tp_conc(~row_bad);
    n_conc     = numel(tide_conc);

    fprintf('Alignment: %d concurrent (tide, surge, Hs) observations (%.1f per year)\n', ...
            n_conc, n_conc / n_years_waves);


    %% ── STEP B: TRIVARIATE PSEUDO-OBSERVATIONS ───────────────────────────────
    % Rank-based transform to uniform marginals [0,1] for all three variables.
    % tiedrank() / (n+1) is O(N log N) and is the standard Weibull plotting
    % position used in copula literature.

    U_tide  = tiedrank(tide_conc)  / (n_conc + 1);
    U_surge = tiedrank(surge_conc) / (n_conc + 1);
    U_Hs    = tiedrank(Hs_conc)   / (n_conc + 1);

    U3 = [U_tide, U_surge, U_Hs];   % [n_conc × 3]

    % Report Spearman correlations as a pre-fit sanity check
    rho_sp = corr(U3, 'type', 'Pearson');   % Pearson on ranks = Spearman
    fprintf('Spearman rank correlations (pre-fit):\n');
    fprintf('  tide–surge : %+.4f\n', rho_sp(1,2));
    fprintf('  tide–Hs    : %+.4f\n', rho_sp(1,3));
    fprintf('  surge–Hs   : %+.4f\n', rho_sp(2,3));

    % Physical guard: surge and Hs are both storm-driven → expect rho_sH >= 0
    if rho_sp(2,3) < 0
        warning(['Negative surge–Hs Spearman correlation (rho = %.4f).\n' ...
                 'Both are storm-driven; negative dependence is unexpected.\n' ...
                 'Check data alignment and site conditions before proceeding.'], ...
                 rho_sp(2,3));
    end


    %% ── STEP C: FIT TRIVARIATE T-COPULA ─────────────────────────────────────

    fprintf('Fitting trivariate t-copula to %d concurrent observations...\n', n_conc);
    [rho_3d, nu_3d] = copulafit('t', U3, 'Method', 'ApproximateML');

    fprintf('  Degrees of freedom: nu = %.2f\n', nu_3d);
    fprintf('  Correlation matrix  [tide, surge, Hs]:\n');
    fprintf('    tide   %+.4f  %+.4f  %+.4f\n', rho_3d(1,:));
    fprintf('    surge  %+.4f  %+.4f  %+.4f\n', rho_3d(2,:));
    fprintf('    Hs     %+.4f  %+.4f  %+.4f\n', rho_3d(3,:));

    % Pairwise upper tail dependence coefficients
    % lambda_U = 2 * t_{nu+1}(-sqrt((nu+1)(1-rho)/(1+rho)))
    % lambda_U = 0  → asymptotic independence (copula recovers convolution result)
    % lambda_U > 0  → co-occurrence of extremes; copula raises the return level
    lambda_U_fn = @(r) 2 * tcdf(-sqrt((nu_3d+1) * (1-r) / (1+r)), nu_3d+1);
    lambda_U_ts = lambda_U_fn(rho_3d(1,2));   % tide–surge
    lambda_U_tH = lambda_U_fn(rho_3d(1,3));   % tide–Hs
    lambda_U_sH = lambda_U_fn(rho_3d(2,3));   % surge–Hs

    % Classify each pairwise dependence for use in output and figures:
    %   near-indep.  |rho| < 0.10 and lambda_U < 0.05
    %                → copula recovers independence; result ≈ convolution
    %   weak dep.    |rho| < 0.25 or lambda_U < 0.15
    %                → copula has a modest but real effect on the return level
    %   IMPORTANT    |rho| >= 0.25 or lambda_U >= 0.15
    %                → copula significantly raises the 100-year return level
    [dep_lbl_ts, dep_clr_ts] = dep_classify(rho_3d(1,2), lambda_U_ts);
    [dep_lbl_tH, dep_clr_tH] = dep_classify(rho_3d(1,3), lambda_U_tH);
    [dep_lbl_sH, dep_clr_sH] = dep_classify(rho_3d(2,3), lambda_U_sH);

    fprintf('  Upper tail dependence coefficients:\n');
    fprintf('    lambda_U  tide–surge : %.4f  [%s]\n', lambda_U_ts, dep_lbl_ts);
    fprintf('    lambda_U  tide–Hs    : %.4f  [%s]\n', lambda_U_tH, dep_lbl_tH);
    fprintf('    lambda_U  surge–Hs   : %.4f  [%s]\n', lambda_U_sH, dep_lbl_sH);


    %% ── STEP D: TP POWER-LAW REGRESSION ─────────────────────────────────────
    % Fit  log(Tp) = a*log(Hs) + b  to storm-peak pairs, giving
    %   Tp_pred = exp(b) * Hs^a
    % Applied per-sample in the Monte Carlo so that Tp scales with the
    % sampled Hs rather than being fixed at a single extreme-storm value.

    p_tp   = polyfit(log(storm_peak_Hs), log(storm_peak_Tp), 1);
    a_tp   = p_tp(1);   % exponent
    b_tp   = p_tp(2);   % log-intercept

    fprintf('Tp regression (raw fit):  Tp = %.3f * Hs^{%.3f}\n', exp(b_tp), a_tp);

    % A negative exponent means Tp weakly decreases with Hs, which is
    % physically counterintuitive (larger storms carry longer-period swell).
    % When a_tp < 0 the trend is flat or inverted — typically a statistical
    % artefact of limited data and high scatter in the Hs–Tp relationship.
    % Clamping to zero uses a fixed Tp = exp(b_tp), which is more defensible
    % than extrapolating a spurious downward trend into the GPD tail where no
    % storm-peak observations exist.
    if a_tp < 0
        fprintf('  Note: exponent is negative (flat/inverted trend) — clamping to 0.\n');
        fprintf('  Using fixed Tp = %.2f s for all Monte Carlo samples.\n', exp(b_tp));
        a_tp = 0;
    end
    fprintf('Tp regression (applied): Tp = %.3f * Hs^{%.3f}\n', exp(b_tp), a_tp);

    % Update Tp100 to be consistent with the regression rather than the
    % conditional mean computed earlier.  This ensures the diagnostic
    % triangle on the scatter plot sits on the regression line, making
    % the two visually consistent.
    Tp100 = exp(b_tp) * Hs100 .^ a_tp;
    fprintf('  Tp100 (regression at Hs100, diagnostic) = %.1f s\n', Tp100);

    % Compute log-space residual standard deviation from the APPLIED regression.
    % log(Tp_obs) - (a_tp*log(Hs) + b_tp) gives the residual for each storm peak.
    % sigma_tp is the standard deviation of these residuals in log space, which
    % corresponds to a lognormal multiplicative scatter around the regression line.
    % A physical minimum of 2 s is enforced after sampling to prevent unrealistically
    % short periods from propagating into the runup calculation.
    log_resid = log(storm_peak_Tp) - (a_tp * log(storm_peak_Hs) + b_tp);
    sigma_tp  = std(log_resid);
    fprintf('  Tp log-space scatter  sigma = %.3f  (%.0f%% lognormal CV)\n', ...
            sigma_tp, 100 * (exp(sigma_tp) - 1));

    % ── DIAGNOSTIC: 2%-exceedance wave contribution at (Hs100, Tp100) ────────
    % Deterministic, single-point evaluation of the chosen shoreline wave
    % model at the diagnostic 100-yr Hs/Tp pair. This is NOT how the
    % Monte-Carlo TWL is built (that convolves the full joint tide-surge-Hs
    % distribution) — it is a representative "how much does wave runup add"
    % number for reporting, e.g. the \WRU macro in the LaTeX summary.
    switch upper(WAVE_MODEL)
        case 'STOCKDON'
            R2_100 = Stockdon2006(Hs100, Tp100, beta);
        case 'BECKER'
            % Representative co-occurring tide/surge for this diagnostic only
            % (the Monte Carlo itself samples the full joint distribution):
            %   tide_rep  = 0   -> MHHW, a typical high-tide assumption
            %   surge_rep = 99th percentile of the concurrent storm surge
            %               used to fit the copula, i.e. a plausible storm
            %               surge co-occurring with an extreme wave event.
            tide_rep  = 0;
            surge_rep = prctile(surge_conc, 99);
            [R2_100, eta_bar_100, eta_prime_100, hr_100] = ...
                Becker2026(Hs100, Tp100, tide_rep, surge_rep, D0);
    end
    fprintf('  R2_100 (%s, diagnostic wave contribution at Hs100/Tp100) = %.2f m\n', ...
            WAVE_MODEL, R2_100);

    % Diagnostic plot: observed storm-peak Hs–Tp with fitted regression
    % and ±1σ / ±2σ lognormal scatter envelopes.
    Hs_plt    = linspace(min(storm_peak_Hs), max(storm_peak_Hs), 200);
    Tp_plt    = exp(b_tp           ) * Hs_plt .^ a_tp;   % regression mean
    Tp_p1sig  = exp(b_tp + sigma_tp) * Hs_plt .^ a_tp;   % +1σ envelope
    Tp_m1sig  = exp(b_tp - sigma_tp) * Hs_plt .^ a_tp;   % −1σ envelope
    Tp_p2sig  = exp(b_tp + 2*sigma_tp) * Hs_plt .^ a_tp; % +2σ envelope
    Tp_m2sig  = exp(b_tp - 2*sigma_tp) * Hs_plt .^ a_tp; % −2σ envelope

    figure;
    % ±2σ shaded band (patch drawn first so it sits behind everything)
    fill([Hs_plt, fliplr(Hs_plt)], [Tp_p2sig, fliplr(Tp_m2sig)], ...
         [1 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.20, ...
         'DisplayName', '\pm2\sigma envelope');
    hold on;
    % ±1σ shaded band
    fill([Hs_plt, fliplr(Hs_plt)], [Tp_p1sig, fliplr(Tp_m1sig)], ...
         [1 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.35, ...
         'DisplayName', '\pm1\sigma envelope');
    scatter(storm_peak_Hs, storm_peak_Tp, 20, 'k', 'filled', ...
            'DisplayName', 'Storm peaks');
    if a_tp == 0
        reg_label = sprintf('Fixed T_p = %.2f s  (exponent clamped to 0)', exp(b_tp));
    else
        reg_label = sprintf('Power law  T_p = %.2f H_s^{%.2f}', exp(b_tp), a_tp);
    end
    plot(Hs_plt, Tp_plt,   'r-',  'LineWidth', 2,   'DisplayName', reg_label);
    plot(Hs_plt, Tp_p1sig, 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    plot(Hs_plt, Tp_m1sig, 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    plot(Hs_plt, Tp_p2sig, 'r:',  'LineWidth', 1.2, 'HandleVisibility', 'off');
    plot(Hs_plt, Tp_m2sig, 'r:',  'LineWidth', 1.2, 'HandleVisibility', 'off');
    scatter(Hs100, Tp100, 120, 'r', 'filled', '^', ...
            'DisplayName', 'H_{s100}, T_{p100}  (diagnostic)');
    xlabel('H_s [m]', 'fontsize', 12);
    ylabel('T_p [s]', 'fontsize', 12);
    title('Storm-peak T_p–H_s regression  (shading: \pm1\sigma, \pm2\sigma)', ...
          'fontsize', 12);
    legend('Location', 'northwest', 'fontsize', 11);
    grid on; hold off;


    %% ── STEP E: HYBRID Hs MARGINAL ───────────────────────────────────────────
    % Below the POT threshold → empirical quantile (interp1 on sorted concurrent Hs).
    % Above the POT threshold → GPD tail (pd_GPD from the POT fit above).
    %
    % This lets the Monte Carlo sample Hs beyond the observed record maximum
    % for the rare, high-U draws that dominate the 100-year TWL.
    %
    % p_thresh_Hs is the empirical CDF value of u_waves in the CONCURRENT
    % dataset (the one whose ranks define U_Hs), ensuring consistency between
    % the copula marginals and the back-transform split point.

    Hs_conc_sorted = sort(Hs_conc);
    q_knots_Hs     = ((1:n_conc)' - 0.5) / n_conc;
    p_thresh_Hs    = mean(Hs_conc <= u_waves);   % F(u_waves) in concurrent data

    fprintf('Hybrid Hs marginal: empirical below u = %.2f m (p = %.4f), GPD above.\n', ...
            u_waves, p_thresh_Hs);

    % Pre-sort marginals for the Monte Carlo back-transform.
    %
    % TIDE: use the FULL hourly predicted tide record (predi_v, NaN-removed),
    % not the 12-hourly concurrent subset.  By Sklar's theorem, the copula
    % dependence structure and the marginal distributions are separable — the
    % copula is correctly fitted to concurrent 12-hourly ranks, while the
    % tidal marginal CDF is better estimated from all ~650,000 hourly values
    % (~75 yr × 8760 h).  The full record gives a much smoother empirical
    % quantile function with accurate tidal extremes, which directly controls
    % the high-TWL tail.  The 12-hourly subsample (~21,000 values) is simply
    % a coarser estimate of the same underlying tidal distribution.
    %
    % SURGE: use the concurrent 12-hourly series (surge at ERA5 timestamps).
    % Surge is stochastic, so its distribution at ERA5 times — the population
    % on which the copula was trained — is the correct conditioning set.
    % Using the full hourly surge would import many calm-weather near-zero
    % values that are not representative of the joint (tide, surge, Hs)
    % population and would bias the surge tail downward.
    tide_full_sorted  = sort(predi_v);                             % full hourly tide marginal
    q_knots_full_tide = ((1:numel(predi_v))' - 0.5) / numel(predi_v);
    surge_conc_sorted = sort(surge_conc);                          % concurrent surge marginal
    q_knots_conc_surge = ((1:n_conc)' - 0.5) / n_conc;

    fprintf('Tide marginal: %d full hourly values (vs %d concurrent 12-h samples)\n', ...
            numel(predi_v), n_conc);


    %% ── STEP F: TRIVARIATE MONTE CARLO → TWL ────────────────────────────────
    % For each repetition, draw N_y ~ Poisson(lambda_waves) storm events for
    % each of 100 simulated years, back-transform to (tide, surge, Hs),
    % compute Tp and runup, sum to get TWL, and record the annual maximum.

    % Storm arrivals follow a Poisson process with rate lambda_waves [storms/yr].
    % Fixing every simulated year to exactly round(lambda_waves) events ignores
    % the year-to-year variability in storm count, which distorts the tail of
    % the annual-maximum distribution: years with anomalously many storms can
    % produce extreme TWL values, but a fixed-n model never generates them.
    % The compound-Poisson annual maximum CDF is:
    %   P(M_year ≤ x) = exp(-lambda_waves * (1 - F(x)))
    % which is the correct form assumed by the POT framework.  Drawing
    %   N_y ~ Poisson(lambda_waves)  for each simulated year
    % is consistent with this and removes the rounding artifact.
    % Note: P(N_y = 0) = exp(-lambda_waves) ≈ exp(-10) ~ 4.5e-5 for typical
    % storm rates, so zero-storm years are extremely rare and handled as NaN.

    fprintf(['Running trivariate Monte Carlo: %d reps × 100 yr, ' ...
             'N_y ~ Poisson(lambda = %.2f storms/yr)...\n'], repeats, lambda_waves);

    annual_maxima_TWL = NaN(100, repeats);   % NaN-initialised; zero-storm years stay NaN

    for mcarlo = 1:repeats

        % Draw per-year storm counts from Poisson(lambda_waves)
        N_per_year = poissrnd(lambda_waves, 100, 1);   % [100 × 1] integer storm counts
        N_total    = sum(N_per_year);

        if N_total == 0
            continue   % astronomically unlikely; skip this repetition
        end

        if USE_PRODUCT_COPULA
            U_samp = rand(N_total, 3);                      % product copula: independent uniforms
        else
            U_samp = copularnd('t', rho_3d, nu_3d, N_total); % fitted t-copula
        end

        % ── Restrict each event to be a true storm (Hs > u) ─────────────
        % N_y ~ Poisson(λ_s) is drawn from the POT-declustered storm rate,
        % which by definition counts only independent Hs peaks above u.
        % A copula draw with U_Hs uniform on [0,1] however sends ~95% of
        % events through the empirical (sub-threshold) branch of the
        % marginal, dropping the effective storm rate to λ_s * (1-p_thresh)
        % and dramatically undersampling the extreme tail. Rescaling U_Hs
        % onto (p_thresh_Hs, 1] enforces that every simulated event has
        % Hs > u, restoring consistency with the POT-derived λ_s and with
        % the paper's stated compound-Poisson formulation.
        %
        % For copulas that recover near-independence (typical for open-
        % coast sites where the joint tide/surge/wave dependence is weak),
        % this rescaling is equivalent to accept-reject filtering on U_Hs.
        % For strongly dependent copulas (rare at open-coast sites but
        % possible in tropical cyclone basins where surge and Hs are
        % tightly coupled) it is a mild approximation to the conditional
        % joint distribution — accept-reject would be exact but slower.
        U_samp(:, 3) = p_thresh_Hs + (1 - p_thresh_Hs) * U_samp(:, 3);

        % ── Back-transform tide (column 1) — full hourly marginal ────────
        % Uses the full NOAA hourly predicted tide record (not the 12-h
        % concurrent subset) for a smoother, more accurate tidal CDF.
        % Justified by Sklar's theorem: copula dependence and marginals
        % are separable, so the marginal can be replaced independently.
        tide_samp = interp1(q_knots_full_tide, tide_full_sorted, ...
                            U_samp(:,1), 'linear', 'extrap');

        % ── Back-transform surge (column 2) — concurrent marginal ─────────
        % Intentionally uses the 12-hourly concurrent surge, not the full
        % hourly record: surge is stochastic and its distribution at ERA5
        % timestamps is the correct conditioning population for this copula.
        surge_samp = interp1(q_knots_conc_surge, surge_conc_sorted, ...
                             U_samp(:,2), 'linear', 'extrap');

        % ── Back-transform Hs (column 3) — hybrid empirical / GPD ─────────
        u_Hs_col  = U_samp(:,3);
        Hs_samp   = zeros(N_total, 1);

        % Body (below POT threshold): empirical quantile
        below = u_Hs_col <= p_thresh_Hs;
        Hs_samp(below) = interp1(q_knots_Hs, Hs_conc_sorted, ...
                                 u_Hs_col(below), 'linear', 'extrap');

        % Tail (above POT threshold): rescale to conditional probability, invert GPD
        %   P(Hs > x | Hs > u) = (U_Hs - p_thresh) / (1 - p_thresh)
        p_cond = (u_Hs_col(~below) - p_thresh_Hs) / (1 - p_thresh_Hs);
        p_cond = min(max(p_cond, 0), 1 - 1e-10);   % clamp to valid range for icdf
        Hs_samp(~below) = u_waves + icdf(pd_GPD, p_cond);

        % ── Physical Hs cap ──────────────────────────────────────────────
        % When the fitted GPD has ξ > 0 (heavy tail, unbounded), rare deep-
        % tail draws from copularnd can produce unphysically large Hs. Cap
        % at Hs_max (user parameter, see USER_INPUT). This keeps the linear-
        % in-Hs Becker runup — and to a lesser extent Stockdon — from
        % amplifying nonphysical Hs into TWL outliers that break the GEV fit.
        Hs_samp = min(Hs_samp, Hs_max);

        % ── Tp from power-law regression + lognormal scatter (per sample) ───
        % Each sample draws an independent log-space residual ε ~ N(0, sigma_tp²),
        % giving Tp = exp(b_tp + ε) * Hs^a_tp.  This captures the real scatter in
        % the Hs–Tp relationship (different storm types: local wind sea vs. swell)
        % and propagates Tp uncertainty into the runup calculation.
        % Physical minimum of 2 s enforced to prevent degenerate Stockdon inputs.
        eps_tp  = sigma_tp * randn(N_total, 1);
        Tp_samp = max(exp(b_tp + eps_tp) .* Hs_samp .^ a_tp, 2.0);

        % ── Shoreline wave contribution (dispatch on WAVE_MODEL) ───────────
        % R2_samp holds the wave-driven contribution to TWL for each sampled
        % storm event. Under Stockdon this is the 2%-exceedance runup R_2%.
        % Under Becker it is the total 2%-exceedance shoreline water level
        % due to waves (setup + variable SS/LF component). Either way it
        % enters the TWL sum the same way, so the copula machinery and the
        % downstream GEV fit are unchanged.
        switch upper(WAVE_MODEL)
            case 'STOCKDON'
                R2_samp = Stockdon2006(Hs_samp, Tp_samp, beta);
            case 'BECKER'
                R2_samp = Becker2026(Hs_samp, Tp_samp, tide_samp, surge_samp, D0);
            otherwise
                error(['Unknown WAVE_MODEL "%s". Set WAVE_MODEL = ''Stockdon'' ' ...
                       'or ''Becker'' in the USER_INPUT section.'], WAVE_MODEL);
        end

        % ── Total Water Level ──────────────────────────────────────────────
        TWL_samp = tide_samp + surge_samp + R2_samp;

        % ── Annual maxima — one per simulated year ─────────────────────────
        % Use cumsum of per-year storm counts to index into the flat TWL_samp
        % vector without looping over individual events.
        cum_n = [0; cumsum(N_per_year)];
        for y = 1:100
            if N_per_year(y) > 0
                annual_maxima_TWL(y, mcarlo) = max(TWL_samp(cum_n(y)+1 : cum_n(y+1)));
            end
            % N_per_year(y) == 0 leaves annual_maxima_TWL(y, mcarlo) = NaN
        end

    end

    fprintf('Monte Carlo complete.\n');


    %% ── 100-YEAR RETURN LEVEL FROM ANNUAL MAXIMA ─────────────────────────────
    % NaN entries (zero-storm years, P ≈ exp(-lambda) per year) are excluded.
    % With lambda_waves ≳ 5, fewer than 1 in 150 year-samples is expected to
    % be NaN, so excluding them has negligible effect on the return-level
    % estimate.
    %
    % The 100-year TWL is reported as the direct empirical 99th percentile
    % of the pooled annual maxima rather than the 99th percentile of a fitted
    % GEV distribution. With N_b × 100 = 30,000 annual maxima and roughly
    % 300 samples above the 99th percentile, the empirical estimator has a
    % standard error of only a few centimetres and does not depend on the
    % stability of a maximum-likelihood shape parameter. GEV MLE is known
    % to be biased for heavy-tailed data (ξ > 0), a regime common at
    % typhoon- and hurricane-exposed sites, and further biased downward
    % when a physical Hs cap piles mass at the upper end of the annual
    % maxima. The GEV distribution is still fitted for the diagnostic
    % plot below so the shape of the annual-max distribution can be
    % visually inspected, but it does not enter the reported TWL_100.

    valid_maxima = annual_maxima_TWL(~isnan(annual_maxima_TWL(:)));
    value_99     = prctile(valid_maxima, 99);   % 100-year TWL, empirical (rel. to MHHW)
    pd_TWL_gev   = fitdist(valid_maxima, 'gev');   % diagnostic GEV fit for plotting
    value_99_gev = icdf(pd_TWL_gev, 0.99);          % diagnostic GEV-based 99th percentile

    % Sanity check: if the GEV-based and empirical percentiles disagree by
    % more than ~10 cm, print a warning. Larger discrepancies typically mean
    % the GEV MLE has undershot the shape parameter (common with ξ_gpd > 0)
    % and the empirical estimate is the one to trust.
    if abs(value_99_gev - value_99) > 0.10
        fprintf(['Note: GEV-based 99th pct (%.3f m) differs from empirical ' ...
                 '99th pct (%.3f m) by %+.2f m.\n' ...
                 '      GEV shape parameter fitted at ξ = %.3f; expected ~ξ_gpd = %.3f.\n' ...
                 '      Reported TWL_100 uses the empirical percentile.\n'], ...
                value_99_gev, value_99, value_99_gev - value_99, ...
                pd_TWL_gev.k, pd_GPD.k);
    end

    % Determine overall copula importance for summary line
    all_near_indep = strcmp(dep_lbl_ts, 'near-indep.') && ...
                     strcmp(dep_lbl_tH, 'near-indep.') && ...
                     strcmp(dep_lbl_sH, 'near-indep.');
    any_important  = strcmp(dep_lbl_ts, 'IMPORTANT')   || ...
                     strcmp(dep_lbl_tH, 'IMPORTANT')   || ...
                     strcmp(dep_lbl_sH, 'IMPORTANT');

    if any_important
        copula_note = 'Copula captures meaningful dependence — result differs from independence assumption.';
    elseif all_near_indep
        copula_note = 'All pairs near-independent — copula recovers the independence case at this site.';
    else
        copula_note = 'Weak dependence detected — copula has a modest effect on the return level.';
    end


    %% ── TRIVARIATE PLOTS ─────────────────────────────────────────────────────

    % 1. Pairwise empirical copula scatter (3 panels)
    %    Subplot title colour signals dependence importance:
    %      grey  = near-independent (copula recovers independence)
    %      blue  = weak dependence
    %      red   = IMPORTANT (copula materially affects return level)
    figure;
    subplot(1,3,1);
    scatter(U_tide, U_surge, 1, 'k', 'filled', 'MarkerFaceAlpha', 0.05);
    xlabel('Uniform tide', 'fontsize', 11);
    ylabel('Uniform surge', 'fontsize', 11);
    ht = title(sprintf('tide–surge  \\rho=%.3f  [%s]', rho_3d(1,2), dep_lbl_ts), ...
               'fontsize', 11);
    ht.Color = dep_clr_ts;
    axis square; grid on;

    subplot(1,3,2);
    scatter(U_tide, U_Hs, 1, 'k', 'filled', 'MarkerFaceAlpha', 0.05);
    xlabel('Uniform tide', 'fontsize', 11);
    ylabel('Uniform H_s', 'fontsize', 11);
    ht = title(sprintf('tide–H_s  \\rho=%.3f  [%s]', rho_3d(1,3), dep_lbl_tH), ...
               'fontsize', 11);
    ht.Color = dep_clr_tH;
    axis square; grid on;

    subplot(1,3,3);
    scatter(U_surge, U_Hs, 1, 'k', 'filled', 'MarkerFaceAlpha', 0.05);
    xlabel('Uniform surge', 'fontsize', 11);
    ylabel('Uniform H_s', 'fontsize', 11);
    ht = title(sprintf('surge–H_s  \\rho=%.3f  [%s]', rho_3d(2,3), dep_lbl_sH), ...
               'fontsize', 11);
    ht.Color = dep_clr_sH;
    axis square; grid on;

    sgtitle({'Empirical trivariate copula (uniform marginals)', ...
             'Title colour: grey = near-indep.  |  blue = weak dep.  |  red = IMPORTANT'}, ...
            'fontsize', 12);

    % 2. Chi-plots: upper tail dependence for all three pairwise combinations
    %    χ(q) = P(V > q | U > q) converges to λ_U as q → 1.
    %    The independence baseline is 1−q (dashed black).
    %    The red horizontal line is the parametric λ_U from the fitted t-copula.
    %    Note: the empirical curve has high variance near q = 1 (few joint
    %    exceedances), so it rarely reaches the parametric λ_U exactly.
    q_levels  = 0.80:0.01:0.99;
    chi_indep = 1 - q_levels;   % independence baseline (same for all pairs)

    chi_ts = arrayfun(@(q) mean(U_surge > q & U_tide  > q) / ...
                           max(mean(U_tide  > q), eps), q_levels);   % tide–surge
    chi_tH = arrayfun(@(q) mean(U_Hs    > q & U_tide  > q) / ...
                           max(mean(U_tide  > q), eps), q_levels);   % tide–Hs
    chi_sH = arrayfun(@(q) mean(U_Hs    > q & U_surge > q) / ...
                           max(mean(U_surge > q), eps), q_levels);   % surge–Hs

    % Pack into struct arrays for compact subplot loop
    chi_pairs = {chi_ts,          chi_tH,         chi_sH        };
    lam_pairs = {lambda_U_ts,     lambda_U_tH,    lambda_U_sH   };
    lbl_pairs = {dep_lbl_ts,      dep_lbl_tH,     dep_lbl_sH    };
    clr_pairs = {dep_clr_ts,      dep_clr_tH,     dep_clr_sH    };
    ttl_pairs = {'tide vs surge', 'tide vs H_s',  'surge vs H_s'};
    ylb_pairs = {'P(surge > q | tide > q)', ...
                 'P(H_s > q | tide > q)',   ...
                 'P(H_s > q | surge > q)'  };

    figure;
    for kp = 1:3
        subplot(1, 3, kp);
        hold on;
        plot(q_levels, chi_pairs{kp}, 'b-o', 'LineWidth', 2, ...
             'DisplayName', 'Empirical \chi(q)');
        plot(q_levels, chi_indep, 'k--', 'LineWidth', 1.5, ...
             'DisplayName', 'Independence baseline');
        yline(lam_pairs{kp}, 'r-', ...
              sprintf('\\lambda_U = %.3f', lam_pairs{kp}), ...
              'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left', ...
              'DisplayName', sprintf('Fitted \\lambda_U = %.3f  (asymptotic)', lam_pairs{kp}));
        text(0.805, 0.97 * max([chi_pairs{kp}, chi_indep]), ...
             sprintf('%s', lbl_pairs{kp}), ...
             'FontSize', 10, 'Color', clr_pairs{kp}, 'FontWeight', 'bold', ...
             'VerticalAlignment', 'top');
        xlabel('Quantile level  q', 'fontsize', 11);
        ylabel(['\chi(q)  =  ' ylb_pairs{kp}], 'fontsize', 10);
        ht = title(['Upper tail dep.: ' ttl_pairs{kp}], 'fontsize', 11);
        ht.Color = clr_pairs{kp};
        if kp == 1
            legend('Location', 'southwest', 'fontsize', 9);
        end
        grid on; hold off;
    end
    sgtitle({'Pairwise upper tail dependence  \chi(q)', ...
             'Title colour: grey = near-indep.  |  blue = weak dep.  |  red = IMPORTANT'}, ...
            'fontsize', 12);

    % 3. GEV fit to TWL annual maxima
    figure;
    histogram(valid_maxima, 60, 'Normalization', 'pdf', ...
              'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none');
    hold on;
    x_gev = linspace(min(valid_maxima), max(valid_maxima), 300);
    plot(x_gev, pdf(pd_TWL_gev, x_gev), 'r-', 'LineWidth', 2.5);
    xline(value_99, 'k--', 'LineWidth', 2, ...
          'Label', sprintf('100-yr  %.3f m', value_99), ...
          'LabelHorizontalAlignment', 'left', 'fontsize', 12);
    xlabel('TWL  [m above MHHW]', 'fontsize', 14);
    ylabel('Probability Density', 'fontsize', 14);
    legend('Annual maxima', 'Fitted GEV', 'Location', 'northeast', 'fontsize', 12);
    title('TWL annual maxima — trivariate t-copula (tide + surge + runup)', 'fontsize', 14);
    grid on; box on; hold off;

    % Plot WRT NAVD88 (only if the NAVD88 datum is defined at this site).
    % For OCONUS stations (e.g. Apra Harbor, Guam, station 1630000) NAVD88
    % is not defined and NAVD88_to_EGM will be NaN. In that case skip the
    % datum-shifted plot with a warning rather than erroring on fitdist.
    offset = MHHW + NAVD88_to_EGM;
    if isfinite(offset)
        VM_NAVD88 = valid_maxima + offset;
        pd_TWL_gevNAVD88 = fitdist(VM_NAVD88, 'gev');
        figure;
        histogram(VM_NAVD88, 60, 'Normalization', 'pdf', ...
                  'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none');
        hold on;
        x_gev = linspace(min(VM_NAVD88), max(VM_NAVD88), 300);
        plot(x_gev, pdf(pd_TWL_gevNAVD88, x_gev), 'r-', 'LineWidth', 2.5);
        xline(value_99+offset, 'k--', 'LineWidth', 2);
        xlabel('TWL  [m]', 'fontsize', 18);
        ylabel('Probability Density', 'fontsize', 18);
        legend('Annual maxima', 'Fitted GEV', 'Location', 'northeast', 'fontsize', 18);
        set(gca,'fontsize',18)
        grid on; box on; hold off;
    else
        warning(['NAVD88 datum offset is not finite (NAVD88_to_EGM = %g). ' ...
                 'Skipping NAVD88-referenced GEV plot. Set MHHW and ' ...
                 'NAVD88_to_EGM in the USER_INPUT section, or leave both ' ...
                 'as NaN if the site does not use NAVD88.'], NAVD88_to_EGM);
    end

elseif CALC_WAVES == 0

    %% ── BIVARIATE COPULA (tide, surge) — CALC_WAVES = 0 FALLBACK ────────────
    %
    % Steps:
    %   1. Map tide and surge to uniform marginals via rank-based
    %      pseudo-observations (Weibull plotting position).
    %   2. Fit the chosen copula to the uniform pairs.
    %   3. Report lambda_U as a diagnostic.
    %   4. Monte Carlo: sample correlated (tide, surge) pairs, back-transform
    %      to physical units, extract annual maxima, fit GEV.

    % --- Step 1: Transform marginals to uniform [0,1] via pseudo-observations ---
    % tiedrank() + Weibull denominator is O(N log N) and avoids the O(N²)
    % kernel CDF evaluation that causes a hang on long records.
    n_data  = length(predi_v);
    u_tide  = tiedrank(predi_v) / (n_data + 1);
    u_surge = tiedrank(SS)      / (n_data + 1);

    U_data = [u_tide, u_surge];   % n × 2 matrix of uniform pairs

    % --- Pre-compute empirical quantile functions for Monte Carlo back-transform ---
    predi_sorted = sort(predi_v);
    SS_sorted    = sort(SS);
    q_knots      = ((1:n_data)' - 0.5) / n_data;

    % --- Step 2: Fit the copula ---
    fprintf('Fitting %s copula to %d hourly (tide, surge) pairs...\n', ...
            COPULA_FAMILY, length(u_tide));

    if strcmpi(COPULA_FAMILY, 't')
        [rho_cop, nu_cop] = copulafit('t', U_data, 'Method', 'ApproximateML');
        rho_val = rho_cop(1, 2);
        fprintf('  t-copula fit:  rho = %.4f,  nu = %.2f\n', rho_val, nu_cop);
        lambda_U = 2 * tcdf(-sqrt((nu_cop + 1) * (1 - rho_val) / (1 + rho_val)), nu_cop + 1);

    elseif strcmpi(COPULA_FAMILY, 'Gumbel')
        rho_spearman = corr(u_tide, u_surge);
        if rho_spearman < 0
            error(['COPULA_FAMILY = ''Gumbel'' requires non-negative tide–surge dependence,\n' ...
                   'but Spearman''s rho = %.4f at this site indicates negative dependence.\n' ...
                   'Switch to COPULA_FAMILY = ''t''.'], rho_spearman);
        end
        theta_cop = copulafit('Gumbel', U_data);
        fprintf('  Gumbel copula fit:  theta = %.4f\n', theta_cop);
        lambda_U = 2 - 2^(1 / theta_cop);

    else
        error('Unrecognised COPULA_FAMILY ''%s''. Choose ''t'' or ''Gumbel''.', COPULA_FAMILY);
    end

    fprintf('  Upper tail dependence coefficient  lambda_U = %.4f\n', lambda_U);
    fprintf('  (0 = asymptotic independence; approaching 1 = strong tail co-occurrence)\n');

    % Classify tide–surge dependence
    [dep_lbl_ts, dep_clr_ts] = dep_classify(rho_val, lambda_U);
    fprintf('  Dependence classification: [%s]\n', dep_lbl_ts);

    % --- Step 3: Tail dependence diagnostic (chi-plot) ---
    q_levels   = 0.80 : 0.01 : 0.99;
    chi_emp    = zeros(size(q_levels));
    chi_indep  = zeros(size(q_levels));
    for qi = 1:length(q_levels)
        q = q_levels(qi);
        chi_emp(qi)   = mean(u_surge > q & u_tide > q) / max(mean(u_tide > q), eps);
        chi_indep(qi) = 1 - q;
    end

    figure;
    subplot(1, 2, 1);
    scatter(u_tide, u_surge, 1, 'k', 'filled', 'MarkerFaceAlpha', 0.05);
    xlabel('Uniform tide  u', 'fontsize', 12);
    ylabel('Uniform surge  v', 'fontsize', 12);
    ht = title(sprintf('tide–surge  \\rho=%.3f  [%s]', rho_val, dep_lbl_ts), 'fontsize', 12);
    ht.Color = dep_clr_ts;
    grid on; axis square;

    subplot(1, 2, 2);
    hold on;
    plot(q_levels, chi_emp,   'b-o', 'LineWidth', 2, 'DisplayName', 'Empirical \chi(q)');
    plot(q_levels, chi_indep, 'k--',  'LineWidth', 1.5, 'DisplayName', 'Independence baseline');
    yline(lambda_U, 'r-', sprintf('Fitted \\lambda_U = %.3f', lambda_U), ...
          'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left', ...
          'DisplayName', sprintf('Fitted \\lambda_U = %.3f  (asymptotic tail dep.)', lambda_U));
    xlabel('Quantile level  q', 'fontsize', 12);
    ylabel('\chi(q)  =  P(V > q | U > q)', 'fontsize', 12);
    title('Upper tail dependence diagnostic', 'fontsize', 12);
    legend('Location', 'northwest');
    grid on;
    sgtitle(sprintf('%s copula — tide vs surge', COPULA_FAMILY), 'fontsize', 13);

    % --- Step 4: Monte Carlo sampling from the bivariate copula ---
    % Derive n_eff_per_year from a POT declustering of the surge record,
    % exactly mirroring the wave analysis.  This removes tau_hours as a
    % tuning parameter for the return-level result: the number of independent
    % storm events per year is measured from the data rather than assumed.
    surge_thresh     = prctile(SS, 95);           % 95th-percentile surge threshold
    min_sep_surge    = ceil(tau_hours);            % minimum gap between surge events (hours; dt=1h)
    surge_exceed_idx = find(SS > surge_thresh);
    if ~isempty(surge_exceed_idx)
        gaps_surge         = diff(surge_exceed_idx);
        cluster_starts_sg  = [1; find(gaps_surge > min_sep_surge) + 1];
        n_surge_storms     = length(cluster_starts_sg);
    else
        n_surge_storms = 0;
    end
    n_years_surge  = numel(SS) / 8760;            % record length in years
    lambda_surge   = n_surge_storms / n_years_surge;

    fprintf('Surge POT declustering (u = %.3f m, 95th pct):  %d storms / %.1f yr  →  lambda = %.2f storms/yr\n', ...
            surge_thresh, n_surge_storms, n_years_surge, lambda_surge);
    fprintf('  N_y ~ Poisson(%.2f) drawn each simulated year (tau_hours used for separation only)\n', lambda_surge);
    fprintf('Running bivariate copula Monte Carlo (%d repetitions × 100 yr)...\n', repeats);

    annual_maxima = NaN(100, repeats);   % NaN-initialised; zero-storm years stay NaN

    for mcarlo = 1:repeats

        N_per_year = poissrnd(lambda_surge, 100, 1);   % [100 × 1] Poisson storm counts
        N_total    = sum(N_per_year);

        if N_total == 0
            continue
        end

        if strcmpi(COPULA_FAMILY, 't')
            U_samp = copularnd('t', rho_cop, nu_cop, N_total);
        elseif strcmpi(COPULA_FAMILY, 'Gumbel')
            U_samp = copularnd('Gumbel', theta_cop, N_total);
        end

        tide_samp  = interp1(q_knots, predi_sorted, U_samp(:, 1), 'linear', 'extrap');
        surge_samp = interp1(q_knots, SS_sorted,    U_samp(:, 2), 'linear', 'extrap');
        twl_samp   = tide_samp + surge_samp;   % R2 = 0 when CALC_WAVES = 0

        cum_n = [0; cumsum(N_per_year)];
        for y = 1:100
            if N_per_year(y) > 0
                annual_maxima(y, mcarlo) = max(twl_samp(cum_n(y)+1 : cum_n(y+1)));
            end
        end

    end

    % Fit GEV to pooled annual maxima (exclude NaN zero-storm years)
    valid_maxima_biv = annual_maxima(~isnan(annual_maxima(:)));
    pd_yrmax = fitdist(valid_maxima_biv, 'gev');
    value_99  = icdf(pd_yrmax, 0.99);

    one_percent_value = prctile(twl_samp, 99);

    % Copula note (bivariate)
    if strcmp(dep_lbl_ts, 'IMPORTANT')
        copula_note = 'Copula captures meaningful tide–surge dependence — result differs from independence.';
    elseif strcmp(dep_lbl_ts, 'near-indep.')
        copula_note = 'Tide–surge near-independent — copula recovers the independence case at this site.';
    else
        copula_note = 'Weak tide–surge dependence — copula has a modest effect on the return level.';
    end

    % --- Bivariate plots ---
    figure;
    histogram(twl_samp, 200, 'Normalization', 'pdf', 'FaceColor', [0.2 0.4 0.8]);
    hold on;
    xline(one_percent_value, '--r', '1% exceedance', 'LineWidth', 1.5, ...
          'LabelHorizontalAlignment', 'left');
    title(sprintf('%s copula — simulated tide + surge distribution (last rep.)', COPULA_FAMILY));
    xlabel('Water Level [m]');
    ylabel('PDF');
    grid on;
    hold off;

    figure;
    histogram(valid_maxima_biv, 'Normalization', 'pdf');
    hold on;
    yrTWL_pts = linspace(min(valid_maxima_biv), max(valid_maxima_biv), 100);
    y_values  = pdf(pd_yrmax, yrTWL_pts);
    plot(yrTWL_pts, y_values, 'r-', 'LineWidth', 2);
    xline(value_99, 'k--', 'LineWidth', 2);
    xlabel('Water Level above MHHW [m]');
    ylabel('Probability Density');
    legend('Annual maxima histogram', sprintf('Fitted GEV (%s copula)', COPULA_FAMILY), '100-yr value','fontsize',12);
    title(sprintf('Tide and Storm Surge — %s copula', COPULA_FAMILY),'fontsize',14);
    grid on; box on;
    hold off;

    % Set trivariate dependence labels to empty for output section
    dep_lbl_tH = '—'; dep_clr_tH = [0.5 0.5 0.5];
    dep_lbl_sH = '—'; dep_clr_sH = [0.5 0.5 0.5];
    lambda_U_ts = lambda_U;
    lambda_U_tH = NaN;
    lambda_U_sH = NaN;

end  % end CALC_WAVES if/elseif


%% Sea Level Rise time series plot
figure
plot(d_t,LPC,'r','linewidth',2)
hold on
plot(d_t,yfit,'b','linewidth',2)
ylabel('Water Depth [m]')
grid on
legend('Low Pass Water Level','Linear Fit')

% Storm surge scatter
figure
scatter(d_t(mask), SS, 1, 'k', 'filled')
ylabel('Storm Surge [m]')
grid on


%% Print Results

fprintf('\n══════════════════════════════════════════════\n');
fprintf('RESULTS\n');
fprintf('══════════════════════════════════════════════\n');

if CALC_WAVES == 1
    % Representative decomposition of the 100-yr TWL into a tide+surge part
    % and a wave-runup part, for LaTeX reporting (\TSS, \WRU macros below).
    % NOTE: this is a diagnostic decomposition only. value_99 is the 99th
    % percentile of the pooled (tide+surge+R2) annual maxima from the joint
    % Monte Carlo, NOT the sum of separately-computed 99th percentiles (the
    % three variables are correlated). Subtracting the deterministic R2_100
    % (evaluated at Hs100/Tp100, see above) from value_99 gives a
    % representative "how much of the 100-yr TWL is tide+surge vs. wave
    % runup" split that is internally consistent (TSS_100 + R2_100 = value_99)
    % but should not be read as two independent return levels.
    TSS_100 = value_99 - R2_100;

    fprintf('100-yr TWL (trivariate copula, tide + surge + runup):\n');
    fprintf('  %.4f m  relative to MHHW\n', value_99);
    if USE_NOAA == 1
        fprintf('  %.4f m  relative to EGM2008\n', value_99 + SLR + MHHW_el + NAVD88_to_EGM);
    end
    fprintf('  -- decomposition (diagnostic, see note in code) --\n');
    fprintf('  Tide + Storm Surge contribution   : %.4f m\n', TSS_100);
    fprintf('  Wave runup/setup contribution (R2): %.4f m  (%s @ Hs100/Tp100)\n', ...
            R2_100, WAVE_MODEL);
    fprintf('----------------------------------------------\n');
    fprintf('Copula: trivariate t,  nu = %.2f\n', nu_3d);
    fprintf('  lambda_U  tide–surge : %.4f  [%s]\n', lambda_U_ts, dep_lbl_ts);
    fprintf('  lambda_U  tide–Hs    : %.4f  [%s]\n', lambda_U_tH, dep_lbl_tH);
    fprintf('  lambda_U  surge–Hs   : %.4f  [%s]\n', lambda_U_sH, dep_lbl_sH);
    fprintf('  >> %s\n', copula_note);
    fprintf('----------------------------------------------\n');
    fprintf('Hs100 (POT+GPD, diagnostic only):  %.2f m\n', Hs100);
    fprintf('Tp100 (regression at Hs100, diagnostic): %.1f s\n', Tp100);
    fprintf('Tp log-space scatter sigma = %.3f  (lognormal CV %.0f%%)\n', ...
            sigma_tp, 100*(exp(sigma_tp)-1));
    if isfinite(Hs_max)
        fprintf('Physical Hs cap applied in Monte Carlo: Hs_max = %.1f m\n', Hs_max);
    end
    fprintf('Shoreline wave model: %s\n', WAVE_MODEL);
    if strcmpi(WAVE_MODEL, 'Stockdon')
        fprintf('  Foreshore beach slope  β  = %.3f\n', beta);
    elseif strcmpi(WAVE_MODEL, 'Becker')
        fprintf('  Reef-flat bed elevation D0 = %+.2f m (tide gauge datum)\n', D0);
        fprintf('  Diagnostic tide_rep = %.2f m (MHHW), surge_rep = %.2f m (99th pct concurrent surge)\n', ...
                tide_rep, surge_rep);
        fprintf('  eta_bar_100 (setup) = %.2f m, eta_prime_100 (SS+LF) = %.2f m, h_r_100 = %.2f m\n', ...
                eta_bar_100, eta_prime_100, hr_100);
    end
    fprintf('Sea level rise over record: %.2f mm/yr\n', P(1)*24*365*1000);
    fprintf('Storm rate (data-derived): lambda = %.2f storms/yr  (N_y ~ Poisson each simulated year)\n', lambda_waves);
    fprintf('  (tau_hours = %d h used for ERA5 POT declustering only; not a tuning parameter for TWL)\n', tau_hours);
    fprintf('  NaN annual-max entries (zero-storm years): %d of %d  (P ≈ %.2e per year)\n', ...
            sum(isnan(annual_maxima_TWL(:))), numel(annual_maxima_TWL), exp(-lambda_waves));
    fprintf('----------------------------------------------\n');
    fprintf('LaTeX summary macros (paste into report):\n');
    fprintf('  \\newcommand{\\TWE}{%.2f m}\n', value_99 + SLR + MHHW_el + NAVD88_to_EGM);
    fprintf('  \\newcommand{\\SLR}{%.2f mm/yr}\n', P(1)*24*365*1000);
    fprintf('  \\newcommand{\\WRU}{%.2f}\n', R2_100);
    fprintf('  \\newcommand{\\TSS}{%.2f m}\n', TSS_100);
    fprintf('  \\newcommand{\\StormRate}{%.2f~storms/yr}\n', lambda_waves);
else
    % No wave runup in the bivariate fallback (R2 = 0 by construction), so
    % the tide+surge contribution IS the reported 100-yr return level.
    R2_100  = 0;
    TSS_100 = value_99;

    fprintf('100-yr tide + surge return level (bivariate copula, R2 = 0):\n');
    fprintf('  %.4f m  relative to MHHW\n', value_99);
    fprintf('  (For comparison to NOAA annual maxima): %.4f m  (includes SLR)\n', value_99 + SLR);
    if USE_NOAA == 1
        fprintf('  %.4f m  relative to EGM2008\n', value_99 + SLR + MHHW_el + NAVD88_to_EGM);
    end
    fprintf('  Tide + Storm Surge contribution   : %.4f m\n', TSS_100);
    fprintf('  Wave runup/setup contribution (R2): 0.0000 m  (CALC_WAVES = 0)\n');
    fprintf('----------------------------------------------\n');
    fprintf('Copula: %s (bivariate)\n', COPULA_FAMILY);
    fprintf('  lambda_U  tide–surge : %.4f  [%s]\n', lambda_U_ts, dep_lbl_ts);
    fprintf('  >> %s\n', copula_note);
    fprintf('----------------------------------------------\n');
    fprintf('Sea level rise over record: %.2f mm/yr\n', P(1)*24*365*1000);
    fprintf('Storm rate (data-derived): lambda = %.2f storms/yr  (N_y ~ Poisson each simulated year)\n', lambda_surge);
    fprintf('  (tau_hours = %d h used for surge POT declustering only; not a tuning parameter for TWL)\n', tau_hours);
    fprintf('----------------------------------------------\n');
    fprintf('LaTeX summary macros (paste into report):\n');
    fprintf('  \\newcommand{\\TWE}{%.2f m}\n', value_99 + SLR + MHHW_el + NAVD88_to_EGM);
    fprintf('  \\newcommand{\\SLR}{%.2f mm/yr}\n', P(1)*24*365*1000);
    fprintf('  \\newcommand{\\WRU}{0.00}\n');
    fprintf('  \\newcommand{\\TSS}{%.2f m}\n', TSS_100);
    fprintf('  \\newcommand{\\StormRate}{%.2f~storms/yr}\n', lambda_surge);
end

if USE_NOAA == 1
    fprintf('Tide Station = %s\n', station);
    fprintf('MHHW = %.3f m,  NAVD88 to EGM2008 = %.3f m\n', MHHW, NAVD88_to_EGM);
end
if CALC_WAVES == 1
    fprintf('ERA5 location (lat, lon): %.2f, %.2f\n', lat, lon);
    fprintf('Beach slope beta = %.3f\n', beta);
    fprintf('Wave POT threshold = %.0f%%  (u = %.2f m)\n', Hs_thresh_pct, u_waves);
    fprintf('Storm rate lambda_waves = %.2f storms/yr  (%d storms / %.1f yr)\n', ...
            lambda_waves, n_storms, n_years_waves);
end
fprintf('══════════════════════════════════════════════\n');


%% Alternative check: GEV fit to annual maxima of observed water level (NOAA method)

% Extract the year from the datetime vector
years = year(d_t);

% Get unique years
unique_years = unique(years);

% Initialize a vector to store the annual maxima
annual_maxima_NOAA = zeros(length(unique_years), 1);

% Loop over each year and find the maximum value for that year
for i = 1:length(unique_years)
    % Find the indices for the current year
    year_idx = years == unique_years(i);

    % Get the maximum value for the current year
    annual_maxima_NOAA(i) = nanmax(obsv(year_idx));
end

% Save the results in T
T = table(unique_years, annual_maxima_NOAA, 'VariableNames', {'Year', 'AnnualMaxima'});
annualMax = T.AnnualMaxima;

% Find the 99% value
pd_yrmax_obs = fitdist(annualMax,'gev'); % Fit a generalized extreme value distribution
value_99_obs = icdf(pd_yrmax_obs, 0.99);

% Plot the evolution of yearly maxima
figure
scatter(unique_years, annual_maxima_NOAA, 'filled', 'k')
grid on
ylabel('Meters above MHHW')
title('Annual maxima of observed water level (NOAA method)')
xlabel('Year')


%% ── LOCAL FUNCTIONS ──────────────────────────────────────────────────────────

function [lbl, clr] = dep_classify(rho_val, lambda_U_val)
% dep_classify  Classify pairwise copula dependence strength.
%
%   [lbl, clr] = dep_classify(rho_val, lambda_U_val)
%
%   Returns a short label string and an RGB colour for use in terminal
%   output and figure titles.
%
%   Thresholds (Spearman rho / upper tail dependence lambda_U):
%     near-indep.  |rho| < 0.10 AND lambda_U < 0.05
%                  Copula recovers the independence (convolution) result.
%                  The copula is statistically valid but adds no return-level
%                  uplift over a simple independence assumption.
%
%     weak dep.    |rho| < 0.25 OR lambda_U < 0.15
%                  The copula modestly raises the return level.  Worth
%                  retaining in the analysis; consider sensitivity tests.
%
%     IMPORTANT    |rho| >= 0.25 OR lambda_U >= 0.15
%                  Meaningful tail co-occurrence.  The copula significantly
%                  raises the 100-year return level relative to the
%                  independence assumption.  Must not be ignored.

    if abs(rho_val) >= 0.25 || lambda_U_val >= 0.15
        lbl = 'IMPORTANT';
        clr = [0.80 0.10 0.10];   % red
    elseif abs(rho_val) >= 0.10 || lambda_U_val >= 0.05
        lbl = 'weak dep.';
        clr = [0.10 0.30 0.75];   % blue
    else
        lbl = 'near-indep.';
        clr = [0.50 0.50 0.50];   % grey
    end
end
