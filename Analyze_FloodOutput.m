clear all

% Path is relative to the repo root; run this script from there (or adjust
% the path below to point at your local copy of Data/).
load("Data/CACC_CoastlineElevations.mat")

%% Extract columns from table
AlongtrackC = CoastlineElevations.Alongtrack_C;  % Along-track distance (m), CoSMoS
AlongtrackT = CoastlineElevations.Alongtrack_T;  % Along-track distance (m), TRCM
ShorelineC  = CoastlineElevations.Shoreline_Cm;   % Flood elevation (m), CoSMoS
ShorelineT  = CoastlineElevations.Shoreline_Tm;   % Flood elevation (m), TRCM
NorthingC   = CoastlineElevations.Northing_C;    % Northing coordinate (m), CoSMoS
NorthingT   = CoastlineElevations.Northing_T;    % Northing coordinate (m), TRCM

%% 1. Spatially low-pass filter the flood elevation along the coastline
% Filter is designed in the spatial domain (cycles per meter).
% Adjust smooth_length_m to change the smoothing scale.
smooth_length_m = 100;  % Smoothing cutoff wavelength (m)

% --- CoSMoS ---
dx_C  = mean(diff(AlongtrackC), 'omitnan');   % Mean along-track spacing (m)
Fs_C  = 1 / dx_C;                             % Spatial sampling rate (cycles/m)
Fc_C  = 1 / smooth_length_m;                  % Cutoff spatial frequency (cycles/m)
Wn_C  = Fc_C / (Fs_C / 2);                    % Normalized cutoff (0-1, Nyquist = 1)
[z, p, k]    = butter(4, Wn_C, 'low');
[sos_C, g_C] = zp2sos(z, p, k);
% Interpolate over NaN gaps, filter, then restore NaNs
nanC = isnan(ShorelineC);
ShorelineC_filled = ShorelineC;
ShorelineC_filled(nanC) = interp1(find(~nanC), ShorelineC(~nanC), find(nanC), 'linear', 'extrap');
ShorelineC_smooth = filtfilt(sos_C, g_C, ShorelineC_filled);
ShorelineC_smooth(nanC) = NaN;  % Re-mask original gap locations

% --- TRCM ---
dx_T  = mean(diff(AlongtrackT), 'omitnan');
Fs_T  = 1 / dx_T;
Fc_T  = 1 / smooth_length_m;
Wn_T  = Fc_T / (Fs_T / 2);
[z, p, k]    = butter(4, Wn_T, 'low');
[sos_T, g_T] = zp2sos(z, p, k);
% Interpolate over NaN gaps, filter, then restore NaNs
nanT = isnan(ShorelineT);
ShorelineT_filled = ShorelineT;
ShorelineT_filled(nanT) = interp1(find(~nanT), ShorelineT(~nanT), find(nanT), 'linear', 'extrap');
ShorelineT_smooth = filtfilt(sos_T, g_T, ShorelineT_filled);
ShorelineT_smooth(nanT) = NaN;  % Re-mask original gap locations

%% 2. Plot smoothed flood elevations vs Northing
figure;
set(gcf, 'Units', 'inches', 'Position', [1, 1, 10, 2.5]);
hold on;
plot(NorthingT-3603280, ShorelineT_smooth, 'r-', 'LineWidth', 1.5, 'DisplayName', 'TRIM');
plot(NorthingC-3603280, ShorelineC_smooth, 'b-', 'LineWidth', 1.5, 'DisplayName', 'CoSMoS');
hold off;
xlabel('Alongshore Distance (m)');
ylabel('Flood Elevation (m)');
% title(sprintf('Smoothed Shoreline Flood Elevation vs Northing  [\\lambda_c = %d m]', smooth_length_m));
legend('Location', 'best');
grid on;
box on;
axis tight;
set(gca, 'fontname', 'times', 'fontsize', 12, 'xdir', 'reverse')