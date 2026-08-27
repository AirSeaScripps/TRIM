function [eta2, eta_bar, eta_prime, hr] = Becker2026(H0, Tp, tide, surge, D0)
    % Becker2026  2%-exceedance shoreline water level due to waves over a
    %             dissipative fringing reef, following Becker & Merrifield (2026).
    %
    % Bulk parameterization of the wave-driven shoreline water level from
    % back-refracted deep-water wave statistics and reef-flat water depth.
    % Intended as a drop-in replacement for Stockdon2006 at sites fronted by
    % a dissipative fringing reef (δ > 1 under typical conditions), where
    % Stockdon's sandy-beach runup assumption is not physically appropriate.
    %
    %   Setup:                 η̄       = p₀(h′) × H₀                (Eqs 13-14)
    %                          p₀(h′)  = 0.194 − 0.063 × h′
    %
    %   Variable component:    η′      = 0.009 × h_r × √(H₀ × L₀)   (Eq 15)
    %
    %   2% exceedance:         η₂      = η̄ + η′                     (Eq 16)
    %
    %   Reef-flat depth:       h_r     = h′ + h_ntr + η̄ − D₀
    %
    % where L₀ = g Tp² / (2π) is the deep-water wavelength. All water-level
    % quantities are referenced to the same vertical datum as the tide gauge.
    %
    % Inputs:
    %   H0    - Deep-water significant wave height [m]. For ERA5 offshore Hs
    %           taken at a nearest grid point in deep water (>~200 m), Hs ≈ H0
    %           and no back-refraction is required. For sites where the source
    %           record is on the continental shelf, back-refract per Eq 4 of
    %           the paper before calling this function.
    %   Tp    - Peak wave period [s].
    %   tide  - Astronomical tidal elevation h′ [m, relative to tide gauge datum].
    %   surge - Non-tidal residual water level h_ntr [m, same datum].
    %   D0    - Mean reef-flat BED elevation relative to the tide gauge datum [m].
    %           NEGATIVE below datum (typical for submerged reef).
    %           Example: Hagåtña Bay, Guam (144.764°E, 13.481°N) — extracted
    %                    from the 2020 NOAA NGS Topobathy Lidar over a 100 m
    %                    window on the reef flat: D0 = -0.63 m (std 0.19 m).
    %
    %           IMPORTANT: D0 must be referenced to the same vertical datum
    %           as `tide` and `surge`. Topobathymetric lidar is often
    %           delivered in ellipsoidal heights (NAD83(2011), etc.) or Mean
    %           High Water; NOAA tide gauges are typically referenced to a
    %           tidal datum (MLLW, MSL) or NAVD88 equivalent. Convert using
    %           VDatum or an equivalent tool before using D0 here.
    %
    % Outputs:
    %   eta2      - 2%-exceedance shoreline water level due to waves [m].
    %               This is the quantity that plays the role of R2 (Stockdon)
    %               in the TWL sum:  TWL = tide + surge + eta2.
    %   eta_bar   - Breaking-wave setup at the shoreline [m].
    %   eta_prime - Variable (SS + LF wave) contribution at the shoreline [m].
    %   hr        - Reef-flat water depth including setup [m], returned so the
    %               calling script can diagnose flooded-reef vs exposed-reef
    %               conditions.
    %
    % Reference:
    %   Becker, J. M. & Merrifield, M. A. (2026). Estimates of extreme and
    %   low-frequency reef-flat water levels at Ipan, Guam from offshore waves.
    %   Journal of Geophysical Research: Oceans.
    %
    % Applicability caveats (from the paper):
    %   - Validated at Ipan reef, Guam (shore-attached dissipative fringing
    %     reef, δ > 1 under typical wave conditions).
    %   - Applicable to reefs of similar character.
    %   - NOT appropriate for wave-like reefs (δ < 1, e.g. Majuro, Roi-Namur
    %     RMI) where the shoreline SS wave contribution must be added
    %     separately (see Becker & Merrifield 2026 §5).
    %   - NOT validated on non-reef configurations — use Stockdon2006 for
    %     sandy beaches or open-coast configurations.
    %
    % Note on the LF coefficient 0.009:
    %   Derived from Ipan reef parameters (damping rate D = 0.008 s⁻¹,
    %   cross-reef length L = 420 m) via 0.009 = 0.0105 × g / (D² L²).
    %   Applying this coefficient at a substantially different reef geometry
    %   (D or L different by more than a factor of ~2) will require
    %   recalibration; see Becker et al. (2016) for measurements of D at
    %   other Pacific reefs (e.g. D = 0.011 s⁻¹ for Saipan, 0.004 s⁻¹ for
    %   Roi-Namur RMI).

    % --- Physical constants ---
    g = 9.81;                              % gravitational acceleration [m/s^2]

    % --- Deep-water wavelength ---
    L0 = (g / (2 * pi)) * Tp .^ 2;

    % --- Breaking-wave setup (Eqs 13-14) ---
    % Tidally-dependent proportionality between H0 and setup.
    % At h′ = 0: p0 ≈ 0.194 → η̄ ≈ 0.19 H0 (matches Fig 9 black line, Ipan).
    % The tidal correction shifts p0 by ~±10% over Guam's tidal range,
    % capturing the observed increase of setup at low tide (Becker et al. 2014).
    p0 = 0.194 - 0.063 .* tide;
    eta_bar = p0 .* H0;

    % --- Reef-flat water depth including setup (Becker's h_r) ---
    % Sign convention: D0 is the reef-bed elevation vs the same datum as
    % tide and surge. D0 < 0 for a submerged reef. h_r ends up positive at
    % Ipan because η̄ typically dominates during extreme events.
    hr = tide + surge + eta_bar - D0;

    % Physical floor: when the reef flat would be exposed (h_r ≤ 0), no LF
    % dynamics on the reef are possible and the variable component collapses
    % to zero. Below h_r ≈ 0.3 m the paper's regression was not fitted
    % (Becker & Merrifield 2026 omit data with h_r < 0.3 m from Eq 15's
    % regression), so results with small h_r should be treated with caution.
    hr = max(hr, 0);

    % --- Variable (SS + LF) shoreline component (Eq 15) ---
    eta_prime = 0.009 .* hr .* sqrt(H0 .* L0);

    % --- Total 2%-exceedance shoreline water level due to waves (Eq 16) ---
    eta2 = eta_bar + eta_prime;

end
