function [R2, setup, swash, sinc, sig] = Stockdon2006(Hs, Tp, beta)
    % Stockdon2006 Runup model from Stockdon et al. (2006)
    % 
    % Inputs:
    %   Hs   - Significant wave height (m)
    %   Tp   - Peak wave period (s)
    %   beta - Beach slope
    % 
    % Outputs:
    %   R2    - 2% exceedence runup level
    %   setup - Setup component
    %   swash - Total swash
    %   sinc  - Incident swash component
    %   sig   - Infragravity swash component
    %
    % Reference:
    %   Stockdon et al. (2006), Coastal Engineering, 53(7), 573–588.
    %   DOI: 10.1016/j.coastaleng.2005.12.005

    % Compute deep-water wavelength (Lp)
    g = 9.81; % gravitational acceleration (m/s^2)
    Lp = (g / (2 * pi)) * Tp.^2; % deep-water wavelength

    % Compute zeta parameter
    zeta = beta ./ sqrt(Hs ./ Lp);

    % Compute setup
    setup = 0.35 * beta .* sqrt(Hs .* Lp);

    % Compute swash components
    sinc = 0.75 * beta .* sqrt(Hs .* Lp);  % Incident swash
    sig  = 0.06 * sqrt(Hs .* Lp);          % Infragravity swash

    % Compute total swash
    swash = sqrt(sinc.^2 + sig.^2);

    % Compute R2 (2% exceedance runup level)
    R2 = 1.1 * (0.35 * beta .* sqrt(Hs .* Lp) + (sqrt(Hs .* Lp .* (0.563 * beta.^2 + 0.004)) / 2));

    % Apply dissipative condition (if zeta < 0.3, use alternative equation)
    dissipative_mask = zeta < 0.3;
    R2(dissipative_mask) = 0.043 * sqrt(Hs(dissipative_mask) .* Lp(dissipative_mask));

end