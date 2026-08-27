function dd = ddm2dd(coordStr)
%DDM2DD Convert degrees-decimal-minutes coordinate strings to decimal degrees.
%
%   dd = DDM2DD(coordStr)
%
% Parses strings of the form "48° 6.7 N" / "122° 45.6 W" into signed
% decimal degrees. Designed to be tolerant of formatting variations
% (degree symbol present or not, minute mark present or not, leading or
% trailing hemisphere letter, extra spaces).
%
% INPUT
%   coordStr : one of
%                - a single char vector or string scalar, e.g. "48° 6.7 N"
%                - a cell array of strings
%                - a string array
%
% OUTPUT
%   dd : decimal degrees, same size/shape as the input. North and East
%        are positive; South and West are negative.
%
% SIGN CONVENTION
%   If the string ends (or starts) with N/S/E/W (case-insensitive), that
%   letter sets the sign: N/E -> positive, S/W -> negative.
%   If there is no hemisphere letter at all, a leading '-' on the
%   degrees value is respected instead (e.g. "-122 45.6" -> -122.76).
%
% EXAMPLES
%   ddm2dd("48° 6.7 N")                    -> 48.1117
%   ddm2dd("122° 45.6 W")                  -> -122.76
%   ddm2dd({"48° 6.7 N", "122° 45.6 W"})   -> [48.1117  -122.76]
%
%   % A lat/lon pair, as you'd typically use it:
%   lat = ddm2dd("48° 6.7 N");
%   lon = ddm2dd("122° 45.6 W");

if iscell(coordStr)
    dd = cellfun(@parseOneDDM, coordStr);
    return
end

if isstring(coordStr) && numel(coordStr) > 1
    dd = arrayfun(@(s) parseOneDDM(char(s)), coordStr);
    return
end

dd = parseOneDDM(char(coordStr));

end

function val = parseOneDDM(s)
% Parse a single degrees-decimal-minutes string into decimal degrees.

s = strtrim(s);
if isempty(s)
    error('ddm2dd:emptyInput', 'Received an empty coordinate string.');
end

% Hemisphere letter, checked at the start and end of the string.
hemiLead  = regexp(s, '^\s*([NSEWnsew])', 'tokens', 'once');
hemiTrail = regexp(s, '([NSEWnsew])\s*$', 'tokens', 'once');

hemiLetter = '';
if ~isempty(hemiLead)
    hemiLetter = hemiLead{1};
elseif ~isempty(hemiTrail)
    hemiLetter = hemiTrail{1};
end

% Numeric tokens in order of appearance: degrees, then minutes (if any).
% This naturally skips over °, ', spaces, commas, hemisphere letters, etc.
nums = regexp(s, '-?\d+\.?\d*', 'match');
if isempty(nums)
    error('ddm2dd:parseError', 'Could not find a numeric degrees value in "%s".', s);
end

degTok = nums{1};
degVal = abs(str2double(degTok));

if numel(nums) >= 2
    minVal = str2double(nums{2});
else
    minVal = 0;
end

val = degVal + minVal / 60;

if ~isempty(hemiLetter)
    if any(upper(hemiLetter) == 'SW')
        val = -val;
    end
elseif degTok(1) == '-'
    val = -val;
end

end