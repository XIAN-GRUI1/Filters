function s = skFmtSI(x, unit)
%SKFMTSI  Format a numeric value with an SI prefix.
%   S = SKFMTSI(3.3e3)   ->  '3.3k'
%   S = SKFMTSI(47e-9)   ->  '47n'
%   S = SKFMTSI(1e3,'Hz')->  '1kHz'
%
%   See also SKESERIES.

if nargin < 2, unit = ''; end
p = [1e-12 1e-9 1e-6 1e-3 1 1e3 1e6 1e9];
l = {'p', 'n', 'u', 'm', '', 'k', 'M', 'G'};
ax = abs(x);
if ~isfinite(ax)
    s = sprintf('%g%s', x, unit);
    return;
end
i = find(ax >= p, 1, 'last');
if isempty(i), i = 4; end
if i > 8, i = 8; end
v = x / p(i);
s = sprintf('%.4g%s%s', v, l{i}, unit);
end
