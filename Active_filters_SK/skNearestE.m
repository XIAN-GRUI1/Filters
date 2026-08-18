function r = skNearestE(ser, target, dir)
%SKNEARESTE  Nearest nominal E-series value (with decade handling).
%   R = SKNEARESTE('E24'|'E96', TARGET)
%   R = SKNEARESTE(..., 'up')    smallest nominal value >= TARGET
%   R = SKNEARESTE(..., 'down')  largest  nominal value <= TARGET
%
%   TARGET may be a vector; R has the same size.
%
%   Example: skNearestE('E96', 47300) -> 47500
%
%   See also SKESERIES, SKSERIESGRID.

if nargin < 3 || isempty(dir), dir = 'nearest'; end
v = skESeries(ser);
sz = size(target);
target = target(:);
r = nan(size(target));
if any(~isfinite(target))
    fi = isfinite(target) & target > 0;
else
    fi = true(size(target));
end
if ~any(fi), r = reshape(r, sz); return; end
t = target(fi);
t = t(:).';
d = 10.^floor(log10(t));
d = d(:).';
tn = t ./ d;
switch lower(dir)
    case 'nearest'
        [~, i] = min(abs(v.' - tn), [], 1);
        r(fi) = v(i) .* d;
    case 'up'
        i = zeros(size(tn));
        for k = 1:numel(tn)
            j = find(v >= tn(k) - 1e-12, 1, 'first');
            if isempty(j), j = 1; d(k) = d(k)*10; end
            i(k) = j;
        end
        r(fi) = v(i) .* d;
    case 'down'
        i = zeros(size(tn));
        for k = 1:numel(tn)
            j = find(v <= tn(k) + 1e-12, 1, 'last');
            if isempty(j), j = numel(v); d(k) = d(k)/10; end
            i(k) = j;
        end
        r(fi) = v(i) .* d;
    otherwise
        error('skNearestE:BadDir', 'dir must be nearest|up|down');
end
r = reshape(r, sz);
end
