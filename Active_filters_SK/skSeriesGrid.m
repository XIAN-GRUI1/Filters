function v = skSeriesGrid(ser, vmin, vmax)
%SKSERIESGRID  All nominal E-series values within [VMIN, VMAX].
%   V = SKSERIESGRID('E24'|'E96', VMIN, VMAX)
%
%   Returns a sorted row vector of every preferred value (all decades)
%   that lies inside the requested range. Used to build the component
%   search grid of the Sallen-Key design.
%
%   See also SKESERIES, SKNEARESTE.

base = skESeries(ser);
d0 = floor(log10(vmin));
d1 = ceil(log10(vmax));
v = [];
for d = d0:d1
    v = [v, base * 10^d]; %#ok<AGROW>
end
v = v(v >= vmin - 1e-15 & v <= vmax + 1e-15);
v = unique(v);
end
