function F = skDesign(varargin)
%SKDESIGN  Design an analog active filter with Sallen-Key sections.
%
%   F = SKDESIGN('Type', 'butter', 'Order', 6, 'Fc', 1e3)
%   F = SKDESIGN('Type', 'cheby1', 'Order', 5, 'Fc', 1e3, 'Rp', 1, ...
%                'PassType', 'highpass', 'GBW', 1e6, 'Plot', true)
%
%   Name-value options:
%     Type       'butter' (default) | 'cheby1' | 'cheby2' | 'bessel' | 'cauer'
%     Order      filter order (>= 1)
%     Fc         cutoff / passband-edge frequency in Hz (default 1000)
%     PassType   'lowpass' (default) | 'highpass'
%     Rp         passband ripple in dB for cheby1 / cauer (default 1)
%     Rs         stopband attenuation in dB for cheby2 / cauer (default 40)
%     Fstop      stopband edge frequency in Hz for cauer (optional; if
%                given, the achievable Rs is computed instead)
%     Rser       resistor series: 'E96' (default) | 'E24'
%     Cser       capacitor series: 'E24' (default)
%     GBW        op-amp gain-bandwidth product in Hz (default Inf = ideal).
%                With a finite GBW the component values are pre-distorted
%                so that the realized response (with the single-pole
%                op-amp model) matches the ideal one - but only for the
%                sections that actually need it (see Predist).
%     Predist    'auto' (default) | 'always' | 'never'
%                'auto' pre-distorts a section only when its GBW margin
%                m = GBW/(f0*Q) (= wt/(w0*Q)) is below PredistThresh.
%                The finite op-amp bandwidth raises the section Q and
%                shifts f0 by ~1/m (verified by exact nodal analysis),
%                so m >= 100 keeps the uncompensated Q error below ~1 %
%                and pre-distortion is not needed there.
%     PredistThresh  margin threshold for 'auto' (default 100)
%     PredistTol  relative (w0,Q) tolerance used to decide, per section,
%                 whether the finite op-amp bandwidth actually distorts
%                 the realized pole pair (default 0.01 = 1 %)
%     Rmin, Rmax  resistor search range in ohm (default 100 ... 1e6)
%     Cmin, Cmax  capacitor search range in farad (default 1e-12 ... 1e-6)
%     SpreadMax   maximum allowed component ratio warning threshold
%                 (default 1000)
%     Plot        plot the response (default false)
%     Verbose     print the design report (default true)
%
%   Output F is a struct with
%     .type, .order, .fc, .passtype, .gbw (rad/s), .options
%     .prototype   prototype struct (see SKPROTO)
%     .sections    array of section structs (components + achieved values;
%                  every Sallen-Key section is unity gain; sections are
%                  ordered by ascending Q so the low-Q stages come first
%                  and the high-Q stages last - best dynamic range)
%     .specs       measured passband ripple / -3 dB / stopband attenuation
%     .report      formatted text report
%
%   Topology notes:
%     * Butterworth / Chebyshev-I / Bessel are all-pole filters and are
%       realized entirely with unity-gain Sallen-Key sections.
%     * Cauer / Chebyshev-II have finite transmission zeros. A unity-gain
%       Sallen-Key cannot realize jw-axis zeros, so the sections that
%       contain a zero use a 3-op-amp state-variable biquad (the standard
%       textbook solution); the real-pole section still uses Sallen-Key.
%     * Every op-amp is modelled as a single integrator A(s) = GBW/s;
%       the finite-GBW pre-distortion is performed automatically by the
%       component search (see GBW option).
%
%   Example:
%     F = skDesign('Type','cauer','Order',5,'Fc',1e3,'Rp',1,'Rs',40,...
%                  'GBW',1e6,'Plot',true);
%
%   See also SKPROTO, SKOPT, SKBIQUAD, SKRESPONSE.

% ======================================================================
% options
opt.Type      = 'butter';
opt.Order     = 6;
opt.Fc        = 1e3;
opt.PassType  = 'lowpass';
opt.Rp        = 1;
opt.Rs        = 40;
opt.Fstop     = [];
opt.Rser      = 'E96';
opt.Cser      = 'E24';
opt.GBW       = Inf;
opt.Predist   = 'auto';
opt.PredistThresh = 100;
opt.PredistTol = 0.01;
opt.Rmin      = 100;
opt.Rmax      = 1e6;
opt.Cmin      = 1e-12;
opt.Cmax      = 1e-6;
opt.SpreadMax = 1000;
opt.Plot      = false;
opt.Verbose   = true;
opt = skParseOpt(opt, varargin{:});

type = lower(opt.Type);
wc = 2*pi*opt.Fc;
wt = inf;
if isfinite(opt.GBW)
    wt = 2*pi*opt.GBW;
end

% High-pass filters with finite transmission zeros (Cauer / Chebyshev-II)
% need a high-pass biquad with the zero below the pole.  The unity-gain
% Sallen-Key cannot realize jw-axis zeros at all, and the state-variable
% biquad included here is only valid for w_z > w_0 (low-pass).  A correct
% high-pass zero section needs a dedicated high-pass biquad (e.g. the
% 4-op-amp state-variable filter), which is not part of this package.
if strcmpi(opt.PassType, 'highpass') && ...
        (strcmp(type, 'cauer') || strcmp(type, 'cheby2'))
    error('skDesign:HPZeros', ...
        ['High-pass %s filters with transmission zeros are not ' ...
         'supported by this package: the unity-gain Sallen-Key cannot ' ...
         'realize jw-axis zeros, and the included state-variable biquad ' ...
         'requires the zero above the pole (low-pass). Use the low-pass ' ...
         'version, an all-pole type (butter/cheby1/bessel), or a ' ...
         'dedicated high-pass biquad.'], upper(type));
end

% ======================================================================
% prototype
if isempty(opt.Fstop)
    P = skProto(type, opt.Order, 'Rp', opt.Rp, 'Rs', opt.Rs);
else
    P = skProto(type, opt.Order, 'Rp', opt.Rp, 'Rs', opt.Rs, ...
        'Fs', opt.Fstop / opt.Fc);   % normalized stopband edge
end

% ======================================================================
% split the prototype into sections
[p, z] = localSplit(P);
% p : struct array of poles with fields w0p, Qp (pairs) or w0p (single)
% z : zero frequencies (rad/s of the prototype) paired with the sections
nsec = numel(p);
sections = struct([]);
for i = 1:nsec
    if ~isempty(p(i).Qp)
        % biquad section
        w0p = p(i).w0p; Qp = p(i).Qp;
        if strcmpi(opt.PassType, 'lowpass')
            w0 = w0p * wc;
        else
            w0 = wc / w0p;
        end
        wzp = z{i};
        if isempty(wzp)
            kind = 'sk';
            sec = skOpt('sk', opt.PassType, w0, Qp, [], ...
                skOptOpt(opt, wt));
        else
            kind = 'sv';
            if strcmpi(opt.PassType, 'lowpass')
                wz = wzp * wc;
            else
                wz = wc / wzp;
            end
            sec = skOpt('sv', opt.PassType, w0, Qp, wz, skOptOpt(opt, wt));
        end
        sec.kind = kind;
        sec.wz = wzp;            % prototype zero frequency
        sec.w0p = w0p;
    else
        % first-order section
        w0p = p(i).w0p;
        if strcmpi(opt.PassType, 'lowpass')
            w0 = w0p * wc;
        else
            w0 = wc / w0p;
        end
        sec = skOpt('rc', opt.PassType, w0, 0.5, [], skOptOpt(opt, wt));
        sec.kind = 'rc';
    end
    % unify the field set so the section structs can be concatenated
    allf = {'C','R','Rx','w3','C1','C2','R1','R2','R1x','R2x','w0a','Qa', ...
            'comp','comp0','w0t','Qt','wzt','err','wt','kind','topology','wz','w0p', ...
            'pred','margin'};
    for f = allf
        if ~isfield(sec, f{1})
            sec.(f{1}) = [];
        end
    end
    sections = [sections, sec]; %#ok<AGROW>
end

% order sections by ascending Q (real-pole sections first)
qs = zeros(1, numel(sections));
for i = 1:numel(sections)
    if strcmp(sections(i).kind, 'rc')
        qs(i) = 0;
    else
        qs(i) = sections(i).Qt;
    end
end
[~, ord] = sort(qs);
sections = sections(ord);

% ======================================================================
% assemble the filter struct
F.type = type;
F.order = opt.Order;
F.fc = opt.Fc;
F.passtype = lower(opt.PassType);
F.gbw = wt;
F.options = opt;
F.prototype = P;
F.sections = sections;

% response + specs
F.specs = skSpecs(F);

% report
F.report = skReport(F);

if opt.Verbose
    fprintf('%s\n', F.report);
end

if opt.Plot
    skPlot(F);
end
end

% ======================================================================
function o = skOptOpt(opt, wt)
% build the skOpt option struct
o = struct('Rser', opt.Rser, 'Cser', opt.Cser, ...
    'Rmin', opt.Rmin, 'Rmax', opt.Rmax, ...
    'Cmin', opt.Cmin, 'Cmax', opt.Cmax, ...
    'wt', wt, 'pred', opt.Predist, 'pthresh', opt.PredistThresh, ...
    'rtol', opt.PredistTol);
end

% ======================================================================
function [p, z] = localSplit(P)
% split prototype poles into pairs + single real pole; pair the
% transmission zeros with the pole pairs (sorted by frequency).
pp = P.p(:);
zz = P.z(:);
pairs_w0 = []; pairs_Q = [];
for i = 1:numel(pp)
    if imag(pp(i)) > 1e-9 * abs(pp(i))     % count each conjugate pair once
        w0 = abs(pp(i));
        pairs_w0(end+1) = w0;          %#ok<AGROW>
        pairs_Q(end+1) = w0 / (2*abs(real(pp(i)))); %#ok<AGROW>
    end
end
p = struct([]);
npair = numel(pairs_w0);
for i = 1:npair
    p(end+1).w0p = pairs_w0(i); %#ok<AGROW>
    p(end).Qp = pairs_Q(i);
end
for i = 1:numel(pp)
    if abs(imag(pp(i))) <= 1e-9 * abs(pp(i))
        p(end+1).w0p = abs(pp(i)); %#ok<AGROW>
    end
end
% zeros (positive frequencies), paired by ascending frequency
zF = sort(abs(zz(1:2:end)));
z = cell(1, numel(p));
nz = numel(zF);
if nz > 0
    [~, ord] = sort(pairs_w0);
    for i = 1:min(nz, npair)
        z{ord(i)} = zF(i);
    end
end
end

% ======================================================================
function opt = skParseOpt(opt, varargin)
if mod(numel(varargin), 2) ~= 0
    error('skDesign:OddArgs', 'Options must be name/value pairs.');
end
for i = 1:2:numel(varargin)
    nm = varargin{i};
    if ~ischar(nm) && ~isstring(nm)
        error('skDesign:BadName', 'Option names must be strings.');
    end
    switch nm
        case {'Type','Order','Fc','PassType','Rp','Rs','Fstop','Rser','Cser', ...
              'GBW','Predist','PredistThresh','PredistTol','Rmin','Rmax','Cmin','Cmax', ...
              'SpreadMax','Plot','Verbose'}
            opt.(nm) = varargin{i+1};
        otherwise
            error('skDesign:UnknownOpt', 'Unknown option "%s".', nm);
    end
end
% normalize
opt.Type = lower(opt.Type);
opt.PassType = lower(opt.PassType);
opt.Rser = upper(opt.Rser);
opt.Cser = upper(opt.Cser);
opt.Predist = lower(opt.Predist);
if ~ismember(opt.Predist, {'auto','always','never'})
    error('skDesign:BadPredist', ...
        'Predist must be ''auto'', ''always'' or ''never''.');
end
end
