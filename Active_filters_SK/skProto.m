function P = skProto(type, n, varargin)
%SKPROTO  Analog low-pass prototype (poles/zeros/gain) of a filter type.
%   P = SKPROTO(TYPE, N) returns the normalized low-pass prototype of
%   order N for TYPE in {'butter','cheby1','cheby2','bessel','cauer'}.
%
%   Name-value options:
%     'Rp'  passband ripple in dB        (cheby1, cauer)  default 1
%     'Rs'  stopband attenuation in dB   (cheby2, cauer)  default 40
%     'Fs'  stopband edge frequency (rad/s, normalized to the passband
%           edge = 1). If given for 'cauer', the achievable Rs is
%           returned instead of the requested one.
%
%   Output P is a struct with fields
%     .p    column vector of poles           (LHP)
%     .z    column vector of zeros           (LP prototype, jw-axis zeros)
%     .k    gain so that |H(0)| = 1 (DC gain 1 for the LP prototype)
%     .wstop   stopband edge frequency (rad/s) for cauer/cheby2
%     .Rp, .Rs achieved ripple/attenuation
%     .info   short description
%
%   Normalization used:
%     butter/bessel : -3 dB point at w = 1 rad/s
%     cheby1        : passband edge (ripple Rp) at w = 1
%     cheby2        : stopband edge (attenuation Rs) at w = 1
%     cauer         : passband edge (ripple Rp) at w = 1, stopband edge at
%                     P.wstop with attenuation P.Rs
%
%   All prototypes are computed self-contained (no Signal Processing
%   Toolbox needed). The elliptic (Cauer) design uses the classical
%   Chebyshev rational function and the base-MATLAB functions ELLIPKE /
%   ELLIPJ only.
%
%   Example:  P = skProto('cauer', 5, 'Rp', 1, 'Rs', 40);
%
%   See also SKDESIGN.

% ----------------------------------------------------------------------
% options
opt.Rp = 1; opt.Rs = 40; opt.Fs = [];
opt = skGetOpt(opt, varargin{:});

if n < 1
    error('skProto:BadOrder', 'Order must be >= 1.');
end
type = lower(type);

switch type
    % ================= Butterworth =================
    case 'butter'
        p = exp(1i * pi * (2*(1:n)' + n - 1) / (2*n));
        z = [];
        k = 1;

    % ================= Chebyshev I =================
    case 'cheby1'
        if opt.Rp <= 0, error('skProto:Rp', 'Rp must be > 0 dB.'); end
        eps = sqrt(10^(opt.Rp/10) - 1);
        phi = asinh(1/eps) / n;
        th = pi * (2*(1:n)' - 1) / (2*n);
        p = -sinh(phi)*sin(th) + 1i*cosh(phi)*cos(th);
        z = [];
        k = 1;                       % |H(0)| = 1
        opt.Rs = nan;

    % ================= Chebyshev II =================
    case 'cheby2'
        if opt.Rs <= 0, error('skProto:Rs', 'Rs must be > 0 dB.'); end
        % poles = reciprocal of the Chebyshev-I poles for ripple eps(Rs)
        eps = 1 / sqrt(10^(opt.Rs/10) - 1);
        phi = asinh(1/eps) / n;
        th = pi * (2*(1:n)' - 1) / (2*n);
        p1 = -sinh(phi)*sin(th) + 1i*cosh(phi)*cos(th);
        p = 1 ./ p1;
        % transmission zeros on the jw axis
        m = floor(n/2);
        z = zeros(2*m, 1);
        for i = 1:m
            z(2*i-1) =  1i * sec(pi*(2*i-1)/(2*n));
            z(2*i)   = -1i * sec(pi*(2*i-1)/(2*n));
        end
        k = 1;                       % |H(0)| = 1
        opt.Rp = nan;

    % ================= Bessel =================
    case 'bessel'
        % reverse Bessel polynomials: th0=1, th1=s+1,
        % th_n = (2n-1)*th_{n-1} + s^2*th_{n-2}
        c0 = 1;
        c1 = [1 1];
        for kk = 2:n
            c2 = (2*kk-1)*[0 c1] + [c0 0 0];
            c0 = c1; c1 = c2;
        end
        if n == 1, c = [1 1]; else, c = c1; end
        p = roots(c);
        % normalize so the -3 dB point is at w = 1 rad/s
        H0 = c(end);
        f3 = @(w) abs(H0 / polyval(c, 1i*w));
        w = logspace(-1, 1.5, 20001);
        Hw = arrayfun(f3, w);
        w3 = w(find(Hw <= 1/sqrt(2), 1, 'first'));
        if isempty(w3) || ~isfinite(w3)
            w3 = 1;
        end
        p = p / w3;                % -3 dB point now at w = 1 rad/s
        z = [];
        k = abs(prod(p));          % |H(0)| = k / prod(-p) = 1

    % ================= Cauer (elliptic) =================
    case 'cauer'
        if opt.Rp <= 0, error('skProto:Rp', 'Rp must be > 0 dB.'); end
        if opt.Rs <= 0, error('skProto:Rs', 'Rs must be > 0 dB.'); end
        [p, z, k, ws, RsOut] = skEllipticProto(n, opt.Rp, opt.Rs, opt.Fs);
        opt.Rs = RsOut;
        P.wstop = ws;

    otherwise
        error('skProto:BadType', ...
            'Unknown type "%s". Use butter|cheby1|cheby2|bessel|cauer.', type);
end

% ---- gain: |H(0)| = 1 ----
if ~strcmp(type, 'cauer')
    % k already set; verify/complete for zero types
    if ~isempty(z)
        k = abs(prod(p)) / abs(prod(z));   % H(0) = k*prod(-z)/prod(-p) = 1
    end
end

% ---- first-order pole goes first, keep pairs ----
P.p = p(:);
P.z = z(:);
P.k = k;
P.type = type;
P.order = n;
P.Rp = opt.Rp;
P.Rs = opt.Rs;
if ~isfield(P, 'wstop')
    switch type
        case {'cheby2'}, P.wstop = 1;
        case {'cauer'},  % set inside skEllipticProto
        otherwise,       P.wstop = inf;
    end
end
P.info = sprintf('%s order %d prototype', type, n);
end

% ======================================================================
function [p, z, k, ws, RsOut] = skEllipticProto(n, Rp, Rs, Fs)
%SK ELLIPTIC PROTO  Elliptic (Cauer) LP prototype, passband edge = 1.
%
%   k1 = eps/es (discrimination); order equation n = K(k)K1'/(K'(k)K1)
%   gives the selectivity k; transmission zeros at 1/(k*cd((2i-1)K/n,k)).
%   The Chebyshev rational function is
%       R_n(w) = C * w^a * prod(w^2-z_i^2) / prod(w^2-p_i^2)
%   with C fixed by R_n(1) = 1; the filter poles are the LHP roots of
%   R_n(s/j) = +-j/eps.

eps = sqrt(10^(Rp/10) - 1);
es  = sqrt(10^(Rs/10) - 1);
k1  = eps / es;
K1  = ellipke(k1^2);
K1p = ellipke(1 - k1^2);

if isempty(Fs)
    % selectivity from the order equation
    f = @(kk) n*ellipke(1-kk.^2)*K1 - ellipke(kk.^2)*K1p;
    k = fzero(f, [1e-8, 1-1e-8]);
    ws = 1 / k;
    RsOut = Rs;                 % requested Rs is met exactly
else
    % fixed stopband edge -> compute achievable Rs
    ws = Fs;
    k = 1 / ws;
    % order eq: n = K(k) K1'/(K'(k) K1)  ->  n K'(k) K1(x) - K(k) K1'(x) = 0
    fk1 = @(x) n*ellipke(1-k^2)*ellipke(x^2) - ellipke(k^2)*ellipke(1-x^2);
    k1 = fzero(fk1, [1e-6, 1-1e-6]);
    K1 = ellipke(k1^2);
    es = eps / k1;
    RsOut = 10*log10(es^2 + 1); % achievable stopband attenuation
end

K = ellipke(k^2);
m = floor(n/2);

% zeros of the rational function R_n (cd values) and its poles
% (the transmission-zero frequencies of the filter)
cdv = zeros(1, m); pp = zeros(1, m);
for i = 1:m
    th = (2*i-1) * K / n;
    [~, cnv, dnv] = ellipj(th, k^2);
    cdv(i) = cnv / dnv;             % cd(th), in (0, 1]
    pp(i)  = 1 / (k * cdv(i));      % transmission-zero frequency, > 1
end

% filter transmission zeros on the jw axis
z = zeros(2*m, 1);
for i = 1:m
    z(2*i-1) =  1i * pp(i);
    z(2*i)   = -1i * pp(i);
end

% poles: LHP roots of R_n(s/j) = +-j/eps
% R_n(w) = C * w^a * prod(w^2 - cdv^2) / prod(w^2 - pp^2)
C = prod(1 - pp.^2) / prod(1 - cdv.^2);    % R_n(1) = 1
if mod(n, 2) == 0
    A = poly([1i*cdv, -1i*cdv]);
    B = poly([1i*pp,  -1i*pp]);
    P1 = C*A - (1i/eps)*B;
    P2 = C*A + (1i/eps)*B;
else
    A = poly([1i*cdv, -1i*cdv, 0]);
    B = poly([1i*pp,  -1i*pp]);
    B = [0, B];
    P1 = (C/1i)*A - (1i/eps)*B;
    P2 = (C/1i)*A + (1i/eps)*B;
end
r = [roots(P1); roots(P2)];
p = sort(r(real(r) < 0));

% gain: |H(0)| = 1  ->  k = prod(-p)/prod(wz^2)  [even]; odd adds |p0|
wz2 = pp.^2;
if mod(n, 2) == 0
    k = abs(prod(p)) / prod(wz2);
else
    [~, im] = min(abs(imag(p)));
    p0 = p(im);
    others = p(imag(p) ~= 0);
    k = abs(p0) * abs(prod(others)) / prod(wz2);
end
end

% ======================================================================
function opt = skGetOpt(opt, varargin)
%SK GET OPT  Parse name/value option pairs into a struct.
if mod(numel(varargin), 2) ~= 0
    error('skGetOpt:OddArgs', 'Options must be name/value pairs.');
end
for i = 1:2:numel(varargin)
    nm = varargin{i};
    if ~ischar(nm) && ~isstring(nm)
        error('skGetOpt:BadName', 'Option names must be strings.');
    end
    if ~isfield(opt, nm)
        error('skGetOpt:UnknownOpt', 'Unknown option "%s".', nm);
    end
    opt.(nm) = varargin{i+1};
end
end
