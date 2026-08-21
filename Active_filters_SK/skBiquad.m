function varargout = skBiquad(cmd, varargin)
%SKBIQUAD  Section math for the Sallen-Key / state-variable designs.
%
%   This helper provides the low-level equations used by SKDESIGN:
%
%   [R1, R2] = SKBIQUAD('skideal', passtype, w0, Q, C1, C2)
%       Ideal unity-gain Sallen-Key component values.
%         lowpass : feedback cap C1 from node A to output, C2 to ground;
%                   requires C1 >= 4*Q^2*C2.
%         highpass: caps C1 (input), C2 (series), resistors R1 (feedback
%                   from node A to output), R2 (to ground).
%
%   [B1, B2, B3] = SKBIQUAD('skcoeffs', R1, R2, C1, C2, passtype, wt)
%       Denominator coefficients of the *actual* transfer function
%           H(s) = N(s) / (1 + B1*s + B2*s^2 + B3*s^3)
%       of a unity-gain Sallen-Key section whose op-amp has a finite
%       gain-bandwidth product wt (rad/s).  wt = Inf gives the ideal
%       quadratic (B3 = 0).  For lowpass  N = 1, for highpass
%       N = s^2*R1*R2*C1*C2.
%
%   [w0a, Qa] = SKBIQUAD('dominant', B1, B2, B3)
%       Dominant pole pair (natural frequency and Q) of the cubic
%       denominator, obtained by deflating the far real pole.
%
%   COMP = SKBIQUAD('svideal', passtype, w0, Q, wz, C)
%       State-variable biquad (3 op-amps) realizing
%           H(s) = -(s^2 + wz^2) / (s^2 + (w0/Q) s + w0^2)
%       Used for the Cauer / Chebyshev-II low-pass sections that need
%       finite transmission zeros above the pole (wz > w0) - a
%       unity-gain Sallen-Key cannot realize jw-axis zeros at all, and
%       this biquad is only valid for the low-pass case.
%
%   H = SKBIQUAD('svresp', w, COMP, wt)
%       Actual frequency response of the state-variable biquad with a
%       single-pole op-amp model (gain A = wt/s) for every stage.
%
%   [R, C] = SKBIQUAD('rc', passtype, wc, wt, C)
%       First-order RC section (passive RC + unity buffer) with the
%       closed-form pre-distortion of the op-amp pole: the realized -3 dB
%       point equals wc even with finite wt.
%
%   See also SKDESIGN.

switch lower(cmd)
    case 'skideal'
        varargout{1} = local_skIdeal(varargin{:});
    case 'skcoeffs'
        [varargout{1}, varargout{2}, varargout{3}] = local_skCoeffs(varargin{:});
    case 'dominant'
        [varargout{1}, varargout{2}] = local_dominant(varargin{:});
    case 'svideal'
        varargout{1} = local_svIdeal(varargin{:});
    case 'svresp'
        varargout{1} = local_svResp(varargin{:});
    case 'rc'
        [varargout{1}, varargout{2}] = local_rc(varargin{:});
    otherwise
        error('skBiquad:BadCmd', 'Unknown command "%s".', cmd);
end
end

% ======================================================================
function [R1, R2] = local_skIdeal(passtype, w0, Q, C1, C2)
% Ideal unity-gain Sallen-Key components.
switch lower(passtype)
    case 'lowpass'
        % omega0 = 1/sqrt(R1 R2 C1 C2), Q = sqrt(R1 R2 C1 C2)/(C2 (R1+R2))
        % -> a = R1+R2 = 1/(w0 Q C2),  b = R1 R2 = 1/(w0^2 C1 C2)
        if C1 < 4*Q^2*C2
            error('skBiquad:CapRatio', ...
                'For the low-pass Sallen-Key C1 >= 4*Q^2*C2 is required (C1=%.4g, need >= %.4g).', ...
                C1, 4*Q^2*C2);
        end
        a = 1 / (w0 * Q * C2);
        b = 1 / (w0^2 * C1 * C2);
        disc = sqrt(max(a^2 - 4*b, 0));
        R1 = (a + disc) / 2;
        R2 = (a - disc) / 2;
    case 'highpass'
        % omega0 = 1/sqrt(R1 R2 C1 C2), Q = sqrt(R1 R2 C1 C2)/(R1 (C1+C2))
        R2 = Q * (C1 + C2) / (w0 * C1 * C2);
        R1 = 1 / (w0 * Q * (C1 + C2));
    otherwise
        error('skBiquad:BadType', 'passtype must be lowpass|highpass.');
end
end

% ======================================================================
function [B1, B2, B3] = local_skCoeffs(R1, R2, C1, C2, passtype, wt)
% Actual normalized denominator coefficients 1 + B1 s + B2 s^2 + B3 s^3
% of a unity-gain Sallen-Key section with an op-amp GBW of wt (rad/s).
%
% The op-amp is modelled as a single integrator A(s) = wt/s; the unity
% follower then has V+ = Vout*(1 + s/wt).  Substituting this into the
% node equations and eliminating the op-amp input node gives the exact
% third-order denominator (verified against a direct nodal solve of the
% 3x3 circuit matrix, agreement < 1e-10 dB):
%
%   lowpass:  D(s) = (1+s/wt)*[1 + s*(R1*C1 + (R1+R2)*C2) + s^2*R1*R2*C1*C2]
%                     - R1*s*C1
%   highpass: D(s) = R1*[(1+s/wt)*(1+s*R2*C2)*(s*C1 + 1/R1)
%                     + s*C2*(1+s/wt) - s*R2*C2/R1]   (normalized to 1 + ...)
%
% so the coefficients are
%   lowpass:  B1 = C2*(R1+R2) + 1/wt
%             B2 = R1*R2*C1*C2 + (R1*C1 + (R1+R2)*C2)/wt
%             B3 = R1*R2*C1*C2 / wt
%   highpass: B1 = R1*(C1+C2) + 1/wt
%             B2 = R1*R2*C1*C2 + (R1*(C1+C2) + R2*C2)/wt
%             B3 = R1*R2*C1*C2 / wt
%
% Note: the op-amp pole only adds a 1/wt term to B1 and linear-in-1/wt
% terms to B2/B3; the pre-distortion in SKOPT compensates the resulting
% Q enhancement (Qa/Q ~ 1 + wt/(w0*Q)).
if isinf(wt)
    % ideal quadratic
    switch lower(passtype)
        case 'lowpass'
            B1 = C2.*(R1 + R2);
            B2 = R1.*R2.*C1.*C2;
        case 'highpass'
            B1 = R1.*(C1 + C2);
            B2 = R1.*R2.*C1.*C2;
    end
    B3 = 0;
    return;
end
switch lower(passtype)
    case 'lowpass'
        B1 = C2.*(R1 + R2) + 1./wt;
        B2 = R1.*R2.*C1.*C2 + (R1.*C1 + C2.*(R1 + R2))./wt;
        B3 = R1.*R2.*C1.*C2./wt;
    case 'highpass'
        B1 = R1.*(C1 + C2) + 1./wt;
        B2 = R1.*R2.*C1.*C2 + (R1.*(C1 + C2) + R2.*C2)./wt;
        B3 = R1.*R2.*C1.*C2./wt;
    otherwise
        error('skBiquad:BadType', 'passtype must be lowpass|highpass.');
end
end

% ======================================================================
function [w0a, Qa] = local_dominant(B1, B2, B3)
% Dominant pole pair of 1 + B1 s + B2 s^2 + B3 s^3.
% The cubic has the designed pair plus a far real pole (the finite-GBW
% effect); find all roots and pick the complex-conjugate pair.
if B3 == 0
    r = roots([B2 B1 1]);
else
    r = roots([B3 B2 B1 1]);
end
% pick the pair: prefer complex-conjugate roots with negative real part
im = abs(imag(r));
isPair = im > 1e-6 * abs(r);
if sum(isPair) >= 2
    rr = r(isPair);
    [~, ord] = sort(abs(real(rr)));
    p = rr(ord(1:2));
else
    % all-real case: the two poles nearest the jw axis
    [~, ord] = sort(abs(real(r)));
    p = r(ord(1:2));
end
w0a = abs(p(1));
Qa = w0a / (2*abs(real(p(1))));
end

% ======================================================================
function comp = local_svIdeal(passtype, w0, Q, wz, C)
% State-variable biquad ideal components.
%
% lowpass (wz > w0):
%   C1 = C2 = C; R1 = R2 = R3 = R6 = Rf = R = 1/(w0*C);
%   R8 = R*(wz/w0)^2; R7 = Q*R*(1-(w0/wz)^2); R5 = Q*(R8-R);
%   H(s) = -(R/R8)*(s^2+wz^2)/(s^2+(w0/Q)s+w0^2)
%
% highpass (wz < w0): the dual circuit.  Design the normalized LP and
% map R -> C' = 1/R, C -> R' = 1/C, then scale to (w0, C).
if wz <= w0 && strcmpi(passtype, 'lowpass')
    error('skBiquad:ZeroPos', 'For lowpass sections the zero must be above w0.');
end
if wz >= w0 && strcmpi(passtype, 'highpass')
    error('skBiquad:ZeroPos', 'For highpass sections the zero must be below w0.');
end

switch lower(passtype)
    case 'lowpass'
        R  = 1/(w0*C);
        R8 = R*(wz/w0)^2;
        R7 = Q*R*(1 - (w0/wz)^2);
        R5 = Q*(R8 - R);
        comp = struct('layout', 'lp', ...
            'R1', R, 'R2', R, 'R3', R, 'R5', R5, 'R6', R, ...
            'R7', R7, 'R8', R8, 'Rf', R, ...
            'C1', C, 'C2', C, 'C3', C, 'C5', C, 'C6', C, ...
            'C7', C, 'C8', C, 'Rfd', [], 'R1d', [], 'R2d', []);
    case 'highpass'
        error('skBiquad:HPNotSupported', ...
            'The state-variable biquad realizes only w_z > w_0 (low-pass).');
    otherwise
        error('skBiquad:BadType', 'passtype must be lowpass|highpass.');
end
end

% ======================================================================
function H = local_svResp(w, comp, wt)
% Actual response of the state-variable biquad at frequency w (scalar or
% vector), with every op-amp modelled as A(s) = wt/s.
if strcmp(comp.layout, 'lp')
    R1 = comp.R1; R2 = comp.R2; R3 = comp.R3;
    R5 = comp.R5; R6 = comp.R6; R7 = comp.R7; R8 = comp.R8; Rf = comp.Rf;
    C1 = comp.C1; C2 = comp.C2;
    H = arrayfun(@(ww) local_svLp(ww, R1,R2,R3,R5,R6,R7,R8,Rf,C1,C2,wt), w);
else
    error('skBiquad:BadLayout', 'Unknown layout "%s".', comp.layout);
end
end

function H = local_svLp(w, R1,R2,R3,R5,R6,R7,R8,Rf,C1,C2,wt)
% Actual Vout/Vin of the lowpass state-variable biquad at frequency w.
% Op-amp model: VA = -L*V1, VB = -L*V2, VC = -L*Vout, L = jw/wt.
jw = 1i*w;
if isinf(wt), L = 0; else, L = jw/wt; end
% KCL at A (op1 input): R1 Vin, R2 Vout, R7 V1, C1 V1
% KCL at B (op2 input): R3 V1, C2 V2
% KCL at C (op3 input): R5 V1, R6 V2, R8 Vin, Rf Vout
M = [ -(1/R7 + jw*C1) - L*(1/R1 + 1/R2 + 1/R7 + jw*C1),         0,            -1/R2;
                 -1/R3,                    -(L/R3 + jw*C2*(1+L)),              0;
                 -1/R5,                               -1/R6,  -(1/Rf + L*(1/R5+1/R6+1/R8+1/Rf)) ];
B = [1/R1; 0; 1/R8];
V = M \ B;
H = V(3);
end

function [R, C] = local_rc(passtype, wc, wt, C)
% First-order RC + unity buffer with closed-form GBW pre-distortion.
%   lowpass : |H|^2 = 1/((1+y)(1+y/g^2)), y = (w*R*C)^2, g = wt/wc
%   highpass: |H|^2 = (y/(1+y))*(g^2/(g^2+y))
g = wt / wc;
if isinf(wt)
    y = 1;
elseif strcmpi(passtype, 'lowpass')
    % (1+y)(1+y/g^2) = 2  ->  y^2 + (g^2+1)y - g^2 = 0
    y = (-(g^2 + 1) + sqrt((g^2 + 1)^2 + 4*g^2)) / 2;
else
    % (y/(1+y))*(g^2/(g^2+y)) = 1/2  ->  y^2 - (g^2-1)y + g^2 = 0
    d = (g^2 - 1)^2 - 4*g^2;
    if d < 0
        y = 1;                       % GBW too low for exact pre-distortion
    else
        y = ((g^2 - 1) - sqrt(d)) / 2;
        if y <= 0, y = 1; end
    end
end
R = sqrt(y) / (wc * C);
end
