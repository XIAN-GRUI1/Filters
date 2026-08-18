function sec = skOpt(kind, passtype, w0, Q, wz, opt)
%SKOPT  E-series component search for one filter section.
%   SEC = SKOPT(KIND, PASSTYPE, W0, Q, WZ, OPT)
%
%   KIND      'sk' : unity-gain Sallen-Key biquad (no zeros)
%             'sv' : state-variable biquad (needs WZ, transmission zero)
%             'rc' : first-order RC + buffer (W0 is the cutoff)
%   PASSTYPE  'lowpass' | 'highpass'
%   OPT       struct with fields:
%             .Rser, .Cser : 'E96'/'E24' etc.
%             .Rmin, .Rmax, .Cmin, .Cmax : component search ranges
%             .wt          : op-amp GBW in rad/s (Inf = ideal)
%             .wQ          : relative weight of the Q error (default 3)
%             .nref        : number of candidates refined (default 25)
%
%   With a finite .wt the search implicitly performs the GBW pre-
%   distortion: the components are chosen so that the *actual* response
%   (computed with the single-pole op-amp model) matches the ideal one.
%
%   SEC is a struct with the rounded components, the achieved (w0,Q)
%   and the relative error.
%
%   See also SKDESIGN, SKBIQUAD.

% ---- defaults ----
if ~isfield(opt, 'wQ'),    opt.wQ = 3;    end
if ~isfield(opt, 'nref'),  opt.nref = 25; end
if ~isfield(opt, 'wt'),    opt.wt = Inf;  end
if ~isfield(opt, 'Rmin'),  opt.Rmin = 100;  end
if ~isfield(opt, 'Rmax'),  opt.Rmax = 1e6;  end
if ~isfield(opt, 'Cmin'),  opt.Cmin = 1e-12; end
if ~isfield(opt, 'Cmax'),  opt.Cmax = 1e-6;  end
if ~isfield(opt, 'Cser'),  opt.Cser = 'E24'; end
if ~isfield(opt, 'Rser'),  opt.Rser = 'E96'; end

switch kind
    case 'sk'
        sec = local_sk(passtype, w0, Q, opt);
    case 'sv'
        sec = local_sv(passtype, w0, Q, wz, opt);
    case 'rc'
        sec = local_rc(w0, passtype, opt);
    otherwise
        error('skOpt:BadKind', 'Unknown kind "%s".', kind);
end
end

% ======================================================================
function sec = local_sk(passtype, w0, Q, opt)
% ---- grids ----
Cg = skSeriesGrid(opt.Cser, opt.Cmin, opt.Cmax);
wt = opt.wt;
wq = opt.wQ;

% ---- pass 1: nearest rounding ----
bestErr = inf; best = [];
switch lower(passtype)
    case 'lowpass'
        C2list = Cg;
        for C2 = C2list
            C1min = max(4*Q^2*C2, opt.Cmin);
            C1 = Cg(Cg >= C1min - 1e-15);
            if isempty(C1), continue; end
            a = 1/(w0*Q*C2);
            b = 1./(w0^2*C1*C2);
            disc = sqrt(max(a^2 - 4*b, 0));
            R1x = (a + disc)/2;  R2x = (a - disc)/2;
            ok = R1x >= opt.Rmin & R1x <= opt.Rmax & R2x >= opt.Rmin & R2x <= opt.Rmax;
            if ~any(ok), continue; end
            [R1r, R2r] = local_round(R1x, R2x, ok, opt.Rser, 0);
            [e, wa, qa] = local_errSk('lowpass', R1r, R2r, C1, C2, w0, Q, wt, wq);
            [em, im] = min(e);
            if em < bestErr
                bestErr = em;
                best = struct('C1', C1(im), 'C2', C2, 'R1', R1r(im), 'R2', R2r(im), ...
                              'R1x', R1x(im), 'R2x', R2x(im), 'w0a', wa(im), 'Qa', qa(im));
            end
        end
    case 'highpass'
        for C1 = Cg
            for C2 = Cg
                R2x = Q*(C1 + C2)/(w0*C1*C2);
                R1x = 1/(w0*Q*(C1 + C2));
                if R1x < opt.Rmin || R1x > opt.Rmax || R2x < opt.Rmin || R2x > opt.Rmax
                    continue;
                end
                [R1r, R2r] = local_round(R1x, R2x, true, opt.Rser, 0);
                [e, wa, qa] = local_errSk('highpass', R1r, R2r, C1, C2, w0, Q, wt, wq);
                if e < bestErr
                    bestErr = e;
                    best = struct('C1', C1, 'C2', C2, 'R1', R1r, 'R2', R2r, ...
                                  'R1x', R1x, 'R2x', R2x, 'w0a', wa, 'Qa', qa);
                end
            end
        end
end

% ---- pass 2: refine the best candidates with +-1 E96 steps ----
if ~isempty(best)
    cands = local_refineCands(passtype, w0, Q, opt, Cg, best, 2*opt.nref);
    for i = 1:numel(cands)
        c = cands(i);
        for d1 = -1:1
            for d2 = -1:1
                [R1r, R2r] = local_round(c.R1x, c.R2x, true, opt.Rser, [d1 d2]);
                [e, wa, qa] = local_errSk(passtype, R1r, R2r, c.C1, c.C2, w0, Q, wt, wq);
                if e < bestErr
                    bestErr = e;
                    best = struct('C1', c.C1, 'C2', c.C2, 'R1', R1r, 'R2', R2r, ...
                                  'R1x', c.R1x, 'R2x', c.R2x, 'w0a', wa, 'Qa', qa);
                end
            end
        end
    end
end

sec = best;
sec.kind = 'sk';
sec.w0t = w0; sec.Qt = Q;
sec.err = bestErr;
sec.wt = wt;
sec.topology = sprintf('Sallen-Key (unity gain) %s', passtype);
end

% ======================================================================
function [e, wa, qa] = local_errSk(passtype, R1, R2, C1, C2, w0, Q, wt, wq)
% error of a set of (possibly vector) components; returns vectors.
R1 = R1(:); R2 = R2(:); C1 = C1(:); C2 = C2(:);
[B1, B2, B3] = skBiquad('skcoeffs', R1, R2, C1, C2, passtype, wt);
n = numel(B1);
B1 = B1(:); B2 = B2(:);
if isscalar(B3), B3 = B3*ones(n,1); else, B3 = B3(:); end
wa = nan(n,1); qa = nan(n,1);
ok = isfinite(B1) & isfinite(B2) & isfinite(B3) & (B1 > 0) & (B2 > 0);
for i = find(ok)'
    [wa(i), qa(i)] = skBiquad('dominant', B1(i), B2(i), B3(i));
end
e = inf(n,1);
ok = ok & isfinite(wa) & isfinite(qa) & (wa > 0) & (qa > 0);
e(ok) = log(wa(ok)/w0).^2 + wq*log(qa(ok)/Q).^2;
end

% ======================================================================
function sec = local_sv(passtype, w0, Q, wz, opt)
Cg = skSeriesGrid(opt.Cser, opt.Cmin, opt.Cmax);
wt = opt.wt;
% ideal (unrounded) components -> target response
comp0 = skBiquad('svideal', passtype, w0, Q, wz, 1e-9);
% search frequencies around w0 and wz
wf = logspace(-0.7, 0.7, 9) * w0;
if isfinite(wz), wf = sort(unique([wf, wz*[0.5 0.8 1 1.2 2]])); end
Ht = skBiquad('svresp', wf, comp0, Inf);
bestErr = inf; best = [];
for C = Cg
    comp = skBiquad('svideal', passtype, w0, Q, wz, C);
    % round all resistors
    fn = fieldnames(comp);
    compr = comp;
    for i = 1:numel(fn)
        v = comp.(fn{i});
        if ~isempty(v) && isnumeric(v) && ~any(isnan(v)) && ~strcmp(fn{i},'layout')
            compr.(fn{i}) = skNearestE(opt.Rser, v);
        end
    end
    Ha = skBiquad('svresp', wf, compr, wt);
    e = sum(abs(Ha - Ht).^2) / sum(abs(Ht).^2);
    if e < bestErr
        bestErr = e;
        best = struct('comp', compr, 'comp0', comp);
    end
end
% refine: +-1 steps on the resistors of the best few cap choices
Cord = local_svOrder(Cg, passtype, w0, Q, wz, opt, wf, Ht, opt.nref);
for j = 1:numel(Cord)
    C = Cord(j);
    comp = skBiquad('svideal', passtype, w0, Q, wz, C);
    fn = fieldnames(comp);
    Rn = zeros(numel(fn),1); ok = false(numel(fn),1);
    for i = 1:numel(fn)
        if ~isempty(comp.(fn{i})) && isnumeric(comp.(fn{i})) && ~any(isnan(comp.(fn{i}))) && ~strcmp(fn{i},'layout')
            ok(i) = true; Rn(i) = comp.(fn{i});
        end
    end
    % enumerate +-1 on each resistor (2^n too many; do one-at-a-time greedy)
    compr = comp;
    for i = 1:numel(fn)
        if ok(i)
            v0 = skNearestE(opt.Rser, Rn(i));
            bestv = v0; beste = inf;
            for d = -1:1
                compr.(fn{i}) = local_stepE(opt.Rser, v0, d);
                Ha = skBiquad('svresp', wf, compr, wt);
                e = sum(abs(Ha - Ht).^2) / sum(abs(Ht).^2);
                if e < beste, beste = e; bestv = compr.(fn{i}); end
            end
            compr.(fn{i}) = bestv;
        end
    end
    Ha = skBiquad('svresp', wf, compr, wt);
    e = sum(abs(Ha - Ht).^2) / sum(abs(Ht).^2);
    if e < bestErr
        bestErr = e;
        best = struct('comp', compr, 'comp0', comp);
    end
end

sec = best;
sec.kind = 'sv';
sec.w0t = w0; sec.Qt = Q; sec.wzt = wz;
sec.err = bestErr;
sec.wt = wt;
sec.topology = sprintf('State-variable biquad %s (3 op-amps, transmission zero)', passtype);
end

function Cord = local_svOrder(Cg, passtype, w0, Q, wz, opt, wf, Ht, nkeep)
% rank the cap choices by nearest-rounding error
errs = zeros(size(Cg));
for i = 1:numel(Cg)
    comp = skBiquad('svideal', passtype, w0, Q, wz, Cg(i));
    fn = fieldnames(comp);
    compr = comp;
    for j = 1:numel(fn)
        if ~isempty(comp.(fn{j})) && isnumeric(comp.(fn{j})) && ~any(isnan(comp.(fn{j}))) && ~strcmp(fn{j},'layout')
            compr.(fn{j}) = skNearestE(opt.Rser, comp.(fn{j}));
        end
    end
    Ha = skBiquad('svresp', wf, compr, opt.wt);
    errs(i) = sum(abs(Ha - Ht).^2) / sum(abs(Ht).^2);
end
[~, ord] = sort(errs);
Cord = Cg(ord(1:min(nkeep, numel(ord))));
end

% ======================================================================
function sec = local_rc(wc, passtype, opt)
Cg = skSeriesGrid(opt.Cser, opt.Cmin, opt.Cmax);
wt = opt.wt;
bestErr = inf; best = [];
for C = Cg
    Rx = skBiquad('rc', passtype, wc, wt, C);
    if Rx < opt.Rmin || Rx > opt.Rmax, continue; end
    R = skNearestE(opt.Rser, Rx);
    for d = -1:1
        Rr = local_stepE(opt.Rser, R, d);
        if Rr < opt.Rmin || Rr > opt.Rmax, continue; end
        % actual -3 dB frequency of the RC + buffer pole
        w3 = local_w3(passtype, Rr, C, wt);
        e = log(w3/wc)^2;
        if e < bestErr
            bestErr = e;
            best = struct('C', C, 'R', Rr, 'Rx', Rx, 'w3', w3);
        end
    end
end
sec = best;
sec.kind = 'rc';
sec.w0t = wc; sec.Qt = 0.5;
sec.err = bestErr;
sec.wt = wt;
sec.topology = sprintf('First-order RC + buffer (%s)', passtype);
end

function w3 = local_w3(passtype, R, C, wt)
% -3 dB frequency of the RC section including the buffer pole (if any)
if isinf(wt)
    w3 = 1/(R*C);
    return;
end
tau = R*C;
f = @(w) abs(local_H1(passtype, w, tau, wt))^2 - 0.5;
if strcmpi(passtype, 'lowpass')
    w3 = fzero(f, [1e-6/tau, 1e6/tau]);
else
    % highpass: the response rises through -3 dB at ~1/tau and falls
    % again beyond the buffer pole; bracket only the rising edge
    w3 = fzero(f, [1e-6/tau, 4/tau]);
end
end

function H = local_H1(passtype, w, tau, wt)
s = 1i*w;
if strcmpi(passtype, 'lowpass')
    H = (1/(1 + s*tau)) * (1/(1 + s/wt));
else
    H = (s*tau/(1 + s*tau)) * (1/(1 + s/wt));
end
end

% ======================================================================
function [R1, R2] = local_round(R1x, R2x, ok, ser, step)
% round vectors to the E-series; step is an index offset (scalar or 2-vec)
if numel(step) == 1, step = [step step]; end
R1 = R1x; R2 = R2x;
R1(~ok) = NaN; R2(~ok) = NaN;
for i = 1:numel(R1x)
    if ok(i)
        R1(i) = local_stepE(ser, skNearestE(ser, R1x(i)), step(1));
        R2(i) = local_stepE(ser, skNearestE(ser, R2x(i)), step(2));
    end
end
end

function r = local_stepE(ser, v, d)
% move d steps (in index) away from the nominal value v
base = skESeries(ser);
[~, i] = min(abs(base - v/10^floor(log10(v))));
i = i + d;
d0 = 10^floor(log10(v));
if i < 1, i = 1; d0 = d0/10; end
if i > numel(base), i = numel(base); d0 = d0*10; end
r = base(i) * d0;
end

function cands = local_refineCands(passtype, w0, Q, opt, Cg, best, n)
% collect the best (C1,C2) candidates for pass-2 refinement
cands = struct('C1', {}, 'C2', {}, 'R1x', {}, 'R2x', {});
switch lower(passtype)
    case 'lowpass'
        list = [];
        for C2 = Cg
            C1min = max(4*Q^2*C2, opt.Cmin);
            C1 = Cg(Cg >= C1min - 1e-15);
            if isempty(C1), continue; end
            a = 1/(w0*Q*C2);
            b = 1./(w0^2*C1*C2);
            disc = sqrt(max(a^2 - 4*b, 0));
            R1x = (a + disc)/2;  R2x = (a - disc)/2;
            ok = R1x >= opt.Rmin & R1x <= opt.Rmax & R2x >= opt.Rmin & R2x <= opt.Rmax;
            if ~any(ok), continue; end
            e = local_errSk('lowpass', skNearestE(opt.Rser, R1x), skNearestE(opt.Rser, R2x), ...
                            C1, C2, w0, Q, opt.wt, opt.wQ);
            [em, im] = min(e);
            list(end+1, :) = [C2, C1(im), R1x(im), R2x(im), em]; %#ok<AGROW>
        end
        [~, ord] = sort(list(:,5));
        for i = 1:min(n, size(list,1))
            cands(end+1) = struct('C1', list(ord(i),2), 'C2', list(ord(i),1), ...
                                  'R1x', list(ord(i),3), 'R2x', list(ord(i),4)); %#ok<AGROW>
        end
    case 'highpass'
        list = [];
        for C1 = Cg
            for C2 = Cg
                R2x = Q*(C1 + C2)/(w0*C1*C2);
                R1x = 1/(w0*Q*(C1 + C2));
                if R1x < opt.Rmin || R1x > opt.Rmax || R2x < opt.Rmin || R2x > opt.Rmax
                    continue;
                end
                e = local_errSk('highpass', skNearestE(opt.Rser, R1x), skNearestE(opt.Rser, R2x), ...
                                C1, C2, w0, Q, opt.wt, opt.wQ);
                list(end+1, :) = [C1, C2, R1x, R2x, e]; %#ok<AGROW>
            end
        end
        [~, ord] = sort(list(:,5));
        for i = 1:min(n, size(list,1))
            cands(end+1) = struct('C1', list(ord(i),1), 'C2', list(ord(i),2), ...
                                  'R1x', list(ord(i),3), 'R2x', list(ord(i),4)); %#ok<AGROW>
        end
end
end
