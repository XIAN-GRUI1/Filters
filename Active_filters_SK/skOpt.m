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
%   With a finite .wt the search can optionally perform the GBW pre-
%   distortion: the components are chosen so that the *actual* response
%   (computed with the single-pole op-amp model) matches the ideal one.
%   Per section this is only done when the section actually needs it
%   (see .pred / .pthresh below): the finite op-amp bandwidth perturbs a
%   unity-gain Sallen-Key section by a Q enhancement and a pole shift of
%   about 1/m, m = wt/(w0*Q) (verified by exact nodal analysis).  With
%   m >= 100 the uncompensated Q error is below ~1 % and the ideal
%   components are kept (pre-distortion would only make the design
%   sensitive to the op-amp GBW and to E-series rounding).
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
if ~isfield(opt, 'pred'),  opt.pred = 'auto';   end   % 'auto'|'always'|'never'
if ~isfield(opt, 'pthresh'), opt.pthresh = 100; end   % margin below which pre-distortion runs
if ~isfield(opt, 'rtol'),  opt.rtol = 0.01;     end   % realized-error threshold (relative)
if ~isfield(opt, 'x0'),    opt.x0 = [];         end   % GBW-optimal C1/C2 ratio
if ~isfield(opt, 'wr'),    opt.wr = 0;          end   % ratio-preference weight

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
% ---- conditional GBW pre-distortion ----
% The finite op-amp bandwidth perturbs the section pole pair by ~1/m
% (m = wt/(w0*Q), verified by exact nodal analysis), but the exact amount
% also depends on the chosen C1/C2 ratio: far from the GBW-optimal split
% (equal R's, i.e. C1/C2 = 4*Q^2 for LP; equal C's for HP) a section can
% become orders of magnitude more sensitive to the op-amp GBW.  Therefore
% with a finite GBW the search keeps the ratio near the optimum, and
% pre-distortion is applied ONLY to the sections that measurably need it:
%   1. fast rule (empirical formula): if m = wt/(w0*Q) < pthresh (100),
%      the section is clearly distorted -> pre-distort directly;
%   2. otherwise search the ideal (component-only) solution first, then
%      MEASURE the realized (w0a, Qa) with the real GBW: if they deviate
%      from the target by more than rtol (1 %) the section is
%      pre-distorted after all.
wt = opt.wt;
wq = opt.wQ;
margin = wt/(w0*Q);
if isinf(wt), margin = Inf; end

% ratio preference: GBW-optimal component split (only for GBW-aware
% designs; 'never' keeps the pure ideal component search)
if isfinite(wt) && ~strcmpi(opt.pred, 'never')
    if strcmpi(passtype, 'lowpass')
        opt.x0 = 4*Q^2;          % equal-R: C1/C2 = 4*Q^2
    else
        opt.x0 = 1;              % equal-C (dual): C1 = C2
    end
    opt.wr = 0.1;
else
    opt.x0 = []; opt.wr = 0;
end

pred = local_needPred(opt, wt, w0, Q, margin);
if pred
    [best, bestErr] = local_searchSk(passtype, w0, Q, opt, wt);
else
    [best, bestErr] = local_searchSk(passtype, w0, Q, opt, Inf);
    if ~isempty(best) && isfinite(wt)
        % measure the REALIZED section with the actual op-amp GBW
        [~, wa, qa] = local_errSk(passtype, best.R1, best.R2, best.C1, best.C2, ...
                                  w0, Q, wt, wq, opt.x0, opt.wr);
        if strcmpi(opt.pred, 'auto') && ...
                max(abs(wa-w0)/w0, abs(qa-Q)/Q) > opt.rtol
            % the ideal parts are distorted by the finite GBW: compensate
            [best, bestErr] = local_searchSk(passtype, w0, Q, opt, wt);
            pred = true;
        else
            best.w0a = wa; best.Qa = qa;   % realized values for the report
        end
    end
end

sec = best;
sec.kind = 'sk';
sec.w0t = w0; sec.Qt = Q;
sec.err = bestErr;
sec.wt = wt;
sec.pred = pred;
sec.margin = margin;
sec.topology = sprintf('Sallen-Key (unity gain) %s', passtype);
end

% ======================================================================
function pred = local_needPred(opt, wt, w0, Q, margin)
% fast rule deciding whether a section needs GBW pre-distortion
if isinf(wt)
    pred = false;
else
    switch lower(opt.pred)
        case 'always', pred = true;
        case 'never',  pred = false;
        otherwise,     pred = margin < opt.pthresh;
    end
end
end

% ======================================================================
function [R1, R2, ok] = local_secR(passtype, w0, Q, wt, C1, C2, Rmin, Rmax)
% Ideal or pre-distorted resistor pair for the given caps.
%   wt = Inf : classical ideal design values realizing (w0, Q);
%   wt finite: exact GBW pre-distortion - the actual dominant pole pair
%              of the section (op-amp A(s) = wt/s) is forced to (w0, Q)
%              by matching the cubic denominator to
%              (1+s/wf)*(1 + s/(Q*w0) + s^2/w0^2).
if isinf(wt)
    switch lower(passtype)
        case 'lowpass'
            a = 1/(w0*Q*C2);
            b = 1/(w0^2*C1*C2);
            disc = sqrt(max(a^2 - 4*b, 0));
            R1 = (a + disc)/2;  R2 = (a - disc)/2;
        case 'highpass'
            R2 = Q*(C1 + C2)/(w0*C1*C2);
            R1 = 1/(w0*Q*(C1 + C2));
    end
    ok = R1 >= Rmin && R1 <= Rmax && R2 >= Rmin && R2 <= Rmax;
    return;
end
ok = false; R1 = []; R2 = [];
Pmin = max(Rmin^2, 1e-20); Pmax = Rmax^2;
switch lower(passtype)
    case 'lowpass'
        % B1: C2*(R1+R2) + 1/wt = 1/(Q*w0) + 1/wf,   B3: 1/wf = w0^2*P*C1*C2/wt
        % -> S = R1+R2 = S0 + k*P,  P = R1*R2
        S0 = (1/(Q*w0) - 1/wt)/C2;
        k  = w0^2*C1/wt;
        g  = @(P) local_gLP(P, C1, C2, w0, Q, wt, S0, k);
        Rfun = @(P) local_RLP(P, S0, k);
    case 'highpass'
        % B1: R1*(C1+C2) + 1/wt = 1/(Q*w0) + 1/wf
        % -> R1 = (1/(Q*w0) - 1/wt + w0^2*P*C1*C2/wt)/(C1+C2),  R2 = P/R1
        g  = @(P) local_gHP(P, C1, C2, w0, Q, wt);
        Rfun = @(P) local_RHP(P, C1, C2, w0, Q, wt);
    otherwise
        return;
end
% scan for sign changes of g over P (log grid) and refine each with fzero
Pgrid = logspace(log10(Pmin), log10(Pmax), 60);
gv = arrayfun(g, Pgrid);
rts = [];
for i = 1:numel(Pgrid)-1
    if isfinite(gv(i)) && isfinite(gv(i+1)) && (gv(i)*gv(i+1) < 0 || gv(i) == 0)
        try
            r = fzero(g, [Pgrid(i), Pgrid(i+1)]);
        catch %#ok<NASGU>
            continue;
        end
        rts(end+1) = r; %#ok<AGROW>
    end
end
if isempty(rts), return; end
% accept the first root whose resistors land inside [Rmin, Rmax]
Pideal = 1/(w0^2*C1*C2);
[~, ord] = sort(abs(rts - Pideal));
for j = ord
    [R1c, R2c] = Rfun(rts(j));
    if isfinite(R1c) && isfinite(R2c) && R1c > 0 && R2c > 0 && ...
            R1c >= Rmin && R1c <= Rmax && R2c >= Rmin && R2c <= Rmax
        R1 = R1c; R2 = R2c; ok = true;
        return;
    end
end
end

function g = local_gLP(P, C1, C2, w0, Q, wt, S0, k)
% residual of the B2 match:  P*C1*C2 + (R1*C1 + S*C2)/wt = 1/w0^2 + (w0/Q)*P*C1*C2/wt
S  = S0 + k*P;
R1 = (S + sqrt(max(S^2 - 4*P, 0)))/2;
g  = P*C1*C2 + (R1*C1 + S*C2)/wt - 1/w0^2 - (w0/Q)*P*C1*C2/wt;
end

function [R1, R2] = local_RLP(P, S0, k)
S = S0 + k*P;
R1 = (S + sqrt(max(S^2 - 4*P, 0)))/2;
R2 = (S - sqrt(max(S^2 - 4*P, 0)))/2;
end

function g = local_gHP(P, C1, C2, w0, Q, wt)
R1 = (1/(Q*w0) - 1/wt + w0^2*P*C1*C2/wt)/(C1 + C2);
R2 = P/R1;
g  = P*C1*C2 + (R1*(C1 + C2) + R2*C2)/wt - 1/w0^2 - (w0/Q)*P*C1*C2/wt;
end

function [R1, R2] = local_RHP(P, C1, C2, w0, Q, wt)
R1 = (1/(Q*w0) - 1/wt + w0^2*P*C1*C2/wt)/(C1 + C2);
R2 = P/R1;
end

% ======================================================================
function [best, bestErr] = local_searchSk(passtype, w0, Q, opt, wtEff)
% full E-series search of one Sallen-Key section; wtEff = Inf searches
% the ideal (component-only) solution, a finite wtEff searches the
% pre-distorted (actual-response) solution.
Cg = skSeriesGrid(opt.Cser, opt.Cmin, opt.Cmax);
wq = opt.wQ;
opt.wt = wtEff;      % downstream helpers (refineCands) use the same wt
bestErr = inf; best = [];
switch lower(passtype)
    case 'lowpass'
        for C2 = Cg
            C1min = max(4*Q^2*C2, opt.Cmin);
            C1 = Cg(Cg >= C1min - 1e-15);
            if isempty(C1), continue; end
            if isinf(wtEff)
                % ideal design: vectorized over C1
                a = 1/(w0*Q*C2);
                b = 1./(w0^2*C1*C2);
                disc = sqrt(max(a^2 - 4*b, 0));
                R1x = (a + disc)/2;  R2x = (a - disc)/2;
                ok = R1x >= opt.Rmin & R1x <= opt.Rmax & R2x >= opt.Rmin & R2x <= opt.Rmax;
                if ~any(ok), continue; end
                [R1r, R2r] = local_round(R1x, R2x, ok, opt.Rser, 0);
                [e, wa, qa] = local_errSk('lowpass', R1r, R2r, C1, C2, w0, Q, wtEff, wq, ...
                                          opt.x0, opt.wr);
                [em, im] = min(e);
                if em < bestErr
                    bestErr = em;
                    best = struct('C1', C1(im), 'C2', C2, 'R1', R1r(im), 'R2', R2r(im), ...
                                  'R1x', R1x(im), 'R2x', R2x(im), 'w0a', wa(im), 'Qa', qa(im));
                end
            else
                % pre-distortion: exact pole compensation per (C1,C2)
                for ii = 1:numel(C1)
                    [R1x, R2x, okk] = local_secR('lowpass', w0, Q, wtEff, C1(ii), C2, ...
                                                  opt.Rmin, opt.Rmax);
                    if ~okk, continue; end
                    [R1r, R2r] = local_round(R1x, R2x, true, opt.Rser, 0);
                    [e, wa, qa] = local_errSk('lowpass', R1r, R2r, C1(ii), C2, w0, Q, wtEff, wq, ...
                                              opt.x0, opt.wr);
                    if e < bestErr
                        bestErr = e;
                        best = struct('C1', C1(ii), 'C2', C2, 'R1', R1r, 'R2', R2r, ...
                                      'R1x', R1x, 'R2x', R2x, 'w0a', wa, 'Qa', qa);
                    end
                end
            end
        end
    case 'highpass'
        for C1 = Cg
            for C2 = Cg
                [R1x, R2x, okk] = local_secR('highpass', w0, Q, wtEff, C1, C2, ...
                                              opt.Rmin, opt.Rmax);
                if ~okk, continue; end
                [R1r, R2r] = local_round(R1x, R2x, true, opt.Rser, 0);
                [e, wa, qa] = local_errSk('highpass', R1r, R2r, C1, C2, w0, Q, wtEff, wq, ...
                                          opt.x0, opt.wr);
                if e < bestErr
                    bestErr = e;
                    best = struct('C1', C1, 'C2', C2, 'R1', R1r, 'R2', R2r, ...
                                  'R1x', R1x, 'R2x', R2x, 'w0a', wa, 'Qa', qa);
                end
            end
        end
end

% ---- pass 2: refine the best candidates with +-1 E-series steps ----
if ~isempty(best)
    cands = local_refineCands(passtype, w0, Q, opt, Cg, best, 2*opt.nref);
    for i = 1:numel(cands)
        c = cands(i);
        for d1 = -1:1
            for d2 = -1:1
                [R1r, R2r] = local_round(c.R1x, c.R2x, true, opt.Rser, [d1 d2]);
                [e, wa, qa] = local_errSk(passtype, R1r, R2r, c.C1, c.C2, w0, Q, wtEff, wq, ...
                                          opt.x0, opt.wr);
                if e < bestErr
                    bestErr = e;
                    best = struct('C1', c.C1, 'C2', c.C2, 'R1', R1r, 'R2', R2r, ...
                                  'R1x', c.R1x, 'R2x', c.R2x, 'w0a', wa, 'Qa', qa);
                end
            end
        end
    end
end
end

% ======================================================================
function [e, wa, qa] = local_errSk(passtype, R1, R2, C1, C2, w0, Q, wt, wq, x0, wr)
% error of a set of (possibly vector) components; returns vectors.
% With a finite wt this is the pre-distortion criterion (match the ACTUAL
% dominant pole pair to the target); wt = Inf gives the ideal criterion.
% x0/wr (optional) add a soft preference for the GBW-optimal C1/C2 ratio.
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
if ~isempty(x0) && wr > 0
    x = C1./C2;
    e(ok) = e(ok) + wr*log(x(ok)/x0).^2;
end
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
% conditional pre-distortion (Q = 1 for a first-order section):
% search the ideal solution first, then MEASURE the realized -3 dB
% frequency with the actual op-amp GBW; pre-distort only if it is off by
% more than rtol (the fast rule m = wt/wc < pthresh skips straight to
% pre-distortion).
margin = wt/wc;
if isinf(wt), margin = Inf; end
pred = local_needPred(opt, wt, wc, 1, margin);
wts = wt; if ~pred, wts = Inf; end
bestErr = inf; best = [];
for C = Cg
    Rx = skBiquad('rc', passtype, wc, wts, C);
    if Rx < opt.Rmin || Rx > opt.Rmax, continue; end
    R = skNearestE(opt.Rser, Rx);
    for d = -1:1
        Rr = local_stepE(opt.Rser, R, d);
        if Rr < opt.Rmin || Rr > opt.Rmax, continue; end
        % search criterion uses the effective wt, reporting uses the
        % actual -3 dB frequency of the RC + buffer pole
        w3s = local_w3(passtype, Rr, C, wts);
        e = log(w3s/wc)^2;
        w3 = local_w3(passtype, Rr, C, wt);
        if e < bestErr
            bestErr = e;
            best = struct('C', C, 'R', Rr, 'Rx', Rx, 'w3', w3);
        end
    end
end
if ~pred && isfinite(wt) && ~isempty(best) && strcmpi(opt.pred, 'auto')
    % measure the realized -3 dB point with the real op-amp GBW
    w3a = local_w3(passtype, best.R, best.C, wt);
    if abs(w3a - wc)/wc > opt.rtol
        wts = wt; pred = true;
        bestErr = inf; best = [];
        for C = Cg
            Rx = skBiquad('rc', passtype, wc, wts, C);
            if Rx < opt.Rmin || Rx > opt.Rmax, continue; end
            R = skNearestE(opt.Rser, Rx);
            for d = -1:1
                Rr = local_stepE(opt.Rser, R, d);
                if Rr < opt.Rmin || Rr > opt.Rmax, continue; end
                w3s = local_w3(passtype, Rr, C, wts);
                e = log(w3s/wc)^2;
                w3 = local_w3(passtype, Rr, C, wt);
                if e < bestErr
                    bestErr = e;
                    best = struct('C', C, 'R', Rr, 'Rx', Rx, 'w3', w3);
                end
            end
        end
    end
end
sec = best;
sec.kind = 'rc';
sec.w0t = wc; sec.Qt = 0.5;
sec.err = bestErr;
sec.wt = wt;
sec.pred = pred;
sec.margin = margin;
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
            if isinf(opt.wt)
                a = 1/(w0*Q*C2);
                b = 1./(w0^2*C1*C2);
                disc = sqrt(max(a^2 - 4*b, 0));
                R1x = (a + disc)/2;  R2x = (a - disc)/2;
                ok = R1x >= opt.Rmin & R1x <= opt.Rmax & R2x >= opt.Rmin & R2x <= opt.Rmax;
                if ~any(ok), continue; end
                e = local_errSk('lowpass', skNearestE(opt.Rser, R1x), skNearestE(opt.Rser, R2x), ...
                                C1, C2, w0, Q, opt.wt, opt.wQ, opt.x0, opt.wr);
                [em, im] = min(e);
                list(end+1, :) = [C2, C1(im), R1x(im), R2x(im), em]; %#ok<AGROW>
            else
                for ii = 1:numel(C1)
                    [R1x, R2x, okk] = local_secR('lowpass', w0, Q, opt.wt, C1(ii), C2, ...
                                                  opt.Rmin, opt.Rmax);
                    if ~okk, continue; end
                    e = local_errSk('lowpass', skNearestE(opt.Rser, R1x), skNearestE(opt.Rser, R2x), ...
                                    C1(ii), C2, w0, Q, opt.wt, opt.wQ, opt.x0, opt.wr);
                    list(end+1, :) = [C2, C1(ii), R1x, R2x, e]; %#ok<AGROW>
                end
            end
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
                [R1x, R2x, okk] = local_secR('highpass', w0, Q, opt.wt, C1, C2, ...
                                              opt.Rmin, opt.Rmax);
                if ~okk, continue; end
                e = local_errSk('highpass', skNearestE(opt.Rser, R1x), skNearestE(opt.Rser, R2x), ...
                                C1, C2, w0, Q, opt.wt, opt.wQ, opt.x0, opt.wr);
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
