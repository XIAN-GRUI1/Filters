function S = skSpecs(F)
%SKSPECS  Measure the achieved filter specifications.
%   S = SKSPECS(F)  returns a struct with
%     .f3dB     -3 dB frequency (Hz) for butter / bessel
%     .ripple   passband ripple (dB) for cheby1 / cauer
%     .astop    stopband attenuation (dB) at the stopband edge
%     .fstop    stopband edge frequency (Hz)
%     .maxerr   max |H| deviation from the ideal prototype (dB)
%
%   See also SKDESIGN.

f = F.fc;
w = 2*pi*logspace(log10(f*1e-4), log10(f*1e4), 8001);
[H, Hi] = skResponse(F, w);
HdB = 20*log10(max(abs(H), 1e-12));
HidB = 20*log10(max(abs(Hi), 1e-12));

S = struct();
S.f3dB = nan;
S.ripple = nan;
S.astop = nan;
S.fstop = nan;
% deviation metric restricted to where both responses are meaningful
% and the single-pole op-amp model is valid (well below the GBW)
idx = HidB > -60 & HdB > -60;
if isfinite(F.gbw)
    idx = idx & (w < F.gbw/10);
end
S.maxerr = max(abs(HdB(idx) - HidB(idx)));
if isempty(S.maxerr), S.maxerr = 0; end

switch lower(F.type)
    case {'butter', 'bessel'}
        % -3 dB crossing
        if strcmpi(F.passtype, 'lowpass')
            i = find(HdB <= -3 + 1e-9, 1, 'first');
        else
            i = find(HdB >= -3 - 1e-9, 1, 'first');
        end
        if ~isempty(i)
            S.f3dB = w(i) / (2*pi);
        end
    case {'cheby1', 'cauer'}
        if strcmpi(F.passtype, 'lowpass')
            pb = w <= 2*pi*f * 1.0001;
        else
            pb = w >= 2*pi*f * 0.9999;
        end
        % for high-pass the passband extends upward without bound; with a
        % finite op-amp GBW the response rolls off above GBW, so keep the
        % ripple measurement in the region where the op-amp model is valid
        if isfinite(F.gbw)
            pb = pb & (w < F.gbw/10);
        end
        if any(pb)
            S.ripple = max(HdB(pb)) - min(HdB(pb));
        end
    otherwise
end

if isfinite(F.prototype.wstop)
    ws = F.prototype.wstop;
    if strcmpi(F.passtype, 'lowpass')
        wstop = ws * 2*pi*f;
    else
        wstop = 2*pi*f / ws;
    end
    [~, i] = min(abs(w - wstop));
    S.astop = -HdB(i);
    S.fstop = w(i) / (2*pi);
end
end
