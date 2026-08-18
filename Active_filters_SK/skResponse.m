function [H, Hideal] = skResponse(F, w)
%SKRESPONSE  Frequency response of a designed filter.
%   [H, HIDEAL] = SKRESPONSE(F, W)
%     F      : filter struct from SKDESIGN
%     W      : frequency vector (rad/s)
%     H      : actual response (with the finite-GBW op-amp model)
%     HIDEAL : ideal prototype response (infinite GBW, exact components)
%
%   See also SKDESIGN.

w = w(:).';
n = numel(w);
H = ones(1, n);
for i = 1:numel(F.sections)
    s = F.sections(i);
    switch s.kind
        case 'rc'
            tau = s.R * s.C;
            wt = F.gbw;
            if isinf(wt)
                H1 = 1 ./ (1 + 1i*w*tau);
                if strcmpi(F.passtype, 'highpass'), H1 = 1 - H1; end
            else
                H1 = (1 ./ (1 + 1i*w*tau)) .* (1 ./ (1 + 1i*w/wt));
                if strcmpi(F.passtype, 'highpass')
                    H1 = (1i*w*tau ./ (1 + 1i*w*tau)) .* (1 ./ (1 + 1i*w/wt));
                end
            end
            H = H .* H1;
        case 'sk'
            [B1, B2, B3] = skBiquad('skcoeffs', s.R1, s.R2, s.C1, s.C2, F.passtype, F.gbw);
            jw = 1i*w;
            D = 1 + B1*jw + B2*jw.^2 + B3*jw.^3;
            if strcmpi(F.passtype, 'lowpass')
                N = ones(size(jw));
            else
                N = jw.^2 * (s.R1*s.R2*s.C1*s.C2);
            end
            H = H .* (N ./ D);
        case 'sv'
            Hs = skBiquad('svresp', w, s.comp, F.gbw);
            H = H .* Hs;
    end
end

% ---- ideal prototype response ----
% H_proto(z) = k * prod(z^2 + wz^2) / prod(z - p),  z = jw/wc (LP) or wc/(jw) (HP)
P = F.prototype;
k = P.k; p = P.p; z = P.z;
wc = 2*pi*F.fc;
if strcmpi(F.passtype, 'lowpass')
    zv = 1i*w/wc;
else
    zv = wc ./ (1i*w);
end
num = ones(size(w));
if ~isempty(z)
    for i = 1:2:numel(z) - 1
        num = num .* (zv.^2 + abs(z(i))^2);
    end
end
den = ones(size(w));
for i = 1:numel(p)
    den = den .* (zv - p(i));
end
Hideal = k * num ./ den;
end
