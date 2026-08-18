function txt = skReport(F)
%SKREPORT  Format the design report text.
%   TXT = SKREPORT(F)  returns the multi-line report string.
%
%   See also SKDESIGN.

L = {};
L{end+1} = '==========================================================';
L{end+1} = sprintf(' Sallen-Key active filter design');
L{end+1} = sprintf('   type     : %s (order %d)', upper(F.type), F.order);
L{end+1} = sprintf('   response : %s', upper(F.passtype));
L{end+1} = sprintf('   fc       : %s Hz', skFmtSI(F.fc, ''));
if isfinite(F.gbw)
    L{end+1} = sprintf('   op-amp GBW: %s Hz (with pre-distortion)', skFmtSI(F.gbw/(2*pi)));
else
    L{end+1} = sprintf('   op-amp GBW: ideal (infinite)');
end
L{end+1} = '----------------------------------------------------------';

for i = 1:numel(F.sections)
    s = F.sections(i);
    L{end+1} = sprintf(' Section %d: %s', i, s.topology);
    switch s.kind
        case 'rc'
            L{end+1} = sprintf('   R = %s ohm,  C = %s F', skFmtSI(s.R), skFmtSI(s.C));
            L{end+1} = sprintf('   target f0 = %s Hz, realized -3 dB = %s Hz (err %.1e)', ...
                skFmtSI(s.w0t/(2*pi)), skFmtSI(s.w3/(2*pi)), sqrt(s.err));
        case 'sk'
            L{end+1} = sprintf('   R1 = %s ohm,  R2 = %s ohm,  C1 = %s F,  C2 = %s F', ...
                skFmtSI(s.R1), skFmtSI(s.R2), skFmtSI(s.C1), skFmtSI(s.C2));
            L{end+1} = sprintf('   f0 target = %s Hz -> realized %s Hz (%.2f%%)', ...
                skFmtSI(s.w0t/(2*pi)), skFmtSI(s.w0a/(2*pi)), ...
                100*abs(s.w0a-s.w0t)/s.w0t);
            L{end+1} = sprintf('   Q  target = %.3f -> realized %.3f (%.2f%%)', ...
                s.Qt, s.Qa, 100*abs(s.Qa-s.Qt)/s.Qt);
            L{end+1} = sprintf('   GBW margin wt/(w0*Q) = %.1f', F.gbw/(s.w0a*s.Qa));
        case 'sv'
            c = s.comp;
            L{end+1} = sprintf('   R1=%s R2=%s R3=%s Rf=%s', ...
                skFmtSI(c.R1), skFmtSI(c.R2), skFmtSI(c.R3), skFmtSI(c.Rf));
            L{end+1} = sprintf('   R5=%s R6=%s R7=%s R8=%s', ...
                skFmtSI(c.R5), skFmtSI(c.R6), skFmtSI(c.R7), skFmtSI(c.R8));
            if strcmp(c.layout, 'hp')
                L{end+1} = sprintf('   (dual) C1=%s C2=%s C3=%s C5=%s C6=%s C7=%s C8=%s', ...
                    skFmtSI(c.C1), skFmtSI(c.C2), skFmtSI(c.C3), skFmtSI(c.C5), ...
                    skFmtSI(c.C6), skFmtSI(c.C7), skFmtSI(c.C8));
            else
                L{end+1} = sprintf('   C1=C2=%s', skFmtSI(c.C1));
            end
            L{end+1} = sprintf('   f0=%s Hz  Q=%.3f  fz=%s Hz  (rel. resp. err %.1e)', ...
                skFmtSI(s.w0t/(2*pi)), s.Qt, skFmtSI(s.wzt/(2*pi)), sqrt(s.err));
    end
end

L{end+1} = '----------------------------------------------------------';
L{end+1} = ' Measured specifications';
if ~isnan(F.specs.f3dB)
    L{end+1} = sprintf('   -3 dB frequency : %s Hz (target %s Hz)', ...
        skFmtSI(F.specs.f3dB), skFmtSI(F.fc));
end
if ~isnan(F.specs.ripple)
    L{end+1} = sprintf('   passband ripple: %.3f dB (target %.3f dB)', ...
        F.specs.ripple, F.prototype.Rp);
end
if ~isnan(F.specs.astop)
    L{end+1} = sprintf('   stopband atten : %.2f dB at %s Hz (target %.2f dB)', ...
        F.specs.astop, skFmtSI(F.specs.fstop), F.prototype.Rs);
end
L{end+1} = sprintf('   max |H| deviation from ideal prototype: %.3f dB', F.specs.maxerr);
L{end+1} = '==========================================================';

txt = strjoin(L, newline);
end
