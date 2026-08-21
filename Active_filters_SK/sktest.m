% sktest.m  comprehensive test of the skDesign library
clear
fprintf('================ TEST 1: LP butter 6th (ideal GBW) ================\n');
F1 = skDesign('Type','butter','Order',6,'Fc',1e3,'Verbose',true);
S1 = F1.specs;
assert(abs(S1.f3dB - 1e3)/1e3 < 0.01, 'butter f3dB off: %.3f', S1.f3dB);
assert(S1.maxerr < 0.2, 'butter deviation too large: %.3f dB', S1.maxerr);
fprintf('PASS: f3dB = %.2f Hz, maxerr = %.3f dB\n', S1.f3dB, S1.maxerr);

fprintf('\n================ TEST 2: LP cheby1 5th Rp=1 ================\n');
F2 = skDesign('Type','cheby1','Order',5,'Fc',1e3,'Rp',1,'Verbose',true);
assert(abs(F2.specs.ripple - 1) < 0.1, 'cheby1 ripple off: %.3f', F2.specs.ripple);
fprintf('PASS: ripple = %.3f dB\n', F2.specs.ripple);

fprintf('\n================ TEST 3: LP bessel 4th ================\n');
F3 = skDesign('Type','bessel','Order',4,'Fc',1e3,'Verbose',true);
assert(abs(F3.specs.f3dB - 1e3)/1e3 < 0.01, 'bessel f3dB off: %.3f', F3.specs.f3dB);
fprintf('PASS: f3dB = %.2f Hz\n', F3.specs.f3dB);

fprintf('\n================ TEST 4: HP cheby1 5th ================\n');
F4 = skDesign('Type','cheby1','Order',5,'Fc',1e3,'Rp',1,'PassType','highpass','Verbose',true);
fprintf('PASS: ripple = %.3f dB\n', F4.specs.ripple);

fprintf('\n================ TEST 5: LP cauer 5th Rp=1 Rs=40 ================\n');
F5 = skDesign('Type','cauer','Order',5,'Fc',1e3,'Rp',1,'Rs',40,'Verbose',true);
assert(abs(F5.specs.ripple - 1) < 0.15, 'cauer ripple off: %.3f', F5.specs.ripple);
assert(F5.specs.astop > 38, 'cauer stopband off: %.2f dB', F5.specs.astop);
fprintf('PASS: ripple = %.3f dB, stopband = %.2f dB\n', F5.specs.ripple, F5.specs.astop);

fprintf('\n================ TEST 6: HP cauer -> clear error ================\n');
try
    F6 = skDesign('Type','cauer','Order',4,'Fc',1e3,'Rp',1,'Rs',40,'PassType','highpass');
    error('should have raised an error');
catch e
    assert(strcmp(e.identifier, 'skDesign:HPZeros'), 'unexpected error: %s', e.message);
    fprintf('PASS: %s\n', e.message(1:min(end,80)));
end

fprintf('\n================ TEST 7: LP cheby2 4th Rs=40 ================\n');
F7 = skDesign('Type','cheby2','Order',4,'Fc',1e3,'Rs',40,'Verbose',true);
assert(F7.specs.astop > 38, 'cheby2 stopband off: %.2f dB', F7.specs.astop);
fprintf('PASS: stopband = %.2f dB\n', F7.specs.astop);

fprintf('\n================ TEST 8: GBW predistortion (LP cheby1 4th, GBW=100kHz) ===\n');
F8 = skDesign('Type','cheby1','Order',4,'Fc',1e3,'Rp',1,'GBW',100e3,'Verbose',true);
% compare with the ideal-GBW design: response should stay close
F8i = skDesign('Type','cheby1','Order',4,'Fc',1e3,'Rp',1,'Verbose',false);
w = 2*pi*logspace(1, 4, 2000);
[H8, ~] = skResponse(F8, w);
[H8i, ~] = skResponse(F8i, w);
err = max(abs(20*log10(abs(H8)) - 20*log10(abs(H8i))));
fprintf('max response deviation ideal vs predistorted (finite GBW): %.3f dB\n', err);
assert(err < 2, 'predistortion failed: %.3f dB deviation', err);
fprintf('PASS\n');

fprintf('\n================ TEST 9: GBW predistortion (LP butter 6th, GBW=50kHz) ===\n');
F9 = skDesign('Type','butter','Order',6,'Fc',1e3,'GBW',50e3,'Verbose',true);
assert(abs(F9.specs.f3dB - 1e3)/1e3 < 0.02, 'butter f3dB with GBW off: %.3f', F9.specs.f3dB);
fprintf('PASS: f3dB = %.2f Hz\n', F9.specs.f3dB);

fprintf('\n================ TEST 10: E24 resistors ================\n');
F10 = skDesign('Type','butter','Order',4,'Fc',1e3,'Rser','E24','Verbose',true);
assert(abs(F10.specs.f3dB - 1e3)/1e3 < 0.03, 'E24 butter off: %.3f', F10.specs.f3dB);
fprintf('PASS: f3dB = %.2f Hz (E24)\n', F10.specs.f3dB);

fprintf('\n================ TEST 11: skcoeffs vs exact nodal solve ================\n');
% regression test: the finite-GBW denominator coefficients must match a
% direct 3x3 nodal solve of the unity-gain Sallen-Key circuit
rng(42);
worstLP = 0; worstHP = 0;
for trial = 1:40
    R1 = 10^(rand*5 + 2); R2 = 10^(rand*5 + 2);
    C1 = 10^(-9 - rand*4); C2 = 10^(-9 - rand*4);
    wt = 2*pi*(10^(rand*6 + 2));
    w  = 2*pi*logspace(0, log10(wt/2), 60);
    s  = 1i*w;
    % lowpass: nodal matrix
    M = zeros(3,3,numel(w));
    M(1,1,:) = 1/R1 + 1/R2 + s*C1; M(1,2,:) = -1/R2;      M(1,3,:) = -s*C1;
    M(2,1,:) = -1/R2;              M(2,2,:) = 1/R2 + s*C2; M(2,3,:) = 0;
    M(3,1,:) = 0;                  M(3,2,:) = 1;           M(3,3,:) = -(1 + s/wt);
    V = zeros(3,numel(w));
    for k = 1:numel(w), V(:,k) = M(:,:,k) \ [1/R1; 0; 0]; end
    [B1,B2,B3] = skBiquad('skcoeffs', R1, R2, C1, C2, 'lowpass', wt);
    Hc = 1./(1 + B1*s + B2*s.^2 + B3*s.^3);
    worstLP = max(worstLP, max(abs(20*log10(abs(Hc./V(3,:))))));
    % highpass
    M = zeros(3,3,numel(w));
    M(1,1,:) = s*C1 + 1/R1 + s*C2; M(1,2,:) = -s*C2; M(1,3,:) = -1/R1;
    M(2,1,:) = -s*C2;              M(2,2,:) = s*C2 + 1/R2; M(2,3,:) = 0;
    M(3,1,:) = 0;                  M(3,2,:) = 1;       M(3,3,:) = -(1 + s/wt);
    V = zeros(3,numel(w));
    for k = 1:numel(w), V(:,k) = M(:,:,k) \ [s(k)*C1; 0; 0]; end
    [B1,B2,B3] = skBiquad('skcoeffs', R1, R2, C1, C2, 'highpass', wt);
    Hc = (s.^2*R1*R2*C1*C2)./(1 + B1*s + B2*s.^2 + B3*s.^3);
    worstHP = max(worstHP, max(abs(20*log10(abs(Hc./V(3,:))))));
end
fprintf('worst LP = %.3e dB, worst HP = %.3e dB\n', worstLP, worstHP);
assert(worstLP < 1e-6 && worstHP < 1e-6, ...
    'skcoeffs deviates from the nodal solution (LP %.3e, HP %.3e dB)', worstLP, worstHP);
fprintf('PASS\n');

fprintf('\n================ TEST 12: conditional pre-distortion ================\n');
% butter-6 fc=1k, GBW=100k: margins m = 193, 141, 51.8 (Q = 0.518, 0.707, 1.932)
% -> only the high-Q section must be pre-distorted (m < PredistThresh = 100)
F12 = skDesign('Type','butter','Order',6,'Fc',1e3,'GBW',100e3,'Verbose',true);
fprintf('margins: %.1f %.1f %.1f, pred flags: %d %d %d\n', ...
    F12.sections(1).margin, F12.sections(2).margin, F12.sections(3).margin, ...
    F12.sections(1).pred, F12.sections(2).pred, F12.sections(3).pred);
assert(F12.sections(1).pred == false && F12.sections(2).pred == false, ...
    'low-Q sections must not be pre-distorted');
assert(F12.sections(3).pred == true, 'high-Q section must be pre-distorted');
assert(abs(F12.specs.f3dB - 1e3)/1e3 < 0.02, 'butter f3dB with GBW off: %.3f', F12.specs.f3dB);
% Predist=always / never
F12a = skDesign('Type','butter','Order',6,'Fc',1e3,'GBW',100e3, ...
    'Predist','always','Verbose',false);
F12n = skDesign('Type','butter','Order',6,'Fc',1e3,'GBW',100e3, ...
    'Predist','never','Verbose',false);
assert(all([F12a.sections.pred]) && ~any([F12n.sections.pred]), ...
    'Predist always/never flags wrong');
fprintf('PASS\n');

fprintf('\n================ TEST 13: odd order + finite GBW (RC section) ================\n');
F13 = skDesign('Type','cheby1','Order',5,'Fc',1e3,'Rp',1,'GBW',50e3,'Verbose',true);
assert(abs(F13.specs.ripple - 1) < 0.2, 'cheby1 ripple with GBW off: %.3f', F13.specs.ripple);
fprintf('PASS: ripple = %.3f dB\n', F13.specs.ripple);

fprintf('\n================ TEST 14: HP + finite GBW (exact pre-distortion) ================\n');
F14 = skDesign('Type','cheby1','Order',5,'Fc',1e3,'Rp',1, ...
    'PassType','highpass','GBW',100e3,'Verbose',true);
assert(abs(F14.specs.ripple - 1) < 0.2, 'HP cheby1 ripple with GBW off: %.3f', F14.specs.ripple);
for i = 1:numel(F14.sections)
    s = F14.sections(i);
    if strcmp(s.kind, 'sk')
        fprintf('sec %d: margin=%.1f pred=%d f0err=%.2f%% Qerr=%.2f%%\n', ...
            i, s.margin, s.pred, 100*abs(s.w0a-s.w0t)/s.w0t, 100*abs(s.Qa-s.Qt)/s.Qt);
        assert(100*abs(s.Qa-s.Qt)/s.Qt < 5, 'HP section %d Q err too large', i);
    end
end
fprintf('PASS: ripple = %.3f dB\n', F14.specs.ripple);

fprintf('\n================ ALL TESTS PASSED ================\n');
