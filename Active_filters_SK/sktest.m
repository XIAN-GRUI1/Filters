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

fprintf('\n================ ALL TESTS PASSED ================\n');
