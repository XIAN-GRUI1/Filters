%SKDEMO  Demo designs for the Sallen-Key filter library.
%
%   Run the whole script to design and plot several filters:
%     1. Butterworth LP, 6th order, E96/E24 components (ideal op-amps)
%     2. Butterworth LP, 6th order, 800 kHz with GBW = 12 MHz (per-section
%        conditional pre-distortion; check the margin / status lines)
%     3. Cauer (elliptic) LP, 5th order, 1 dB / 40 dB (transmission-zero
%        sections realized with 3-op-amp state-variable biquads)
%     4. Bessel LP, 4th order
%     5. Chebyshev-II LP, 4th order, 40 dB stopband
%
%   See also SKDESIGN.

clc
fprintf('==========================================================\n');
fprintf(' SK filter design demo\n');
fprintf('==========================================================\n');

% 1 ---------------------------------------------------------------
fprintf('\n[1] Butterworth 6th-order low-pass, 1 kHz, ideal op-amps\n');
F1 = skDesign('Type', 'butter', 'Order', 6, 'Fc', 1e3, ...
    'Plot', true);

% 2 ---------------------------------------------------------------
fprintf('\n[2] Butterworth 6th-order low-pass, 800 kHz,\n');
fprintf('    GBW = 12 MHz (only the sections that need it are\n');
fprintf('    pre-distorted; check the margin / status lines)\n');
F2 = skDesign('Type', 'butter', 'Order', 9, 'Fc', 800*1e3, ...
    'GBW', 12e6, 'Plot', true);

% 3 ---------------------------------------------------------------
fprintf('\n[3] Cauer (elliptic) 5th-order low-pass, 1 dB ripple, 40 dB\n');
fprintf('    stopband (zero sections use 3-op-amp biquads)\n');
F3 = skDesign('Type', 'cauer', 'Order', 5, 'Fc', 1e3, 'Rp', 1, 'Rs', 40, ...
    'Plot', true);

% 4 ---------------------------------------------------------------
fprintf('\n[4] Bessel 4th-order low-pass (maximally flat group delay)\n');
F4 = skDesign('Type', 'bessel', 'Order', 4, 'Fc', 1e3, 'Plot', true);

% 5 ---------------------------------------------------------------
fprintf('\n[5] Chebyshev-II 4th-order low-pass, 40 dB stopband\n');
F5 = skDesign('Type', 'cheby2', 'Order', 4, 'Fc', 1e3, 'Rs', 40, ...
    'Plot', true);

fprintf('\n==========================================================\n');
fprintf(' All demo designs completed.\n');
fprintf(' Inspect the F1..F5 structs for components and responses.\n');
fprintf('==========================================================\n');
