% test_butterworth.m — Butterworth.m 功能验证脚本 (v2: 多标称系列)
% 验证: 1)原型g值 vs 闭式解  2)理想-3dB==fc(LP/HP)  3)E24/E96等系列合法性
%       4)优化不劣化性能(E24/E96)  5)大阶数与边界情形  6)结果展示
fprintf('==== Test 1: 原型 g 值 vs 闭式解 (n=1..10) ====\n');
for n = 1:10
    [~, info] = Butterworth(n, 1e6, 'lowpass', 'Series', 'none', 'Verbose', false);
    gClosed = [2*sin((2*(1:n)-1)*pi/(2*n)), 1];
    err = max(abs(info.g - gClosed));
    assert(err < 1e-8, 'g mismatch n=%d: %g', n, err);
end
fprintf('PASS\n');

fprintf('==== Test 2: 理想 -3dB == fc, 低通 n=1..8 ====\n');
for n = 1:8
    [~, info] = Butterworth(n, 1e6, 'lowpass', 'Series', 'none', 'Verbose', false);
    rel = abs(info.f3dB_ideal - 1e6)/1e6;
    assert(rel < 1e-8, 'LP f3dB mismatch n=%d: %g', n, rel);
end
fprintf('PASS\n');

fprintf('==== Test 3: 理想 -3dB == fc, 高通 n=1..8 ====\n');
for n = 1:8
    [~, info] = Butterworth(n, 1e6, 'highpass', 'Series', 'none', 'Verbose', false);
    rel = abs(info.f3dB_ideal - 1e6)/1e6;
    assert(rel < 1e-8, 'HP f3dB mismatch n=%d: %g', n, rel);
end
fprintf('PASS\n');

fprintf('==== Test 4: 标称值均属于所选 IEC 系列 (E24 默认 / E12 / E48 / E96) ====\n');
E12 = [1.0 1.2 1.5 1.8 2.2 2.7 3.3 3.9 4.7 5.6 6.8 8.2];
E24 = [1.0 1.1 1.2 1.3 1.5 1.6 1.8 2.0 2.2 2.4 2.7 3.0 ...
       3.3 3.6 3.9 4.3 4.7 5.1 5.6 6.2 6.8 7.5 8.2 9.1];
E48 = [1.00 1.05 1.10 1.15 1.21 1.27 1.33 1.40 1.47 1.54 ...
       1.62 1.69 1.78 1.87 1.96 2.05 2.15 2.26 2.37 2.49 ...
       2.61 2.74 2.87 3.01 3.16 3.32 3.48 3.65 3.83 4.02 ...
       4.22 4.42 4.64 4.87 5.11 5.36 5.62 5.90 6.19 6.49 ...
       6.81 7.15 7.50 7.87 8.25 8.66 9.09 9.53];
E96 = [1.00 1.02 1.05 1.07 1.10 1.13 1.15 1.18 1.21 1.24 ...
     1.27 1.30 1.33 1.37 1.40 1.43 1.47 1.50 1.54 1.58 ...
     1.62 1.65 1.69 1.74 1.78 1.82 1.87 1.91 1.96 2.00 ...
     2.05 2.10 2.15 2.21 2.26 2.32 2.37 2.43 2.49 2.55 ...
     2.61 2.67 2.74 2.80 2.87 2.94 3.01 3.09 3.16 3.24 ...
     3.32 3.40 3.48 3.57 3.65 3.74 3.83 3.92 4.02 4.12 ...
     4.22 4.32 4.42 4.53 4.64 4.75 4.87 4.99 5.11 5.23 ...
     5.36 5.49 5.62 5.76 5.90 6.04 6.19 6.34 6.49 6.65 ...
     6.81 6.98 7.15 7.32 7.50 7.68 7.87 8.06 8.25 8.45 ...
     8.66 8.87 9.09 9.31 9.53 9.76];
for n = [1 3 5 7 9]
    % E24 (默认)
    [el, info] = Butterworth(n, 10e6, 'lowpass', 'Verbose', false);
    assert(strcmpi(info.series, 'E24'), 'default series should be E24');
    for k = 1:n
        m = info.nominalValues(k) / 10^floor(log10(info.nominalValues(k)));
        assert(min(abs(E24 - m)) < 0.05, 'not E24 n=%d k=%d m=%g', n, k, m);
    end
    % E12 / E48 / E96
    for s = {'E12','E48','E96'}
        [~, info2] = Butterworth(n, 10e6, 'lowpass', 'Series', s{1}, 'Verbose', false);
        assert(strcmpi(info2.series, s{1}));
        Eref = eval(s{1});
        for k = 1:n
            v = info2.nominalValues(k);
            m = v / 10^floor(log10(v));
            assert(min(abs(Eref - m)) < 0.05, 'not %s n=%d k=%d m=%g', s{1}, n, k, m);
        end
    end
end
fprintf('PASS\n');

fprintf('==== Test 5: 优化不劣化性能(就近 vs 优化), E24 与 E96 ====\n');
for s = {'E24','E96'}
    for n = [3 5 7]
        [~, in1] = Butterworth(n, 10e6, 'lowpass', 'Series', s{1}, 'Optimize', false, 'Verbose', false);
        [~, in2] = Butterworth(n, 10e6, 'lowpass', 'Series', s{1}, 'Optimize', true, 'Verbose', false);
        fprintf('%s n=%d: 就近 rms=%.5f dB, -3dB偏差=%+.4f%% | 优化 rms=%.5f dB, -3dB偏差=%+.4f%%\n', ...
            s{1}, n, in1.responseError_rms_dB, in1.fcDeviationPct, in2.responseError_rms_dB, in2.fcDeviationPct);
        assert(in2.responseError_rms_dB <= in1.responseError_rms_dB + 1e-12, '%s optimization worsened error', s{1});
    end
end
fprintf('PASS\n');

fprintf('==== Test 6: 大阶数 n=20, 32 与 Steps=2, Metric=max ====\n');
[el, info] = Butterworth(20, 1e9, 'lowpass', 'Verbose', false);
assert(numel(info.nominalValues) == 20);
fprintf('n=20 低通(E24): RMS=%.4f dB, 最大=%.4f dB, -3dB偏差=%+.4f%%\n', ...
    info.responseError_rms_dB, info.responseError_max_dB, info.fcDeviationPct);
[el, info] = Butterworth(32, 1e9, 'highpass', 'Verbose', false);
fprintf('n=32 高通(E24): RMS=%.4f dB, 最大=%.4f dB, -3dB偏差=%+.4f%%\n', ...
    info.responseError_rms_dB, info.responseError_max_dB, info.fcDeviationPct);
[el, info] = Butterworth(9, 10e6, 'lowpass', 'Series', 'E24', 'Steps', 2, 'Metric', 'max', 'Verbose', false);
fprintf('n=9 E24 Steps=2 Metric=max: 搜索=%s, 最大误差=%.4f dB\n', info.searchMethod, info.searchError);
fprintf('PASS\n');

fprintf('==== Test 7: 边界情形 n=1, 自定义 Z0/Range, 兼容旧 E96 选项 ====\n');
el = Butterworth(1, 1e6, 'lowpass', 'Verbose', false);
assert(abs(el.Ideal(1) - 2*50/(2*pi*1e6)) < 1e-12);
el = Butterworth(6, 10e6, 'lowpass', 'Z0', 75, 'Range', [1e-9, 1e-3], 'Verbose', false);
[~, in] = Butterworth(5, 10e6, 'lowpass', 'E96', true, 'Verbose', false);  % 旧接口
assert(strcmpi(in.series, 'E96'));
[~, in] = Butterworth(5, 10e6, 'lowpass', 'E96', false, 'Verbose', false);
assert(strcmpi(in.series, 'none'));
fprintf('PASS\n');

fprintf('==== Test 8: 结果展示 (E24 默认) ====\n');
disp(Butterworth(5, 10e6, 'lowpass'));
disp(Butterworth(3, 100e3, 'highpass', 'Z0', 75));

fprintf('\nALL TESTS PASSED\n');
