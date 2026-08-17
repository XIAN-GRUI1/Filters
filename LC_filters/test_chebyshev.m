% test_chebyshev.m — Chebyshev.m 功能验证脚本
% 验证: 1)g值 vs 闭式解  2)功率响应匹配  3)-3dB公式  4)系列合法性
%       5)优化不劣化性能  6)大阶数  7)边界  8)展示
fprintf('==== Test 1: 原型 g 值 vs 闭式解 (n=1..10, 纹波 0.1 dB) ====\n');
for n = 1:10
    [~, info] = Chebyshev(n, 1e6, 'lowpass', 'Ripple', 0.1, 'Series', 'none', 'Verbose', false);
    beta = log(coth(0.1 * log(10) / 40));
    gam = sinh(beta/(2*n));
    a = sin((2*(1:n)-1)*pi/(2*n));
    b = gam^2 + sin((1:n)*pi/n).^2;
    gCl = zeros(1, n+1);
    gCl(1) = 2*a(1)/gam;
    for k = 2:n, gCl(k) = 4*a(k-1)*a(k)/(b(k-1)*gCl(k-1)); end
    if mod(n,2)==0, gCl(n+1) = coth(beta/4)^2; else, gCl(n+1) = 1; end
    rel = max(abs(info.g - gCl))/max(abs(gCl));
    assert(rel < 1e-8, 'g mismatch n=%d: %g', n, rel);
end
fprintf('PASS\n');

fprintf('==== Test 2: 功率传输 |S21|^2 = 1/(1+eps^2*T_n^2) (LP/HP, n=1..8) ====\n');
worst = 0;
for n = 1:8
    for tt = {'lowpass','highpass'}
        [~, info] = Chebyshev(n, 1e6, tt{1}, 'Ripple', 0.1, 'Series', 'none', 'Verbose', false);
        f = 1e6 * linspace(0.05, 3, 30);
        H = ladderT(f, tt{1}, info.idealValues, info.Rs, info.RL);
        p = 4*info.Rs*abs(H).^2/info.RL;
        if strcmpi(tt{1},'lowpass'), x = f/info.fc; else, x = info.fc./f; end
        pIdeal = 1./(1 + info.epsilon^2*chebT(x, n).^2);
        worst = max(worst, max(abs(p - pIdeal)));
    end
end
fprintf('最大功率误差: %.3g\n', worst);
assert(worst < 1e-9, 'power response mismatch: %g', worst);
fprintf('PASS\n');

fprintf('==== Test 3: 理想 -3dB 频率公式 (LP: fc*cosh, HP: fc/cosh) ====\n');
for n = 1:8
    for tt = {'lowpass','highpass'}
        [~, info] = Chebyshev(n, 1e6, tt{1}, 'Ripple', 0.1, 'Series', 'none', 'Verbose', false);
        f3 = 1e6*cosh(acosh(1/info.epsilon)/n);
        if strcmpi(tt{1}, 'highpass'), f3 = 1e12/f3; end
        rel = abs(info.f3dB_ideal - f3)/f3;
        assert(rel < 1e-12, 'f3dB formula mismatch %s n=%d: %g', tt{1}, n, rel);
    end
end
fprintf('PASS\n');

fprintf('==== Test 4: 标称值均属于所选 IEC 系列 (E24 默认 / E96) ====\n');
E24 = [1.0 1.1 1.2 1.3 1.5 1.6 1.8 2.0 2.2 2.4 2.7 3.0 ...
       3.3 3.6 3.9 4.3 4.7 5.1 5.6 6.2 6.8 7.5 8.2 9.1];
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
    [~, info] = Chebyshev(n, 10e6, 'lowpass', 'Verbose', false);
    assert(strcmpi(info.series, 'E24'));
    for k = 1:n
        m = info.nominalValues(k)/10^floor(log10(info.nominalValues(k)));
        assert(min(abs(E24 - m)) < 0.05, 'not E24 n=%d k=%d m=%g', n, k, m);
    end
    [~, info2] = Chebyshev(n, 10e6, 'lowpass', 'Series', 'E96', 'Verbose', false);
    for k = 1:n
        v = info2.nominalValues(k);
        m = v/10^floor(log10(v));
        assert(min(abs(E96 - m)) < 0.005, 'not E96 n=%d k=%d m=%g', n, k, m);
    end
end
fprintf('PASS\n');

fprintf('==== Test 5: 优化不劣化性能 (E24 与 E96) ====\n');
for s = {'E24','E96'}
    for n = [3 5 7]
        [~, in1] = Chebyshev(n, 10e6, 'lowpass', 'Series', s{1}, 'Optimize', false, 'Verbose', false);
        [~, in2] = Chebyshev(n, 10e6, 'lowpass', 'Series', s{1}, 'Optimize', true, 'Verbose', false);
        fprintf('%s n=%d: 就近 rms=%.5f dB | 优化 rms=%.5f dB, -3dB偏差=%+.4f%%\n', ...
            s{1}, n, in1.responseError_rms_dB, in2.responseError_rms_dB, in2.f3dBDeviationPct);
        assert(in2.responseError_rms_dB <= in1.responseError_rms_dB + 1e-12, '%s optimization worsened', s{1});
    end
end
fprintf('PASS\n');

fprintf('==== Test 6: 大阶数 n=20, 纹波 0.01/1 dB ====\n');
[~, info] = Chebyshev(20, 1e9, 'lowpass', 'Ripple', 0.5, 'Verbose', false);
assert(numel(info.nominalValues) == 20);
fprintf('n=20 纹波0.5dB: RMS=%.4f dB, -3dB偏差=%+.4f%%\n', info.responseError_rms_dB, info.f3dBDeviationPct);
[~, info] = Chebyshev(20, 1e9, 'lowpass', 'Ripple', 0.01, 'Verbose', false);
fprintf('n=20 纹波0.01dB: RMS=%.4f dB, -3dB偏差=%+.4f%%\n', info.responseError_rms_dB, info.f3dBDeviationPct);
fprintf('PASS\n');

fprintf('==== Test 7: 边界情形 n=1, 自定义纹波/Z0 ====\n');
el = Chebyshev(1, 1e6, 'lowpass', 'Verbose', false);
beta1 = log(coth(0.1 * log(10) / 40));              % n=1: gamma = sinh(beta/2)
g1 = 2 / sinh(beta1 / 2);
assert(abs(el.Ideal(1) - g1*50/(2*pi*1e6)) < 1e-9, 'n=1 L1 mismatch');
el = Chebyshev(6, 10e6, 'lowpass', 'Z0', 75, 'Ripple', 1.0, 'Verbose', false);
fprintf('PASS\n');

fprintf('==== Test 8: 结果展示 ====\n');
disp(Chebyshev(5, 10e6, 'lowpass'));
disp(Chebyshev(4, 100e3, 'highpass', 'Z0', 75, 'Ripple', 0.5));

fprintf('\nALL TESTS PASSED\n');

% ---------- 辅助(复现局部函数, 便于独立验证) ----------
function H = ladderT(f, type, vals, Rs, RL)
w = 2*pi*f(:).';
nf = numel(w);
A = ones(1,nf); Bm = zeros(1,nf); Cm = zeros(1,nf); Dm = ones(1,nf);
for k = 1:numel(vals)
    if mod(k,2)==1
        if strcmpi(type,'lowpass'), Z = 1j*w*vals(k); else, Z = 1./(1j*w*vals(k)); end
        Bm = A.*Z + Bm;
        Dm = Cm.*Z + Dm;
    else
        if strcmpi(type,'lowpass'), Y = 1j*w*vals(k); else, Y = 1./(1j*w*vals(k)); end
        A = A + Bm.*Y;
        Cm = Cm + Dm.*Y;
    end
end
H = RL./(A.*RL + Bm + Cm.*Rs.*RL + Dm.*Rs);
H = H(:);
end

function T = chebT(x, n)
x = x(:);
T = zeros(size(x));
a = abs(x);
k = a <= 1;
T(k) = cos(n*acos(x(k)));
T(~k) = cosh(n*acosh(a(~k)));
T(~k & x<0) = (-1)^n * T(~k & x<0);
end
