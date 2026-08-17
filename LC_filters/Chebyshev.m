function [elements, info] = Chebyshev(n, fc, type, varargin)
%CHEBYSHEV 基于现代网络函数综合法的切比雪夫(I型)LC滤波器设计
%
%   [elements, info] = Chebyshev(n, fc, type)
%   [elements, info] = Chebyshev(n, fc, type, 'Name', Value)
%   Chebyshev('demo')          % 运行演示(含自检)
%   Chebyshev()                % 显示本帮助
%
% 功能
% ----
%   1. 基于"现代网络函数综合法"(Darlington 综合: 功率传输函数 ->
%      特征函数 -> 反射系数 -> 输入阻抗 -> 连分式展开), 求解理论上
%      任意阶数 n、任意截止频率 fc 的切比雪夫(等纹波)低通/高通滤波器的
%      LC 元件精确值(仅使用 MATLAB 基础功能, 不依赖任何工具箱);
%   2. 元件值就近匹配 IEC 60063 标称系列(E12/E24/E48/E96, 默认 E24),
%      并默认在标称值邻域内做穷举/坐标下降搜索, 使幅频响应与理想响应
%      偏差最小, 即"在尽可能保证性能的前提下完成标称值匹配"。
%
% 综合原理
% --------
%   归一化原型(|S21|^2 = 1/(1 + eps^2*T_n^2), eps 由纹波决定):
%     Tn(s)   : 第一类切比雪夫多项式, 特征函数 F(s) = eps*(-1j)^n*T_n(s/1j);
%     poles   : 椭圆分布极点 s_k = -sinh(phi)*sin(th_k) + j*cosh(phi)*cos(th_k),
%               phi = asinh(1/eps)/n, th_k = (2k-1)*pi/(2n);
%     E(s)    : 极点多项式(首一); D(s) = E(s)*lead(F), 使 D-D(-) 与 F 首项
%               系数一致, 从而 D(s)D(-s) = 1 + eps^2*T_n^2(s/1j);
%     Zin(s)  : 输入阻抗 = (D + F) / (D - F);
%     g1..gn  : 对 Zin(s) 作连分式展开(交替提取 s=inf 处的串联感抗与
%               并联容抗); g(n+1) 为负载电阻(奇阶 = 1, 偶阶 = coth^2(beta/4),
%               即偶阶为不等端接设计)。
%   结果与经典闭式解一致, 程序自动交叉校验; n > 10 时双精度多项式
%   连分式数值病态, 自动采用等效闭式解(数值精确, 阶数理论上任意)。
%   经阻抗标度(×Z0)与频率标度(×wc = 2*pi*fc)去归一化得到实际 L/C;
%   高通由低通原型经 s -> 1/s 变换得到(串联电容 / 并联电感)。
%
% 输入
% ----
%   n    正整数, 滤波器阶数(理论上任意; 连分式展开精确到 n=10,
%         更大阶数自动采用等效闭式解, 数值精确)
%   fc   通带边缘(纹波边缘)频率, 单位 Hz; 注意切比雪夫 -3dB 频率
%         > fc, 由 info.f3dB_* 给出
%   type 'lowpass'|'lp' 或 'highpass'|'hp'(不区分大小写)
%
% 可选名值对
% ----------
%   'Z0'        源端接阻抗(欧姆), 默认 50; 偶阶时负载 = Z0*g(n+1) (不等端接)
%   'Ripple'    通带纹波(dB), 默认 0.1; 要求 0 < Ripple < 3.01
%               (纹波 > 3.01 dB 时通带凹点低于 -3dB, "相对通带峰值 -3dB"
%                的定义失效, 仍可计算但含义需注意)
%   'Series'    标称系列: 'E12'|'E24'|'E48'|'E96'|'none'(精确值),
%               默认 'E24'(贴片常用)
%   'E96'       兼容别名: true -> Series='E96', false -> Series='none'
%   'Range'     标称值允许范围 [min, max](SI 单位), 默认 [1e-15, 1e6]
%   'Optimize'  true/false, 是否在标称值邻域内优化响应, 默认 true
%   'Search'    'auto'|'exhaustive'|'greedy'|'none'
%               默认 'auto': n<=10 穷举全部组合, 更大阶数用坐标下降
%   'Steps'     每个元件搜索的标称邻域档数(±Steps), 默认 1
%   'Metric'    'rms'|'max', 优化目标函数, 默认 'rms'
%   'Plot'      true/false, 绘制理想与标称幅频响应对比, 默认 false
%   'Verbose'   true/false, 是否在命令行打印结果汇总, 默认 true
%
% 输出
% ----
%   elements  table, 每行一个元件: 名称/类型/位置/理想值/标称值/偏差%
%   info      struct, 原型 g 值、极点、纹波、-3dB 频率、响应误差等
%
% 示例
% ----
%   Chebyshev('demo');
%   el = Chebyshev(5, 10e6, 'lowpass');                  % 5 阶 10MHz, 0.1dB纹波
%   el = Chebyshev(4, 100e3, 'highpass', 'Z0', 75, 'Ripple', 0.5);
%   el = Chebyshev(7, 1e9, 'lp', 'Series', 'E96');       % E96 精密系列
%   [el, in] = Chebyshev(3, 1e8, 'lp', 'Plot', true);
%
% 版本: 1.0   适用 MATLAB R2016b 及以上(无工具箱依赖)

% ======================================================================
% 0. 演示 / 帮助
% ======================================================================
if (ischar(n) || isstring(n)) && strcmpi(char(n), 'demo')
    demo();
    if nargout > 0, elements = []; info = []; end
    return;
end
if nargin == 0
    help Chebyshev;
    if nargout > 0, elements = []; info = []; end
    return;
end

% ======================================================================
% 1. 输入校验
% ======================================================================
validateattributes(n, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'n', 1);
validateattributes(fc, {'numeric'}, {'scalar','positive'}, mfilename, 'fc', 2);
if ischar(type) || isstring(type)
    t = lower(char(type));
    switch t
        case {'lowpass','lp','low'},  type = 'lowpass';
        case {'highpass','hp','high'}, type = 'highpass';
        otherwise, error('Chebyshev:badType', 'type 必须是 lowpass 或 highpass。');
    end
else
    error('Chebyshev:badType', 'type 必须是 lowpass 或 highpass。');
end

% ======================================================================
% 2. 可选参数
% ======================================================================
opts = struct('Z0', 50, 'Ripple', 0.1, 'Series', 'E24', 'Range', [1e-15, 1e6], ...
              'Optimize', true, 'Search', 'auto', 'Steps', 1, ...
              'Metric', 'rms', 'Plot', false, 'Verbose', true);
optNames = {'Z0','Ripple','Series','E96','Range','Optimize','Search','Steps','Metric','Plot','Verbose'};
k = 1;
while k <= numel(varargin)
    nm = validatestring(varargin{k}, optNames, mfilename);
    k = k + 1;
    if k > numel(varargin)
        error('Chebyshev:opt', '名值对 %s 缺少取值。', nm);
    end
    v = varargin{k};
    k = k + 1;
    switch nm
        case 'Z0'
            validateattributes(v, {'numeric'}, {'scalar','positive'}, mfilename, 'Z0');
            opts.Z0 = v;
        case 'Ripple'
            validateattributes(v, {'numeric'}, {'scalar','positive'}, mfilename, 'Ripple');
            opts.Ripple = v;
            if v >= 3.01
                warning('Chebyshev:ripple', '纹波 >= 3.01 dB 时, 通带凹点低于 -3dB, "-3dB 相对通带峰值"定义失效。');
            end
        case 'Series'
            opts.Series = validatestring(lower(char(v)), {'e12','e24','e48','e96','none'}, ...
                                         mfilename, 'Series');
        case 'E96'
            opts.Series = 'none';
            if tf(v), opts.Series = 'E96'; end
        case 'Range'
            validateattributes(v, {'numeric'}, {'numel',2,'increasing','positive'}, mfilename, 'Range');
            opts.Range = v(:).';
        case 'Optimize'
            opts.Optimize = tf(v);
        case 'Search'
            opts.Search = validatestring(lower(char(v)), {'auto','exhaustive','greedy','none'}, mfilename, 'Search');
        case 'Steps'
            validateattributes(v, {'numeric'}, {'scalar','integer','nonnegative'}, mfilename, 'Steps');
            opts.Steps = v;
        case 'Metric'
            opts.Metric = validatestring(lower(char(v)), {'rms','max'}, mfilename, 'Metric');
        case 'Plot'
            opts.Plot = tf(v);
        case 'Verbose'
            opts.Verbose = tf(v);
    end
end

% ======================================================================
% 3. 原型综合: 切比雪夫多项式/极点 -> 反射系数 -> 输入阻抗 -> 连分式展开
% ======================================================================
[g, poles, D, F, method] = chebyshevPrototype(n, opts.Ripple);
gL = real(g(end));                                % 负载原型值(偶阶 > 1)
g0 = 1;                                           % 归一化源电阻
eps_ = sqrt(10^(opts.Ripple/10) - 1);

% ======================================================================
% 4. 阻抗/频率去归一化 -> 实际 L/C 元件值与端接
% ======================================================================
[idealVals, elemType, elemPos, units] = scaleElements(g, type, fc, opts.Z0);
Rs = opts.Z0;                                     % 源阻抗
RL = opts.Z0 * gL;                                % 负载阻抗(偶阶不等端接)

% ======================================================================
% 5. 性能评估网格与理想响应(过渡带网格, 优化与报告一致)
% ======================================================================
fgrid = fc * logspace(-0.6, 0.6, 41);
HidB  = idealChebyshevDB(fgrid, fc, n, eps_, gL, type);

% ======================================================================
% 6. 标称系列匹配(就近取值 + 邻域搜索保性能)
% ======================================================================
nominalVals = idealVals;
searchInfo  = struct('method', 'none', 'error', 0);
if ~strcmpi(opts.Series, 'none')
    nominalVals = arrayfun(@(x) seriesNearest(x, opts.Series, opts.Range), idealVals);
    if any(isnan(nominalVals))
        warning('Chebyshev:range', '部分元件值超出标称值范围, 相关元件保留理想值。');
        nominalVals(isnan(nominalVals)) = idealVals(isnan(nominalVals));
    end
    if opts.Optimize && ~strcmpi(opts.Search, 'none')
        [nominalVals, searchInfo] = optimizeSeries(nominalVals, idealVals, fgrid, HidB, type, Rs, RL, opts);
    end
end

% ======================================================================
% 7. -3dB 截止频率(相对设计通带峰值)与响应误差评估
% ======================================================================
f3ideal = fc * cosh(acosh(1/eps_) / n);           % 理想 -3dB 频率(低通)
if strcmpi(type, 'highpass')
    f3ideal = fc^2 / f3ideal;                     % 高通: 取倒数 fc/cosh(...)
end
f3nom   = cutoffFreq(fc, type, nominalVals, Rs, RL, f3ideal, gL);
rmsErr  = responseError(nominalVals, fgrid, HidB, type, Rs, RL, 'rms');
maxErr  = responseError(nominalVals, fgrid, HidB, type, Rs, RL, 'max');

% ======================================================================
% 8. 输出表格 elements
% ======================================================================
dev    = (nominalVals - idealVals) ./ idealVals * 100;
nElem  = numel(idealVals);
names  = cell(1, nElem);
for k = 1:nElem
    if strcmpi(elemPos{k}, 'series')
        if strcmpi(type, 'lowpass'), names{k} = sprintf('L%d', k);
        else,                        names{k} = sprintf('C%d', k); end
    else
        if strcmpi(type, 'lowpass'), names{k} = sprintf('C%d', k);
        else,                        names{k} = sprintf('L%d', k); end
    end
end
txtI = cell(nElem, 1); txtN = cell(nElem, 1);
for k = 1:nElem
    txtI{k} = engstr(idealVals(k),    units{k});
    txtN{k} = engstr(nominalVals(k),  units{k});
end
elements = table(names.', elemType.', elemPos.', idealVals.', units.', ...
                 nominalVals.', units.', dev.', ...
                 'VariableNames', {'Element','Type','Position','Ideal','Unit', ...
                                   'Nominal','UnitN','DevPct'});
elements.IdealText   = txtI;
elements.NominalText = txtN;
Rrows = table({'Rs';'RL'}, {'R';'R'}, {'source';'load'}, ...
              [Rs; RL], {'ohm';'ohm'}, [Rs; RL], {'ohm';'ohm'}, [0; 0], ...
              'VariableNames', {'Element','Type','Position','Ideal','Unit', ...
                                'Nominal','UnitN','DevPct'});
Rrows.IdealText   = {engstr(Rs,'ohm'); engstr(RL,'ohm')};
Rrows.NominalText = Rrows.IdealText;
elements = [elements; Rrows];

% ======================================================================
% 9. 输出结构体 info
% ======================================================================
info = struct();
info.n = n;
info.fc = fc;
info.type = type;
info.Z0 = opts.Z0;
info.Ripple = opts.Ripple;
info.epsilon = eps_;
info.series = opts.Series;
info.g0 = g0;
info.g  = g;                                   % g(1..n+1), g(n+1) 为负载
info.gL = gL;
info.Rs = Rs;
info.RL = RL;
info.poles = poles;
info.D = D;
info.F = F;
info.method = method;
info.idealValues   = idealVals;
info.nominalValues = nominalVals;
info.deviationPct  = dev;
info.f3dB_ideal    = f3ideal;
info.f3dB_nominal  = f3nom;
info.f3dBDeviationPct = (f3nom - f3ideal) / f3ideal * 100;
info.responseError_rms_dB = rmsErr;
info.responseError_max_dB = maxErr;
info.searchMethod = searchInfo.method;
info.searchError  = searchInfo.error;

% ======================================================================
% 10. 绘图(可选)
% ======================================================================
if opts.Plot
    plotResponse(fc, n, type, eps_, gL, Rs, RL, idealVals, nominalVals);
end

% ======================================================================
% 11. 命令行汇总(可选)
% ======================================================================
if opts.Verbose
    fprintf('\n===== 切比雪夫 %d 阶%s, fc = %s (纹波 %g dB), Z0 = %g ohm =====\n', ...
            n, type, engstr(fc, 'Hz'), opts.Ripple, opts.Z0);
    disp(elements);
    fprintf('综合方法 : %s\n', method);
    if strcmpi(opts.Series, 'none')
        fprintf('标称匹配 : false (精确值)\n');
    else
        fprintf('标称匹配 : true  (系列: %s, 搜索: %s, 目标: %s)\n', ...
                upper(opts.Series), searchInfo.method, opts.Metric);
    end
    fprintf('端接     : Rs = %s, RL = %s %s\n', engstr(Rs,'ohm'), engstr(RL,'ohm'), ...
            ternaryStr(mod(n,2)==0, '(偶阶不等端接)', ''));
    fprintf('理想 -3dB 频率 : %s\n', engstr(f3ideal, 'Hz'));
    fprintf('标称 -3dB 频率 : %s (相对理想偏差 %+.3f%%)\n', engstr(f3nom, 'Hz'), info.f3dBDeviationPct);
    fprintf('响应 RMS 误差  : %.4f dB, 最大误差 %.4f dB\n', rmsErr, maxErr);
    fprintf('============================================================\n');
end
end

% ======================================================================
% 以下为局部函数
% ======================================================================

function demo()
% 演示与自检
fprintf('\n================ 切比雪夫 LC 滤波器设计演示 ================\n');
fprintf('\n--- 示例 1: 5 阶低通, fc = 10 MHz, 纹波 0.1 dB, E24 ---\n');
[el1, in1] = Chebyshev(5, 10e6, 'lowpass', 'Plot', true);
disp(el1);
fprintf('综合方法 : %s\n', in1.method);
fprintf('-3dB     : %s | RMS 误差: %.4f dB\n', engstr(in1.f3dB_nominal, 'Hz'), in1.responseError_rms_dB);

fprintf('\n--- 示例 2: 4 阶高通, fc = 100 kHz, 纹波 0.5 dB, Z0 = 75 ---\n');
[el2, in2] = Chebyshev(4, 100e3, 'highpass', 'Z0', 75, 'Ripple', 0.5, 'Plot', true);
disp(el2);
fprintf('综合方法 : %s\n', in2.method);
fprintf('-3dB     : %s | RMS 误差: %.4f dB\n', engstr(in2.f3dB_nominal, 'Hz'), in2.responseError_rms_dB);

fprintf('\n--- 对比: 同一 5 阶 10 MHz 低通, 不同纹波 ---\n');
for r = [0.01 0.1 0.5 1.0]
    [~, in] = Chebyshev(5, 10e6, 'lowpass', 'Ripple', r, 'Verbose', false);
    fprintf('  纹波 %4.2f dB: -3dB = %s (fc 的 %5.2f 倍), RMS 误差 %.4f dB\n', r, ...
            engstr(in.f3dB_nominal, 'Hz'), in.f3dB_nominal/in.fc, in.responseError_rms_dB);
end

fprintf('\n--- 自检 1: 原型 g 值 vs 闭式解 (n=1..8, 纹波 0.1 dB) ---\n');
worst = 0;
for nn = 1:8
    [~, in] = Chebyshev(nn, 1e6, 'lowpass', 'Series', 'none', 'Verbose', false);
    gCl = closedG(nn, 0.1);
    e = max(abs(in.g - gCl)) / max(abs(gCl));
    worst = max(worst, e);
end
fprintf('  最大相对误差 %.3g (应远小于 1e-6)\n', worst);

fprintf('\n--- 自检 2: 梯形网络功率传输 |S21|^2 = 1/(1+eps^2*T_n^2) ---\n');
worst = 0;
for nn = 1:6
    for tt = {'lowpass','highpass'}
        [~, in] = Chebyshev(nn, 1e6, tt{1}, 'Series', 'none', 'Verbose', false);
        f = 1e6 * linspace(0.08, 3, 25);
        H = ladderResponse(f, tt{1}, in.idealValues, in.Rs, in.RL);
        p = 4 * in.Rs * abs(H).^2 / in.RL;                  % 功率传输
        if strcmpi(tt{1}, 'lowpass'), x = f/in.fc; else, x = in.fc./f; end
        pIdeal = 1 ./ (1 + in.epsilon^2 * chebyshevT(x, nn).^2);
        e = max(abs(p - pIdeal));
        worst = max(worst, e);
    end
end
fprintf('  最大功率误差 %.3g (应远小于 1e-9)\n', worst);

fprintf('\n--- 自检 3: 理想 -3dB 频率公式 fc*cosh(acosh(1/eps)/n) ---\n');
worst = 0;
for nn = 1:8
    [~, in] = Chebyshev(nn, 1e6, 'lowpass', 'Series', 'none', 'Verbose', false);
    f3 = 1e6 * cosh(acosh(1/in.epsilon)/nn);
    e = abs(in.f3dB_ideal - f3)/f3;
    worst = max(worst, e);
end
fprintf('  最大相对误差 %.3g (应远小于 1e-12)\n', worst);
fprintf('============================================================\n');
end

% ----------------------------------------------------------------------
function [g, poles, D, F, method] = chebyshevPrototype(n, rippleDb)
% 归一化切比雪夫原型综合(网络函数综合法: 连分式展开)
method = '现代网络函数综合法 (特征函数->反射系数->输入阻抗->连分式展开)';
eps_ = sqrt(10^(rippleDb/10) - 1);

% 椭圆极点: s_k = -sinh(phi)sin(th_k) + j*cosh(phi)cos(th_k)
phi = asinh(1/eps_) / n;
th  = (2*(1:n).' - 1) * pi / (2*n);
poles = -sinh(phi)*sin(th) + 1j*cosh(phi)*cos(th);

E = poly(poles);                                  % 首一极点多项式
E = E(:).';

% 特征函数 F(s) = eps*(-1j)^n*T_n(s/1j) (实系数, 首项 = eps*2^(n-1)*(-1)^n)
F = chebyshevF(n, eps_);
fLead = F(1);                                     % F 的首项系数
D = E * fLead;                                    % 使 D-D(-) 首项系数为 0

% 闭式解(综合结果的解析形式), 用于交叉校验与大阶数回退
gClosed = closedG(n, rippleDb);

if n <= 10
    Znum = D + F;
    Zden = D - F;
    g = ladderExtract(Znum, Zden, n);
    rel = max(abs(g - gClosed)) / max(abs(gClosed));
    if rel > 1e-6 * max(1, n / 10)
        warning('Chebyshev:cf', '连分式展开数值误差较大 (rel=%.2e), 改用闭式解。', rel);
        g = gClosed;
        method = [method, ' (数值误差回退: 采用闭式解)'];
    end
else
    g = gClosed;
    method = [method, ' (n>10 双精度连分式病态, 采用等效闭式解)'];
end
g = real(g);                                      % 消除复数舍入噪声
end

function F = chebyshevF(n, eps_)
% F(s) = eps * (-1j)^n * T_n(s/1j), 实系数多项式(降幂)
c = chebyshevPoly(n);
m = numel(c);
F = zeros(1, m);
for k = 1:m
    deg = m - k;                                  % x 的次数
    switch mod(n + deg, 4)                        % (-1j)^(n+deg)
        case 0, coef = 1;
        case 1, coef = -1j;
        case 2, coef = -1;
        case 3, coef = 1j;
    end
    F(k) = c(k) * coef * eps_;
end
F = real(F);
end

function c = chebyshevPoly(n)
% 第一类切比雪夫多项式 T_n(x) 系数(降幂, 实系数)
if n == 0, c = 1; return; end
t0 = 1; t1 = [1 0];                               % T_0 = 1, T_1 = x
for kk = 2:n
    t2 = 2*[t1 0] - [0 0 t0];                     % T_{k+1} = 2x*T_k - T_{k-1}
    t0 = t1; t1 = t2;
end
c = t1;
end

function g = closedG(n, rippleDb)
% 切比雪夫原型闭式解(经典公式)
% 注意: 文献常用 17.37, 精确值为 40/log(10) = 17.37155, 此处用精确值
% 以保证与连分式综合完全一致
beta = log(coth(rippleDb * log(10) / 40));
gam  = sinh(beta / (2*n));
a = sin((2*(1:n).' - 1) * pi / (2*n));
b = gam^2 + sin((1:n).' * pi / n).^2;
g = zeros(1, n+1);
g(1) = 2*a(1) / gam;
for kk = 2:n
    g(kk) = 4*a(kk-1)*a(kk) / (b(kk-1)*g(kk-1));
end
if mod(n, 2) == 0
    g(n+1) = coth(beta/4)^2;
else
    g(n+1) = 1;
end
g = g(:).';
end

% ----------------------------------------------------------------------
function g = ladderExtract(Znum, Zden, n)
% 对输入阻抗 Zin(s) = Znum/Zden 作连分式展开:
% 交替提取 s->inf 处的极点(串联元件 g奇 / 并联元件 g偶), 余项为负载电阻。
% 数值稳定策略: 多项式阶数完全由结构确定, 每次减法后余式阶数精确下降 2
% (前两项系数数学上为 0), 末步下降 1 —— 直接按阶数截断, 无需容差剥离。
g = zeros(1, n + 1);
N = Znum(:).';                        % deg = n
D = Zden(:).';                        % 首项已精确为 0
D = D(2:end);                         % 去掉该零项, deg = n-1
for kk = 1:n
    g(kk) = N(1) / D(1);
    Ds = [D, 0];                      % s * D
    N = N - g(kk) * Ds;
    if numel(D) >= 2
        N = N(3:end);                 % 余式阶数下降 2: 前两项精确为 0
    else
        N = N(2:end);                 % 末步: 仅首项为 0
    end
    if isempty(N), N = 0; end
    tmp = N; N = D; D = tmp;          % 取倒数进入下一轮
end
g(n + 1) = N(1) / D(1);               % 负载电阻
end

% ----------------------------------------------------------------------
function [vals, elemType, elemPos, units] = scaleElements(g, type, fc, Z0)
% 原型 g(1..n) -> 实际 L/C (SI 单位), 源阻抗 = Z0
wc = 2 * pi * fc;
n  = numel(g) - 1;
vals = zeros(1, n);
elemType = cell(1, n);
elemPos  = cell(1, n);
units    = cell(1, n);
for k = 1:n
    if mod(k, 2) == 1                                     % 串联元件
        elemPos{k} = 'series';
        if strcmpi(type, 'lowpass')
            vals(k) = g(k) * Z0 / wc;                     % 串联电感
            elemType{k} = 'L'; units{k} = 'H';
        else
            vals(k) = 1 / (g(k) * Z0 * wc);               % 串联电容
            elemType{k} = 'C'; units{k} = 'F';
        end
    else                                                  % 并联元件
        elemPos{k} = 'shunt';
        if strcmpi(type, 'lowpass')
            vals(k) = g(k) / (Z0 * wc);                   % 并联电容
            elemType{k} = 'C'; units{k} = 'F';
        else
            vals(k) = Z0 / (g(k) * wc);                   % 并联电感
            elemType{k} = 'L'; units{k} = 'H';
        end
    end
end
end

% ----------------------------------------------------------------------
function E = seriesMantissas(series)
% IEC 60063 标称系列: E12/E24/E48/E96 (容差分别约 ±10%/±5%/±2%/±1%)
switch lower(series)
    case 'e12'
        E = [1.0 1.2 1.5 1.8 2.2 2.7 3.3 3.9 4.7 5.6 6.8 8.2];
    case 'e24'
        E = [1.0 1.1 1.2 1.3 1.5 1.6 1.8 2.0 2.2 2.4 2.7 3.0 ...
             3.3 3.6 3.9 4.3 4.7 5.1 5.6 6.2 6.8 7.5 8.2 9.1];
    case 'e48'
        E = [1.00 1.05 1.10 1.15 1.21 1.27 1.33 1.40 1.47 1.54 ...
             1.62 1.69 1.78 1.87 1.96 2.05 2.15 2.26 2.37 2.49 ...
             2.61 2.74 2.87 3.01 3.16 3.32 3.48 3.65 3.83 4.02 ...
             4.22 4.42 4.64 4.87 5.11 5.36 5.62 5.90 6.19 6.49 ...
             6.81 7.15 7.50 7.87 8.25 8.66 9.09 9.53];
    otherwise                                            % 'e96'
        E = [1.00 1.02 1.05 1.07 1.10 1.13 1.15 1.18 1.21 1.24 ...
             1.27 1.30 1.33 1.37 1.40 1.43 1.47 1.50 1.54 1.58 ...
             1.62 1.65 1.69 1.74 1.78 1.82 1.87 1.91 1.96 2.00 ...
             2.05 2.10 2.15 2.21 2.26 2.32 2.37 2.43 2.49 2.55 ...
             2.61 2.67 2.74 2.80 2.87 2.94 3.01 3.09 3.16 3.24 ...
             3.32 3.40 3.48 3.57 3.65 3.74 3.83 3.92 4.02 4.12 ...
             4.22 4.32 4.42 4.53 4.64 4.75 4.87 4.99 5.11 5.23 ...
             5.36 5.49 5.62 5.76 5.90 6.04 6.19 6.34 6.49 6.65 ...
             6.81 6.98 7.15 7.32 7.50 7.68 7.87 8.06 8.25 8.45 ...
             8.66 8.87 9.09 9.31 9.53 9.76];
end
end

function v = seriesNearest(x, series, range)
% 就近标称值; 超出 range 返回 NaN
E = seriesMantissas(series);
e = floor(log10(x));
m = x / 10^e;
[~, ci] = min(abs(E - m));
v = E(ci) * 10^e;
if v < range(1) || v > range(2), v = NaN; end
end

function vals = seriesCandidates(x, series, steps, range)
% 目标值附近的标称候选集(±steps 档), 按与理想值距离升序排列
E = seriesMantissas(series);
nE = numel(E);
e = floor(log10(x));
m = x / 10^e;
[~, ci] = min(abs(E - m));
idxs = (ci - steps) : (ci + steps);
mm = zeros(size(idxs));
dd = zeros(size(idxs));
for t = 1:numel(idxs)
    j = idxs(t);
    if j < 1
        mm(t) = E(j + nE);  dd(t) = -1;
    elseif j > nE
        mm(t) = E(j - nE);  dd(t) = 1;
    else
        mm(t) = E(j);       dd(t) = 0;
    end
end
vals = mm .* (10 .^ (e + dd));
vals = unique(vals);
vals = vals(vals >= range(1) & vals <= range(2));
[~, order] = sort(abs(vals - x));
vals = vals(order);
end

% ----------------------------------------------------------------------
function H = ladderResponse(f, type, vals, Rs, RL)
% 逐级 ABCD 矩阵级联, 计算梯形网络电压传输函数 H = Vload/Vsource
w = 2 * pi * f(:).';
nf = numel(w);
A  = ones(1, nf);
Bm = zeros(1, nf);
Cm = zeros(1, nf);
Dm = ones(1, nf);
n = numel(vals);
for k = 1:n
    if mod(k, 2) == 1                                     % 串联元件
        if strcmpi(type, 'lowpass')
            Z = 1j * w * vals(k);                         % 串联电感
        else
            Z = 1 ./ (1j * w * vals(k));                  % 串联电容
        end
        Bm = A .* Z + Bm;
        Dm = Cm .* Z + Dm;
    else                                                  % 并联元件
        if strcmpi(type, 'lowpass')
            Y = 1j * w * vals(k);                         % 并联电容
        else
            Y = 1 ./ (1j * w * vals(k));                  % 并联电感
        end
        A  = A + Bm .* Y;                                 % 新A = A + B*Y
        Cm = Cm + Dm .* Y;                                % 新C = C + D*Y
    end
end
H = RL ./ (A .* RL + Bm + Cm .* Rs .* RL + Dm .* Rs);
H = H(:);
end

% ----------------------------------------------------------------------
function HdB = idealChebyshevDB(f, fc, n, eps_, gL, type)
% 理想切比雪夫幅频响应(含 6.02dB 端接损耗与偶阶负载因子)
% |S21|^2 = 1/(1+eps^2*T_n^2), |H| = |S21|/2 * sqrt(gL)
if strcmpi(type, 'lowpass')
    x = f(:) / fc;
else
    x = fc ./ f(:);
end
T = chebyshevT(x, n);
H = 0.5 * sqrt(gL) ./ sqrt(1 + eps_^2 * T .^ 2);
HdB = 20 * log10(max(abs(H), 1e-300));
end

function T = chebyshevT(x, n)
% 第一类切比雪夫多项式 T_n(x) 数值求值(实轴, 数值稳定)
x = x(:);
T = zeros(size(x));
a = abs(x);
k = a <= 1;
T(k)  = cos(n * acos(x(k)));
T(~k) = cosh(n * acosh(a(~k)));
T(~k & x < 0) = (-1)^n * T(~k & x < 0);
end

% ----------------------------------------------------------------------
function err = responseError(vals, fgrid, HidB, type, Rs, RL, metricName)
% 标称元件响应与理想响应的偏差(dB 域)
H = ladderResponse(fgrid, type, vals, Rs, RL);
HdB = 20 * log10(max(abs(H), 1e-300));
d = HdB - HidB;
switch lower(metricName)
    case 'max'
        err = max(abs(d));
    otherwise
        err = sqrt(mean(d .^ 2));
end
end

% ----------------------------------------------------------------------
function [bestVals, out] = optimizeSeries(startVals, idealVals, fgrid, HidB, type, Rs, RL, opts)
% 在标称系列候选集内搜索使响应误差最小的元件组合(保性能标称值匹配)
n = numel(idealVals);
cands = cell(1, n);
nCand = zeros(1, n);
for k = 1:n
    c = seriesCandidates(idealVals(k), opts.Series, opts.Steps, opts.Range);
    if isempty(c)
        c = startVals(k);
        warning('Chebyshev:noCand', '元件 %d 无可用标称候选, 保留就近值。', k);
    end
    cands{k} = c;
    nCand(k) = numel(c);
end
metric = @(v) responseError(v, fgrid, HidB, type, Rs, RL, opts.Metric);
useEx  = strcmpi(opts.Search, 'exhaustive') || ...
         (strcmpi(opts.Search, 'auto') && n <= 10);
total  = prod(nCand);
if useEx && total <= 1e6
    [bestVals, bestErr] = exhaustiveSearch(cands, metric);
    out.method = 'exhaustive';
else
    [bestVals, bestErr] = greedySearch(cands, metric);
    out.method = 'coordinate-descent';
    if useEx && total > 1e6
        warning('Chebyshev:combos', '组合数 %g 过大, 改用坐标下降搜索。', total);
    end
end
out.error = bestErr;
end

function [best, bestErr] = exhaustiveSearch(cands, metric)
% 全组合穷举(组合数受限时使用)
n = numel(cands);
sizes = cellfun(@numel, cands);
total = prod(sizes);
idx = ones(1, n);
best = arrayfun(@(k) cands{k}(idx(k)), 1:n);
bestErr = metric(best);
for c = 2:total
    for k = n:-1:1
        idx(k) = idx(k) + 1;
        if idx(k) <= sizes(k), break; end
        idx(k) = 1;
    end
    vals = arrayfun(@(k) cands{k}(idx(k)), 1:n);
    e = metric(vals);
    if e < bestErr, bestErr = e; best = vals; end
end
end

function [best, bestErr] = greedySearch(cands, metric)
% 坐标下降( Gauss-Seidel ) + 多起点重启, 降低局部最优风险
% 对粗系列(E12/E24, 档距 ~10%)尤为重要。
% 起点: A 就近 / B 次近 / C 邻域远端 / D 奇偶交错 / E,F 固定种子随机
n = numel(cands);
rs = RandStream('mt19937ar', 'Seed', 20240517);   % 局部随机流(确定性, 不动全局RNG)
nStarts = 6;
best = arrayfun(@(k) cands{k}(1), 1:n);       % 默认起点: 就近标称值
bestErr = metric(best);
for st = 1:nStarts
    switch st
        case 1                                    % 起点 A: 就近标称值
            start = arrayfun(@(k) cands{k}(1), 1:n);
        case 2                                    % 起点 B: 次近值(交替取整)
            start = arrayfun(@(k) cands{k}(min(2, numel(cands{k}))), 1:n);
        case 3                                    % 起点 C: 邻域远端
            start = arrayfun(@(k) cands{k}(end), 1:n);
        case 4                                    % 起点 D: 奇偶交错取整
            start = arrayfun(@(k) cands{k}(1 + mod(k, 2)), 1:n);
        otherwise                                 % 起点 E/F: 随机扰动
            start = arrayfun(@(k) cands{k}(randi(rs, numel(cands{k}))), 1:n);
    end
    [cur, curErr] = coordinateDescent(cands, metric, start);
    if curErr < bestErr
        best = cur; bestErr = curErr;
    end
end
end

function [best, bestErr] = coordinateDescent(cands, metric, start)
% 单起点坐标下降: 每次固定其余元件, 单元件取最优候选
n = numel(cands);
best = start;
bestErr = metric(best);
for it = 1:8
    improved = false;
    for k = 1:n
        for j = 1:numel(cands{k})
            trial = best;
            trial(k) = cands{k}(j);
            e = metric(trial);
            if e < bestErr - 1e-12
                best = trial; bestErr = e; improved = true;
            end
        end
    end
    if ~improved, break; end
end
end

% ----------------------------------------------------------------------
function f3 = cutoffFreq(fc, type, vals, Rs, RL, f3hint, gL)
% 求 |H(f)| 相对"设计通带峰值"下降 3.01 dB 处的频率。
% 参考电平取理想通带峰值 0.5*sqrt(gL) 而非实测峰值: 标称元件会改变
% 纹波峰值, 若以实测峰值为参考, -3dB 电平随之漂移, 偏差数字失真。
% 扫描全部过零点并选取最接近理想 -3dB(f3hint) 者, 避免误捕通带凹点。
lvl = 0.5 * sqrt(gL) / sqrt(2);                   % 设计 -3dB 电平
if strcmpi(type, 'lowpass')
    fscan = logspace(log10(fc) - 0.5, log10(fc) + 2, 601);% 过渡带扫描
else
    fscan = logspace(log10(fc) + 2, log10(fc) - 2, 601);  % 过渡带扫描
                                                          % (n=1 时 -3dB 可低至 eps*fc)
end
H  = abs(ladderResponse(fscan, type, vals, Rs, RL));
d  = H - lvl;
idx = find(d(1:end-1) .* d(2:end) < 0);
if isempty(idx)
    [~, i2] = min(abs(d));
    f3 = fscan(i2);
    warning('Chebyshev:cutoff', '未能稳定定位 -3dB 点, 返回近似值。');
    return;
end
f3c = zeros(size(idx));
dfun = @(lf) abs(ladderResponse(10^lf, type, vals, Rs, RL)) - lvl;
for t = 1:numel(idx)
    i = idx(t);
    try
        f3c(t) = 10^fzero(dfun, [log10(fscan(i)), log10(fscan(i + 1))]);
    catch
        f3c(t) = sqrt(fscan(i) * fscan(i + 1));           % 中点近似
    end
end
[~, bi] = min(abs(f3c - f3hint));
f3 = f3c(bi);
end

% ----------------------------------------------------------------------
function s = engstr(v, unit)
% 工程单位格式化, 如 7.958e-9 H -> '7.96 nH'
if ~isfinite(v), s = '---'; return; end
if v == 0, s = sprintf('0 %s', unit); return; end
pref = {'a','f','p','n','u','m','','k','M','G','T','P','E'};
e  = floor(log10(abs(v)));
ie = round(e / 3) + 7;
ie = min(max(ie, 1), numel(pref));
sc = 10^(3 * (ie - 7));
while abs(v) / sc >= 1000 && ie < numel(pref)
    ie = ie + 1; sc = 10^(3 * (ie - 7));
end
while abs(v) / sc < 1 && ie > 1
    ie = ie - 1; sc = 10^(3 * (ie - 7));
end
s = sprintf('%.3g %s%s', v / sc, pref{ie}, unit);
end

% ----------------------------------------------------------------------
function plotResponse(fc, n, type, eps_, gL, Rs, RL, idealVals, nominalVals)
% 绘制理想与标称幅频响应对比
f = fc * logspace(-1.5, 1.5, 801);
HdBnom   = 20 * log10(max(abs(ladderResponse(f, type, nominalVals, Rs, RL)), 1e-300));
HdBideal = idealChebyshevDB(f, fc, n, eps_, gL, type);
if fc >= 1e6
    u = 1e6; ul = 'MHz';
elseif fc >= 1e3
    u = 1e3; ul = 'kHz';
else
    u = 1;   ul = 'Hz';
end
figure('Name', sprintf('Chebyshev %d-order %s (fc=%s, ripple=%g dB)', n, type, engstr(fc, 'Hz'), 10*log10(1+eps_^2)), ...
       'NumberTitle', 'off');
semilogx(f/u, HdBideal, 'k--', 'LineWidth', 1.1);
hold on; grid on;
semilogx(f/u, HdBnom, 'b-', 'LineWidth', 1.6);
yl = ylim;
line([fc fc]/u, yl, 'Color', [0.7 0.3 0.3], 'LineStyle', ':', 'LineWidth', 1.2);
xlabel(sprintf('频率 (%s)', ul));
ylabel('|H| (dB)');
title(sprintf('切比雪夫 %d 阶%s 幅频响应 (纹波 %g dB)', n, type, 10*log10(1+eps_^2)));
legend('理想', '标称', 'fc', 'Location', 'southwest');
hold off;
end

% ----------------------------------------------------------------------
function b = tf(v)
% 逻辑型选项解析: 接受 logical / 数值 / 'true' 'on' 等字符串
if islogical(v)
    b = v;
elseif isnumeric(v)
    b = v ~= 0;
elseif ischar(v) || isstring(v)
    s = lower(char(v));
    b = strcmp(s, 'true') || strcmp(s, 'on') || strcmp(s, 'yes') || strcmp(s, '1');
else
    error('Chebyshev:tf', '逻辑型选项取值无效。');
end
end

function s = ternaryStr(cond, ifTrue, ifFalse)
% 简单三元字符串辅助
if cond, s = ifTrue; else, s = ifFalse; end
end
