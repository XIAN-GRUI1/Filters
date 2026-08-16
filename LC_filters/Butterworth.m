function [elements, info] = Butterworth(n, fc, type, varargin)
%BUTTERWORTH 基于现代网络函数综合法的巴特沃斯 LC 滤波器设计
%
%   [elements, info] = Butterworth(n, fc, type)
%   [elements, info] = Butterworth(n, fc, type, 'Name', Value)
%   Butterworth('demo')          % 运行演示(含自检)
%   Butterworth()                % 显示本帮助
%
% 功能
% ----
%   1. 基于"现代网络函数综合法"(Darlington 综合: 功率传输函数 ->
%      反射系数 -> 输入阻抗 -> 连分式展开), 求解理论上任意阶数 n、
%      任意截止频率 fc 的低通/高通巴特沃斯滤波器的 LC 元件精确值
%      (仅使用 MATLAB 基础功能, 不依赖任何工具箱);
%   2. 元件值就近匹配 IEC 60063 E96 标称值(±1% 精密系列), 并默认在
%      E96 邻域内做穷举/坐标下降搜索, 使幅频响应与理想巴特沃斯响应
%      偏差最小, 即"在尽可能保证性能的前提下完成标称值匹配"。
%
% 综合原理
% --------
%   归一化原型(|S21|^2 = 1/(1 + w^(2n)), 源/负载各 1 欧姆):
%     Bn(s)   : 巴特沃斯多项式(由左半平面极点构造, Bn(0) = 1);
%     S11(s)  : 反射系数 = s^n / Bn(s);
%     Zin(s)  : 输入阻抗 = (1 + S11) / (1 - S11);
%     g1..gn  : 对 Zin(s) 作连分式展开(交替提取 s=inf 处的串联感抗与
%               并联容抗), 得到原型元件值; g(n+1) 为负载电阻(1 欧姆)。
%   结果与经典闭式解 gk = 2*sin((2k-1)*pi/(2n)), g(n+1)=1 一致,
%   程序自动交叉校验; 由于双精度多项式连分式在 n>10 时数值病态,
%   n>10 自动采用同一综合结果的闭式解(数值精确, 阶数理论上任意)。
%   经阻抗标度(×Z0)与频率标度(×wc = 2*pi*fc)去归一化得到实际 L/C;
%   高通由低通原型经 s -> 1/s 变换得到(串联电容 / 并联电感)。
%
% 输入
% ----
%   n    正整数, 滤波器阶数(理论上任意; 连分式展开精确到 n=10,
%         更大阶数自动采用等效闭式解, 数值精确)
%   fc   截止频率, 单位 Hz(|H(fc)| 相对通带下降 3.01 dB)
%   type 'lowpass'|'lp' 或 'highpass'|'hp'(不区分大小写)
%
% 可选名值对
% ----------
%   'Z0'        源与负载端接阻抗(欧姆), 默认 50
%   'E96'       true/false, 是否进行 E96 标称值匹配, 默认 true
%   'Range'     标称值允许范围 [min, max](SI 单位), 默认 [1e-15, 1e6]
%   'Optimize'  true/false, 是否在标称值邻域内优化响应, 默认 true
%   'Search'    'auto'|'exhaustive'|'greedy'|'none'
%               默认 'auto': n<=10 穷举全部组合, 更大阶数用坐标下降
%   'Steps'     每个元件搜索的 E96 邻域档数(±Steps), 默认 1
%   'Metric'    'rms'|'max', 优化目标函数, 默认 'rms'
%   'Plot'      true/false, 绘制理想与标称幅频响应对比, 默认 false
%   'Verbose'   true/false, 是否在命令行打印结果汇总, 默认 true
%
% 输出
% ----
%   elements  table, 每行一个元件: 名称/类型/位置/理想值/标称值/偏差%
%   info      struct, 原型 g 值、极点、-3dB 频率、响应误差等详细结果
%
% 示例
% ----
%   Butterworth('demo');
%   el = Butterworth(5, 10e6, 'lowpass');            % 5 阶 10MHz 低通
%   el = Butterworth(3, 100e3, 'highpass', 'Z0', 75);% 3 阶 100kHz 高通
%   [el, in] = Butterworth(7, 1e9, 'lp', 'Plot', true);
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
    help Butterworth;
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
        otherwise, error('Butterworth:badType', 'type 必须是 lowpass 或 highpass。');
    end
else
    error('Butterworth:badType', 'type 必须是 lowpass 或 highpass。');
end

% ======================================================================
% 2. 可选参数
% ======================================================================
opts = struct('Z0', 50, 'E96', true, 'Range', [1e-15, 1e6], ...
              'Optimize', true, 'Search', 'auto', 'Steps', 1, ...
              'Metric', 'rms', 'Plot', false, 'Verbose', true);
optNames = {'Z0','E96','Range','Optimize','Search','Steps','Metric','Plot','Verbose'};
k = 1;
while k <= numel(varargin)
    nm = validatestring(varargin{k}, optNames, mfilename);
    k = k + 1;
    if k > numel(varargin)
        error('Butterworth:opt', '名值对 %s 缺少取值。', nm);
    end
    v = varargin{k};
    k = k + 1;
    switch nm
        case 'Z0'
            validateattributes(v, {'numeric'}, {'scalar','positive'}, mfilename, 'Z0');
            opts.Z0 = v;
        case 'E96'
            opts.E96 = tf(v);
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
% 3. 原型综合: 巴特沃斯多项式 -> 反射系数 -> 输入阻抗 -> 连分式展开
% ======================================================================
[g, poles, B, method] = butterworthPrototype(n);
g0 = 1;                                       % 归一化源电阻

% ======================================================================
% 4. 阻抗/频率去归一化 -> 实际 L/C 元件值
% ======================================================================
[idealVals, elemType, elemPos, units] = scaleElements(g, type, fc, opts.Z0);

% ======================================================================
% 5. 优化用频率网格与理想响应(用于 E96 匹配的保性能搜索)
% ======================================================================
fgrid = fc * logspace(-0.6, 0.6, 41);
HidB  = idealButterworthDB(fgrid, fc, n, type);

% ======================================================================
% 6. E96 标称值匹配(就近取值 + 邻域搜索保性能)
% ======================================================================
nominalVals = idealVals;
searchInfo  = struct('method', 'none', 'error', 0);
if opts.E96
    nominalVals = arrayfun(@(x) e96Nearest(x, opts.Range), idealVals);
    if any(isnan(nominalVals))
        warning('Butterworth:range', '部分元件值超出标称值范围, 相关元件保留理想值。');
        nominalVals(isnan(nominalVals)) = idealVals(isnan(nominalVals));
    end
    if opts.Optimize && ~strcmpi(opts.Search, 'none')
        [nominalVals, searchInfo] = optimizeE96(nominalVals, idealVals, fgrid, HidB, type, opts.Z0, opts);
    end
end

% ======================================================================
% 7. -3dB 截止频率与响应误差评估
% ======================================================================
f3ideal = cutoffFreq(fc, type, idealVals, opts.Z0);
f3nom   = cutoffFreq(fc, type, nominalVals, opts.Z0);
fchk    = fc * logspace(-1.2, 1.2, 201);
HidBchk = idealButterworthDB(fchk, fc, n, type);
rmsErr  = responseError(nominalVals, fchk, HidBchk, type, opts.Z0, 'rms');
maxErr  = responseError(nominalVals, fchk, HidBchk, type, opts.Z0, 'max');

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
              [opts.Z0; opts.Z0], {'ohm';'ohm'}, [opts.Z0; opts.Z0], ...
              {'ohm';'ohm'}, [0; 0], ...
              'VariableNames', {'Element','Type','Position','Ideal','Unit', ...
                                'Nominal','UnitN','DevPct'});
Rrows.IdealText   = {engstr(opts.Z0,'ohm'); engstr(opts.Z0,'ohm')};
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
info.g0 = g0;
info.g  = g;                                   % g(1..n+1), g(n+1) 为负载
info.poles = poles;
info.B  = B;
info.method = method;
info.idealValues   = idealVals;
info.nominalValues = nominalVals;
info.deviationPct  = dev;
info.f3dB_ideal    = f3ideal;
info.f3dB_nominal  = f3nom;
info.fcDeviationPct = (f3nom - fc) / fc * 100;
info.responseError_rms_dB = rmsErr;
info.responseError_max_dB = maxErr;
info.searchMethod = searchInfo.method;
info.searchError  = searchInfo.error;

% ======================================================================
% 10. 绘图(可选)
% ======================================================================
if opts.Plot
    plotResponse(fc, n, type, opts.Z0, idealVals, nominalVals);
end

% ======================================================================
% 11. 命令行汇总(可选)
% ======================================================================
if opts.Verbose
    fprintf('\n===== 巴特沃斯 %d 阶%s, fc = %s, Z0 = %g ohm =====\n', ...
            n, type, engstr(fc, 'Hz'), opts.Z0);
    disp(elements);
    fprintf('综合方法 : %s\n', method);
    fprintf('E96 匹配 : %s', mat2str(opts.E96));
    if opts.E96
        fprintf('  (搜索: %s, 目标: %s)', searchInfo.method, opts.Metric);
    end
    fprintf('\n');
    fprintf('理想 -3dB 频率 : %s\n', engstr(f3ideal, 'Hz'));
    fprintf('标称 -3dB 频率 : %s (相对 fc 偏差 %+.3f%%)\n', engstr(f3nom, 'Hz'), info.fcDeviationPct);
    fprintf('响应 RMS 误差  : %.4f dB, 最大误差 %.4f dB\n', rmsErr, maxErr);
    fprintf('============================================================\n');
end
end

% ======================================================================
% 以下为局部函数
% ======================================================================

function demo()
% 演示与自检
fprintf('\n================ 巴特沃斯 LC 滤波器设计演示 ================\n');
fprintf('\n--- 示例 1: 5 阶低通, fc = 10 MHz, Z0 = 50 ohm ---\n');
[el1, in1] = Butterworth(5, 10e6, 'lowpass', 'Plot', true);
disp(el1);
fprintf('综合方法 : %s\n', in1.method);
fprintf('-3dB     : %s | RMS 误差: %.4f dB\n', engstr(in1.f3dB_nominal, 'Hz'), in1.responseError_rms_dB);

fprintf('\n--- 示例 2: 4 阶高通, fc = 100 kHz, Z0 = 75 ohm ---\n');
[el2, in2] = Butterworth(4, 100e3, 'highpass', 'Z0', 75, 'Plot', true);
disp(el2);
fprintf('综合方法 : %s\n', in2.method);
fprintf('-3dB     : %s | RMS 误差: %.4f dB\n', engstr(in2.f3dB_nominal, 'Hz'), in2.responseError_rms_dB);

fprintf('\n--- 自检: 理想元件值下, 各阶 -3dB 频率应精确等于 fc ---\n');
worst = 0;
for nn = 1:8
    [~, in] = Butterworth(nn, 1e6, 'lowpass', 'E96', false, 'Verbose', false);
    e = abs(in.f3dB_ideal - 1e6) / 1e6;
    worst = max(worst, e);
    fprintf('  n = %d: f3dB = %.12g Hz, 相对误差 %.3g\n', nn, in.f3dB_ideal, e);
end
fprintf('  自检结论: 最大相对误差 %.3g (应远小于 1e-6)\n', worst);
fprintf('============================================================\n');
end

% ----------------------------------------------------------------------
function [g, poles, B, method] = butterworthPrototype(n)
% 归一化巴特沃斯原型综合(网络函数综合法: 连分式展开)
method = '现代网络函数综合法 (反射系数->输入阻抗->连分式展开)';
kk = (1:n).';
poles = exp(1j * pi * (2*kk - 1) / (2*n) + 1j * pi / 2);  % 左半平面极点
B = poly(poles);
B = B / B(end);                                           % 首一且 B(0)=1

% 闭式解(综合结果的解析形式), 用于交叉校验与大阶数回退
gClosed = [2 * sin((2*kk - 1) * pi / (2*n)); 1].';        % g1..g(n+1)

if n <= 10
    Znum = B;  Znum(1) = 2 * Znum(1);                     % B + s^n
    Zden = B;  Zden(1) = 0;                               % B - s^n (首项消去)
    g = ladderExtract(Znum, Zden, n);
    rel = max(abs(g - gClosed)) / max(abs(gClosed));
    if rel > 1e-6 * max(1, n / 10)
        warning('Butterworth:cf', '连分式展开数值误差较大 (rel=%.2e), 改用闭式解。', rel);
        g = gClosed;
        method = [method, ' (数值误差回退: 采用闭式解)'];
    end
else
    g = gClosed;
    method = [method, ' (n>10 双精度连分式病态, 采用等效闭式解)'];
end
g = real(g);                                              % 消除复数舍入噪声
end

% ----------------------------------------------------------------------
function g = ladderExtract(Znum, Zden, n)
% 对输入阻抗 Zin(s) = Znum/Zden 作连分式展开:
% 交替提取 s->inf 处的极点(串联元件 g奇 / 并联元件 g偶), 余项为负载电阻。
% 数值稳定策略: 多项式阶数完全由结构确定, 每次减法后余式阶数精确下降 2
% (前两项系数数学上为 0), 末步下降 1 —— 直接按阶数截断, 无需容差剥离,
% 舍入噪声因此永远无法进入主导系数。
g = zeros(1, n + 1);
N = Znum(:).';                        % deg = n
D = Zden(:).';                        % 首项已显式置零
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
% 原型 g(1..n) -> 实际 L/C (SI 单位), g(n+1) 为负载电阻
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
function E = e96Mantissas()
% IEC 60063 E96 系列(1.00 ~ 9.76, 三有效数字, 容差 ±1%)
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

function v = e96Nearest(x, range)
% 就近 E96 标称值; 超出 range 返回 NaN
E = e96Mantissas();
e = floor(log10(x));
m = x / 10^e;
[~, ci] = min(abs(E - m));
v = E(ci) * 10^e;
if v < range(1) || v > range(2), v = NaN; end
end

function vals = e96Candidates(x, steps, range)
% 目标值附近的 E96 候选集(±steps 档), 按与理想值距离升序排列
E = e96Mantissas();
e = floor(log10(x));
m = x / 10^e;
[~, ci] = min(abs(E - m));
idxs = (ci - steps) : (ci + steps);
mm = zeros(size(idxs));
dd = zeros(size(idxs));
for t = 1:numel(idxs)
    j = idxs(t);
    if j < 1
        mm(t) = E(j + 96);  dd(t) = -1;
    elseif j > 96
        mm(t) = E(j - 96);  dd(t) = 1;
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
function H = ladderResponse(f, type, vals, Z0)
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
H = Z0 ./ (A .* Z0 + Bm + Cm .* Z0 .* Z0 + Dm .* Z0);
H = H(:);
end

% ----------------------------------------------------------------------
function HdB = idealButterworthDB(f, fc, n, type)
% 理想巴特沃斯幅频响应(两端等阻抗端接, 通带损耗 6.02 dB)
x = (f(:) / fc) .^ (2 * n);
if strcmpi(type, 'lowpass')
    H = 0.5 ./ sqrt(1 + x);
else
    H = 0.5 .* sqrt(x ./ (1 + x));
end
HdB = 20 * log10(max(abs(H), 1e-300));
end

% ----------------------------------------------------------------------
function err = responseError(vals, fgrid, HidB, type, Z0, metricName)
% 标称元件响应与理想响应的偏差(dB 域)
H = ladderResponse(fgrid, type, vals, Z0);
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
function [bestVals, out] = optimizeE96(startVals, idealVals, fgrid, HidB, type, Z0, opts)
% 在 E96 候选集内搜索使响应误差最小的元件组合
n = numel(idealVals);
cands = cell(1, n);
nCand = zeros(1, n);
for k = 1:n
    c = e96Candidates(idealVals(k), opts.Steps, opts.Range);
    if isempty(c)
        c = startVals(k);
        warning('Butterworth:noCand', '元件 %d 无可用 E96 候选, 保留就近值。', k);
    end
    cands{k} = c;
    nCand(k) = numel(c);
end
metric = @(v) responseError(v, fgrid, HidB, type, Z0, opts.Metric);
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
        warning('Butterworth:combos', '组合数 %g 过大, 改用坐标下降搜索。', total);
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
% 坐标下降( Gauss-Seidel ): 每次固定其余元件, 单元件取最优候选
n = numel(cands);
best = arrayfun(@(k) cands{k}(1), 1:n);       % 起点: 就近标称值
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
function f3 = cutoffFreq(fc, type, vals, Z0)
% 求 |H(f)| 相对通带电平下降 3.01 dB 处的频率
if strcmpi(type, 'lowpass')
    fpb   = fc * 1e-6;                                    % 低频通带采样
    fscan = logspace(log10(fc) - 3, log10(fc) + 3, 401);  % 通带->阻带
else
    fpb   = fc * 1e6;                                     % 高频通带采样
    fscan = logspace(log10(fc) + 3, log10(fc) - 3, 401);  % 通带->阻带
end
Hpb = abs(ladderResponse(fpb, type, vals, Z0));
lvl = Hpb / sqrt(2);
H  = abs(ladderResponse(fscan, type, vals, Z0));
d  = H - lvl;
idx = find(d(1:end-1) > 0 & d(2:end) <= 0, 1, 'first');
if isempty(idx)
    [~, i2] = min(abs(d));
    f3 = fscan(i2);
    warning('Butterworth:cutoff', '未能稳定定位 -3dB 点, 返回近似值。');
    return;
end
try
    dfun = @(lf) abs(ladderResponse(10^lf, type, vals, Z0)) - lvl;
    f3 = 10^fzero(dfun, [log10(fscan(idx)), log10(fscan(idx + 1))]);
catch
    f3 = sqrt(fscan(idx) * fscan(idx + 1));               % 中点近似
end
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
function plotResponse(fc, n, type, Z0, idealVals, nominalVals)
% 绘制理想与标称幅频响应对比
f = fc * logspace(-1.5, 1.5, 801);
HdBnom   = 20 * log10(max(abs(ladderResponse(f, type, nominalVals, Z0)), 1e-300));
HdBideal = idealButterworthDB(f, fc, n, type);
if fc >= 1e6
    u = 1e6; ul = 'MHz';
elseif fc >= 1e3
    u = 1e3; ul = 'kHz';
else
    u = 1;   ul = 'Hz';
end
figure('Name', sprintf('Butterworth %d-order %s (fc=%s)', n, type, engstr(fc, 'Hz')), ...
       'NumberTitle', 'off');
semilogx(f/u, HdBideal, 'k--', 'LineWidth', 1.1);
hold on; grid on;
semilogx(f/u, HdBnom, 'b-', 'LineWidth', 1.6);
yl = ylim;
line([fc fc]/u, yl, 'Color', [0.7 0.3 0.3], 'LineStyle', ':', 'LineWidth', 1.2);
xlabel(sprintf('频率 (%s)', ul));
ylabel('|H| (dB)');
title(sprintf('巴特沃斯 %d 阶%s 幅频响应 (Z0 = %g ohm)', n, type, Z0));
legend('理想', 'E96 标称', 'fc', 'Location', 'southwest');
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
    error('Butterworth:tf', '逻辑型选项取值无效。');
end
end
