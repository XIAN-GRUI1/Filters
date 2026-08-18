# Sallen-Key 有源滤波器自动设计工具箱 (MATLAB)

基于 **Sallen-Key（SK）拓扑** 的模拟有源滤波器设计程序：
理论上支持**任意阶数**、**任意截止频率**、**多种滤波器类型**与
**低通/高通**，电阻按 **E24/E96 标称值**匹配，电容按 **E24** 标称值匹配，
并可针对运放有限的 **增益带宽积 (GBW)** 进行**预失真补偿**，
在给定 GBW 的条件下直接求得经补偿的元件参数。

## 文件清单

| 文件 | 说明 |
|---|---|
| `skDesign.m` | **主函数**：一键设计滤波器（见下方用法） |
| `skProto.m` | 归一化低通原型：Butterworth / Chebyshev-I / Chebyshev-II / Bessel / Cauer（自含实现，不依赖信号处理工具箱） |
| `skBiquad.m` | 节级数学：SK 双二阶元件解、含有限 GBW 的实际传递系数、主导极点提取、状态变量双二阶（传输零点节）、一阶 RC 节 |
| `skOpt.m` | 每节 E 系列元件搜索（含 GBW 预失真优化） |
| `skResponse.m` | 整滤波器频响计算（实际 + 理想） |
| `skSpecs.m` | 实测指标：−3 dB 点 / 通带纹波 / 阻带衰减 |
| `skReport.m` | 文本设计报告 |
| `skPlot.m` | 幅频/相频绘图 |
| `skESeries.m`, `skNearestE.m`, `skSeriesGrid.m`, `skFmtSI.m` | E 系列工具 |
| `skDemo.m` | 演示脚本 |
| `sktest.m` | 自动测试（10 组用例） |

## 快速开始

```matlab
% 6 阶 Butterworth 低通，fc = 1 kHz，理想运放
F = skDesign('Type','butter','Order',6,'Fc',1e3, 'Plot',true);

% 5 阶 Chebyshev-I 高通，1 dB 纹波，运放 GBW = 200 kHz（自动预失真）
F = skDesign('Type','cheby1','Order',5,'Fc',1e3,'Rp',1, ...
             'PassType','highpass','GBW',200e3, 'Plot',true);

% 5 阶 Cauer（椭圆）低通，1 dB / 40 dB
F = skDesign('Type','cauer','Order',5,'Fc',1e3,'Rp',1,'Rs',40,'Plot',true);

% 4 阶 Bessel 低通（群延迟最平坦）
F = skDesign('Type','bessel','Order',4,'Fc',1e3);
```

运行后终端输出完整设计报告（每节元件值、实现精度、GBW 裕量、实测指标），
返回结构体 `F` 内含 `sections`（各节元件与实现值）、`specs`（实测指标）、
`prototype`（原型）等。

## 主要选项（名称-值对）

| 选项 | 默认 | 说明 |
|---|---|---|
| `Type` | `'butter'` | `butter` / `cheby1` / `cheby2` / `bessel` / `cauer` |
| `Order` | 6 | 阶数（≥1，理论上任意） |
| `Fc` | 1000 | 截止/通带边沿频率（Hz）。butter/bessel 为 −3 dB 点；cheby1/cauer 为通带边沿；cheby2 为阻带边沿 |
| `PassType` | `'lowpass'` | `lowpass` / `highpass` |
| `Rp` | 1 | cheby1/cauer 通带纹波（dB） |
| `Rs` | 40 | cheby2/cauer 阻带衰减（dB） |
| `Fstop` | — | cauer 的阻带边沿频率（Hz）；给出时改求可实现的最大 Rs |
| `Rser` | `'E96'` | 电阻系列：`E96` / `E24` |
| `Cser` | `'E24'` | 电容系列 |
| `GBW` | `Inf` | 运放增益带宽积（Hz）；有限值时自动预失真 |
| `Rmin/Rmax` | 100 / 1e6 | 电阻搜索范围（Ω） |
| `Cmin/Cmax` | 1e-12 / 1e-6 | 电容搜索范围（F） |
| `Plot` | false | 绘制幅频/相频曲线 |
| `Verbose` | true | 打印报告 |

## 拓扑说明（重要）

* **全极点类型**（Butterworth、Chebyshev-I、Bessel）：
  全部采用**单位增益 Sallen-Key 二阶节**（奇数阶含一个一阶 RC + 缓冲节）。
  低通节取 R1=R2、C1/C2 = 4Q² 的经典等电阻形式（搜索中自动放宽为
  任意 E 系列组合）；高通为其对偶。
* **带传输零点类型**（Cauer、Chebyshev-II）：
  单位增益 Sallen-Key 的传递函数分子为常数，**无法实现有限频率传输零点**，
  因此含零点的二阶节改用教科书标准的 **3 运放状态变量双二阶**
  （Tow-Thomas 型，可精确实现 `(s²+ωz²)/(s²+(ω0/Q)s+ω0²)`，要求
  **ωz > ω0，即低通情形**）；奇数阶的实极点节仍用 Sallen-Key。
  所有节统一做 E 系列匹配与 GBW 预失真。
  **限制**：高通 Cauer/Chebyshev-II 需要零点低于极点的专用高通双二阶
  （例如 4 运放状态变量滤波器），本工具箱未包含，`skDesign` 会给出明确
  错误提示。这类滤波器请使用低通版本、全极点类型，或专用高通双二阶电路。

## GBW 预失真原理

运放按单极点积分器建模 `A(s) = GBW/s`（跟随器即 `V+ = Vout(1+s/GBW)`）。
对 SK 节，代入电路方程后实际分母为三阶
`1 + B1 s + B2 s² + B3 s³`（有限 GBW 引入第三个极点）；
元件搜索以"**实际主导极点对 (ω0', Q') 尽量等于目标 (ω0, Q)**"为准则，
在 E 系列网格上自动选出补偿后的元件。有限 GBW 会使 SK 节的 Q 升高
（Q 增强效应），预失真会自动降低设计 Q 以抵消。
对状态变量双二阶，则在目标频点集上最小化实际响应与理想响应的偏差。

一阶节给出预失真的解析解：设 `g = GBW/(2πfc)`，低通取
`y = [−(g²+1)+√((g²+1)²+4g²)]/2`，`RC = √y/ωc`，使含缓冲极点后的
−3 dB 点精确落在 fc。

## 设计报告中的裕量提示

* `GBW margin wt/(w0·Q)`：经验上 ≥ 10 时未预失真的 Q 误差约 <1%；
  预失真后对裕量的要求放宽，但极低的裕量（<3）仍会警告精度受限。
* 低通 SK 节要求 `C1 ≥ 4Q²C2`（等电阻形式）；Q 很高时电容比值大，
  元件搜索会自动挑选最小展开的 E 系列组合。

## 备注

* 全部原型公式自含实现且经数值验证（与 MATLAB `ellipap`/`cheby1` 等
  对比到机器精度）；不需要 Signal Processing Toolbox。
  Cauer 原型仅使用基础 MATLAB 的 `ellipke`/`ellipj`。
* 每节按 Q 升序级联（低 Q 在前）以获得较好的动态范围。
* 若要进一步提高有限 GBW 下的精度，可改用"增益 2"型 SK 节
  （本工具箱按单位增益实现，增益 2 版本留作扩展）。

## 在 DSH 沙箱中运行 MATLAB（重要）

本会话的 DSH 文件沙箱为 `workspace-write` 模式：只允许在会话工作区与
平台临时区写入。MATLAB 启动时要向 `%USERPROFILE%`（AppData、偏好目录、
临时目录）等处写入文件，因此**直接运行 `matlab` 会以
"File system inconsistency" / "CreateFile failed: Access is denied" 失败**。

解决办法：把 MATLAB 的全部可写位置重定向到工作区内的可写目录。
仓库里已备好两个启动包装脚本，直接使用即可（无需任何权限提升）：

```bat
:: 运行批处理命令（推荐）
run_matlab.cmd -batch "sktest"
run_matlab.cmd -batch "skDemo"
run_matlab.cmd -batch "F = skDesign('Type','cauer','Order',5,'Fc',1e3,'Rp',1,'Rs',40)"

:: 启动交互式 MATLAB（桌面界面）
run_matlab.cmd

:: PowerShell 版（若执行策略拦截脚本，可用
::   powershell -ExecutionPolicy Bypass -File run_matlab.ps1 -Batch "sktest"）
.\run_matlab.ps1 -Batch "sktest"
```

原理：包装脚本将 `USERPROFILE / HOMEDRIVE / HOMEPATH / APPDATA /
LOCALAPPDATA / MATLAB_PREFDIR / TEMP / TMP` 全部指向本目录下的
`.matlab_sandbox\`（位于可写的工作区内），然后调用 `matlab.exe` 并把
其余参数原样传递。脚本结束时自动恢复原环境变量，不影响其他程序。

> 注意：`.matlab_sandbox` 存放 MATLAB 的偏好与临时文件，可随时删除
> （删除后下次运行会重建，MATLAB 偏好会重置）。

### 在"其他文件夹"里设计滤波器

不需要把库文件复制到每个文件夹。分两种情况：

**1. 其他文件夹位于当前会话工作区之内（例如本目录下的子文件夹）**
直接引用本目录的包装脚本，并用 `-sd`（cmd 版）或 `-WorkDir`（ps1 版）
指定 MATLAB 的工作目录即可，无需复制任何文件：

```bat
:: 在子文件夹 my_filter_proj 里运行你的脚本
run_matlab.cmd -sd "D:\文档\MATLAB\Filters\Source_filters\my_filter_proj" -batch "my_design_script"
```

```powershell
.\run_matlab.ps1 -Batch "my_design_script" -WorkDir "D:\文档\MATLAB\Filters\Source_filters\my_filter_proj"
```

MATLAB 启动后会进入该文件夹（`-sd` 是 MATLAB 自带参数，包装脚本原样透传）。
若该文件夹里的脚本要用本库的 `skDesign` 等函数，在脚本开头加一句
`addpath('D:\文档\MATLAB\Filters\Source_filters')` 即可。

**2. 换了一个全新的 DSH 会话 / 工作区（文件夹不再是本工作区）**
此时才需要复制文件。只需把 `run_matlab.cmd` 拷到新工作区文件夹（与你的
滤波器代码放一起），再按需拷入 `sk*.m` 库文件。规则只有一条：
**`run_matlab.cmd` 所在文件夹必须是当前会话可写的工作区**（因为它要在
自己旁边创建 `.matlab_sandbox`）。复制后直接使用：

```bat
run_matlab.cmd -batch "你的脚本"
```

**3. 只想跑 MATLAB、不需要本滤波器库**
只需一个 `run_matlab.cmd`（拷到可写工作区即可），运行任意 `.m` 代码
都用它启动，无需 `sk*.m`。
