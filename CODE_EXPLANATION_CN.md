# porous_ehl_final_v1 代码逐段解读与控制方程（中文）

> 说明：你提供的程序本身已经带有大量模块级注释。这个文档进一步把“每个模块在做什么、对应哪条控制方程、离散成什么格式”串起来，便于你后续逐行排查与改模。

---

## 1. 物理问题与无量纲化

程序求解的是 **刚体球-多孔弹性层等温点接触** 的耦合问题，耦合链条是：

\[
P \rightarrow P^* \rightarrow U_z \rightarrow H \rightarrow P
\]

采用的核心无量纲变量：

\[
X=\frac{x}{b},\quad Y=\frac{y}{b},\quad Z=\frac{z}{b},\quad
P=\frac{p}{p_H},\quad P^*=\frac{p_p}{p_H},\quad
H=\frac{hR}{b^2},\quad U_z=\frac{u_zR}{b^2}
\]

程序里 `closeDimensionlessGroups_local` 闭合了 4 个核心群：

\[
\Lambda = \frac{12\eta_0 U_e R^2}{p_H b^3},\quad
\Pi_p = \frac{k_0 p_H R}{\eta_0 U_e b^2},\quad
\Gamma = \frac{G b}{R p_H},\quad
\Lambda_s = \frac{\lambda b}{R p_H}
\]

---

## 2. 主控方程（与代码一一对应）

## 2.1 二维界面 Reynolds 方程（含渗流源项）

代码在 `buildReynoldsSource_local` 和 `solveReynolds2D_local` 中使用：

\[
\nabla_{XY}\cdot\left(\Xi\nabla_{XY}P\right)=S
\]

其中

\[
\Xi = \frac{\bar\rho H^3}{\Lambda\bar\eta}
\]

源项

\[
S = \frac{\partial(\bar\rho H)}{\partial X}
- \Pi_p\,\bar\rho\,\kappa\,\frac{\partial P^*}{\partial Z}
\]

- 第一项是经典楔形（卷吸）项；
- 第二项是多孔界面法向渗流对膜内质量守恒的修正。

离散方式：
- 空间二阶中心型通量离散（面系数 `Ae Aw An As`）。
- 压力迭代：GS/SOR 平滑（`omegaP`）+ 多重网格 V-cycle。
- 空化：负压截断 `P<0 -> 0`。
- 边界：四周 `P=0`。

---

## 2.2 多孔层稳态渗流方程

`solveSeepage3D_local` 对应：

\[
\nabla\cdot(\kappa\nabla P^*)=0
\]

渗透-黏性比系数：

\[
\kappa = \frac{K}{\bar\eta_p},\quad K=\exp(M\hat{}\,E_v)
\]

边界条件：
- 顶面（界面）Dirichlet：\(P^*(X,Y,0)=P(X,Y)\)
- 底面 Dirichlet：\(P^*(X,Y,-H_p)=P_0\)
- 侧面 Neumann：\(\partial P^*/\partial n=0\)（代码里用镜像赋值实现）

离散方式：
- 3D 七点模板（x/y/z 六邻居）。
- 面系数采用 **谐均值**（`harmonicMean_local`）增强非均匀介质稳定性。
- 内层 SOR 松弛（`omegaSeep`）。

---

## 2.3 骨架法向位移简化方程（标量 Uz）

`solveSkeleton3D_local` 使用了仅 Uz 的标量形式（不是完整 3D Navier 位移向量）：

\[
\Gamma\left(\frac{\partial^2 U_z}{\partial X^2}+\frac{\partial^2 U_z}{\partial Y^2}\right)
+ (\Lambda_s+2\Gamma)\frac{\partial^2 U_z}{\partial Z^2}
= \alpha\frac{\partial P^*}{\partial Z}
\]

顶面牵引边界（以二阶单边差分离散）：

\[
(\Lambda_s+2\Gamma)\frac{\partial U_z}{\partial Z}\Big|_{Z=0}
= -\text{skelTopTractionScale}\,(1-\alpha)P
\]

代码离散成：

\[
U_z^{k=1}=\frac{4U_z^{k=2}-U_z^{k=3}+2\Delta z\,T_n/(\Lambda_s+2\Gamma)}{3}
\]

底面固支：\(U_z=0\)。
侧面 Neumann：\(\partial U_z/\partial n=0\)。

体应变（这里简化为竖向梯度主导）在 `computeEvFromUz_local` 中：

\[
E_v\approx\frac{\partial U_z}{\partial Z}
\]

---

## 2.4 膜厚方程（几何 + 弹性 + 多孔附加）

`updateFilmAndProps_local` 中：

\[
H = H_{00}+\frac{X^2+Y^2}{2}+\delta_{EHL}+\delta_{porous}
\]

其中：
- \(\delta_{EHL}\)：非局部弹性卷积主变形（A2）
- \(\delta_{porous}=-w_p\,U_z|_{Z=0}\)：多孔附加修正项（A1）

并施加下限：\(H\ge H_{min}\)。

---

## 2.5 非局部表面弹性变形（A2）

`prepareElasticInfluence_local` + `solveElasticSurfaceDeflection_local`：

\[
\delta_{EHL}(X,Y)=\iint \frac{C_e\,P(\xi,\eta)}{\sqrt{(X-\xi)^2+(Y-\eta)^2+\varepsilon^2}}\,d\xi d\eta
\]

- 用 FFT 卷积实现：`ifft2( fft2(P).*fft2(Kernel) )`
- 系数 `C_e` 可选：
  - `physical`: \(C_e=\frac{2p_HR}{\pi E' b}\), \(E'=E/(1-\nu^2)\)
  - `legacy`: \(2/\pi^2\)

---

## 3. 外层耦合算法（主程序 for itOuter）

每轮执行：
1. 解渗流 `P*`
2. 解骨架 `Uz`
3. 更新 `H, rho, eta, Xi, dP*/dZ`
4. 解 Reynolds 得到 `P`
5. 载荷平衡修正 `H00`（带阻尼 Secant）
6. 再次更新并重解 Reynolds，保证 `P-H-H00` 同步

收敛判据同时监控：
- `errP, errH, errPp, errUz0, errH00, errDelta`
- 以及载荷误差 `errW = |W-Wtarget|`

这相当于“强耦合外迭代 + 子方程内迭代”的分裂策略。

---

## 4. 载荷闭合（Secant）

`adjustH00_local` 中定义

\[
F(H_{00}) = W(H_{00}) - W_{target}
\]

若历史足够，采用阻尼割线：

\[
\Delta H_{00} = -\omega_s F_k\frac{H_{00,k}-H_{00,k-1}}{F_k-F_{k-1}}
\]

否则回退比例控制：

\[
\Delta H_{00} = \text{loadGain}\cdot F_k
\]

并做步长限制与区间裁剪（`maxH00Step`, `H00min/max`）。

另外，程序还有“冻结形态后的最终载荷闭合”：
- 固定 \(\delta_{EHL},\delta_{porous},P^*,Uz\) 的形状，
- 仅调 `H00` + 重解 Reynolds，
- 把终态载荷误差再压低。

---

## 5. 你最关心的几个数值细节（非常关键）

1. **`dPp_dZ` 的符号**：
   - 代码 `Z` 从 0 到负值（向下），`dZs<0`。
   - 二阶公式 `(-3P1+4P2-P3)/(2*dZs)` 自动带入几何方向，后续符号要与 `seepTerm`、顶面牵引符号整体核对。

2. **Uz 正负与膜厚耦合**：
   - 当前实现 `deltaPorous = -Uz_surf * weight`。
   - 即在当前符号约定下，受压区 `Uz<0` 会抬高 `deltaPorous`。

3. **对称化处理**：
   - Reynolds 后、渗流后、骨架后都做了 `Y` 镜像平均，抑制扫描顺序导致的横向偏置。

4. **物性更新顺序**：
   - 每次 `P/H/P*` 变化后都重算 `eta, rho, Xi, kappa`，这点非常重要，否则会出现“伪收敛”。

---

## 6. 我建议你下一步怎么做

- 先用 `baseline_ehl_shape` 验证“无孔 + A2 主弹性”是否产生合理中心抬压与出口收敛膜厚。
- 再切 `current_reference`，确认 `deltaPorous` 仍明显小于 `deltaEHL`（程序诊断里也在检查）。
- 最后开 `porous_diag_strong` 做敏感性扫描：`k0_dim, E_dim, alpha, skelTopTractionScale, porousDispWeight`。

---

如果你愿意，我下一条消息可以直接给你：
1) **“逐函数逐行”解释清单（按你这份代码原顺序，一段一段拆）**；
2) 一版更“论文可写”的**控制方程+边界条件+离散格式**整理（可直接放方法章节）。

---

## 7. 按你原始代码顺序的“逐函数近似逐行”讲解

> 约定：这里严格按你贴出的函数出现顺序讲解；每个函数都给出
> 1) 输入/输出；2) 循环里每个关键系数在物理上的含义；3) 数值稳定性关注点。

### 7.1 `porous_ehl_final_v1`（总控主程序）

**主流程对应你代码中的 Step A~E：**
1. 初始化参数、网格、影响核、场变量
2. 以 Hertz 压力作初值
3. 外层耦合：渗流 → 骨架 → 膜厚物性更新 → Reynolds → 载荷闭合
4. 额外“冻结形态后”只调 `H00` 做最终载荷闭合
5. 输出诊断 + 作图

**逐段要点：**
- `P_old/H_old/Pp_old/Uz0_old/H00_old/deltaEHL_old`：保存上一外层步，用于多指标收敛判断。
- Step A（渗流）若关闭，则强制 `Pp(:,:,1)=P`、底部为 `P0_bottom`，等价于“无多孔传递”。
- Step B（骨架）若关闭，`Ux/Uy/Uz/Ev` 全清零，避免残留历史场污染。
- Step C 更新 `H, eta, rho, Xi, dPp_dZ`，是耦合中最关键的数据同步点。
- Step D 解 Reynolds 获得新膜压。
- Step E 用载荷误差修正 `H00`，然后**再**更新一次膜厚并重解 Reynolds，确保 `P-H00` 一致。

**稳定性意义：**
- 你做了“Reynolds 二次同步”，这对避免 `H00` 更新后压力场滞后非常关键。
- 收敛判据不是只看压力，而是把 `H/Pp/Uz/H00/deltaEHL` 都纳入，这是强耦合系统应有做法。

---

### 7.2 `initParameters_local`

这是“数值实验总开关”。

- 网格：`Nx,Ny,Nz` 与计算域边界。
- 物理参数：`pH, eta0, rho0, Ue_dim, R_dim, b_dim, k0_dim, E_dim, nu_dim`。
- 由 `closeDimensionlessGroups_local` 闭合 `Lambda, Pi_p, Gamma, Lambda_s`。
- 数值保护：`Hmin/Kmin/Kmax/kappaFloor/kappaMax/XiMin/XiMax/sourceMax/...`。
- 迭代与松弛：`maxOuter/maxRe/maxSeep/maxSkel` 与 `omegaP/omegaSeep/skelRelax/deflRelax`。
- 多重网格参数：`mgLevel, mgPre, mgPost`。
- 工况预设：`applyCasePreset_local`。

**稳定性关注：**
- `Hmin`、`XiMin`、`kappaFloor` 是防止系数退化（0 或 NaN）的“底座”。
- `secantDamping`、`maxH00Step` 决定载荷闭合是否振荡。

---

### 7.3 `closeDimensionlessGroups_local`

先做参数合法性检查，再计算：
- `G_dim = E/(2(1+nu))`
- `lambda_dim = E*nu/((1+nu)(1-2nu))`
- 再映射为 `Gamma/Lambda_s` 等无量纲群。

**稳定性关注：**
- 对 `nu` 的检查非常必要；`nu→0.5` 时 Lamé 常数会发散。

---

### 7.4 `applyCasePreset_local`

三套工况通过同一个入口切换：
- `baseline_ehl_shape`
- `current_reference`
- `porous_diag_strong`

每个工况同时改：
1) 物理强度（`Ue_dim, k0_dim, E_dim, alpha`）
2) 耦合权重（`porousDispWeight`, `skelTopTractionScale`）
3) 数值策略（`maxOuter`, `secantDamping`, `maxH00Step`）

**稳定性关注：**
- 工况切换后你又调用了一次 `closeDimensionlessGroups_local`，避免参数改了但无量纲群没同步。

---

### 7.5 `applyDefaultSkeletonParams_local`

这是容错层：`isfield` 检查确保缺省参数可运行。

**工程意义：**
- 便于以后拆分配置文件或分支脚本，不会因缺字段崩溃。

---

### 7.6 `buildGrid2D_local` / `buildGrid3D_local`

- 2D：`X,Y` + `XX,YY`。
- 3D：`Z` 从 `0` 到 `-Hp`，所以 `dZs<0`，`dZ=abs(dZs)`。

**稳定性关注：**
- 你后续所有 `dP/dZ` 公式都必须和这个坐标方向一致（这是最常见符号坑）。

---

### 7.7 `prepareElasticInfluence_local`

- 先给出核系数 `elasticInfluenceCoeff`（physical 或 legacy）。
- 构造卷积核 `K ~ 1/sqrt(r^2+eps^2)`。
- `epsReg` 去奇异化。
- 预存 `fft2(K)`，并记录裁剪索引。

**稳定性关注：**
- `epsReg` 太小会放大中心奇异，太大则过度平滑；你现在取 `0.5*sqrt(dX^2+dY^2)` 较稳妥。

---

### 7.8 `initField2D_local` / `initField3D_local` / `initPressureHertz_local`

- 初始化所有场到可用状态，避免空数组参与运算。
- 初值压力采用截断 Hertz：`P = sqrt(1-X^2-Y^2)`（接触圆内）。

**稳定性关注：**
- 这个初值比全零更容易收敛到物理态，尤其对强耦合工况。

---

### 7.9 `solveReynolds2D_local`

每次内迭代流程：
1. `S = buildReynoldsSource_local(...)`
2. MG V-cycle 或 GS 平滑
3. 非法值清理、负压截断、边界置零
4. `err = max|P-Pold|` 判断收敛

结束后再做一次 `Y` 对称化。

**循环系数物理意义：**
- `Xi` 是“局部流动能力”（膜厚三次方 + 压力黏度密度效应）
- `S` 里 `convTerm` 给卷吸供给，`seepTerm` 给渗流吸/注通量。

**稳定性关注：**
- 你在每 sweep 后都做空化与边界重置，这能避免负压和边界漂移进入下一 sweep。

---

### 7.10 `buildReynoldsSource_local`

核心是：
- `convTerm = (rho*H)_i - (rho*H)_{i-1}` 的一阶迎风式写法
- `seepTerm = Pi_p * rho * kappa * dPp_dZ`
- `S = convTerm - seepTerm`
- 再做 `sourceMax` 截断

**稳定性关注：**
- `sourceMax` 对强渗流工况是必要的“保险丝”，防止单步爆源。

---

### 7.11 `reynoldsVcycle_local`

标准 V-cycle 结构：
1. 预平滑
2. 算残差 `R = S - A(U)`
3. 限制到粗网格（`restrictFW_local`）
4. 粗网格误差方程递归
5. 双线性延拓回细网格修正
6. 后平滑

**稳定性关注：**
- `Xic = max(restrictFW(Xi), XiMin)` 避免粗网格系数退化。

---

### 7.12 `reynoldsSmoothPressureGS_local` / `reynoldsSmoothErrorGS_local`

每个内点系数：
- `Ae,Aw,An,As`：四个方向“导通能力”
- `Ap = Ae+Aw+An+As`
- `rhs = 邻点贡献 - S`
- `pNew = rhs/Ap`
- SOR 更新 `P = (1-omega)P + omega*pNew`

**物理解释：**
- `Ae` 大意味着该方向压力扩散更强（本质受 `Xi` 支配）。

**稳定性关注：**
- `Ap<=0` 时兜底到 1.0，防止除零。

---

### 7.13 `reynoldsApplyA_local`

计算离散算子 `A(P)`，用于多重网格残差。

---

### 7.14 `restrictFW_local` / `prolongBilinear_local`

- 限制：full-weighting（中心 1/4，十字 1/8，角点 1/16）
- 延拓：先注入再线性插值

**稳定性关注：**
- 延拓后边界清零，和 Reynolds 边界条件保持一致。

---

### 7.15 `computeLoad2D_local`

双重 `trapz` 做 2D 载荷积分：
\[
W=\iint P\,dX\,dY
\]

---

### 7.16 `adjustH00_local`

- 构造 `F = W-Wtarget`
- 若割线条件成立，用阻尼 Secant
- 否则比例反馈 `dH = loadGain*F`
- 步长限制 + 区间裁剪
- 更新历史对 `(H_prev,F_prev)`

**稳定性关注：**
- `secantTolDen` 防止分母接近 0 导致大跳步。

---

### 7.17 `plotResults_local`

输出 3D 面图、等值线、中心线剖面，并打印 `(0,0)` 附近的关键量。

**工程意义：**
- 这一步是“肉眼验算”：看压力峰位、膜厚谷位、对称性与中心线形态是否合理。

---

### 7.18 `updateFilmAndProps_local`

这是全耦合中最重要的“状态刷新器”。

1. 主弹性变形：
   - `deltaEHL_raw = weight * elasticConvolution(P+)`
   - 截断到 `[0, deltaEHLMax]`
   - 用 `deflRelax` 松弛更新
2. 多孔附加项：`deltaPorous = porousDispWeight * (-Uz_surf)`
3. 膜厚：`H = H00 + Hgeom + deltaEHL + deltaPorous`
4. 物性：`eta(P), rho(P), etaP(Pp)`
5. 渗流系数：`K=exp(Mhat*Ev)`，`kappa=K/etaP`
6. Reynolds 系数：`Xi=rho*H^3/(Lambda*eta)`
7. 界面法向梯度：二阶或一阶公式算 `dPp_dZ`

**稳定性关注：**
- 多处 `~isfinite` 清理 + 上下限夹逼，是这个程序能跑强耦合的关键。

---

### 7.19 `updateFilmAndPropsFrozenShape_local`

和上一个函数区别：
- `deltaEHL/deltaPorous/Pp/Uz/Ev/K` 全冻结；
- 只让 `H00` 与当前 `P` 驱动 `H,eta,rho,Xi` 更新。

**用途：**
- 用于最终载荷精闭合，避免形态漂移掩盖 `H00` 调节效果。

---

### 7.20 `solveElasticSurfaceDeflection_local`

- 使用预存 FFT 核做卷积得到 `deltaEHL`
- 裁剪回原网格
- 去 NaN/Inf，负值截断为 0

**稳定性关注：**
- 负值截断使弹性附加项保持“开膜”贡献，避免与 A1 项符号冲突。

---

### 7.21 `solveSeepage3D_local`

核心迭代：
1. 顶面/底面 Dirichlet，侧面 Neumann
2. 六个方向系数 `aE..aB`（由谐均值 `kappa` 给出）
3. `aP = sum(a*)`，`pNew = rhs/aP`
4. SOR 更新
5. 收敛后做 Y 对称化，再强制边界

**每个离散系数的物理意义：**
- `aE=ke/dX^2`：沿 +X 方向渗流导通能力
- `aT=kt/dZ^2`：沿上方（靠界面）方向导通能力
- 系数越大，邻点对本点孔压影响越强

**稳定性关注：**
- `seepUpdateKappa` 允许迭代中周期重算 `kappa(Pp)`，使非线性渗流更一致。

---

### 7.22 `solveSkeleton3D_local`

核心迭代：
1. 侧面 Neumann
2. 顶面牵引边界离散更新 `Uz(:,:,1)`
3. 底面 `Uz(:,:,Nz)=0`
4. 内点 7 点格式解 `Uz`
5. 迭代收敛后 Y 对称化
6. 由 `Uz` 反算 `Ev`

**离散系数含义：**
- `aE,aW,aN,aS = G/dx^2...` 对应横向剪切刚度贡献
- `aT,aB = (Lambda_s+2G)/dz^2` 对应竖向法向刚度贡献
- 右端 `-alpha*dPdz` 是孔压梯度耦合体力项

**稳定性关注：**
- `skelRelax` 太大会振荡，太小会慢；强耦合工况通常 0.5~0.8 较稳。

---

### 7.23 `etaBarFromP_local` / `rhoBarFromP_local` / `rebuildSeepageKappa_local`

- 黏度：指数压黏模型并做指数限幅（`etaExpMin/Max`）
- 密度：有理式压密模型
- `kappa = K/etaP` 后限幅

**稳定性关注：**
- 指数限幅是防止高压区 `exp()` 溢出的关键措施。

---

### 7.24 `harmonicMean_local`

给界面系数用谐均值：在高对比介质中，较算术均值更保守、更稳定。

---

### 7.25 `enforceYSymmetry2D_local` / `enforceYSymmetry3D_local`

`A = 0.5*(A + 镜像A)`：消除扫描顺序导致的微弱偏侧。

---

### 7.26 `computeEvFromUz_local`

中心差分/单边差分计算 `Ev≈dUz/dZ`。

---

### 7.27 `safeContour_local`

防止全常数场时 `contourf` 报错，自动退化为 `surf(...);view(2)`。

---

## 8. 你后续调参时可直接套用的“稳定性检查清单”

1. **先看 `errW` 是否单调下降**，若震荡，先减小 `secantDamping` 与 `maxH00Step`。
2. **看 `cavFracRe` 是否异常增大**，若大，检查 `H00` 是否过小导致大面积空化。
3. **看 `deltaEHL` 与 `deltaPorous` 的量级比**，避免 A1 抢占 A2 主导权。
4. **看 Y 对称误差**，若变大，优先检查边界更新与对称化调用顺序。
5. **看 `kappa` 是否频繁触及上下限**，若是，说明参数过激或模型需要再标度。


---

## 9. “真正逐行”版（按源码从头到尾，每 20~40 行一组）

> 说明：你在对话里给的是整段 MATLAB 源码文本而不是仓库内 `.m` 文件，因此下面采用“**连续逻辑行号**”方式逐组讲解（从函数第一行开始递增理解）。如果你把 `.m` 文件放到仓库里，我下一版可以给你“精确文件行号版（L1/L2/...）”。

### G01（约 1~30 行）：文件头与目标定义
- 这段是总览注释：列出 13 个版本特性（A1/A2、预设工况、Y 对称化、冻结闭合等）。
- 变量无量纲定义明确了 `X,Y,Z,P,P*,H,Uz` 的尺度，后续所有公式都基于这一套。
- **潜在 bug 点**：若后续新增模块未遵守该无量纲定义，会出现“看似收敛但量纲错位”的隐蔽错误。

### G02（约 31~70 行）：主函数开场与参数初始化
- `clc/clear/close all` 清理环境，防止残留变量污染。
- 初始化 `par` 后调用 `applyDefaultSkeletonParams_local` 做参数兜底。
- `Nx/Ny/Nz >= 3` 的检查保障二阶差分合法。
- **潜在 bug 点**：若网格过粗（例如 Nz=3 临界），边界二阶式对噪声更敏感。

### G03（约 71~110 行）：网格与影响核准备
- 建 2D/3D 网格，再预计算 A2 弹性影响核（FFT）。
- 初始化 `field2D/field3D`，用 Hertz 压力作 `P` 初值。
- **潜在 bug 点**：若初值全零，强耦合工况下可能进入“弱解分支”导致中心压力上不去。

### G04（约 111~150 行）：3D 初始场设定
- 顶面孔压 `Pp(:,:,1)=P`，其余深度层给底压。
- 位移场初值全 0。
- 首次调用 `updateFilmAndProps_local` 同步 `H/eta/rho/Xi/kappa/dPp_dZ`。
- **潜在 bug 点**：若此处不先更新，第一次 Reynolds 的 `Xi/S` 会不一致。

### G05（约 151~220 行）：外层耦合迭代框架
- 保存旧场用于收敛判据。
- Step A：渗流；Step B：骨架；Step C：属性更新；Step D：Reynolds；Step E：载荷闭合。
- 每轮 `H00` 更新后再解一次 Reynolds，防止载荷与膜厚不同步。
- **潜在 bug 点**：若删掉“二次 Reynolds 同步”，常见表现是 `errW` 卡住不降。

### G06（约 221~280 行）：收敛判据与迭代退出
- 同时监控 `errP/errH/errPp/errUz0/errH00/errDelta/errW`。
- 采用双阈值：场变量阈值 + 载荷阈值。
- **潜在 bug 点**：仅盯 `errP` 会误判；`P` 稳了但 `H00` 可能还在漂。

### G07（约 281~340 行）：冻结形态后的最终载荷闭合
- 冻结 `deltaEHL/deltaPorous/Pp/Uz/Ev/K`。
- 清空 Secant 历史，避免把外层历史曲线带入新阶段。
- 仅调 `H00` + 重解 Reynolds，压紧终态 `errW`。
- **潜在 bug 点**：若不清空 Secant 历史，可能出现“跨阶段割线误导”导致大跳步。

### G08（约 341~430 行）：结果输出与快速诊断
- 打印最终 `Load/H00/Pmax/Hmin/Uzmax/Uzmin`。
- 对 `Uz` 正负做自检，提示膜厚耦合符号是否一致。
- **潜在 bug 点**：`Uz` 符号与 `deltaPorous=-Uz` 若配置相反，会出现反物理“受压变厚/变薄”问题。

### G09（约 431~520 行）：A1+A2 形态诊断 + 对称诊断 + 绘图入口
- 看中心点 `deltaEHL` 与 `deltaPorous` 的主次关系。
- 诊断 `P(0,0)/Pmax` 是否抬升。
- 做 Y 对称误差统计并输出。
- **潜在 bug 点**：若 `symErr` 大，优先查边界更新顺序与对称化调用时机。

### G10（约 521~640 行）：`initParameters_local` 参数总表
- 定义网格、物理量、数值阈值、迭代控制、MG 参数、骨架参数等。
- 调 `applyCasePreset_local` 后再做基线覆盖逻辑。
- **潜在 bug 点**：预设修改有量纲参数后必须重新闭合无量纲群（你已做对）。

### G11（约 641~700 行）：`closeDimensionlessGroups_local`
- 做参数合法性检查（正值、泊松比奇异点）。
- 计算 `G_dim/lambda_dim/Lambda/Pi_p/Gamma/Lambda_s`。
- **潜在 bug 点**：`nu_dim` 接近 0.5 时 `lambda_dim` 巨大，骨架方程强刚性。

### G12（约 701~820 行）：`applyCasePreset_local`
- 三个 preset 分别控制“基线验证/参考工况/强多孔诊断”。
- 每个 preset 同时调整物理与数值参数（不是只改物理）。
- **潜在 bug 点**：强工况里若 `secantDamping` 太大，`H00` 闭合常振荡。

### G13（约 821~870 行）：`applyDefaultSkeletonParams_local`
- `isfield` 防御式赋默认值，保证参数缺省可运行。
- **潜在 bug 点**：多人协作改参数结构时，这层能避免字段缺失崩溃。

### G14（约 871~910 行）：`buildGrid2D_local` / `buildGrid3D_local`
- `Z` 方向从 `0` 到负值，所以 `dZs<0`。
- **潜在 bug 点**：`dPp_dZ` 的离散公式必须跟 `dZs` 符号匹配。

### G15（约 911~990 行）：`prepareElasticInfluence_local`
- 计算弹性系数 `C_e`，构造核 `K=const/sqrt(r^2+eps^2)`。
- 预存 FFT 与索引窗口供快速卷积裁剪。
- **潜在 bug 点**：`epsReg` 太小导致核中心过尖、噪声放大。

### G16（约 991~1040 行）：初始化 2D/3D 场 + Hertz 初值
- `field2D` 包含 `P/H/rho/eta/Xi/S/deltaEHL/deltaPorous/dPp_dZ`。
- `field3D` 包含 `Pp/Ux/Uy/Uz/Ev/K/etaP/kappa`。
- **潜在 bug 点**：字段若缺失，后续 `update` 会触发尺寸不匹配。

### G17（约 1041~1130 行）：`solveReynolds2D_local` 主体
- 每轮先 `build source` 再 `MG/GS`。
- 非法值、负压、边界反复强制，最后再对称化一次。
- **离散意义**：`A(P)=S` 的变系数椭圆方程，`Xi` 在面上取平均。
- **潜在 bug 点**：`tolRe` 太紧 + `omegaP` 偏大时会慢振荡。

### G18（约 1131~1175 行）：`buildReynoldsSource_local`
- `convTerm` 表示卷吸供给梯度；`seepTerm` 表示渗流交换量。
- `S = conv - seep` 并做 `sourceMax` 限幅。
- **潜在 bug 点**：`dPp_dZ` 符号错会让 `seepTerm` 方向反了。

### G19（约 1176~1255 行）：`reynoldsVcycle_local`
- 预平滑 → 残差限制 → 粗网格误差 → 延拓修正 → 后平滑。
- **潜在 bug 点**：粗网格 `Xi` 若不设下限，可能出现 `Ap≈0`。

### G20（约 1256~1370 行）：`reynoldsSmoothPressureGS_local` 与 `ErrorGS`
- `Ae/Aw/An/As` 是四向导通系数，`Ap` 为中心系数总和。
- SOR 更新用于加速收敛。
- **潜在 bug 点**：`omegaP>1` 且源项强时更易振荡。

### G21（约 1371~1460 行）：`reynoldsApplyA_local` + 限制/延拓
- `A(P)` 用于残差。
- `restrictFW` full-weighting；`prolongBilinear` 双线性插值。
- **潜在 bug 点**：延拓后边界若不清零，会破坏 Dirichlet 条件。

### G22（约 1461~1530 行）：`computeLoad2D_local` + `adjustH00_local`
- 载荷积分 `W=∬P dX dY`。
- H00 更新优先 Secant，退化时比例反馈。
- **潜在 bug 点**：`Fcur-Fprev` 太小会造成割线跳步（你已用 `secantTolDen` 防护）。

### G23（约 1531~1700 行）：`plotResults_local`
- 输出压力/膜厚/位移/孔压梯度等曲面与中心线。
- **潜在 bug 点**：仅看曲面不够，必须结合中心线与定量输出判断是否“形似而神不似”。

### G24（约 1701~1840 行）：`updateFilmAndProps_local`
- A2：`deltaEHL`（卷积 + 松弛 + 限幅）。
- A1：`deltaPorous = -Uz_surf * weight`。
- 更新 `H/eta/rho/etaP/K/kappa/Xi/dPp_dZ`。
- **潜在 bug 点**：任何一步未清 NaN/Inf 都会在 MG 或 SOR 中扩散成全场异常。

### G25（约 1841~1940 行）：`updateFilmAndPropsFrozenShape_local`
- 固定形态量，只让 `H00` 驱动膜厚与 Reynolds 系数变化。
- **潜在 bug 点**：若此阶段仍更新 `deltaEHL`，冻结闭合将失去意义。

### G26（约 1941~1985 行）：`solveElasticSurfaceDeflection_local`
- FFT 卷积 + 裁剪回原尺寸。
- **潜在 bug 点**：padding 尺寸不足会出现循环卷积污染（你已用较大 pad）。

### G27（约 1986~2145 行）：`solveSeepage3D_local`
- 内迭代中可周期性重算 `kappa(Pp)`，处理非线性。
- 侧面 Neumann 通过“复制邻点”实现。
- 6 个方向系数 `aE..aB` + `aP` 组装局部离散方程。
- **潜在 bug 点**：`omegaSeep` 太大 + 强非线性会在深层出现锯齿震荡。

### G28（约 2146~2295 行）：`solveSkeleton3D_local`
- 顶面牵引二阶单边离散，底面固支，内点 7 点更新。
- 右端耦合项 `-alpha*dPdz` 把孔压梯度传给骨架。
- **潜在 bug 点**：`dZs<0` 下 `dPdz` 正负与牵引正负必须整体一致。

### G29（约 2296~2360 行）：物性与渗流重建函数
- `etaBarFromP_local`：指数模型并限制指数范围。
- `rhoBarFromP_local`：压密有理式。
- `rebuildSeepageKappa_local`：`kappa=K/etaP` + 限幅。
- **潜在 bug 点**：不限制指数会导致 `exp` 溢出。

### G30（约 2361~结尾）：工具函数
- `harmonicMean_local`：高反差介质更稳。
- `enforceYSymmetry2D/3D`：抑制横向偏置。
- `computeEvFromUz_local`：差分计算 `Ev`。
- `safeContour_local`：常数场容错绘图。
- **潜在 bug 点**：`Ev` 差分若使用 `dZ` 而非 `dZs` 会改变量纲符号。

---

## 10. 给你一版“可直接对照调试”的检查顺序（从第 1 行跑到最后）

1. **先固定 `baseline_ehl_shape`**，只看 `P/H` 是否合理（不看多孔）。
2. **再开 `current_reference`**，检查 `deltaPorous` 相对 `deltaEHL` 的比例是否在可控范围。
3. **最后 `porous_diag_strong`**，逐步扫描 `k0_dim/E_dim/alpha/skelTopTractionScale/porousDispWeight`。
4. 每次改参数后都记录：`errW、cavFracRe、Pmax、Hmin、symErrP/H/Uz`。
5. 若收敛差：优先减小 `secantDamping`、`maxH00Step`，其次减 `omegaSeep/skelRelax`。

