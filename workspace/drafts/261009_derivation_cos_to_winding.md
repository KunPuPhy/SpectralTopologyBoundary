# 从强度的余弦形式到绕数：推导链（统一符号，2026-10-09）

> **目的**：一条线走完"两类探针的强度都是余弦 → 数据的镜像对称 → 绕数之模可读、符号不可读 → δ 弦给出 no-go 的精确边界"，每一步给出推导，全文符号一次约定。
> **数值核验**：每条恒等式由 `261009_verify_cos_to_winding.py`（同目录）验证到机器精度（24 项，§7）。
> **符号沿**：`spectral_function_vs_bragg.md`（球坐标规范）、`delta_string_2d_bragg.md`（δ 弦定理）；绕数统一记 $W$（不用 $\nu$）。

---

## 0. 约定与符号表（全文唯一来源）

| 符号 | 含义 |
|---|---|
| $k\in[-\pi,\pi]$ | BZ（晶格常数 $a=1$；恢复 $a$：$k\to ka$、$\delta\to\delta/a$）|
| $H(k)=d_0(k)\mathbb 1+\mathbf d(k)\cdot\boldsymbol\sigma$ | 二能级；$\mathbf d=(d_x,d_y,d_z)$ 实函数；能隙 $\|\mathbf d\|>0$ |
| $\hat{\mathbf d}(k)=\mathbf d/\|\mathbf d\|\in S^2$ | Bloch 矢量；轨迹 $\mathcal C=\{\hat{\mathbf d}(k)\}$ |
| $\hat z$ | 被分辨内态轴（子格/Nambu/自旋基）；$\lvert A\rangle=(1,0)^T$，$\lvert B\rangle=(0,1)^T$ |
| $\theta_k,\ \phi_k$ | $\hat{\mathbf d}(k)$ 绕 $\hat z$ 的球坐标：$\cos\theta_k=\hat d_z(k)$，$\phi_k=\arg(\hat d_x+i\hat d_y)$ |
| $\Delta\varphi(k,k_0)\equiv\phi_k-\phi_{k_0}$ | 方位差 |
| $\Theta(k_0,k)$ | $\hat{\mathbf d}(k_0)$ 与 $\hat{\mathbf d}(k)$ 的夹角：$\cos\Theta=\hat{\mathbf d}(k_0)\cdot\hat{\mathbf d}(k)$ |
| $W$ | 绕数（§2 定义）|
| $\Gamma=\hat{\mathbf n}_\Gamma\cdot\boldsymbol\sigma$ | 手性算符（若 $\{H,\Gamma\}=0$）；$W$ 的保护轴 $\hat n_W=\hat n_\Gamma$ |
| $\psi_k$ | 手性大圆上的带符号角（§4 起；构型 A 中 $\psi=\phi$）|
| $\beta=\angle(\hat z,\hat n_\Gamma)$，$\rho=\sin\beta$ | $\hat z$ 与手性圆平面的几何（§4–5）|
| $w^A_\pm(k)$；$w^S_\pm(k_0,k)$ | SPSF 权重；Bragg 权重（$+$ 带间、$-$ 带内）|
| $P_{\hat{\mathbf a}}=\tfrac12(\mathbb 1+\hat{\mathbf a}\cdot\boldsymbol\sigma)$ | 沿 $\hat{\mathbf a}$ 的自旋投影算符（$P_\pm(k)=P_{\pm\hat{\mathbf d}(k)}$）|
| $F(q)=\mathrm{diag}(1,e^{iq\delta})$ | δ 弦（§6）|
| $M_P$ | 过 $\hat z$ 平面 $P$ 的镜面反射；具体实现（$P=x$–$z$ 面）：$d_y\to-d_y$ |

---

## 1. 主公式：一切强度是两个自旋投影的重叠

**引理 P（带投影算符）。** $P_\pm(k)=\lvert u_\pm(k)\rangle\langle u_\pm(k)\rvert=\tfrac12(\mathbb 1\pm\hat{\mathbf d}(k)\cdot\boldsymbol\sigma)$。
*证明*：Pauli 代数 $(\mathbf d\cdot\boldsymbol\sigma)^2=\|\mathbf d\|^2\mathbb 1$，故 $(\hat{\mathbf d}\cdot\boldsymbol\sigma)^2=\mathbb 1$，谱为 $\{\pm1\}$，谱投影即上式。$\blacksquare$

**引理 O（重叠公式）。** 对任意单位矢量 $\hat{\mathbf a},\hat{\mathbf b}$：
$$\mathrm{Tr}\,P_{\hat{\mathbf a}}P_{\hat{\mathbf b}}=\tfrac12\big(1+\hat{\mathbf a}\cdot\hat{\mathbf b}\big).$$
*证明*：$P_{\hat a}P_{\hat b}=\tfrac14\big[\mathbb 1+(\hat a+\hat b)\cdot\boldsymbol\sigma+(\hat a\cdot\boldsymbol\sigma)(\hat b\cdot\boldsymbol\sigma)\big]$，用 $(\hat a\cdot\sigma)(\hat b\cdot\sigma)=(\hat a\cdot\hat b)\mathbb 1+i(\hat a\times\hat b)\cdot\sigma$ 与 $\mathrm{Tr}\,\sigma_i=0$。$\blacksquare$

两个探针都是"初态投影算符 × 末态投影算符"的迹：

- **SPSF**（子格 $A$、带 $\pm$）：
  $$w^A_\pm(k)=\big|\langle A|u_\pm(k)\rangle\big|^2=\mathrm{Tr}\,P_A\,P_\pm(k),\qquad P_A=P_{\hat z}
  \;\Longrightarrow\; w^A_\pm(k)=\tfrac12\big(1\pm\cos\theta_k\big).$$
- **Bragg**（下带 $k_0$ 初态 → 带 $\pm$ @ $k$）：
  $$w^S_\pm(k_0,k)=\big|\langle u_\pm(k)|u_-(k_0)\rangle\big|^2=\mathrm{Tr}\,P_{-\hat{\mathbf d}(k_0)}\,P_{\pm\hat{\mathbf d}(k)}
  \;\Longrightarrow\; w^S_\pm(k_0,k)=\tfrac12\big(1\mp\cos\Theta(k_0,k)\big).$$

统一表（"强度的余弦形式"的唯一来源）：

| 跃迁 | 初态自旋 $\hat a$ | 末态自旋 $\hat b$ | 权重 $\tfrac12(1+\hat a\cdot\hat b)$ |
|---|---|---|---|
| SPSF $A\to$ 带 $\pm$ | $\hat z$（**晶格钉死**）| $\pm\hat{\mathbf d}(k)$ | $\tfrac12(1\pm\cos\theta_k)$ |
| Bragg 带间 | $-\hat{\mathbf d}(k_0)$（**可动，轨迹上**）| $+\hat{\mathbf d}(k)$ | $\tfrac12(1-\cos\Theta)$ |
| Bragg 带内 | $-\hat{\mathbf d}(k_0)$ | $-\hat{\mathbf d}(k)$ | $\tfrac12(1+\cos\Theta)$ |

一行结论：**一切强度 = $\tfrac12(1\pm\cos$ 夹角$)$；两类探针的差别只是参考 $\hat a$ 被钉死（SPSF）还是可动（Bragg）。** 求和规则 $w^S_++w^S_-=1$ 顺带成立（$F$ 幺正）。

**球坐标展开（$\varphi$ 的唯一入口）。**
$$\cos\Theta(k_0,k)=\cos\theta_{k_0}\cos\theta_k+\sin\theta_{k_0}\sin\theta_k\cos\Delta\varphi(k,k_0).$$
方位差只以 $\cos\Delta\varphi$ 进入强度，通道增益 $\sin\theta_{k_0}\sin\theta_k$。SPSF 是增益为零的极端：参考在极轴上，点积只剩单项 $\hat d_z$。

---

## 2. 绕数的定义（与数据同一坐标系）

$$W(\hat n_W)\;\equiv\;\frac{1}{2\pi}\oint_{\mathrm{BZ}} dk\,\partial_k\varphi^{(\hat n_W)}(k),$$
$\varphi^{(\hat n_W)}$ = 绕轴 $\hat n_W$ 的方位角。**良定义条件**：轨迹不与 $\hat n_W$ 轴相交（方位角沿闭回路连续回绕）。等价表述：$\perp\hat n_W$ 的两个分量 $d^{(1)}+i\,d^{(2)}$ 绕复平面原点的圈数。

$W$ 的物理地位（两条性质，决定"符号值不值得读"）：

1. **一般（非手性）二能级**：对任选轴是良定义整数，但**不受保护**——形变可让轨迹穿过该轴（能隙不闭合），$W$ 跳 $\pm1$。
2. **手性系统**（$\{H,\Gamma\}=0$，$\Gamma=\hat n_\Gamma\cdot\sigma$）：$\hat{\mathbf d}(k)\cdot\hat n_\Gamma\equiv0$ ⟹ 轨迹钉死在平面大圆 $\mathcal C_\Gamma$ 上 ⟹ 穿 $\hat n_\Gamma$ 轴 = 穿极点 = $\|\mathbf d\|\to0$ 能隙闭合 ⟹ $W(\hat n_\Gamma)$ **拓扑保护**。此时取 $\hat n_W=\hat n_\Gamma$，且 $W=\frac{1}{2\pi}\oint d\psi_k$（$\psi$：沿圆带符号角，§4）。

两个标准构型（§5 展开）：SSH：$\hat n_\Gamma=\hat z$（圆 = 赤道，$\psi=\phi$）；Kitaev：$\hat n_\Gamma=\hat x$（圆 = $y$–$z$ 大圆，$d_z+i\,d_y=\|\mathbf d\|e^{i\psi}$）。

---

## 3. 镜像分离定理（核心推导）

**数据的对称群**：$O(2)_{\hat z}\triangleq\{R\in O(3):R\hat z=\hat z\}$ —— 绕 $\hat z$ 的旋转 $R_{\hat z}(\alpha)$ 与过 $\hat z$ 平面 $P$ 的反射 $M_P$。

**断言 1（数据不变）。** $O(2)_{\hat z}$ 的每个元都保持全部 $\delta=0$ 强度数据。
*证明*（一行一组）：旋转保 $z$ 分量、反射 $P\ni\hat z$ 的每个点固定故保 $z$ 分量 ⟹ $\hat d_z(k)$ 全同 ⟹ SPSF 全同；$R$ 正交 ⟹ $(R\hat{\mathbf d}(k))\cdot(R\hat{\mathbf d}(k'))=\hat{\mathbf d}(k)\cdot\hat{\mathbf d}(k')$ ⟹ Bragg 全同。$\blacksquare$

**断言 2（$W$ 的变换）。** $M_P\,R(\hat n,\alpha)\,M_P^{-1}=R(M_P\hat n,-\alpha)$（$\det M_P=-1$）⟹
- $\hat n_W\in P$：方位反演，$W\to-W$；
- $\hat n_W\perp P$：$W\to W$；
- 旋转 $R_{\hat z}(\alpha)$：$W$ 不变（缠绕对参考平移/旋转不变；系统幺正等价）。

手性系统取 $P=\mathrm{span}(\hat z,\hat n_\Gamma)$（$\hat n_\Gamma\parallel\hat z$ 时任取）⟹ $\hat n_W=\hat n_\Gamma\in P$ ⟹ $W\to-W$，且 $M_P\hat n_\Gamma=\hat n_\Gamma$ ⟹ 镜像系统仍以同一 $\Gamma$ 手性。

**分离定理。** $\delta=0$ 强度数据在 $O(2)_{\hat z}$ 下不变；$W$ 在其旋转分量下不变、在反射分量下变号。故：

- **(a) no-go**：数据决定的是 $O(2)_{\hat z}$ 轨道 $\{RS\}$，轨道上 $W$ 恰取 $\pm$ 两值 ⟹ 任何 $\delta=0$ 强度泛函（任意探针、参考、组合）只能给出**镜像偶**的 $W$ 泛函，**不能确定 $\mathrm{sign}\,W$**。这是存在性论证而非难度论证：镜像对 $(S,M_PS)$ 是两个真实物理系统——对实 $\mathbf d(k)$，$H'=(d_x,-d_y,d_z)\cdot\boldsymbol\sigma+d_0\mathbb 1$ 仍是合法哈密顿量（反幺正 $K$ 实现，能谱全同、带序不变），数据逐项全同而 $W$ 相反。
- **(b) $\lvert W\rvert$ 镜像偶**：原则上在数据可及范围内——但"原则上可定"≠"免费可读"，可操作判据需 §4 的条件。
- **(c) 适用范围**：对一切带隙二能级系统成立（含非手性——核验 c2 的 random 模型：非手性、$W_{\hat z}=1$、数据镜像全同、$W\to-1$）。分离是**运动学事实**；$W$ 值不值得读由 §2 的保护性决定。

---

## 4. $\lvert W\rvert$ 的操作判据（手性大圆上）

设置：轨迹在圆 $\mathcal C_\Gamma$ 上，$\psi_k$ = 沿圆带符号角，$W=\frac{1}{2\pi}\oint d\psi_k$。两个信号的圆上形式（全部从 §1 主公式代入）：

- **Bragg**：圆上两点的夹角余弦 $=\cos(\psi_k-\psi_{k_0})$ ⟹
  $w^S_\pm(k_0,k)=\tfrac12\big(1\mp\cos(\psi_k-\psi_{k_0})\big)$。
- **SPSF**：分解 $\hat z=\hat z_\parallel+\hat z_\perp$（$\perp$ 圆平面部分幅值 $\rho=\sin\beta$，$\beta=\angle(\hat z,\hat n_\Gamma)$；记其圆上位置 $\psi_z$）⟹ $\hat d_z(k)=\rho\cos(\psi_k-\psi_z)$ ⟹
  $w^A_\pm(k)=\tfrac12\big(1\pm\rho\cos(\psi_k-\psi_z)\big)$ —— **移位余弦，振幅 $\rho$**。

**判据一（满转移 ⟺ 对跖在像内）**：$w^S_+(k_0,k)=\tfrac12\big(1-\cos(\psi_k-\psi_{k_0})\big)$ 达 $1$ ⟺ $\cos(\psi_k-\psi_{k_0})=-1$ ⟺ $\psi_k=\psi_{k_0}+\pi$。一行，直接从余弦形式。记 $W_{\max}(k_0)\triangleq\max_k w^S_+(k_0,k)=1\iff-\hat{\mathbf d}(k_0)\in\mathrm{image}\,\hat{\mathbf d}$。

**判据二（手性满射定理）**：$W\neq0$ ⟹ $\psi:\mathrm{BZ}\to S^1$ 度数非零 ⟹ 满射（否则像含于真闭弧、可收缩，度数只能为 $0$）⟹ 像盖满整个圆 ⟹ 对**一切** $k_0$：对跖在像内 ⟹ $W_{\max}(k_0)\equiv1$。

**判据三（计数）**：$\psi(k)$ 单调（净绕 $\lvert W\rvert$ 圈）时，方程 $\psi_k=\psi_{k_0}+\pi\ (\mathrm{mod}\ 2\pi)$ 的解数 $=\lvert W\rvert$。操作形式：
$$R(k)\equiv\frac{w^S_+(k_0,k)}{w^S_-(k_0,k)}=\tan^2\frac{\psi_k-\psi_{k_0}}{2},$$
恰在对跖点发散 ⟹ **发散次数 $=\lvert W\rvert$**（无需绝对校准）。

**判据四（SPSF 侧，$\rho\neq0$ 时）**：$\hat d_z(k)=0\iff\psi_k=\psi_z\pm\tfrac\pi2$；单调绕 $\lvert W\rvert$ 圈每圈穿越 $2$ 次 ⟹ **$\pi$ 穿越数 $=2\lvert W\rvert$**（Estake–Roy 判据）。

**反向失效的三个洞（诚实清单，均有构造反例）**：
1. **满转移 $\not\Rightarrow W\neq0$**：DC 双程整圆（$\psi$ 一去一回：度数 $0$、像盖满整圆）⟹ $W_{\max}\equiv1$、每带 $\pi$ 穿越 $4$ 次（读作 $\lvert W\rvert=2$）——**像 $\neq$ 绕数**（核验 c6）。
2. **平庸弧可"擦过"对跖**（touch-and-go，测度为零的精细调节）⟹ 满转移 $\not\Rightarrow$ 非平庸（`generalization_to_general_systems.md` §5）。
3. **通道条件**：判据三、四要求 $\Delta\psi$ 通道开启——参考离开极点（§5 构型 B）、$\delta\neq a/2$（§6）；$\rho=0$ 时判据四整条失明（§5 构型 A）。

---

## 5. 两个构型（同一套公式，两个 $\beta$）

| | 构型 A（SSH 型）| 构型 B（Kitaev 型）|
|---|---|---|
| 手性轴 $\hat n_\Gamma$ | $\hat z$（$=$ 被分辨轴）| $\hat x$（$\perp$ 被分辨轴）|
| 手性圆 | 赤道 | $y$–$z$ 大圆（过两极）|
| $\psi$ | $=\phi_k$ | 圆内角：$d_z+i\,d_y=\|\mathbf d\|e^{i\psi}$ |
| $\beta=\angle(\hat z,\hat n_\Gamma)$ | $0$ | $\tfrac\pi2$ |
| $\rho=\sin\beta$ | $0$ | $1$ |
| SPSF | $w^A_\pm\equiv\tfrac12$ **严格平——全盲** | $\tfrac12(1\pm\cos\psi_k)$ 满对比 |
| Bragg | $\tfrac12(1\mp\cos\Delta\varphi)$ —— $\varphi$ 的唯一入口 | $\tfrac12(1\mp\cos\Delta\psi)$ |
| 退化点 | 无（极点不在赤道上）| $k_0$ 在极点（$\psi_{k_0}=0$ 或 $\pi$）|
| 天然 BEC 点 | —— | $k=0,\pi$ 都在极点上（$d_x\equiv0$、$d_y|_{0,\pi}=0$ ⟹ $\hat{\mathbf d}=\pm\hat z$）；能带底取其一 |

**极点退化（构型 B）**：参考在极点（$\psi_{k_0}=\pi$ ⟺ $\hat{\mathbf d}(k_0)=-\hat z$，此时 $\lvert u_-(k_0)\rangle=\lvert A\rangle$ 字面成立）：
$$w^S_\pm(k_0,k)=\tfrac12\big(1\mp\cos(\psi_k-\pi)\big)=\tfrac12\big(1\pm\cos\psi_k\big)=w^A_\pm(k)\quad\text{逐点}$$
—— **Bragg 退化为 SPSF：可动参考撞上被钉死的参考**（核验 c7：$5.6\times10^{-17}$）。后果：极点参考的 Bragg 继承 SPSF 的一切假阳性（$\Delta\psi$ 通道在参考端关闭，$\sin\theta_{k_0}=0$）。

**统一读法**：$\rho=\sin\beta$ 把两个构型缝成同一件事——SPSF 是"参考被钉死、振幅折减的 Bragg"。$\beta=0$：$\hat z\perp$ 圆平面 ⟹ 到圆上各点等距 ⟹ SPSF 平（判据四失明），Bragg 独担全部 $\Delta\psi$ 通道；$\beta=\pi/2$：$\hat z$ 在圆平面内 ⟹ 极点在圆上 ⟹ SPSF 满对比，且可动参考走到极点时两者字面相撞。判据一至四与 §3 分离定理在两个构型中逐字成立，无需新机制。

---

## 6. δ 弦：no-go 的精确边界

**来源**：子格 $B$ 相对 $A$ 错位 $\delta$ ⟹ 总密度带形状因子
$$\hat\rho_q=\sum_k\hat{\boldsymbol\psi}_{k+q}^\dagger\,F(q)\,\hat{\boldsymbol\psi}_k,\qquad
F(q)=\mathrm{diag}\big(1,e^{iq\delta}\big)=e^{i\beta/2}\,e^{-i(\beta/2)\sigma_z},\quad\beta=q\delta.$$
整体相位不进强度；剩下的是**绕 $\hat z$ 转 $\beta$ 的旋转算符**——一个外部手征对象。

**修正的主公式**（δ 弦定理，沿 `delta_string_2d_bragg.md` §2；核验 c1）：
$$w^S_+(k_0,k;\delta)=\tfrac12\Big[1-\cos\theta_{k_0}\cos\theta_k-\sin\theta_{k_0}\sin\theta_k\cos\big(\Delta\varphi(k,k_0)-q\delta\big)\Big],$$
等效几何：参考被弦移到 $R_{\hat z}(q\delta)\hat{\mathbf d}(k_0)$（方位移 $q\delta$、纬度不动）——弦走的就是 §1 的 $\Delta\varphi$ 通道，增益同 $\sin\theta_{k_0}\sin\theta_k$。手性圆上：$w^S_+=\tfrac12\big(1-\cos(\psi_k-\psi_{k_0}-q\delta)\big)$。

**精确对称性（镜像把 $\delta$ 映到 $-\delta$；核验 c4：$0.00$）**：
$$w^S_\pm(M_PS;\,\delta)=w^S_\pm(S;\,-\delta)\qquad\text{逐 }(k_0,k).$$
*一行推导*：镜像本征态 $=K\lvert u\rangle$（复共轭规范），$\langle Ku_+(k)\rvert F(\delta)\rvert Ku_-(k_0)\rangle=\langle u_+(k)\rvert F(\delta)^*\rvert u_-(k_0)\rangle^*=\langle u_+(k)\rvert F(-\delta)\rvert u_-(k_0)\rangle^*$，模方即 $w(S;-\delta)$。$\blacksquare$

**三个推论（no-go 的精确边界）**：

1. **$\delta=0$**：镜像是数据的精确对称 ⟹ §3 的 no-go 成立。**no-go 的精确范围 = 镜像对称探针类**（$F=\mathbb 1$ 纯密度，或 $\pm\delta$ 对称协议）。
2. **$\delta\neq0$**：弦携带外部手征（绕 $\hat z$ 的旋转有定向）⟹ 破坏镜像简并 ⟹ 原则上 $\mathrm{sign}\,W$ 可从强度读出。定量（单调缠绕）：$w\propto\cos\big((W-\delta)(k-k_0)\big)$——数据由 $\lvert W-\delta\rvert$ 决定，$\delta$ 扫描的 V 形顶点在 $\delta=W$（含符号）。**物理约束（诚实）**：$\delta$ 是样品几何参数（设计期固定）；自然对称晶格 $\delta=a/2$ 恰为擦除值——$\mathbb Z_2$ 对 $(W,W\oplus1)$：$\lvert W-\tfrac12\rvert=\lvert W\oplus1-\tfrac12\rvert$ ⟹ 数据全同（核验 c8：均匀对 $3.3\times10^{-16}$；真实 SSH 两相 $\max w$ 均 $0.500000$）；二维极点惯性（极点是 $\hat z$ 旋转不动点 ⟹ 弦在极点参考严格惰性——$\delta$ 旋钮与 SPSF 轴在同一退化点失明，`delta_string_2d_bragg.md` §3–4）。
3. **机制判别**（已验证）：$\delta$ 扫描区分缠绕驱动（真拓扑：$\delta\to a/2$ 被擦除，$1.0\to0.5$）与幅度驱动（touch-and-go：全 $\delta$ 存留 $1.0$）——方向与直觉相反。

---

## 7. 数值核验（`261009_verify_cos_to_winding.py`，24 项全过）

纯 numpy，本征分解无规范假设（模量规范不变）。模型库：random（非手性 $W_{\hat z}=1$）、A1/A2（赤道 $W=1/2$）、A0（常数 $W=0$）、B1$=$T（大圆 $W=1$）、R（回折 $W=0$）、DC（双程整圆 $W=0$）、fold（折叠 Kitaev $W=0$）、sshT/ssh0（真实 SSH 两相）。

| 检验 | 内容 | 结果 |
|---|---|---|
| c1 | 主公式 + δ 弦公式（6 模型含非手性）| $\le6.5\times10^{-15}$ |
| c2/c3 | 镜像 $d_y\to-d_y$：完整 $\delta=0$ 数据逐项相同；$W\to-W$ | $0.00$；$(1,-1),(2,-2),(1,-1)$ |
| c4 | $w(M_PS;\delta)=w(S;-\delta)$（3 模型 × $\delta\in\{0.3,0.7\}$）| $0.00$ |
| c5 | 手性 $W\neq0$：$W_{\max}(k_0)\equiv1\ \forall k_0$；对跖簇数 $=\lvert W\rvert$（A1:1、A2:2、B1:1）| $1.000000$ |
| c6 | DC（$W=0$）：$W_{\max}\equiv1$、簇数 $2$（假阳性，读作 $\lvert W\rvert=2$）| $1.000000$ |
| c7 | 极点参考 $k_0=\pi$：Bragg ≡ SPSF 逐点 | $5.6\times10^{-17}$ |
| c8 | $\delta=\tfrac12$ 擦除：均匀对 $(W{=}1,W{=}0)$ 全同；真实 SSH 两相 $\max w$ | $3.3\times10^{-16}$；$0.500000/0.500000$ |
| c9 | SPSF $\pi$ 穿越数：$W=1\to2$；R/fold（$W=0$）$\to2$（假阳性）| $2$ |
| c10 | fold（$W=0$，SPSF 读 $\nu=1$）：$W_{\max}(k_0)<1\ \forall k_0$ | $[0.517,0.934]$ |
| c11 | T vs R 最小对：SPSF 全同、极点 Bragg 全同、一般 $k_0$ 分离 | $0$；$1.1\times10^{-16}$；$1.000$ vs $0.500$ |
| c12 | 旋转 $R_{\hat z}(0.7)$：数据不变、$W$ 不变 | $\sim10^{-15}$（eigh 舍入；镜像为精确符号翻转故 $0.00$）|

**c11 端到端**：T（$\psi=k$，$W=1$）vs R（$k>0$ 半支镜像回折，$W=0$）——SPSF 逐点全同（$\cos$ 偶）、$\pi$ 穿越 $2$ 次同读 $\nu=1$（假阳性）、极点参考 Bragg ≡ SPSF 同骗；一般参考 $k_0=-\pi/2$：T 满转移 $1.000$、R 仅 $0.500$（对跖不在半圆弧像内）——被揭穿。整条推导链（余弦 → 镜像 → 判据 → 揭穿）在最小对上一次走完。

---

## 8. 对两句定稿的含义（推导的直接产出）

1. **研究问题句不动**——"确定到什么程度"正是 §3 分离定理 + §4 判据 + §6 边界的总纲。
2. **核心命题句的修复从句需两处精化**：
   (a) 缺失比特的相位载体要分翼：回路相位（Zak）$=e^{i\pi W}$ 给的是 **$W\bmod2$**——它自身是镜像偶的（$e^{i\pi(\pm1)}=-1$，分不开镜像对），其独特价值在**与 $k_0$ 无关、极点照常工作**；**符号那一比特**住在奇宇称 sector 的相位里（1D：非对角谱函数的幅角 $\propto\mathrm{sign}(\sin\psi)$；2D：Berry 曲率/Chern）。
   (b) "只能由相位补上"的**精确范围 = 镜像对称探针类**（$\delta=0$）；δ 弦是强度层内唯一重新注入取向的通道——代价是探针自带外部手征参考，且受几何约束（$\delta=a/2$ 擦除、极点惯性）。
3. 两处定稿已并入命题句（2026-10-15 三层版：共性 → 特性 → 限定，见 `261009_research_question_and_thesis.md`）。
