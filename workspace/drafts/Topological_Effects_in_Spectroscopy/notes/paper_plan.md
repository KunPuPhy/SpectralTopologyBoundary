# PRB 论文计划：最小机制版（v2，2026-09-11）

- **v2 替代 v1**。v1 的"完备刻画（complete characterization）"定位是 clutter 的来源
  （完备 → 分类学 → 等级表 → 几何/信息论词汇）；v2 改为**机制文**定位，装置自动消失。
- **定位声明**：机制文，**不是协议文**。
- **目的一句话**（自我约束，写在稿纸开头）：
  *本文回答强度谱探测能带拓扑为什么会错、何时错、错在哪一比特，
  并说明为什么修复必然来自相位测量。*

---

## 标题与摘要（初版）

### 标题

**主选**：
> **Why spectroscopic intensities mistake trivial bands for topological ones: an even-angle mechanism with phase repair**

**备选**：
1. *Even functions of an angle: structural false positives in spectral-function and Bragg detection of band topology*
2. *Spectroscopic false positives in band topology: blind condensation points, an off-axis bound, and the one-bit phase repair*

对照 PRB 要求：sentence case；主选 16 词（≤20）；关键词齐
（spectroscopic intensities / topological bands / false positives / phase repair）；
问题式主句给记忆点，副题给机制与检索。

### 摘要（完整版，~210 词）

> The single-particle spectral function (measured in photoemission and tunneling)
> and the dynamical structure factor (measured in Bragg spectroscopy) are the two
> standard intensity probes of topological bands. We show that their false
> positives in detecting band topology stem from one mechanism. In any two-band
> system, both signals depend on the Bloch trajectory only through the cosine of
> an angle---the angle to a fixed internal axis for the spectral function, and
> the angle to the ground-state Bloch vector for Bragg. Because the cosine is
> even, a mirror deformation of the trajectory about the resolved axis leaves
> every intensity unchanged while erasing or flipping its winding number: the
> false positives are structural, not accidental. Three consequences follow.
> (i) When the Bragg reference lies on the resolved axis---as at the natural
> condensation points of the Kitaev chain and the Qi-Wu-Zhang model---the Bragg
> weights reduce exactly to spectral-function weights. (ii) An off-axis
> reference bounds the maximum inter-band Bragg intensity of any folded
> trajectory by $[1+\sqrt{1-d_\perp^2}]/2<1$, where $d_\perp$ is the reference's
> transverse component; this unmasks the known false-positive models and
> quantifies how far the condensate must be moved. (iii) No combination of
> intensity data can determine the sign of the winding number: the missing
> information is a single phase bit, carried by off-diagonal coherences or by
> the loop phase of Bragg matrix elements, which equals the Zak phase. The
> results are derived for noninteracting two-band systems at zero temperature
> and verified numerically for spinless-fermion chains and a two-dimensional
> Chern insulator.

**压缩版（~700 字符，备 APS 短摘要要求）**：

> Single-particle spectroscopy and Bragg spectroscopy are the standard intensity
> probes of topological bands. In any two-band system both signals depend on the
> Bloch trajectory only through the cosine of an angle: to a fixed internal axis
> for the spectral function, to the ground-state Bloch vector for Bragg. This
> even dependence makes false positives structural---mirror deformations of the
> trajectory leave all intensities unchanged while erasing its winding. A Bragg
> reference on the resolved axis (the natural condensation point of the Kitaev
> and Qi-Wu-Zhang models) restores the full blindness; an off-axis reference
> bounds folded trajectories' inter-band intensity below unity. No intensity
> combination reads the winding's sign: the missing bit is a phase, carried by
> coherences or the loop phase of Bragg matrix elements (the Zak phase).

**中文对照（校对故事线用，不进稿）**：
单粒子谱函数（光电子谱/隧道谱）与动力学结构因子（Bragg 谱）是探测拓扑能带的两种
标准强度观测量。我们指出二者探测能带拓扑时的假阳性源于同一机制：任何二能带系统中，
两种信号对 Bloch 轨迹的依赖都只通过某个夹角的余弦——谱函数取轨迹对固定内态轴的夹角，
Bragg 取对基态 Bloch 矢量的夹角。由于余弦是偶函数，把轨迹关于被分辨轴作镜面变形，
一切强度不变而绕数被抹去或翻转：假阳性是结构性的、非偶然的。三条结论：
(i) Bragg 参考点位于被分辨轴上（Kitaev 链与 QWZ 模型的天然凝聚点即是）时，
Bragg 权重严格退化为谱函数权重；(ii) 参考点离轴后，折叠型轨迹的带间 Bragg 强度
被 $[1+\sqrt{1-d_\perp^2}]/2<1$ 压制——既揭穿已知假阳性模型，也量化"凝聚点要挪多远"；
(iii) 任何强度组合都读不出绕数符号：缺失的信息是一个相位比特，由非对角相干或
Bragg 矩阵元的回路相位（等于 Zak 相位）携带。结果在零温、无相互作用二能带框架下导出，
并对无自旋费米子链与二维 Chern 绝缘体数值验证。

### 摘要自检（对照 PRB 要求）

- [x] 自含（定义两种探针，无需正文即可读懂）
- [x] 无引用（PRB 摘要不引文献；Estake–Roy 只在 §I 正文引）
- [x] 含定量结果（严格退化 + 显式界公式 + "below unity"）
- [x] 三条贡献编号呈现，与 §I 贡献四条一一对应（机制句 + i/ii/iii）
- [x] 结尾诚实 scope（两带、无相互作用、零温）+ 数值验证对象
- [x] 关键词覆盖：spectral function / dynamical structure factor / Bragg /
      band topology / winding number / Zak phase / false positive

### 写作注记（摘要→正文的一致性）

- 界公式前提是**折叠型（半空间 confined）轨迹**——摘要用 "folded trajectories" 限定；
  §III/§IV 行文保持一致。cos-k 弧型走另一机制（弧跨度 $<\pi$），**不套此界**；
- "natural condensation points" 承担极点预言的检索功能，正文 §III 推论 (b) 与之对应；
- 摘要中机制句"erasing **or flipping**"对应折叠（抹除）与整轨镜像（翻转）两种变形，
  §III 机制段要写全这两种；
- "off-diagonal coherences" 泛指一阶相干/非对角谱函数，避免无引用地具体引 E-R Eq. (25)。

---

## 0. 为什么这个定位够发 PRB

1. **锚点**：E-R 2025 (PRB 111, 085406) 刚制造了问题（逐模型结果 + 二维假阳性 + 相干补救）。
   一篇"解释其全部结果为特例、证明假阳性是结构性的、证明其补救是必然"的机制文
   是天然可发表续作——机制/解释型论文在 PRB 的生态位条件是
   "改变人们对一个现成测量方案的理解"，此处满足。
2. **闭环**：偶函数机制 → 解释 Bragg 的真好处 → 解释假阳性普适性 →
   **符号推论（3 行）把门关死** → "点到为止"合法。
   没有 no-go，停笔会被追打；有了它，停笔是逻辑完备的。
3. **两个定量物理句**防"纯负面"质疑（均现成）：
   - 极点退化 → "天然 BEC 点恰是最坏测点"（可检验的实验指导）；
   - 偏轴界 $w_{\max}=\frac{1+\sqrt{1-\hat d_{0\perp}^2}}{2}<1$（定量、可证伪）。

---

## 1. 故事闭环（全文四步）

1. 两个信号都是**夹角的余弦**：SPSF = 轨迹点对**固定 $z$ 轴**的夹角；
   Bragg = 轨迹点对**基态点**的夹角。
2. 这解释 Bragg 的真好处：**盲区随参考轴移动**——SPSF 的盲区
   （绕固定轴的镜像折叠）钉死不动；Bragg 的盲区跟着基态走，一般躲得开，
   **除非基态恰好坐在固定轴上（极点）**。
3. 这解释假阳性的普适性：把轨迹**关于被分辨轴作镜面折叠**，
   一切强度不变、绕数归零——假阳性不是模型意外，是机制必然。
4. 符号推论：没有任何强度组合能读出绕数符号 → 修复必然是相位测量 → 停笔。

---

## 2. 保留的定量要素（只有四件，其余删/进抽屉）

| # | 内容 | 为什么留 |
|---|---|---|
| 1 | 两个公式 $\pi(1\pm\hat d_z)$、$\frac12(1\mp\hat d\cdot\hat d_0)$ | 全文语言，纯点积 |
| 2 | 镜像机制 + 符号推论 | 机制本体 + 合法停笔 |
| 3 | 极点退化：天然 BEC 点（Kitaev $k_0=0$、QWZ $(\pi,\pi)$）是最坏测点 | 可检验预言 |
| 4 | 偏轴界 $\frac{1+\sqrt{1-\hat d_{0\perp}^2}}{2}$ | 定量、可证伪、物理 |

---

## 3. 语言处方（执行到成稿）

| 删/改 | 改成 |
|---|---|
| 量子距离 / Fubini–Study / 测地角距 | "Bloch 球上两个矢量的夹角"（公式本来就是点积） |
| ℤ₂ 盲区 / 模糊等级 / 等价类 | 一句话："分不清 $\Phi$ 与 $-\Phi$：两条镜像轨迹的一切强度相同" |
| Wilson 圈 / Berry 联络 | "把一圈带内 Bragg 矩阵元的相位乘起来，等于 Zak 相位"（引用即可） |
| 像 / 满射 / degree | "轨迹盖没盖满球、过没过对跖点" |
| **保留** | **偶函数依赖**（全文唯一机制词）、绕数、泡利矩阵、点积、$\delta$ 线强度 |

---

## 4. 论文结构（~6 页双栏 + SM，4–5 图）

### §I 引言（~1 页）
- 问题来自 E-R 2025；贡献四条：一个机制 / 极点结果 / 定量界 / 一比特修复；
- 一句反差埋线：一维绕数符号是良性的（Kitaev 链物理不变量是 $\nu\bmod 2$），
  二维 Chern 符号 = Hall 符号、不可回避；
- 诚实 scope 句：两带、零温、无相互作用。

### §II 两个公式，一句话机制（~1 页）
- 定义两种探针，写出两个公式（点积语言）；
- Bloch 球图像各一句：SPSF 量"到固定轴的夹角"，Bragg 量"到基态点的夹角"；
- **结尾钩子一句，不展开**："两者都是夹角的余弦——偶函数"。

### §III 机制与三条推论（~1.5 页）
- **机制段**：把轨迹关于被分辨轴作镜面折叠（横向分量 $d_\perp\to-d_\perp$），
  一切强度不变、绕数翻转或归零；
- **推论 (a) 符号不可见**（3 行）：$d_\perp\to-d_\perp$ 保持所有点积，$\nu\to-\nu$；
- **推论 (b) 极点退化**：基态在极点（$d_\perp(k_0)=0$）时
  $\hat d\cdot\hat d_0=\pm\hat d_z$，Bragg 权重逐点化为 SPSF；
  天然 BEC 点恰好是极点（Kitaev $\mathbf d(0)=(0,0,-\mu/2-J)$；QWZ $(\pi,\pi)$）；
- **推论 (c) 偏轴界**：基态横向分量非零时，折叠轨迹的带间线最大强度
  被严格压在 1 之下——$k_0$ 挪出固定轴多远、安全多少，读数即得；
- 诚实一句：双覆盖型轨迹可绕过此界，由推论 (a) 统摄（不做分类学）；
- 可选 remark 一句：形状因子（$\delta\neq0$）旋钮也打不破极点退化。

### §IV 假阳性是一个家族（~1.5 页）
同一机制读四个已知案例，每家一段 + 一张总结小表：

| 案例 | 折叠方式 | $k_0$ | 谁被骗 | 机制角色 |
|---|---|---|---|---|
| Kitaev $\mu<2J$ | —— | —— | 无人（正面参照） | 判据及其单向性 |
| $\sin^2k$ 变体 | $d_y\geq0$ 半平面折叠 | $0$（极点） | SPSF + Bragg 同骗 | 机制最小实例 + 极点 |
| QWZ 假阳性（E-R 自己的） | $d_y=\sin^2k_y$ | $(\pi,\pi)$（极点）→一般点 | 极点骗 / 挪开不骗 | E-R 衔接 + 界起作用 |
| $\cos k$ 变体 | 线段退化 | $0$（**非**极点） | 只有 SPSF | 偏轴修复起效（$0.93<1$） |

### §V 修复是一比特相位（~0.75 页）
- 一段：缺失的信息 = 轨迹侧别/取向 = 一个相位比特；
  非对角谱函数 / 一阶相干测量（引 E-R Eq. 25）；
- 一段：带内 Bragg 矩阵元的回路相位乘积 $=e^{i\pi\nu}$（Zak 相位），
  任意 $k_0$、非手性也给 Zak——原理上可测的陈述 + 现有数值一张图，**不展开协议**；
- 结尾一句：no-go 与修复在此闭合。

### §VI 展望（~0.5 页）
- 双覆盖残余（一句，指向 SM/后续）；
- 二维：符号 = Hall 符号，值得单独处理；
- 平台一句级（光晶格 SSH / 抖动格点 / 合成维度）；相互作用、多带一句级。

---

## 5. 图件（4–5 张）

| 图 | 节 | 内容 | 状态 |
|---|---|---|---|
| F1 | §II | 示意图：Bloch 球、两种夹角（固定轴 vs 基态轴） | 新画（简单） |
| F2 | §III | 极点退化最小演示（可用 even_blindness 的 (c)(d) 重组） | 现成素材 |
| F3 | §IV | 六面板 kitaev 或 falsepositive | 现成 |
| F4 | §V | 回路相位 = Zak | 现成（phase_channel） |
| F5（可选） | §III | 偏轴界曲线 $w_{\max}(\hat d_{0\perp})$（很物理的小图） | 现成数据 |

页数预算：I 1 + II 1 + III 1.5 + IV 1.5 + V 0.75 + VI 0.5 ≈ **6 页** + SM。
格式：`revtex4-2`，`prb,twocolumn`。

---

## 6. 抽屉清单（不进正文，别扔）

| 内容 | 去向 |
|---|---|
| 模糊等级表的形式化（cos_even_blindness §5） | SM 或 reply letter |
| DC 双覆盖分类学、$\Phi=2\pi\sin k$ 光滑实例 | SM / 第二篇 |
| 2D 定理化（$C\neq0\Rightarrow$ 满射、二维 DC、镜像 no-go） | 第二篇或 reply |
| 双参考 $\lvert\nu\rvert$ 重构算法 + 含噪 demo | SM / 第二篇 |
| 干涉 Bragg 协议设计 | 后续协议文 |
| 相互作用 / 多带 / 平台参数估计 | 后续 |

抽屉三个用途：审稿索要时的弹药；第二篇坯子；reply letter。

---

## 7. 任务清单（最小版）

- **必做（唯一）**：P2 —— 把 §III 或 §IV 一张关键图重画为
  Lorentzian 展宽 $\Gamma$ + 有限温细致平衡版本（约 1 天）；
- 顺手（可选，各半天）：F1 示意图、F5 界曲线；
- 其余全部降级进抽屉。

---

## 8. 审稿风险（精简版）

| 风险 | 预防 |
|---|---|
| "纯负面 / 无提议" | 两个定量句（极点 + 界）+ §V 一比特修复；no-go 使停笔合法 |
| "全篇 $T=0$、$\delta$ 线" | P2 图 |
| "相对 E-R 新意" | 四条差异化：解释其结果为特例 / 证明其补救必然 / 极点预言 / 定量界 |
| "数学味重" | §3 词表执行；全文只用点积与偶函数 |
| "只做两带" | scope 句诚实划界 |

---

## 9. 素材 → 章节索引（文件级）

| 章节 | 来源 |
|---|---|
| §II | `notes/main.tex` §1–2（压缩）；`spsf_vs_bragg_universal.py` |
| §III | `notes/cos_even_blindness.tex` T1–T3（**翻译成物理语言**）；`phase_channel_verification.py`（横向界）；`delta_string_verification.py`（极点惰性，可选 remark） |
| §IV | `codes/plot_1d_bloch.py`（六面板）；`kitaev_dy2cosk.py`；`qwz_false_positive_bragg.py`；`verify_even_blindness.py`（选面板） |
| §V | `phase_channel_verification.py` + `notes/phase_information_paradox.md` |
| §VI | `kitaev_dy2cosk.py`（开放链零模，一句）；各展望 |
| SM | `main.tex` 附录 A–C；其余 md 笔记按需 |

---

## 10. 写作顺序

1. **§III 先写**（从 cos_even_blindness.tex 把 T1–T3 翻译成物理语言，删全部等级表词汇）；
2. §II、§IV（公式 + 案例四段 + 小表）；
3. P2 图插入；
4. §V、§VI；
5. §I 引言最后写（贡献四条此时定型）；
6. 全文按 §3 词表过一遍。
