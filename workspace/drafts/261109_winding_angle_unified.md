# Winding number 的参考无关表述：从二能级定义到“任意面内参考走满 0 到 2π”

> **日期**：2026-11-09（重写版 2）  
> **任务**：从 `notes/Note.tex` 中 `Winding number` 一节的二能级定义出发，逐步推导“以 winding 平面内**任意**固定矢量为参考，相对角都取遍 $0$ 到 $2\pi$”，并说明它对 ARPES / Bragg 判据的含义。  
> **符号**：与 `Note.tex` 完全一致——$H$，$\Psi_k$，$h(k)$，$d_0(k)$，$\vec{d}(k)$，$\vec{\sigma}$，$\hat n_\Gamma$，$\Gamma$，$\hat{d}(k)=\vec{d}(k)/|\vec{d}(k)|$，$W$。全文用 $W$（不用 $\nu$）记 winding number。  
> **核验**：同目录 `261109_verify_reference_free_winding.py`，全部断言数值验证通过（结果见 §4）。

---

## 0. 主定理（本文要证明的结论）

设 $\vec{d}(k)$ 满足手性约束 $\vec{d}(k)\cdot\hat n_\Gamma=0$（轨迹在垂直手性轴的平面内），$\phi_k$ 为平面内极角，$W$ 为绕数。取 winding 平面内**任意固定的非零矢量** $\vec{a}$（归一化为 $\hat a$；轨迹点 $\hat a=\hat{d}(k_0)$ 只是其特例），定义相对角
$$
\Delta\phi_{k,a}\equiv\phi_k-\phi_a\pmod{2\pi},
$$
其中 $\phi_a$ 为 $\hat a$ 在平面内的固定方位角。则：

1. **参考无关性**：
   $$
   \deg\bigl(k\mapsto\Delta\phi_{k,a}\bigr)=W
   \qquad\text{对平面内一切固定 }\hat a.
   $$
2. **全覆盖定理**：若 $W\neq0$，则对每一个 $\hat a$，
   $$
   \{\Delta\phi_{k,a}:k\in{\rm BZ}\}=[0,2\pi).
   $$
   若角度提升单调，每个值恰被取到 $|W|$ 次。
3. **$\pi$ 点**：$[0,2\pi)$ 包含 $\pi$，故对每个 $\hat a$ 存在 $k$ 使
   $$
   \Delta\phi_{k,a}=\pi
   \quad\Longleftrightarrow\quad
   \hat{d}(k)=-\hat a
   $$
   （$\hat a$ 的对跖点在像内）。这是拓扑信号能发生的点。
4. **单向性**：原路返回的整圆轨迹对每个 $\hat a$ 也覆盖整个 $[0,2\pi)$，但 $W=0$。故全覆盖是 $W\neq0$ 的**必要非充分**条件；确定 $W$ 必须数带符号的净圈数。

**关键点**：推导中用到的唯一性质是 $\phi_a$ 为**常数**——参考矢量是否落在轨迹上完全无关。轨迹参考与面内任意固定参考给出同一个度数、同一份覆盖、同一个 $\pi$ 点（核验 c6）。

---

## 1. 出发点：`Note.tex` 的二能级定义

### 1.1 二能级模型

每个原胞有两个内态 $\alpha,\beta$，Bloch 哈密顿量为
$$
H=\sum_k \Psi_k^\dagger\,h(k)\,\Psi_k,
\qquad
\Psi_k=\begin{pmatrix}c_{k,\alpha}\\ c_{k,\beta}\end{pmatrix},
$$
$h(k)$ 是 $2\times2$ 厄米矩阵。由 Pauli 基展开，
$$
h(k)=d_0(k)\,\sigma_0+\vec{d}(k)\cdot\vec{\sigma},
$$
$d_0$ 只平移能谱、不动本征态，吸收后
$$
h(k)=\vec{d}(k)\cdot\vec{\sigma},
\qquad
E_\pm(k)=\pm|\vec{d}(k)|,
$$
能隙打开 $\iff$ $|\vec{d}(k)|>0$ 对一切 $k$。归一化 Bloch 矢量
$$
\hat{d}(k)=\frac{\vec{d}(k)}{|\vec{d}(k)|}\in S^2 .
$$

### 1.2 手性对称 $\Rightarrow$ 轨迹在平面内

手性算符 $\Gamma=\hat n_\Gamma\cdot\vec{\sigma}$（$\hat n_\Gamma$ 为 Bloch 空间任意单位矢量），手性对称为 $\{h(k),\Gamma\}=0$。用 Pauli 代数（$(\vec{a}\cdot\vec{\sigma})(\vec{b}\cdot\vec{\sigma})=(\vec{a}\cdot\vec{b})\,\sigma_0+i(\vec{a}\times\vec{b})\cdot\vec{\sigma}$）：
$$
\{h(k),\Gamma\}
=\{\vec{d}(k)\cdot\vec{\sigma},\,\hat n_\Gamma\cdot\vec{\sigma}\}
=2\bigl(\vec{d}(k)\cdot\hat n_\Gamma\bigr)\,\sigma_0
=0
\quad\Longrightarrow\quad
\boxed{\ \vec{d}(k)\cdot\hat n_\Gamma=0\ \ \forall k\ }
$$
即 $\vec{d}(k)$ 整条轨迹被限制在**垂直于手性轴的二维平面**内；能隙打开保证闭轨迹不过原点。这就是“存在手性轴，使矢量只在二维平面变化”的定义级来源。

### 1.3 Winding number：手性轴分量定义（`Note.tex` 式 eq:estake_winding）

在平面内取右手正交基 $\hat e_1,\hat e_2$（$\hat e_1\times\hat e_2=\hat n_\Gamma$）。定义
$$
\boxed{\ 
W=\frac{1}{2\pi}\int_{-\pi}^{\pi}dk\,
\hat n_\Gamma\cdot
\left(\hat{d}(k)\times\frac{d}{dk}\hat{d}(k)\right).
\ }
$$
$W=0$：轨迹不绕原点；$W=\pm1$：绕一圈，符号由相对 $\hat n_\Gamma$ 的取向决定。

---

## 2. 推导链（五步）

### 2.1 第一步：平面内极角参数化

因 $\hat{d}(k)$ 在平面内，可写
$$
\hat{d}(k)=\cos\phi_k\,\hat e_1+\sin\phi_k\,\hat e_2,
\qquad
\phi_k\ \text{定义 mod }2\pi .
$$
沿 BZ 取连续提升 $\tilde\phi_k$。注意 $k=-\pi$ 与 $k=\pi$ 是同一物理点：主值角 $\phi_\pi=\phi_{-\pi}\ (\mathrm{mod}\ 2\pi)$，提升两端只能差 $2\pi$ 的整数倍。

**角度的参考任意性**：换一个平面内参考方向 $\hat a$（方位角 $\phi_a$，**常数**），只是把整个角度表平移一个常数：
$$
\phi^{(a)}_k=\phi_k-\phi_a .
$$
这正是一开始选 $\hat e_1$ 为零方位角的做法——$\hat e_1$ 本身就是“任意平面内参考”的一个选取。

### 2.2 第二步：手性轴分量 $=$ 角度提升速率

直接计算：
$$
\frac{d}{dk}\hat{d}(k)
=\partial_k\tilde\phi_k\,
\bigl(-\sin\phi_k\,\hat e_1+\cos\phi_k\,\hat e_2\bigr),
$$
$$
\hat{d}\times\frac{d}{dk}\hat{d}
=\partial_k\tilde\phi_k\,
(\cos^2\phi_k+\sin^2\phi_k)\,
(\hat e_1\times\hat e_2)
=\partial_k\tilde\phi_k\,\hat n_\Gamma .
$$
因此
$$
\hat n_\Gamma\cdot\Bigl(\hat{d}\times\frac{d}{dk}\hat{d}\Bigr)=\partial_k\tilde\phi_k,
$$
代回定义（数值核验 c1：随机手性轴、随机半径 $r(k)$ 下恒等，说明该式与 $|\vec d|$ 无关）：
$$
\boxed{\ 
W=\frac{1}{2\pi}\int_{-\pi}^{\pi}dk\,\partial_k\tilde\phi_k
=\frac{\tilde\phi_\pi-\tilde\phi_{-\pi}}{2\pi}.
\ }
$$
**这一步是整个推导的枢纽：`Note.tex` 的手性轴分量定义自动等于角度提升的度数。** 两端主值相同，绕数完全来自提升多出的 $2\pi W$。

### 2.3 第三步：任意面内固定参考，度数不变

取平面内任意固定非零矢量 $\hat a$（方位角 $\phi_a$ 为常数；轨迹点 $\hat a=\hat{d}(k_0)$ 是特例），定义相对角及其提升
$$
\Delta\phi_{k,a}=\phi_k-\phi_a\ (\mathrm{mod}\ 2\pi),
\qquad
\Delta\tilde\phi_{k,a}=\tilde\phi_k-\phi_a .
$$
两端提升差：
$$
\Delta\tilde\phi_{\pi,a}-\Delta\tilde\phi_{-\pi,a}
=\bigl(\tilde\phi_\pi-\phi_a\bigr)-\bigl(\tilde\phi_{-\pi}-\phi_a\bigr)
=\tilde\phi_\pi-\tilde\phi_{-\pi}
=2\pi W.
$$
**参考角 $\phi_a$ 完全抵消**——推导用到的唯一性质就是它是常数。故
$$
\boxed{\ 
\deg\bigl(k\mapsto\Delta\phi_{k,a}\bigr)=W
\qquad\text{对平面内一切固定 }\hat a.
\ }
$$
（核验 c2：轨迹参考；核验 c6：随机面内非轨迹参考——两者给出同一度数到机器精度。）这就是“换平面内任何一个矢量作参考，绕数不变”的严格含义。

### 2.4 第四步：全覆盖定理与 $\pi$ 点

**定理**：$W\neq0$ $\Longrightarrow$ 对每个面内固定 $\hat a$，$\{\Delta\phi_{k,a}\}=S^1$。

**证明**：$\{\Delta\phi_{k,a}\}=\{\phi_k\}-\phi_a$，即像只是 $\{\phi_k\}$ 整体平移一个常数，平移不改变是否满射。而 $W\neq0$ 时映射 $k\mapsto\phi_k$ 的度数非零：若像不满，则存在 $p\in S^1$ 不在像中；像为连通紧集，含于 $S^1\setminus\{p\}$，后者同胚于开区间、可缩，于是映射同伦于常值映射，度数为 $0$——矛盾。故满射。$\blacksquare$

**推论**（核验 c3、c6）：

- 对每个 $\hat a$，相对角取遍 $[0,2\pi)$ 每个值（网格分辨率内最大圆间隙 $\sim$ 采样间距）；
- 对每个 $\hat a$ 存在 $k$ 使 $\Delta\phi_{k,a}=\pi$，即
  $$
  \hat{d}(k)=-\hat a
  $$
  ——**$\pi$ 点 $=$ 参考方向的对跖点在像内**，这就是“$0$ 到 $2\pi$ 包含 $\pi$”的数学内容；
- 若 $\tilde\phi_k$ 单调，每个值恰被取到 $|W|$ 次。

### 2.5 第五步：单向性——原路返回反例与单调等价

全覆盖**不**蕴涵 $W\neq0$。反例（核验 c4）：
$$
\tilde\phi_k=
\begin{cases}
2k+\pi, & k\in[-\pi,0],\\[4pt]
\pi-2k, & k\in[0,\pi],
\end{cases}
$$
$\tilde\phi:-\pi\to+\pi\to-\pi$，像为整圆、对每个 $\hat a$ 都含 $\pi$ 点，但
$$
\tilde\phi_\pi-\tilde\phi_{-\pi}=0
\quad\Longrightarrow\quad
W=0.
$$
这就是“轨迹到了 $\pi$ 也可能折回”：度数是带符号积分
$$
W=\frac{1}{2\pi}\int_{-\pi}^{\pi}dk\,\partial_k\tilde\phi_k,
$$
去程与回程贡献相消。

**单调情形**：若 $\tilde\phi_k$ 单调，则 $W=0$ 迫使 $\tilde\phi$ 常值（像为一点），故
$$
\text{单调时}\quad
\text{全覆盖}\iff W\neq0;
\qquad
\text{一般情形}\quad
W\neq0\Rightarrow\text{全覆盖，反向不成立}.
$$

---

## 3. 与 ARPES / Bragg 判据的对应

上面的定理说明：**数学上，面内任意固定参考都等价**。那么两个探针的差别到底在哪？差别不在“参考是否在轨迹上”，而在**参考是否在 winding 平面内**——以及参考是否随系统自动就位。

### 3.1 Bragg：参考自动在平面内

Bragg 参考是初态下带 Bloch 矢量 $\hat{a}=\hat{d}(k_0)$。它**本来就在轨迹上**，而轨迹整条在平面内——所以无论手性轴 $\hat n_\Gamma$ 取什么方向，Bragg 的参考**自动**满足“面内固定参考”条件，不需要任何对准。带间权重
$$
w_B(k_0,k)=\frac12\bigl[1-\hat{d}(k)\cdot\hat{d}(k_0)\bigr]
=\frac12\bigl[1-\cos\Delta\phi_{k,k_0}\bigr].
$$
由全覆盖定理：$W\neq0$ $\Rightarrow$ 对每个 $k_0$ 存在 $\pi$ 点 $\Rightarrow w_B=1$（满转移）。但原路返回时 $\pi$ 点仍在像内（$W=0$ 仍满转移），故

> **Bragg 满转移是 $W\neq0$ 的必要非充分条件。**

Bragg 的实际优势有两点，都不改变数学内容：其一，参考自动在面内（无需知道或对准平面）；其二，$k_0$ 可变，给出一族参考，任一 $k_0$ 不满转移即可一票否决 $W\neq0$。

### 3.2 ARPES：参考是固定外部轴，未必在平面内

ARPES 的参考是固定内态轴 $\hat z$，**不受系统对称性支配**。分解 $\hat z=\rho\,\hat z_\parallel+\sqrt{1-\rho^2}\,\hat n_\Gamma$（$\rho=|\hat z_\parallel|$ 为面内投影长度）。因 $\hat{d}\cdot\hat n_\Gamma=0$：
$$
\hat{d}(k)\cdot\hat z=\rho\cos(\phi_k-\phi_z),
\qquad
w_A(k)=\frac12\bigl[1+\rho\cos(\phi_k-\phi_z)\bigr].
$$
（核验 c5：随机 $\hat z$ 下该式精确成立。）于是：

| $\rho$ | 几何 | ARPES 能否看到 $\pi$ 点（$\hat{d}=-\hat z$） |
|---|---|---|
| $0$ | $\hat z\parallel\hat n_\Gamma$（轴即手性轴） | 不能；权重恒 $\tfrac12$，全盲 |
| $0<\rho<1$ | $\hat z$ 斜出平面 | 不能；$\min_k\hat{d}\cdot\hat z=-\rho>-1$，权重只在 $\tfrac12(1\pm\rho)$ 内 |
| $1$ | $\hat z$ 在 winding 平面内 | 能；此时 $\hat z$ 正是 §2 的“面内固定参考”，$W\neq0$ 时像为整圆，$\pm\hat z$ 都在像内，权重取到 $0$ 与 $1$ |

所以“ARPES 不行是因为 $d_z$ 不在轨迹上”的更准确表述是：**ARPES 的固定轴 $\hat z$ 只有恰落在 winding 平面内（$\hat n_\Gamma\cdot\hat z=0$）时，才成为 §2 定理中的合格参考；一般情形它只有长度 $\rho$ 的面内投影，看到的是被 $\rho$ 折减的影子。** 一个假想的“参考轴可任意调向面内”的探针，在数学上与 Bragg 等价（覆盖与 $\pi$ 点判据完全相同）；标准探针中不存在这种自由度——ARPES 的轴由晶格钉死，Bragg 的参考由初态自动就位。

### 3.3 两者都不行：像到了 $\pi$ 但折回

对应 2.5 节反例：像覆盖整圆、含 $\pi$ 点——Bragg 可满转移、ARPES 可见符号变化——但 $W=0$。数学根源：**满射 $\neq$ 非零度数**。任何只依赖“像是否覆盖 $[0,2\pi)$”的判据都无法单独确定 $W$。

---

## 4. 数值核验（`261109_verify_reference_free_winding.py`，全部通过）

模型：随机手性轴 $\hat n_\Gamma$、随机正半径 $r(k)=1+0.3\cos2k$（半径必须不影响任何结论）；轨迹 T（$\phi=k+0.5\sin k$，$W=1$）、A（$\phi=2k$，$W=2$）、R（原路返回，$W=0$）。

| 断言 | 条目 | 结果 |
|---|---|---|
| $\hat n_\Gamma\cdot(\hat{d}\times\partial_k\hat{d})=\partial_k\tilde\phi$（随机轴、随机半径） | c1 | $\le4.2\times10^{-10}$（差分精度） |
| 积分 $W$ = 提升差 = 各轨迹参考 $k_0$ 的相对角度数 | c2 | T/A/R 全部 $\le5.7\times10^{-11}$ |
| $W\neq0$ ⟹ 每个轨迹参考覆盖整圆、命中 $\pi$ 点 | c3 | 最大圆间隙 $\le3.1\times10^{-5}$（=网格间距）；$\pi$ 点 $\le3.7\times10^{-11}$ |
| 原路返回：$W=0$ 但整圆覆盖、$\pi$ 点命中 | c4 | $W=0$，覆盖与 $\pi$ 点照常成立 |
| $\hat{d}\cdot\hat z=\rho\cos(\phi-\phi_z)$；$\pi$ 点可达 $\iff\rho=1$ | c5 | 随机 $\hat z$：$\min\hat{d}\cdot\hat z=-\rho$；面内 $\hat z$：$=-1.000$ |
| **任意面内固定参考（不在轨迹上）与轨迹参考完全等价** | c6 | 度数偏差 $\le3.3\times10^{-10}$；覆盖与 $\pi$ 点与 c3/c4 逐项相同 |

---

## 5. 小结

| 说法 | 数学表述 |
|---|---|
| 定义（`Note.tex`） | $W=\dfrac{1}{2\pi}\displaystyle\int_{-\pi}^{\pi}dk\,\hat n_\Gamma\cdot(\hat{d}\times\partial_k\hat{d})$ |
| 化为角度提升 | $\hat n_\Gamma\cdot(\hat{d}\times\partial_k\hat{d})=\partial_k\tilde\phi$，$W=\dfrac{\tilde\phi_\pi-\tilde\phi_{-\pi}}{2\pi}$ |
| 换参考不改变绕数 | $\deg\bigl(k\mapsto\phi_k-\phi_a\bigr)=W$，对平面内一切固定 $\hat a$（轨迹点是特例） |
| $W\neq0$ 的参考无关特征 | 对一切面内 $\hat a$，$\{\phi_k-\phi_a\}=S^1$ |
| $\pi$ 点 | $\exists k:\ \phi_k-\phi_a=\pi\iff-\hat a\in C$ |
| Bragg | 参考在轨迹上 ⟹ 自动面内（任意手性轴方向皆然）；$\pi$ 点即满转移；必要非充分 |
| ARPES | 参考是固定 $\hat z$；$\hat{d}\cdot\hat z=\rho\cos(\phi-\phi_z)$，仅 $\rho=1$（$\hat z$ 恰在平面内）时才是合格参考 |
| 两者都失败 | 像可覆盖整圆但 $W=0$（原路返回）；满射 $\neq$ 非零度数 |

一句话：

> **从二能级定义出发：手性对称把 $\vec{d}(k)$ 压进垂直手性轴的平面，绕数的手性轴分量自动等于角度提升的度数；因而 $W\neq0$ 时，以平面内任意固定矢量为参考（轨迹点只是特例），相对角都必走满 $0$ 到 $2\pi$ 并含 $\pi$ 点。探针之间的真正分界不是“参考在不在轨迹上”，而是“参考在不在平面内”：Bragg 的参考自动面内，ARPES 的固定轴 $\hat z$ 一般只有长度 $\rho$ 的面内投影。而“走满”本身只说明像覆盖整圆，确定 $W$ 还须数带符号的净圈数。**

---

*笔记日期：2026-11-09。核验脚本：`261109_verify_reference_free_winding.py`（同目录，c1–c6 全部通过）。*
