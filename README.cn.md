# 🧲 Physics++ 物理计算框架

> **引言：从数字到物理世界的桥梁**

在经典力学中，牛顿第二定律告诉我们：**F = ma**。然而，当我们要用计算机模拟一个力场在空间中的传播时，这个简单的公式瞬间变得复杂起来——我们需要离散化空间、离散化时间、处理无穷远处的边界条件、确保数值稳定性。这正是 **Physics++** 诞生的原因：**为物理量赋予计算生命**。

---

## 📐 数学基础层：`mathematics.h`

### 🛡️ 安全函数：从极限定理到数值鲁棒性

在实数分析中，函数极限的定义要求：对于任意 ε > 0，存在 δ > 0。但计算机中的浮点数遵循 **IEEE 754**，存在两个违反极限定理的特殊值：**NaN**（非数）和 **Inf**（无穷）。

一旦 NaN 或 Inf 进入计算图，根据**错误传播定理**，整个计算结果将不可信。因此，`mathematics.h` 中的每一个函数都在运算前后进行异常检测：

```cpp
若 x = NaN ∨ x = Inf → 抛出 domain_error
若 result = NaN ∨ result = Inf → 抛出 overflow_error
```

> 这一设计基于 **契约式编程**：前置条件保护输入，后置条件保护输出。

---

### 🔁 超越函数：从泰勒展开到安全封装

| 函数 | 数学基础 |
|----|----|
| `exponential(x)` | 泰勒级数：$e^x = \sum_{n=0}^{\infty} \frac{x^n}{n!}$ |
| `sine(x)` / `cosine(x)` | 三角恒等式 + 周期性归约 |

- 指数函数调用 `std::exp` 后执行溢出检查；
- 三角函数利用 $\sin(-x) = -\sin(x)$ 和 $2\pi$ 周期性，将参数归约至收敛最快区间。

---

### 🧰 数值稳定工具

| 工具 | 数学原理 | 作用 |
|----|----|----|
| `exact(x, ε)` | 零判定定理 | 消除 $10^{-16}$ 级浮点残余误差 |
| `limit(x, low, high)` | 夹逼定理 | 防止数值漂移出物理合理范围 |
| `smooth(x, t, s)` | $\tanh$ S 形曲线 | 提供连续可导的平滑过渡 |

$$
\text{smooth}(x) = \frac{1 + \tanh(s \cdot (x - t) / 2)}{2}
$$

---

## 🧮 代数结构层：`quaternion.h`

### 🌀 复数：从代数闭域到信号表示

根据**代数基本定理**，任何 $n$ 次多项式在复数域 $\mathbb{C}$ 中恰好有 $n$ 个根。

- `(ℂ, +)`：加法群  
- `(ℂ\{0\}, ×)`：乘法群  

**欧拉公式**：

$$
e^{j\theta} = \cos\theta + j\sin\theta
$$

`exponential()` 将复数指数转化为实指数与三角函数的乘积：

$$
e^{a+jb} = e^a \cdot (\cos b + j\sin b)
$$

> 这是频域分析与傅里叶变换的理论基石。

---

### 🎯 四元数：从旋转群到三维场表示

哈密顿于 1843 年在都柏林布鲁姆桥上刻下：

$$
i^2 = j^2 = k^2 = ijk = -1
$$

四元数构成**非交换可除代数**，是描述三维旋转的理想结构：

- 不满足交换律：$q_1 q_2 \neq q_2 q_1$
- 满足结合律
- 编码“轴–角”旋转（欧拉旋转定理）

**球面线性插值（SLERP）**：

$$
\text{slerp}(q_1, q_2, t) =
\frac{\sin((1-t)\theta)}{\sin\theta} q_1 +
\frac{\sin(t\theta)}{\sin\theta} q_2
$$

其中 $\theta = \arccos(q_1 \cdot q_2)$。

---

### 🧩 自由函数：类型匹配的多态

`quaternion.h` 提供统一的自由函数接口：

- `module`
- `conjugate`
- `exponential`
- `logarithm`

利用 **C++ 函数重载**，同时支持：

- `Complex<N>`
- `Quaternion<N>`
- `std::complex<N>`

> 高层代码（如 `Wave`、`Field`）无需关心具体类型，只管调用。

---

## ⚖️ 物理量纲层：`quantization.h`

### 📏 量纲分析：从 π 定理到类型安全

**白金汉 π 定理**：任何物理定律都可表示为无量纲参数的函数。

| 运算 | 结果 |
|----|----|
| Time + Time | ✅ Time |
| Time + Resistance | ❌ 编译错误 |
| Voltage = Current × Resistance | ✅ |
| Power = Voltage × Current | ✅ |
| Energy = Power × Time | ✅ |

---

### 🌊 频域物理量

拉普拉斯变量：

$$
s = \sigma + j\omega
$$

- $\sigma$：衰减率（neper/s）
- $\omega$：角频率（rad/s）

品质因数与阻尼比：

$$
Q = \frac{1}{2\zeta}
$$

**相量（Phasor）** 基于欧拉公式：

$$
A\cos(\omega t + \phi) \longleftrightarrow A e^{j\phi}
$$

---

### 🧭 矢量物理量

内部存储为 `Quaternion<double>`（标量部分为 0）：

- 点积：$\mathbf{a} \cdot \mathbf{b}$
- 叉积：$\mathbf{a} \times \mathbf{b}$

物理关系：

- **牛顿第二定律**：$\mathbf{F} = m\mathbf{a}$
- **功的定义**：$W = \mathbf{F} \cdot \mathbf{s}$

---

## 🌊 信号处理层：`wave.h`

### 🎼 时域与频域：傅里叶分析

**傅里叶定理**：任何周期信号都可分解为正弦波的加权和。

连续傅里叶变换：

$$
F(\omega) = \int_{-\infty}^{\infty} f(t) e^{-j\omega t} dt
$$

离散傅里叶变换（DFT）：

$$
X[k] = \sum_{n=0}^{N-1} x[n] \cdot e^{-j2\pi kn/N}
$$

- `discreteFourierTransform()`：直接实现
- `fastFourierTransform()`：Cooley–Tukey 算法，复杂度 $O(N \log N)$

---

### 📊 频谱的两种表示

| 表示 | 插值 | 用途 |
|----|----|----|
| `Spectrum<N>` | ❌ | 原始变换结果 |
| `Domain<T>` | ✅ | 分析与可视化 |

$$
f(x) = f_1 + (f_2 - f_1) \cdot \frac{x - x_1}{x_2 - x_1}
$$

---

### 📐 从离散到连续：分辨率与置信度

`Spectral<T>` 额外携带：

- **分辨率**：基于不确定性原理 $\Delta f \cdot \Delta t \ge 1/(4\pi)$
- **置信度**：谱线估计的可靠性 $[0, 1]$

---

### 🔄 丰富的变换

| 变换 | 公式 | 应用 |
|----|----|----|
| FFT | $X[k] = \sum_n x[n] W_N^{kn}$ | 通用频谱分析 |
| DCT | 余弦基 | 图像 / 音频压缩 |
| DST | 正弦基 | 边界值问题 |
| Hilbert | Cauchy 主值积分 | 包络检测 |
| Chirp-Z | 螺旋采样 | 细粒度频谱放大 |
| Cepstrum | 对数谱反变换 | 回波识别 |

---

### 🧠 类型分派

通过 `enable_if_t` 实现返回类型分派：

- **算术类型** → `Domain<Complex<double>>`
- **物理类型** → `Spectrum<N>`

> 算术类型不携带量纲，适合连续分析；  
> 物理类型携带量纲，适合离散存储。

---

## 🌐 场论层：`field.h`

### 🌊 四元数波：从波动方程到离散表示

$$
\mathbf{Q}(t) = (s(t),\; v_x(t),\; v_y(t),\; v_z(t))
$$

- 标量部分：势（压力、电势）
- 矢量部分：有向梯度（电场、速度）

---

### ⏱️ 力学量的时间演化

一阶导数（中心差分，$O(\Delta t^2)$）：

$$
v(t) \approx \frac{u(t+\Delta t) - u(t-\Delta t)}{2\Delta t}
$$

二阶导数（五点模板，$O(\Delta t^4)$）：

$$
a(t) \approx \frac{-u(t+2\Delta t) + 16u(t+\Delta t) - 30u(t) + 16u(t-\Delta t) - u(t-2\Delta t)}{12\Delta t^2}
$$

---

### 🪨 扰动算法：共振涟漪

波动方程含源项：

$$
\frac{\partial^2 \mathbf{Q}}{\partial t^2} + \gamma \frac{\partial \mathbf{Q}}{\partial t} = c^2 \nabla^2 \mathbf{Q} + \mathbf{S}(t)
$$

**路径一：空场初始扰动（高斯扩散）**

$$
\mathbf{Q}'(t_i) = \mathbf{Q}(t_i) + \sum_j \mathbf{S}(t_j) \cdot \exp\left(-\frac{(i-j)^2}{2\sigma^2}\right)
$$

**路径二：共振调制扰动**

$$
R(t) = \frac{|\mathbf{Q}(t) \cdot \mathbf{S}(t)|}{\max_t |\mathbf{Q}(t) \cdot \mathbf{S}(t)|}
$$

$$
\mathbf{Q}'(t) = \mathbf{Q}(t) + R(t) \cdot \mathbf{S}(t)
$$

---

### 🔋 能量与熵

**能量（哈密顿力学）**：

$$
E = \frac{1}{2} \sum_i \left[ \left(\frac{du_i}{dt}\right)^2 + |\nabla u_i|^2 \right]
$$

**熵（香农信息论）**：

$$
H = -\sum_i p_i \ln p_i,\quad
p_i = \frac{|\mathbf{Q}_i|^2}{\sum_j |\mathbf{Q}_j|^2}
$$

- 熵大 → 能量分布均匀  
- 熵小 → 能量集中于局部（如涟漪峰）

---

### 🧬 物理类型特化

`Field<Voltage>`、`Field<Force>` 等通过**模板特化**实现：

- 外层：物理类型接口
- 内层：数值类型存储（`Complex<double>` / `Quaternion<double>`）

> 这是**桥梁模式**的编译时版本。

---

## 🏗️ 架构总览

```
mathematics.h        ← 安全数学函数（极限定理 + IEEE 754）
        ↓
quaternion.h         ← 复数 + 四元数 + 自由函数（代数结构）
        ↓
quantization.h       ← 物理量纲类型（量纲齐次性）
        ↓
wave.h               ← 频域 / 时域分析（傅里叶定理）
        ↓
field.h              ← 四元数波场（波动方程 + 共振理论）
```

> 从最底层的 `add(a,b)` 到最顶层的 `Field<Force>::perturb()`，  
> **每一步计算都有理论依据，每一个结果都可追溯、可验证。**
