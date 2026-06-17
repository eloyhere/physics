# 🧲 Physics++: A Physical Computation Framework

> **Foreword: Bridging the Gap Between Numbers and Reality**

In classical mechanics, Newton’s Second Law tells us: **F = ma**. Simple enough—until you try to simulate how a force field propagates through space using a computer. Suddenly, you’re juggling spatial discretisation, time-stepping, boundary conditions at infinity, and numerical stability. That’s precisely why **Physics++** exists: to breathe computational life into physical quantities.

---

## 📐 Mathematical Foundations: `mathematics.h`

### 🛡️ Safe Functions: From Limit Theorems to Numerical Robustness

In real analysis, the definition of a function limit states that for every ε > 0, there exists a δ > 0. But computers don’t deal in pure reals—they use **IEEE 754 floating point**, complete with two troublemakers: **NaN** (Not-a-Number) and **Inf** (Infinity).

Once NaN or Inf sneaks into your computation graph, the **error propagation theorem** ensures the entire result becomes untrustworthy. Hence, every function in `mathematics.h` performs rigorous sanity checks:

```cpp
if x == NaN OR x == Inf → throw domain_error
if result == NaN OR result == Inf → throw overflow_error
```

> This follows **Design by Contract**: preconditions guard inputs, postconditions guard outputs.

---

### 🔁 Transcendental Functions: From Taylor Expansions to Safe Wrappers

| Function | Mathematical Basis |
|--------|-------------------|
| `exponential(x)` | Taylor series: $e^x = \sum_{n=0}^{\infty} \frac{x^n}{n!}$ |
| `sine(x)` / `cosine(x)` | Trigonometric identities + periodicity reduction |

- The exponential implementation wraps `std::exp` with overflow detection.
- Trigonometric functions exploit $\sin(-x) = -\sin(x)$ and $2\pi$-periodicity to reduce arguments to their fastest-converging intervals.

---

### 🧰 Numerical Stability Utilities

| Utility | Principle | Purpose |
|------|---------|--------|
| `exact(x, ε)` | Zero-detection theorem | Eliminates floating-point noise (~$10^{-16}$) |
| `limit(x, low, high)` | Squeeze theorem | Prevents numerical drift outside physical bounds |
| `smooth(x, t, s)` | $\tanh$ sigmoid curve | Continuous, differentiable transitions |

$$
\text{smooth}(x) = \frac{1 + \tanh(s \cdot (x - t) / 2)}{2}
$$

---

## 🧮 Algebraic Structures: `quaternion.h`

### 🌀 Complex Numbers: From Algebraically Closed Fields to Signal Representation

The **Fundamental Theorem of Algebra** guarantees that any $n$th-degree polynomial has exactly $n$ roots in the complex plane $\mathbb{C}$.

- `(ℂ, +)`: Additive group  
- `(ℂ\{0\}, ×)`: Multiplicative group  

**Euler’s Formula**:

$$
e^{j\theta} = \cos\theta + j\sin\theta
$$

The `exponential()` method transforms complex exponentiation into real exponentiation multiplied by trigonometry:

$$
e^{a+jb} = e^a \cdot (\cos b + j\sin b)
$$

> This is the bedrock of frequency-domain analysis and Fourier theory.

---

### 🎯 Quaternions: From Rotation Groups to 3D Field Representation

Hamilton carved into Broome Bridge in Dublin in 1843:

$$
i^2 = j^2 = k^2 = ijk = -1
$$

Quaternions form a **non-commutative division algebra**—perfect for describing 3D rotations:

- Non-commutative: $q_1 q_2 \neq q_2 q_1$
- Associative
- Encodes axis-angle rotations (Euler’s rotation theorem)

**Spherical Linear Interpolation (SLERP)**:

$$
\text{slerp}(q_1, q_2, t) =
\frac{\sin((1-t)\theta)}{\sin\theta} q_1 +
\frac{\sin(t\theta)}{\sin\theta} q_2
$$

where $\theta = \arccos(q_1 \cdot q_2)$.

---

### 🧩 Free Functions: Type-Aware Polymorphism

`quaternion.h` exposes free functions such as:

- `module`
- `conjugate`
- `exponential`
- `logarithm`

Using **C++ overload resolution**, these support:

- `Complex<N>`
- `Quaternion<N>`
- `std::complex<N>`

> Higher-level code (`Wave`, `Field`) remains blissfully type-agnostic.

---

## ⚖️ Physical Dimensional Layer: `quantization.h`

### 📏 Dimensional Analysis: From the Buckingham π Theorem to Type Safety

The **Buckingham π Theorem** states that any physical law can be expressed as a function of dimensionless parameters.

| Operation | Result |
|--------|--------|
| Time + Time | ✅ Time |
| Time + Resistance | ❌ Compile-time error |
| Voltage = Current × Resistance | ✅ |
| Power = Voltage × Current | ✅ |
| Energy = Power × Time | ✅ |

---

### 🌊 Frequency-Domain Quantities

Laplace variable:

$$
s = \sigma + j\omega
$$

- $\sigma$: Attenuation rate (neper/s)
- $\omega$: Angular frequency (rad/s)

Quality factor vs damping ratio:

$$
Q = \frac{1}{2\zeta}
$$

**Phasors** represent sinusoids in steady state via Euler’s formula:

$$
A\cos(\omega t + \phi) \longleftrightarrow A e^{j\phi}
$$

---

### 🧭 Vector Quantities

Internally stored as `Quaternion<double>` (scalar part = 0):

- Dot product: $\mathbf{a} \cdot \mathbf{b}$
- Cross product: $\mathbf{a} \times \mathbf{b}$

Physical relationships:

- **Newton’s Second Law**: $\mathbf{F} = m\mathbf{a}$
- **Work-energy principle**: $W = \mathbf{F} \cdot \mathbf{s}$

---

## 🌊 Signal Processing Layer: `wave.h`

### 🎼 Time and Frequency Domains: Fourier Analysis

**Fourier’s Theorem**: Any periodic signal can be decomposed into a weighted sum of sinusoids.

Continuous Fourier Transform:

$$
F(\omega) = \int_{-\infty}^{\infty} f(t) e^{-j\omega t} dt
$$

Discrete Fourier Transform (DFT):

$$
X[k] = \sum_{n=0}^{N-1} x[n] \cdot e^{-j2\pi kn/N}
$$

- `discreteFourierTransform()`: Direct implementation
- `fastFourierTransform()`: Cooley–Tukey algorithm, $O(N \log N)$

---

### 📊 Two Flavours of Spectrum

| Representation | Interpolation | Use Case |
|--------------|-------------|----------|
| `Spectrum<N>` | ❌ | Raw transform output |
| `Domain<T>` | ✅ | Analysis & visualisation |

Linear interpolation based on Lagrange’s Mean Value Theorem:

$$
f(x) = f_1 + (f_2 - f_1) \cdot \frac{x - x_1}{x_2 - x_1}
$$

---

### 📐 Resolution and Confidence

`Spectral<T>` carries two meta-values:

- **Resolution**: Governed by the uncertainty principle $\Delta f \cdot \Delta t \ge 1/(4\pi)$
- **Confidence**: Reliability estimate of spectral lines, ranging from 0 to 1

---

### 🔄 Rich Set of Transforms

| Transform | Formula | Application |
|---------|--------|-------------|
| FFT | $X[k] = \sum_n x[n] W_N^{kn}$ | General spectrum analysis |
| DCT | Cosine basis | Image & audio compression |
| DST | Sine basis | Boundary value problems |
| Hilbert | Cauchy principal value integral | Envelope detection |
| Chirp-Z | Spiral sampling | Zoomed-in spectral analysis |
| Cepstrum | Log-spectrum inverse transform | Echo & pitch detection |

---

### 🧠 Type Dispatch

Return types are dispatched using `enable_if_t`:

- **Arithmetic types** → `Domain<Complex<double>>`
- **Physical types** → `Spectrum<N>`

> Arithmetic types lack dimensions—ideal for continuous analysis.  
> Physical types carry units—ideal for discrete storage.

---

## 🌐 Field Theory Layer: `field.h`

### 🌊 Quaternion Waves: From Wave Equations to Discrete Representations

$$
\mathbf{Q}(t) = (s(t),\; v_x(t),\; v_y(t),\; v_z(t))
$$

- Scalar part: Potential (pressure, electric potential)
- Vector part: Directed gradient (electric field, velocity)

---

### ⏱️ Temporal Evolution of Mechanical Quantities

First derivative (central difference, $O(\Delta t^2)$):

$$
v(t) \approx \frac{u(t+\Delta t) - u(t-\Delta t)}{2\Delta t}
$$

Second derivative (five-point stencil, $O(\Delta t^4)$):

$$
a(t) \approx \frac{-u(t+2\Delta t) + 16u(t+\Delta t) - 30u(t) + 16u(t-\Delta t) - u(t-2\Delta t)}{12\Delta t^2}
$$

---

### 🪨 Perturbation Algorithms: Resonant Ripples

Wave equation with source term:

$$
\frac{\partial^2 \mathbf{Q}}{\partial t^2} + \gamma \frac{\partial \mathbf{Q}}{\partial t} = c^2 \nabla^2 \mathbf{Q} + \mathbf{S}(t)
$$

**Path 1: Initial perturbation into an empty field (Gaussian diffusion)**

$$
\mathbf{Q}'(t_i) = \mathbf{Q}(t_i) + \sum_j \mathbf{S}(t_j) \cdot \exp\left(-\frac{(i-j)^2}{2\sigma^2}\right)
$$

**Path 2: Resonant modulation**

$$
R(t) = \frac{|\mathbf{Q}(t) \cdot \mathbf{S}(t)|}{\max_t |\mathbf{Q}(t) \cdot \mathbf{S}(t)|}
$$

$$
\mathbf{Q}'(t) = \mathbf{Q}(t) + R(t) \cdot \mathbf{S}(t)
$$

---

### 🔋 Energy and Entropy

**Energy (Hamiltonian Mechanics)**:

$$
E = \frac{1}{2} \sum_i \left[ \left(\frac{du_i}{dt}\right)^2 + |\nabla u_i|^2 \right]
$$

**Entropy (Shannon Information Theory)**:

$$
H = -\sum_i p_i \ln p_i,\quad
p_i = \frac{|\mathbf{Q}_i|^2}{\sum_j |\mathbf{Q}_j|^2}
$$

- High entropy → energy evenly spread  
- Low entropy → energy concentrated (e.g., ripple peaks)

---

### 🧬 Physical Type Specialisation

`Field<Voltage>`, `Field<Force>`, etc., are implemented via **template specialisation**:

- Outer layer: Physical-type interface
- Inner layer: Numerical storage (`Complex<double>` / `Quaternion<double>`)

> This is a compile-time rendition of the **Bridge Pattern**.

---

## 🏗️ Architectural Overview

```
mathematics.h        ← Safe mathematical functions (limits + IEEE 754)
        ↓
quaternion.h         ← Complex numbers, quaternions, free functions
        ↓
quantization.h       ← Physical dimensional types (dimensional homogeneity)
        ↓
wave.h               ← Time/frequency domain analysis (Fourier theory)
        ↓
field.h              ← Quaternion wave fields (wave equations + resonance)
```

> From the humble `add(a,b)` all the way up to `Field<Force>::perturb()`,  
> **every computation rests on firm theoretical ground, and every result is traceable and verifiable.**
