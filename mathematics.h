#pragma once
#include <complex>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>
#include <map>
#include <queue>

namespace mathematics
{
const double pi = std::acos(-1.0);
const double tolerance = 1e-12;

template <typename N>
inline bool abnormal(N value)
{
    if constexpr (std::is_arithmetic_v<N>)
    {
        return std::isnan(value) || std::isinf(value);
    }
    throw std::domain_error("abnormal: Must be arithmetic type.");
}

template <typename N>
N absolute(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("absolute: NaN or Inf input");
    }
    if (value < N(0))
    {
        return -value;
    }
    return value;
}

template <typename N>
N negative(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("negative: NaN or Inf input");
    }
    return -value;
}

template <typename N>
N modulus(N value, N divisor)
{
    if (abnormal(value) || abnormal(divisor))
    {
        throw std::domain_error("modulus: NaN or Inf input");
    }
    if (divisor == N(0))
    {
        throw std::domain_error("modulus: division by zero");
    }
    N result = value % divisor;
    if (result < N(0))
    {
        result = result + absolute(divisor);
    }
    return result;
}

template <typename N>
N floor(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("floor: NaN or Inf input");
    }
    return std::floor(value);
}

template <typename N>
N ceil(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("ceil: NaN or Inf input");
    }
    return std::ceil(value);
}

template <typename N>
N round(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("round: NaN or Inf input");
    }
    return std::round(value);
}

template <typename N>
N sign(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("sign: NaN or Inf input");
    }
    if (value > N(0))
    {
        return N(1);
    }
    if (value < N(0))
    {
        return N(-1);
    }
    return N(0);
}

template <typename N>
N squareRoot(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("squareRoot: NaN or Inf input");
    }
    if (value < N(0))
    {
        throw std::domain_error("squareRoot: negative input");
    }
    N result = std::sqrt(value);
    if (abnormal(result))
    {
        throw std::overflow_error("squareRoot: result overflow");
    }
    return result;
}

template <typename N>
N exponential(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("exponential: NaN or Inf input");
    }
    N result = std::exp(value);
    if (abnormal(result))
    {
        throw std::overflow_error("exponential: result overflow");
    }
    return result;
}

template <typename N>
N logarithm(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("logarithm: NaN or Inf input");
    }
    if (value <= N(0))
    {
        throw std::domain_error("logarithm: non-positive input");
    }
    N result = std::log(value);
    if (abnormal(result))
    {
        throw std::overflow_error("logarithm: result overflow");
    }
    return result;
}

template <typename N>
N logarithm(N base, N value)
{
    if (abnormal(base) || abnormal(value))
    {
        throw std::domain_error("logarithm: NaN or Inf input");
    }
    if (base <= N(0) || absolute(base - N(1)) < tolerance)
    {
        throw std::domain_error("logarithm: base must be positive and not 1");
    }
    if (value <= N(0))
    {
        throw std::domain_error("logarithm: value must be positive");
    }
    N result = logarithm(value) / logarithm(base);
    if (abnormal(result))
    {
        throw std::overflow_error("logarithm: result overflow");
    }
    return result;
}

template <typename N>
N power(N base, N exponent)
{
    if (abnormal(base) || abnormal(exponent))
    {
        throw std::domain_error("power: NaN or Inf input");
    }
    if (base < N(0) && std::abs(exponent - std::round(exponent)) > tolerance)
    {
        throw std::domain_error("power: negative base with non-integer exponent");
    }
    if (absolute(base) < tolerance && exponent <= N(0))
    {
        throw std::domain_error("power: zero base with non-positive exponent");
    }
    N result = std::pow(base, exponent);
    if (abnormal(result))
    {
        throw std::overflow_error("power: result overflow");
    }
    return result;
}

template <typename N>
N sine(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("sine: NaN or Inf input");
    }
    N result = std::sin(value);
    if (abnormal(result))
    {
        throw std::overflow_error("sine: result overflow");
    }
    return result;
}

template <typename N>
N cosine(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("cosine: NaN or Inf input");
    }
    N result = std::cos(value);
    if (abnormal(result))
    {
        throw std::overflow_error("cosine: result overflow");
    }
    return result;
}

template <typename N>
N tangent(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("tangent: NaN or Inf input");
    }
    N result = std::tan(value);
    if (abnormal(result))
    {
        throw std::overflow_error("tangent: result overflow");
    }
    return result;
}

template <typename N>
N hyperbolicSine(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("hyperbolicSine: NaN or Inf input");
    }
    N result = std::sinh(value);
    if (abnormal(result))
    {
        throw std::overflow_error("hyperbolicSine: result overflow");
    }
    return result;
}

template <typename N>
N hyperbolicCosine(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("hyperbolicCosine: NaN or Inf input");
    }
    N result = std::cosh(value);
    if (abnormal(result))
    {
        throw std::overflow_error("hyperbolicCosine: result overflow");
    }
    return result;
}

template <typename N>
N hyperbolicTangent(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("hyperbolicTangent: NaN or Inf input");
    }
    N result = std::tanh(value);
    if (abnormal(result))
    {
        throw std::overflow_error("hyperbolicTangent: result overflow");
    }
    return result;
}

template <typename N>
N arcsine(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("arcsine: NaN or Inf input");
    }
    if (value < N(-1) || value > N(1))
    {
        throw std::domain_error("arcsine: argument out of [-1, 1]");
    }
    N result = std::asin(value);
    if (abnormal(result))
    {
        throw std::overflow_error("arcsine: result overflow");
    }
    return result;
}

template <typename N>
N arccosine(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("arccosine: NaN or Inf input");
    }
    if (value < N(-1) || value > N(1))
    {
        throw std::domain_error("arccosine: argument out of [-1, 1]");
    }
    N result = std::acos(value);
    if (abnormal(result))
    {
        throw std::overflow_error("arccosine: result overflow");
    }
    return result;
}

template <typename N>
N arctangent(N value)
{
    if (abnormal(value))
    {
        throw std::domain_error("arctangent: NaN or Inf input");
    }
    N result = std::atan(value);
    if (abnormal(result))
    {
        throw std::overflow_error("arctangent: result overflow");
    }
    return result;
}

template <typename N>
N arctangent2(N y, N x)
{
    if (abnormal(y) || abnormal(x))
    {
        throw std::domain_error("arctangent2: NaN or Inf input");
    }
    if (absolute(y) < tolerance && absolute(x) < tolerance)
    {
        throw std::domain_error("arctangent2: both arguments zero");
    }
    N result = std::atan2(y, x);
    if (abnormal(result))
    {
        throw std::overflow_error("arctangent2: result overflow");
    }
    return result;
}

template <typename N>
N add(N left, N right)
{
    N result = left + right;
    if (abnormal(result))
    {
        throw std::overflow_error("add: result overflow");
    }
    return result;
}

template <typename N>
N subtract(N left, N right)
{
    N result = left - right;
    if (abnormal(result))
    {
        throw std::overflow_error("subtract: result overflow");
    }
    return result;
}

template <typename N>
N multiply(N left, N right)
{
    N result = left * right;
    if (abnormal(result))
    {
        throw std::overflow_error("multiply: result overflow");
    }
    return result;
}

template <typename N>
N divide(N numerator, N denominator)
{
    if (abnormal(numerator) || abnormal(denominator))
    {
        throw std::domain_error("divide: NaN or Inf input");
    }
    if (absolute(denominator) < tolerance)
    {
        throw std::domain_error("divide: division by zero");
    }
    N result = numerator / denominator;
    if (abnormal(result))
    {
        throw std::overflow_error("divide: result overflow");
    }
    return result;
}

template <typename N>
int cycle(N index, int limit)
{
    if (limit < 1)
    {
        throw std::domain_error("cycle: limit must be >= 1");
    }
    int intIndex = static_cast<int>(index);
    return (intIndex % limit + limit) % limit;
}

template <typename N>
N exact(N value, N epsilon = N(tolerance))
{
    if (abnormal(value))
    {
        throw std::domain_error("exact: NaN or Inf input");
    }
    if (absolute(value) < epsilon)
    {
        return N(0);
    }
    return value;
}

template <typename N>
N limit(N value, N start, N end)
{
    if (abnormal(start) || abnormal(end))
    {
        return N(0);
    }
    if (abnormal(value))
    {
        return (start + end) / N(2);
    }
    if (value >= start && value <= end)
    {
        return value;
    }
    N range = end - start;
    if (absolute(range) < tolerance)
    {
        return (start + end) / N(2);
    }
    N middle = (start + end) / N(2);
    if (value < start)
    {
        return start + (middle - start) * arctangent((start - value) / range) * N(2) / N(pi);
    }
    return end - (end - middle) * arctangent((value - end) / range) * N(2) / N(pi);
}

template <typename N>
N smooth(N value, N threshold, N steepness)
{
    if (abnormal(value) || abnormal(threshold) || abnormal(steepness))
    {
        return N(0);
    }
    N argument = steepness * (value - threshold);
    N tanhArgument = std::tanh(argument / N(2));
    if (abnormal(tanhArgument))
    {
        return N(0);
    }
    return (N(1) + tanhArgument) / N(2);
}

inline std::complex<double> euler(double angle)
{
    double cosValue = cosine(angle);
    double sinValue = sine(angle);
    if (abnormal(cosValue) || abnormal(sinValue))
    {
        throw std::overflow_error("euler: result overflow");
    }
    return std::complex<double>(cosValue, sinValue);
}

inline std::complex<double> rootOfUnity(size_t n, size_t k)
{
    double angle = (-2.0 * pi * static_cast<double>(k)) / static_cast<double>(n);
    return euler(angle);
}

inline size_t nextPowerOfTwo(size_t value)
{
    if (value == 0)
    {
        return 1;
    }
    size_t power = 1;
    while (power < value)
    {
        power <<= 1;
        if (power == 0)
        {
            throw std::overflow_error("nextPowerOfTwo: overflow");
        }
    }
    return power;
}

template <typename Container>
void bitReversalPermutation(Container &data)
{
    size_t size = data.size();
    if (size == 0)
    {
        return;
    }
    for (size_t i = 1, j = 0; i < size; ++i)
    {
        size_t bit = size >> 1;
        for (; j & bit; bit >>= 1)
        {
            j ^= bit;
        }
        j ^= bit;
        if (i < j)
        {
            std::swap(data[i], data[j]);
        }
    }
}

template <typename Container>
void fastFourierTransform1D(Container &data, bool inverse = false)
{
    size_t dataSize = data.size();
    if (dataSize <= 1)
    {
        return;
    }
    bitReversalPermutation(data);
    for (size_t length = 2; length <= dataSize; length <<= 1)
    {
        double angle = (2.0 * pi) / static_cast<double>(length);
        if (inverse)
        {
            angle = -angle;
        }
        double cosAngle = cosine(angle);
        double sinAngle = sine(angle);
        std::complex<double> wlen(cosAngle, sinAngle);
        for (size_t i = 0; i < dataSize; i += length)
        {
            std::complex<double> w(1.0, 0.0);
            for (size_t j = 0; j < length / 2; ++j)
            {
                std::complex<double> u = data[i + j];
                std::complex<double> v = data[i + j + length / 2] * w;
                data[i + j] = u + v;
                data[i + j + length / 2] = u - v;
                w = w * wlen;
            }
        }
    }
    if (inverse)
    {
        double scale = 1.0 / static_cast<double>(dataSize);
        for (auto &v : data)
        {
            v = v * scale;
        }
    }
    for (auto &v : data)
    {
        v = exact(v);
        if (abnormal(v))
        {
            throw std::overflow_error("fastFourierTransform1D: result overflow");
        }
    }
}

template <typename InputContainer>
auto discreteFourierTransform(const InputContainer &input) -> std::vector<std::complex<double>>
{
    size_t length = input.size();
    if (length == 0)
    {
        return {};
    }
    std::vector<std::complex<double>> result(length);
    double baseAngle = (-2.0 * pi) / static_cast<double>(length);
    for (size_t k = 0; k < length; ++k)
    {
        std::complex<double> sum(0.0, 0.0);
        for (size_t t = 0; t < length; ++t)
        {
            double angle = baseAngle * static_cast<double>(k) * static_cast<double>(t);
            double cosValue = cosine(angle);
            double sinValue = sine(angle);
            double inputValue = static_cast<double>(input[t]);
            sum.real(sum.real() + inputValue * cosValue);
            sum.imag(sum.imag() + inputValue * sinValue);
        }
        if (abnormal(sum))
        {
            throw std::overflow_error("discreteFourierTransform: result overflow");
        }
        result[k] = sum;
    }
    return result;
}

template <typename InputContainer>
void fastFourierTransform2D(InputContainer &data, bool inverse = false)
{
    size_t rows = data.size();
    if (rows == 0)
    {
        return;
    }
    size_t cols = data[0].size();
    if (cols == 0)
    {
        return;
    }
    for (size_t i = 0; i < rows; ++i)
    {
        fastFourierTransform1D(data[i], inverse);
    }
    std::vector<std::complex<double>> column(rows);
    for (size_t j = 0; j < cols; ++j)
    {
        for (size_t i = 0; i < rows; ++i)
        {
            column[i] = data[i][j];
        }
        fastFourierTransform1D(column, inverse);
        for (size_t i = 0; i < rows; ++i)
        {
            data[i][j] = column[i];
        }
    }
}

template <typename InputContainer>
auto shortTimeFourierTransform(const InputContainer &input, size_t windowSize, size_t hopSize)
    -> std::vector<std::vector<std::complex<double>>>
{
    size_t inputSize = input.size();
    if (windowSize == 0 || hopSize == 0 || inputSize == 0 || windowSize > inputSize)
    {
        return {};
    }
    if (windowSize == 1)
    {
        std::vector<std::vector<std::complex<double>>> spectrogram;
        for (size_t start = 0; start < inputSize; start += hopSize)
        {
            std::vector<std::complex<double>> frame(1);
            frame[0] = std::complex<double>(static_cast<double>(input[start]), 0.0);
            fastFourierTransform1D(frame);
            spectrogram.push_back(std::move(frame));
        }
        return spectrogram;
    }
    std::vector<std::vector<std::complex<double>>> spectrogram;
    std::vector<double> window(windowSize);
    double windowBase = (2.0 * pi) / static_cast<double>(windowSize - 1);
    for (size_t i = 0; i < windowSize; ++i)
    {
        window[i] = 0.5 * (1.0 - std::cos(windowBase * static_cast<double>(i)));
    }
    for (size_t start = 0; start + windowSize <= inputSize; start += hopSize)
    {
        std::vector<std::complex<double>> frame(windowSize);
        for (size_t i = 0; i < windowSize; ++i)
        {
            double value = static_cast<double>(input[start + i]) * window[i];
            frame[i] = std::complex<double>(value, 0.0);
        }
        fastFourierTransform1D(frame);
        spectrogram.push_back(std::move(frame));
    }
    return spectrogram;
}

template <typename InputContainer>
auto chirpZTransform(const InputContainer &input, double startFrequency, double endFrequency,
                     size_t numPoints, double sampleRate) -> std::vector<std::complex<double>>
{
    size_t N = input.size();
    size_t M = numPoints;
    if (M < 2 || N == 0 || sampleRate <= 0.0)
    {
        return {};
    }
    double startOmega = (2.0 * pi * startFrequency) / sampleRate;
    double endOmega = (2.0 * pi * endFrequency) / sampleRate;
    double deltaOmega = (endOmega - startOmega) / static_cast<double>(M - 1);
    std::complex<double> A(cosine(startOmega), sine(startOmega));
    std::complex<double> W(cosine(-deltaOmega), sine(-deltaOmega));
    size_t L = nextPowerOfTwo(N + M - 1);
    std::vector<std::complex<double>> x(L, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> y(L, std::complex<double>(0.0, 0.0));
    for (size_t n = 0; n < N; ++n)
    {
        double powerA = -static_cast<double>(n);
        double powerW = static_cast<double>(n * n) * 0.5;
        x[n] = std::complex<double>(static_cast<double>(input[n]), 0.0) * std::pow(A, powerA) * std::pow(W, powerW);
    }
    for (size_t k = 0; k < M; ++k)
    {
        y[k] = std::pow(W, -static_cast<double>(k * k) * 0.5);
    }
    for (size_t k = 1; k < N; ++k)
    {
        size_t index = L - k;
        y[index] = std::pow(W, -static_cast<double>(index * index) * 0.5);
    }
    fastFourierTransform1D(x, false);
    fastFourierTransform1D(y, false);
    for (size_t i = 0; i < L; ++i)
    {
        x[i] = x[i] * y[i];
    }
    fastFourierTransform1D(x, true);
    std::vector<std::complex<double>> spectrum(M);
    for (size_t k = 0; k < M; ++k)
    {
        spectrum[k] = x[k] * std::pow(W, static_cast<double>(k * k) * 0.5);
        if (abnormal(spectrum[k]))
        {
            throw std::overflow_error("chirpZTransform: result overflow");
        }
    }
    return spectrum;
}

template <typename InputContainer>
auto fractionalFourierTransform(const InputContainer &input, double alpha) -> std::vector<std::complex<double>>
{
    size_t length = input.size();
    if (length == 0)
    {
        return {};
    }
    std::vector<std::complex<double>> result(length);
    double factor = alpha * pi / 2.0;
    double baseAngle = (-2.0 * pi) / static_cast<double>(length);
    for (size_t k = 0; k < length; ++k)
    {
        std::complex<double> sum(0.0, 0.0);
        for (size_t t = 0; t < length; ++t)
        {
            double chirp = factor * static_cast<double>(t * t) / static_cast<double>(length);
            double angle = baseAngle * static_cast<double>(k) * static_cast<double>(t) + chirp;
            double cosValue = cosine(angle);
            double sinValue = sine(angle);
            double inputValue = static_cast<double>(input[t]);
            sum.real(sum.real() + inputValue * cosValue);
            sum.imag(sum.imag() + inputValue * sinValue);
        }
        if (abnormal(sum))
        {
            throw std::overflow_error("fractionalFourierTransform: result overflow");
        }
        result[k] = sum;
    }
    return result;
}

template <typename N>
N fastPower(N base, int exponent)
{
    if (abnormal(base))
    {
        throw std::domain_error("fastPower: NaN or Inf input");
    }
    if (exponent == 0)
    {
        return N(1);
    }
    if (exponent < 0)
    {
        return divide(N(1), fastPower(base, -exponent));
    }
    N result = N(1);
    N currentBase = base;
    for (int bit = 0; bit < 32; ++bit)
    {
        if (exponent & (1 << bit))
        {
            result = result * currentBase;
        }
        currentBase = currentBase * currentBase;
        if ((exponent >> (bit + 1)) == 0)
        {
            break;
        }
    }
    if (abnormal(result))
    {
        throw std::overflow_error("fastPower: result overflow");
    }
    return result;
}

template <typename Container>
auto fastConvolution(const Container &a, const Container &b) -> std::vector<typename Container::value_type>
{
    size_t sizeA = a.size();
    size_t sizeB = b.size();
    if (sizeA == 0 || sizeB == 0)
    {
        return {};
    }
    size_t resultSize = sizeA + sizeB - 1;
    size_t fftSize = nextPowerOfTwo(resultSize);
    std::vector<std::complex<double>> fa(fftSize, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> fb(fftSize, std::complex<double>(0.0, 0.0));
    for (size_t i = 0; i < sizeA; ++i)
    {
        fa[i] = std::complex<double>(static_cast<double>(a[i]), 0.0);
    }
    for (size_t i = 0; i < sizeB; ++i)
    {
        fb[i] = std::complex<double>(static_cast<double>(b[i]), 0.0);
    }
    fastFourierTransform1D(fa, false);
    fastFourierTransform1D(fb, false);
    for (size_t i = 0; i < fftSize; ++i)
    {
        fa[i] = fa[i] * fb[i];
    }
    fastFourierTransform1D(fa, true);
    std::vector<typename Container::value_type> result(resultSize);
    for (size_t i = 0; i < resultSize; ++i)
    {
        result[i] = static_cast<typename Container::value_type>(fa[i].real());
    }
    return result;
}

template <typename Container>
auto fastCrossCorrelation(const Container &a, const Container &b) -> std::vector<typename Container::value_type>
{
    size_t sizeA = a.size();
    size_t sizeB = b.size();
    if (sizeA == 0 || sizeB == 0)
    {
        return {};
    }
    size_t resultSize = sizeA + sizeB - 1;
    size_t fftSize = nextPowerOfTwo(resultSize);
    std::vector<std::complex<double>> fa(fftSize, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> fb(fftSize, std::complex<double>(0.0, 0.0));
    for (size_t i = 0; i < sizeA; ++i)
    {
        fa[i] = std::complex<double>(static_cast<double>(a[i]), 0.0);
    }
    for (size_t i = 0; i < sizeB; ++i)
    {
        fb[i] = std::complex<double>(static_cast<double>(b[sizeB - 1 - i]), 0.0);
    }
    fastFourierTransform1D(fa, false);
    fastFourierTransform1D(fb, false);
    for (size_t i = 0; i < fftSize; ++i)
    {
        fa[i] = fa[i] * fb[i];
    }
    fastFourierTransform1D(fa, true);
    std::vector<typename Container::value_type> result(resultSize);
    for (size_t i = 0; i < resultSize; ++i)
    {
        result[i] = static_cast<typename Container::value_type>(fa[i].real());
    }
    return result;
}

template <typename Container>
auto fastAutoCorrelation(const Container &a) -> std::vector<typename Container::value_type>
{
    return fastCrossCorrelation(a, a);
}

template <typename Container>
typename Container::value_type hornerEvaluate(const Container &coefficients, double x)
{
    size_t degree = coefficients.size();
    if (degree == 0)
    {
        return typename Container::value_type(0);
    }
    if (abnormal(x))
    {
        throw std::domain_error("hornerEvaluate: NaN or Inf input");
    }
    auto result = coefficients[degree - 1];
    for (size_t i = degree - 1; i > 0; --i)
    {
        result = result * x + coefficients[i - 1];
    }
    if (abnormal(result))
    {
        throw std::overflow_error("hornerEvaluate: result overflow");
    }
    return result;
}

template <typename Container>
auto polynomialMultiply(const Container &a, const Container &b) -> std::vector<typename Container::value_type>
{
    return fastConvolution(a, b);
}

template <typename Container>
typename Container::value_type dotProduct(const Container &a, const Container &b)
{
    if (a.size() != b.size())
    {
        throw std::domain_error("dotProduct: container size mismatch");
    }
    using ValueType = typename Container::value_type;
    ValueType sum(0);
    for (size_t i = 0; i < a.size(); ++i)
    {
        sum = sum + a[i] * b[i];
    }
    if (abnormal(sum))
    {
        throw std::overflow_error("dotProduct: result overflow");
    }
    return sum;
}

template <typename Container>
double euclideanNorm(const Container &v)
{
    double sum = 0.0;
    for (const auto &val : v)
    {
        double d = static_cast<double>(val);
        sum += d * d;
    }
    return squareRoot(sum);
}

template <typename Container>
auto convolution(const Container &a, const Container &b) -> std::vector<typename Container::value_type>
{
    return fastConvolution(a, b);
}

template <typename Container>
auto crossCorrelation(const Container &a, const Container &b) -> std::vector<typename Container::value_type>
{
    return fastCrossCorrelation(a, b);
}

template <typename Container>
auto autoCorrelation(const Container &a) -> std::vector<typename Container::value_type>
{
    return fastAutoCorrelation(a);
}

template <typename Container>
double shannonEntropy(const Container &probabilities)
{
    double entropyValue = 0.0;
    for (const auto &p : probabilities)
    {
        double prob = static_cast<double>(p);
        if (prob > tolerance)
        {
            entropyValue -= prob * std::log(prob);
        }
    }
    if (abnormal(entropyValue))
    {
        throw std::overflow_error("shannonEntropy: result overflow");
    }
    return entropyValue;
}

template <typename N>
N gaussian(N x, N mean, N standardDeviation)
{
    if (standardDeviation <= N(0))
    {
        throw std::domain_error("gaussian: standard deviation must be positive");
    }
    N diff = (x - mean) / standardDeviation;
    N exponent = N(-0.5) * diff * diff;
    N result = std::exp(exponent) / (standardDeviation * std::sqrt(N(2) * N(pi)));
    if (abnormal(result))
    {
        throw std::overflow_error("gaussian: result overflow");
    }
    return result;
}

template <typename N>
N uniformDistribution(N x, N a, N b)
{
    if (abnormal(x) || abnormal(a) || abnormal(b))
    {
        throw std::domain_error("uniformDistribution: NaN or Inf input");
    }
    if (a >= b)
    {
        throw std::domain_error("uniformDistribution: a must be less than b");
    }
    if (x < a || x > b)
    {
        return N(0);
    }
    return N(1) / (b - a);
}

template <typename N>
N exponentialDistribution(N x, N lambda)
{
    if (abnormal(x) || abnormal(lambda))
    {
        throw std::domain_error("exponentialDistribution: NaN or Inf input");
    }
    if (lambda <= N(0))
    {
        throw std::domain_error("exponentialDistribution: lambda must be positive");
    }
    if (x < N(0))
    {
        return N(0);
    }
    N result = lambda * std::exp(-lambda * x);
    if (abnormal(result))
    {
        throw std::overflow_error("exponentialDistribution: result overflow");
    }
    return result;
}

template <typename N>
N poissonDistribution(int k, N lambda)
{
    if (abnormal(lambda))
    {
        throw std::domain_error("poissonDistribution: NaN or Inf input");
    }
    if (lambda <= N(0) || k < 0)
    {
        return N(0);
    }
    N logResult = k * std::log(lambda) - lambda - std::lgamma(k + 1);
    N result = std::exp(logResult);
    if (abnormal(result))
    {
        throw std::overflow_error("poissonDistribution: result overflow");
    }
    return result;
}

template <typename N>
N binomialDistribution(int k, int n, N p)
{
    if (abnormal(p))
    {
        throw std::domain_error("binomialDistribution: NaN or Inf input");
    }
    if (k < 0 || k > n || p < N(0) || p > N(1))
    {
        return N(0);
    }
    N logResult = std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1) + k * std::log(p) + (n - k) * std::log(N(1) - p);
    N result = std::exp(logResult);
    if (abnormal(result))
    {
        throw std::overflow_error("binomialDistribution: result overflow");
    }
    return result;
}

inline std::vector<double> hanningWindow(size_t size)
{
    if (size == 0)
    {
        return {};
    }
    if (size == 1)
    {
        return {1.0};
    }
    std::vector<double> window(size);
    double base = (2.0 * pi) / static_cast<double>(size - 1);
    for (size_t i = 0; i < size; ++i)
    {
        window[i] = 0.5 * (1.0 - std::cos(base * static_cast<double>(i)));
    }
    return window;
}

inline std::vector<double> hammingWindow(size_t size)
{
    if (size == 0)
    {
        return {};
    }
    if (size == 1)
    {
        return {1.0};
    }
    std::vector<double> window(size);
    double base = (2.0 * pi) / static_cast<double>(size - 1);
    for (size_t i = 0; i < size; ++i)
    {
        window[i] = 0.54 - 0.46 * std::cos(base * static_cast<double>(i));
    }
    return window;
}

inline std::vector<double> blackmanWindow(size_t size)
{
    if (size == 0)
    {
        return {};
    }
    if (size == 1)
    {
        return {1.0};
    }
    std::vector<double> window(size);
    double base2 = (2.0 * pi) / static_cast<double>(size - 1);
    double base4 = (4.0 * pi) / static_cast<double>(size - 1);
    for (size_t i = 0; i < size; ++i)
    {
        window[i] = 0.42 - 0.5 * std::cos(base2 * static_cast<double>(i)) + 0.08 * std::cos(base4 * static_cast<double>(i));
    }
    return window;
}

template <typename Container>
auto prefixSum(const Container &data) -> std::vector<typename Container::value_type>
{
    size_t dataSize = data.size();
    if (dataSize == 0)
    {
        return {};
    }
    std::vector<typename Container::value_type> result(dataSize);
    result[0] = data[0];
    for (size_t i = 1; i < dataSize; ++i)
    {
        result[i] = result[i - 1] + data[i];
    }
    return result;
}

template <typename Container>
auto differenceArray(const Container &data) -> std::vector<typename Container::value_type>
{
    size_t dataSize = data.size();
    if (dataSize == 0)
    {
        return {};
    }
    std::vector<typename Container::value_type> result(dataSize);
    result[0] = data[0];
    for (size_t i = 1; i < dataSize; ++i)
    {
        result[i] = data[i] - data[i - 1];
    }
    return result;
}

template <typename N>
N greatestCommonDivisor(N a, N b)
{
    if (abnormal(a) || abnormal(b))
    {
        throw std::domain_error("greatestCommonDivisor: NaN or Inf input");
    }
    a = absolute(a);
    b = absolute(b);
    while (b != N(0))
    {
        N temp = b;
        b = modulus(a, b);
        a = temp;
    }
    return a;
}

template <typename N>
N extendedGreatestCommonDivisor(N a, N b, N &x, N &y)
{
    if (abnormal(a) || abnormal(b))
    {
        throw std::domain_error("extendedGreatestCommonDivisor: NaN or Inf input");
    }
    if (b == N(0))
    {
        x = N(1);
        y = N(0);
        return a;
    }
    N x1 = N(0);
    N y1 = N(0);
    N gcd = extendedGreatestCommonDivisor(b, modulus(a, b), x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return gcd;
}

template <typename N>
N modularExponentiation(N base, N exponent, N modulusValue)
{
    if (abnormal(base) || abnormal(exponent) || abnormal(modulusValue))
    {
        throw std::domain_error("modularExponentiation: NaN or Inf input");
    }
    if (modulusValue <= N(0))
    {
        throw std::domain_error("modularExponentiation: modulus must be positive");
    }
    N result = N(1);
    base = modulus(base, modulusValue);
    while (exponent > N(0))
    {
        if (modulus(exponent, N(2)) == N(1))
        {
            result = modulus(result * base, modulusValue);
        }
        base = modulus(base * base, modulusValue);
        exponent = exponent / N(2);
    }
    return result;
}

template <typename N>
std::vector<bool> sieveOfEratosthenes(N limit)
{
    if (limit < N(2))
    {
        return {};
    }
    int intLimit = static_cast<int>(limit);
    std::vector<bool> isPrime(intLimit + 1, true);
    isPrime[0] = false;
    isPrime[1] = false;
    for (int i = 2; i * i <= intLimit; ++i)
    {
        if (isPrime[i])
        {
            for (int j = i * i; j <= intLimit; j += i)
            {
                isPrime[j] = false;
            }
        }
    }
    return isPrime;
}

template <typename Container>
int binarySearch(const Container &data, const typename Container::value_type &target)
{
    int left = 0;
    int right = static_cast<int>(data.size()) - 1;
    while (left <= right)
    {
        int middle = left + (right - left) / 2;
        if (data[middle] == target)
        {
            return middle;
        }
        else if (data[middle] < target)
        {
            left = middle + 1;
        }
        else
        {
            right = middle - 1;
        }
    }
    return -1;
}

template <typename Function, typename N>
N ternarySearch(Function f, N left, N right, int iterations)
{
    if (abnormal(left) || abnormal(right))
    {
        throw std::domain_error("ternarySearch: NaN or Inf input");
    }
    if (left >= right)
    {
        throw std::domain_error("ternarySearch: left must be less than right");
    }
    for (int i = 0; i < iterations; ++i)
    {
        N middleLeft = left + (right - left) / N(3);
        N middleRight = right - (right - left) / N(3);
        N valueLeft = f(middleLeft);
        N valueRight = f(middleRight);
        if (abnormal(valueLeft) || abnormal(valueRight))
        {
            throw std::overflow_error("ternarySearch: function returned NaN or Inf");
        }
        if (valueLeft < valueRight)
        {
            right = middleRight;
        }
        else
        {
            left = middleLeft;
        }
    }
    return (left + right) / N(2);
}

template <typename Container>
typename Container::value_type quickSelect(Container &data, size_t k)
{
    size_t dataSize = data.size();
    if (k >= dataSize)
    {
        throw std::out_of_range("quickSelect: k out of range");
    }
    size_t left = 0;
    size_t right = dataSize - 1;
    while (left < right)
    {
        size_t pivotIndex = left + (right - left) / 2;
        auto pivotValue = data[pivotIndex];
        std::swap(data[pivotIndex], data[right]);
        size_t storeIndex = left;
        for (size_t i = left; i < right; ++i)
        {
            if (data[i] < pivotValue)
            {
                std::swap(data[i], data[storeIndex]);
                ++storeIndex;
            }
        }
        std::swap(data[storeIndex], data[right]);
        if (storeIndex == k)
        {
            return data[k];
        }
        else if (storeIndex < k)
        {
            left = storeIndex + 1;
        }
        else
        {
            right = storeIndex - 1;
        }
    }
    return data[left];
}

template <typename Container>
typename Container::value_type median(Container &data)
{
    if (data.empty())
    {
        throw std::domain_error("median: container is empty");
    }
    size_t middle = data.size() / 2;
    return quickSelect(data, middle);
}

template <typename Container>
size_t longestIncreasingSubsequence(const Container &data)
{
    size_t dataSize = data.size();
    if (dataSize == 0)
    {
        return 0;
    }
    std::vector<typename Container::value_type> tails;
    tails.reserve(dataSize);
    for (const auto &value : data)
    {
        auto it = std::lower_bound(tails.begin(), tails.end(), value);
        if (it == tails.end())
        {
            tails.push_back(value);
        }
        else
        {
            *it = value;
        }
    }
    return tails.size();
}

template <typename WeightContainer, typename ValueContainer>
auto knapsack(const WeightContainer &weights, const ValueContainer &values, int capacity)
    -> typename ValueContainer::value_type
{
    size_t itemCount = weights.size();
    if (itemCount != values.size() || capacity <= 0)
    {
        return typename ValueContainer::value_type(0);
    }
    std::vector<typename ValueContainer::value_type> dp(capacity + 1, typename ValueContainer::value_type(0));
    for (size_t i = 0; i < itemCount; ++i)
    {
        int weight = static_cast<int>(weights[i]);
        auto value = values[i];
        for (int w = capacity; w >= weight; --w)
        {
            auto candidate = dp[w - weight] + value;
            if (candidate > dp[w])
            {
                dp[w] = candidate;
            }
        }
    }
    return dp[capacity];
}

class DisjointSet
{
  private:
    std::vector<int> parent;
    std::vector<int> rank;

  public:
    DisjointSet(size_t size) : parent(size), rank(size, 0)
    {
        for (size_t i = 0; i < size; ++i)
        {
            parent[i] = static_cast<int>(i);
        }
    }

    int find(int x)
    {
        if (parent[x] != x)
        {
            parent[x] = find(parent[x]);
        }
        return parent[x];
    }

    void unite(int x, int y)
    {
        int rootX = find(x);
        int rootY = find(y);
        if (rootX == rootY)
        {
            return;
        }
        if (rank[rootX] < rank[rootY])
        {
            parent[rootX] = rootY;
        }
        else if (rank[rootX] > rank[rootY])
        {
            parent[rootY] = rootX;
        }
        else
        {
            parent[rootY] = rootX;
            ++rank[rootX];
        }
    }

    bool connected(int x, int y)
    {
        return find(x) == find(y);
    }
};

template <typename N>
std::vector<N> dijkstra(const std::vector<std::vector<std::pair<size_t, N>>> &graph, size_t start)
{
    size_t vertexCount = graph.size();
    if (vertexCount == 0 || start >= vertexCount)
    {
        return {};
    }
    std::vector<N> distances(vertexCount, std::numeric_limits<N>::max());
    distances[start] = N(0);
    using Pair = std::pair<N, size_t>;
    std::priority_queue<Pair, std::vector<Pair>, std::greater<Pair>> queue;
    queue.push({N(0), start});
    while (!queue.empty())
    {
        auto [currentDistance, currentVertex] = queue.top();
        queue.pop();
        if (currentDistance > distances[currentVertex])
        {
            continue;
        }
        for (const auto &[neighbor, weight] : graph[currentVertex])
        {
            N newDistance = currentDistance + weight;
            if (newDistance < distances[neighbor])
            {
                distances[neighbor] = newDistance;
                queue.push({newDistance, neighbor});
            }
        }
    }
    return distances;
}

} // namespace mathematics