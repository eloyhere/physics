#pragma once
#include <deque>
#include <vector>
#include <complex>
#include <algorithm>
#include <type_traits>
#include "mathematics.h"
#include "quaternion.h"
#include "quantization.h"

namespace physics
{
template <typename T>
struct is_arithmetic : std::false_type
{
};

template <typename T>
constexpr bool is_arithmetic_v = is_arithmetic<T>::value;

template <>
struct is_arithmetic<bool> : std::true_type
{
};
template <>
struct is_arithmetic<char> : std::true_type
{
};
template <>
struct is_arithmetic<signed char> : std::true_type
{
};
template <>
struct is_arithmetic<unsigned char> : std::true_type
{
};
template <>
struct is_arithmetic<short> : std::true_type
{
};
template <>
struct is_arithmetic<unsigned short> : std::true_type
{
};
template <>
struct is_arithmetic<int> : std::true_type
{
};
template <>
struct is_arithmetic<unsigned int> : std::true_type
{
};
template <>
struct is_arithmetic<long> : std::true_type
{
};
template <>
struct is_arithmetic<unsigned long> : std::true_type
{
};
template <>
struct is_arithmetic<long long> : std::true_type
{
};
template <>
struct is_arithmetic<unsigned long long> : std::true_type
{
};
template <>
struct is_arithmetic<float> : std::true_type
{
};
template <>
struct is_arithmetic<double> : std::true_type
{
};
template <>
struct is_arithmetic<long double> : std::true_type
{
};

template <typename N>
struct is_arithmetic<std::complex<N>> : is_arithmetic<N>
{
};

template <typename N>
struct is_arithmetic<::mathematics::Complex<N>> : is_arithmetic<N>
{
};

template <typename N>
struct is_arithmetic<::mathematics::Quaternion<N>> : is_arithmetic<N>
{
};

template <typename T>
struct is_scalar : std::false_type
{
};

template <typename T>
constexpr bool is_scalar_v = is_scalar<T>::value;

template <>
struct is_scalar<Time> : std::true_type
{
};
template <>
struct is_scalar<Period> : std::true_type
{
};
template <>
struct is_scalar<Resistance> : std::true_type
{
};
template <>
struct is_scalar<Conductance> : std::true_type
{
};
template <>
struct is_scalar<Capacitance> : std::true_type
{
};
template <>
struct is_scalar<Inductance> : std::true_type
{
};
template <>
struct is_scalar<Charge> : std::true_type
{
};
template <>
struct is_scalar<Flux> : std::true_type
{
};
template <>
struct is_scalar<Temperature> : std::true_type
{
};
template <>
struct is_scalar<HeatCapacity> : std::true_type
{
};
template <>
struct is_scalar<Noise> : std::true_type
{
};

template <typename T>
struct is_phasor : std::false_type
{
};

template <typename T>
constexpr bool is_phasor_v = is_phasor<T>::value;

template <typename N>
struct is_phasor<std::complex<N>> : std::true_type
{
};

template <typename N>
struct is_phasor<::mathematics::Complex<N>> : std::true_type
{
};

template <>
struct is_phasor<Phasor> : std::true_type
{
};
template <>
struct is_phasor<Voltage> : std::true_type
{
};
template <>
struct is_phasor<Current> : std::true_type
{
};
template <>
struct is_phasor<Impedance> : std::true_type
{
};
template <>
struct is_phasor<Admittance> : std::true_type
{
};
template <>
struct is_phasor<Power> : std::true_type
{
};
template <>
struct is_phasor<Energy> : std::true_type
{
};
template <>
struct is_phasor<Entropy> : std::true_type
{
};

template <typename T>
struct is_vector : std::false_type
{
};

template <typename T>
constexpr bool is_vector_v = is_vector<T>::value;

template <typename N>
struct is_vector<::mathematics::Quaternion<N>> : std::true_type
{
};

template <>
struct is_vector<Displacement> : std::true_type
{
};
template <>
struct is_vector<Velocity> : std::true_type
{
};
template <>
struct is_vector<Acceleration> : std::true_type
{
};
template <>
struct is_vector<Force> : std::true_type
{
};
template <>
struct is_vector<Momentum> : std::true_type
{
};
template <>
struct is_vector<Torque> : std::true_type
{
};

template <typename T>
constexpr bool is_physical_v = is_arithmetic_v<T> || is_scalar_v<T> || is_phasor_v<T> || is_vector_v<T>;

template <typename N>
class Sample
{
    static_assert(is_physical_v<N>, "Sample<N>: N must be arithmetic, scalar, phasor, or vector type");

  public:
    Frequency frequency;
    N value;

    Sample() : frequency(), value() {}

    Sample(const Frequency &frequencyValue, const N &sampleValue)
        : frequency(frequencyValue), value(sampleValue) {}

    Sample(double hertz, const N &sampleValue)
        : frequency(frequencyFromHertz(hertz)), value(sampleValue) {}

    bool operator<(const Sample &other) const { return frequency < other.frequency; }
    bool operator>(const Sample &other) const { return other < *this; }
    bool operator<=(const Sample &other) const { return !(other < *this); }
    bool operator>=(const Sample &other) const { return !(*this < other); }
    bool operator==(const Sample &other) const { return frequency == other.frequency && value == other.value; }
    bool operator!=(const Sample &other) const { return !(*this == other); }

    Sample operator+(const Sample &other) const
    {
        if (frequency != other.frequency)
        {
            throw std::domain_error("Sample::operator+: frequency mismatch");
        }
        return Sample(frequency, value + other.value);
    }

    Sample operator-(const Sample &other) const
    {
        if (frequency != other.frequency)
        {
            throw std::domain_error("Sample::operator-: frequency mismatch");
        }
        return Sample(frequency, value - other.value);
    }

    Sample operator*(double scalar) const { return Sample(frequency, value * scalar); }

    Sample operator/(double scalar) const
    {
        if (::mathematics::absolute(scalar) < ::mathematics::tolerance)
        {
            throw std::domain_error("Sample::operator/: division by zero");
        }
        return Sample(frequency, value / scalar);
    }

    Sample &operator+=(const Sample &other)
    {
        if (frequency != other.frequency)
        {
            throw std::domain_error("Sample::operator+=: frequency mismatch");
        }
        value = value + other.value;
        return *this;
    }

    Sample &operator-=(const Sample &other)
    {
        if (frequency != other.frequency)
        {
            throw std::domain_error("Sample::operator-=: frequency mismatch");
        }
        value = value - other.value;
        return *this;
    }

    Sample &operator*=(double scalar)
    {
        value = value * scalar;
        return *this;
    }

    Sample &operator/=(double scalar)
    {
        if (::mathematics::absolute(scalar) < ::mathematics::tolerance)
        {
            throw std::domain_error("Sample::operator/=: division by zero");
        }
        value = value / scalar;
        return *this;
    }

    friend Sample operator*(double scalar, const Sample &sample)
    {
        return Sample(sample.frequency, sample.value * scalar);
    }

    friend std::ostream &operator<<(std::ostream &output, const Sample &sample)
    {
        return output << sample.frequency << ": " << sample.value;
    }
};

template <typename T>
class Spectral
{
    static_assert(is_arithmetic_v<T>, "Spectral<T>: T must be double, Complex<N>, or Quaternion<N>");

  private:
    Frequency frequency;
    T phasor;
    double resolution;
    double confidence;

    double magnitude() const
    {
        if constexpr (std::is_arithmetic_v<T>)
        {
            return ::mathematics::absolute(static_cast<double>(phasor));
        }
        else
        {
            return phasor.module();
        }
    }

    ::mathematics::Complex<double> to() const
    {
        if constexpr (std::is_arithmetic_v<T>)
        {
            return ::mathematics::Complex<double>(static_cast<double>(phasor), 0.0);
        }
        else if constexpr (std::is_same_v<T, ::mathematics::Complex<double>>)
        {
            return phasor;
        }
        else if constexpr (std::is_same_v<T, ::mathematics::Quaternion<double>>)
        {
            return ::mathematics::Complex<double>(static_cast<double>(phasor.getReal()), static_cast<double>(phasor.getImagI()));
        }
        else if constexpr (is_phasor_v<T>)
        {
            return ::mathematics::Complex<double>(static_cast<double>(phasor.getReal()), static_cast<double>(phasor.getImag()));
        }
        else
        {
            return ::mathematics::Complex<double>(static_cast<double>(phasor.getReal()), static_cast<double>(phasor.getImagI()));
        }
    }

    T from(const ::mathematics::Complex<double> &complex) const
    {
        if constexpr (std::is_arithmetic_v<T>)
        {
            return static_cast<T>(complex.module());
        }
        else if constexpr (std::is_same_v<T, ::mathematics::Complex<double>>)
        {
            return complex;
        }
        else if constexpr (std::is_same_v<T, ::mathematics::Quaternion<double>>)
        {
            return ::mathematics::Quaternion<double>(0.0, complex.getReal(), complex.getImag(), 0.0);
        }
        else if constexpr (is_phasor_v<T>)
        {
            return T(static_cast<typename T::value_type>(complex.getReal()), static_cast<typename T::value_type>(complex.getImag()));
        }
        else
        {
            return T(static_cast<typename T::value_type>(0), static_cast<typename T::value_type>(complex.getReal()), static_cast<typename T::value_type>(complex.getImag()), static_cast<typename T::value_type>(0));
        }
    }

  public:
    Spectral() : frequency(), phasor(), resolution(0.0), confidence(1.0) {}

    Spectral(const Frequency &frequencyValue, const T &phasorValue, double resolutionValue = 0.0, double confidenceValue = 1.0)
        : frequency(frequencyValue), phasor(phasorValue), resolution(resolutionValue), confidence(confidenceValue) {}

    Spectral(double hertz, const T &phasorValue, double resolutionValue = 0.0, double confidenceValue = 1.0)
        : frequency(frequencyFromHertz(hertz)), phasor(phasorValue), resolution(resolutionValue), confidence(confidenceValue) {}

    const Frequency &getFrequency() const { return frequency; }
    double hertz() const { return frequency.hertz(); }
    double sigma() const { return frequency.sigma(); }
    double omega() const { return frequency.omega(); }
    const T &getPhasor() const { return phasor; }
    double amplitude() const { return magnitude(); }
    double phase() const
    {
        ::mathematics::Complex<double> complex = to();
        return ::mathematics::arctangent2(complex.getImag(), complex.getReal());
    }
    ::mathematics::Complex<double> complexCoefficient() const { return to(); }
    double getResolution() const { return resolution; }
    double getConfidence() const { return confidence; }
    bool valid() const { return amplitude() > ::mathematics::tolerance; }

    double phaseDifference(const Spectral &other) const
    {
        return ::mathematics::arctangent2(std::sin(phase() - other.phase()), std::cos(phase() - other.phase()));
    }

    double frequencyDifference(const Spectral &other) const { return frequency.hertz() - other.frequency.hertz(); }

    double amplitudeRatio(const Spectral &other) const
    {
        if (other.amplitude() > ::mathematics::tolerance)
        {
            return ::mathematics::divide(amplitude(), other.amplitude());
        }
        return std::numeric_limits<double>::infinity();
    }

    Spectral withFrequency(const Frequency &newFrequency) const { return Spectral(newFrequency, phasor, resolution, confidence); }
    Spectral withConfidence(double newConfidence) const { return Spectral(frequency, phasor, resolution, newConfidence); }

    bool operator<(const Spectral &other) const { return frequency < other.frequency; }
    bool operator>(const Spectral &other) const { return other < *this; }
    bool operator<=(const Spectral &other) const { return !(other < *this); }
    bool operator>=(const Spectral &other) const { return !(*this < other); }
    bool operator==(const Spectral &other) const
    {
        return frequency == other.frequency && phasor == other.phasor &&
               ::mathematics::absolute(resolution - other.resolution) < ::mathematics::tolerance &&
               ::mathematics::absolute(confidence - other.confidence) < ::mathematics::tolerance;
    }
    bool operator!=(const Spectral &other) const { return !(*this == other); }

    Spectral operator+(const Spectral &other) const
    {
        return Spectral(frequency, from(to() + other.to()),
                        ::mathematics::limit(resolution, 0.0, other.resolution),
                        ::mathematics::limit(confidence, 0.0, other.confidence));
    }

    Spectral operator-(const Spectral &other) const
    {
        return Spectral(frequency, from(to() - other.to()),
                        ::mathematics::limit(resolution, 0.0, other.resolution),
                        ::mathematics::limit(confidence, 0.0, other.confidence));
    }

    Spectral operator*(double scalar) const { return Spectral(frequency, from(to() * scalar), resolution, confidence); }

    Spectral operator/(double scalar) const
    {
        if (::mathematics::absolute(scalar) < ::mathematics::tolerance)
        {
            throw std::domain_error("Spectral division by zero scalar");
        }
        return Spectral(frequency, from(to() / scalar), resolution, confidence);
    }

    Spectral &operator+=(const Spectral &other)
    {
        *this = *this + other;
        return *this;
    }
    Spectral &operator-=(const Spectral &other)
    {
        *this = *this - other;
        return *this;
    }
    Spectral &operator*=(double scalar)
    {
        *this = *this * scalar;
        return *this;
    }
    Spectral &operator/=(double scalar)
    {
        *this = *this / scalar;
        return *this;
    }

    friend Spectral operator*(double scalar, const Spectral &spectral) { return spectral * scalar; }

    friend std::ostream &operator<<(std::ostream &output, const Spectral &spectral)
    {
        return output << spectral.hertz() << " Hz, " << spectral.amplitude()
                      << " amplitude, " << spectral.phase()
                      << " rad, resolution=" << spectral.resolution
                      << " Hz, confidence=" << spectral.confidence;
    }
};

template <typename N>
class Spectrum
{
    static_assert(is_physical_v<N>, "Spectrum<N>: N must be arithmetic, scalar, phasor, or vector type");

  private:
    std::deque<Sample<N>> lines;

    double magnitude(const N &value) const
    {
        if constexpr (std::is_arithmetic_v<N>)
        {
            return ::mathematics::absolute(static_cast<double>(value));
        }
        else if constexpr (is_phasor_v<N> || is_vector_v<N>)
        {
            return value.module();
        }
        else if constexpr (std::is_same_v<N, Time>)
        {
            return ::mathematics::absolute(value.seconds());
        }
        else if constexpr (std::is_same_v<N, Period>)
        {
            return ::mathematics::absolute(value.getDuration());
        }
        else if constexpr (std::is_same_v<N, Resistance>)
        {
            return value.ohms();
        }
        else if constexpr (std::is_same_v<N, Conductance>)
        {
            return value.siemens();
        }
        else if constexpr (std::is_same_v<N, Capacitance>)
        {
            return value.farads();
        }
        else if constexpr (std::is_same_v<N, Inductance>)
        {
            return value.henries();
        }
        else if constexpr (std::is_same_v<N, Charge>)
        {
            return ::mathematics::absolute(value.coulombs());
        }
        else if constexpr (std::is_same_v<N, Flux>)
        {
            return ::mathematics::absolute(value.webers());
        }
        else if constexpr (std::is_same_v<N, Temperature>)
        {
            return value.kelvin();
        }
        else if constexpr (std::is_same_v<N, HeatCapacity>)
        {
            return value.joulesPerKelvin();
        }
        else if constexpr (std::is_same_v<N, Noise>)
        {
            return value.density();
        }
        else
        {
            return 1.0;
        }
    }

  public:
    Spectrum() = default;

    void addLine(const Frequency &frequency, const N &value)
    {
        Sample<N> sample(frequency, value);
        auto iterator = std::lower_bound(lines.begin(), lines.end(), sample);
        if (iterator != lines.end() && iterator->frequency == frequency)
        {
            iterator->value = iterator->value + value;
        }
        else
        {
            lines.insert(iterator, sample);
        }
    }

    void addLine(double hertz, const N &value) { addLine(frequencyFromHertz(hertz), value); }

    void removeLine(const Frequency &frequency)
    {
        Sample<N> sample(frequency, N());
        auto iterator = std::lower_bound(lines.begin(), lines.end(), sample);
        if (iterator != lines.end() && iterator->frequency == frequency)
        {
            lines.erase(iterator);
        }
    }

    void removeLine(double hertz) { removeLine(frequencyFromHertz(hertz)); }

    N at(const Frequency &frequency) const
    {
        Sample<N> sample(frequency, N());
        auto iterator = std::lower_bound(lines.begin(), lines.end(), sample);
        if (iterator != lines.end() && iterator->frequency == frequency)
        {
            return iterator->value;
        }
        return N();
    }

    N at(double hertz) const { return at(frequencyFromHertz(hertz)); }

    size_t size() const { return lines.size(); }
    bool empty() const { return lines.empty(); }
    const std::deque<Sample<N>> &getLines() const { return lines; }

    double centroid() const
    {
        double weightedSum = 0.0;
        double totalWeight = 0.0;
        for (const auto &line : lines)
        {
            double weight = magnitude(line.value);
            weightedSum = ::mathematics::add(weightedSum, line.frequency.hertz() * weight);
            totalWeight = ::mathematics::add(totalWeight, weight);
        }
        if (totalWeight < ::mathematics::tolerance)
        {
            return 0.0;
        }
        return ::mathematics::divide(weightedSum, totalWeight);
    }

    double spread() const
    {
        double center = centroid();
        double totalWeight = 0.0;
        double variance = 0.0;
        for (const auto &line : lines)
        {
            double weight = magnitude(line.value);
            double difference = ::mathematics::subtract(line.frequency.hertz(), center);
            variance = ::mathematics::add(variance, difference * difference * weight);
            totalWeight = ::mathematics::add(totalWeight, weight);
        }
        if (totalWeight < ::mathematics::tolerance)
        {
            return 0.0;
        }
        return ::mathematics::squareRoot(::mathematics::divide(variance, totalWeight));
    }

    std::vector<double> peaks(double threshold = 0.0) const
    {
        if (lines.size() < 3)
        {
            return {};
        }
        std::vector<double> result;
        for (size_t index = 1; index < lines.size() - 1; ++index)
        {
            double previous = magnitude(lines[index - 1].value);
            double current = magnitude(lines[index].value);
            double next = magnitude(lines[index + 1].value);
            if (current > previous && current > next && current >= threshold)
            {
                result.push_back(lines[index].frequency.hertz());
            }
        }
        return result;
    }

    double correlation(const Spectrum &other) const
    {
        double numerator = 0.0;
        double denomFirst = 0.0;
        double denomSecond = 0.0;
        auto iterator1 = lines.begin();
        auto iterator2 = other.lines.begin();
        while (iterator1 != lines.end() && iterator2 != other.lines.end())
        {
            if (iterator1->frequency < iterator2->frequency)
            {
                double weight = magnitude(iterator1->value);
                denomFirst = ::mathematics::add(denomFirst, weight * weight);
                ++iterator1;
            }
            else if (iterator2->frequency < iterator1->frequency)
            {
                double weight = magnitude(iterator2->value);
                denomSecond = ::mathematics::add(denomSecond, weight * weight);
                ++iterator2;
            }
            else
            {
                double weight1 = magnitude(iterator1->value);
                double weight2 = magnitude(iterator2->value);
                numerator = ::mathematics::add(numerator, weight1 * weight2);
                denomFirst = ::mathematics::add(denomFirst, weight1 * weight1);
                denomSecond = ::mathematics::add(denomSecond, weight2 * weight2);
                ++iterator1;
                ++iterator2;
            }
        }
        while (iterator1 != lines.end())
        {
            double weight = magnitude(iterator1->value);
            denomFirst = ::mathematics::add(denomFirst, weight * weight);
            ++iterator1;
        }
        while (iterator2 != other.lines.end())
        {
            double weight = magnitude(iterator2->value);
            denomSecond = ::mathematics::add(denomSecond, weight * weight);
            ++iterator2;
        }
        if (denomFirst < ::mathematics::tolerance || denomSecond < ::mathematics::tolerance)
        {
            return 0.0;
        }
        return ::mathematics::limit(
            ::mathematics::divide(numerator, ::mathematics::squareRoot(denomFirst * denomSecond)),
            -1.0, 1.0);
    }

    Spectrum operator+(const Spectrum &other) const
    {
        Spectrum result;
        auto iterator1 = lines.begin();
        auto iterator2 = other.lines.begin();
        while (iterator1 != lines.end() && iterator2 != other.lines.end())
        {
            if (iterator1->frequency < iterator2->frequency)
            {
                result.lines.push_back(*iterator1);
                ++iterator1;
            }
            else if (iterator2->frequency < iterator1->frequency)
            {
                result.lines.push_back(*iterator2);
                ++iterator2;
            }
            else
            {
                result.lines.push_back(*iterator1 + *iterator2);
                ++iterator1;
                ++iterator2;
            }
        }
        while (iterator1 != lines.end())
        {
            result.lines.push_back(*iterator1);
            ++iterator1;
        }
        while (iterator2 != other.lines.end())
        {
            result.lines.push_back(*iterator2);
            ++iterator2;
        }
        return result;
    }

    Spectrum operator-(const Spectrum &other) const
    {
        Spectrum result;
        auto iterator1 = lines.begin();
        auto iterator2 = other.lines.begin();
        while (iterator1 != lines.end() && iterator2 != other.lines.end())
        {
            if (iterator1->frequency < iterator2->frequency)
            {
                result.lines.push_back(*iterator1);
                ++iterator1;
            }
            else if (iterator2->frequency < iterator1->frequency)
            {
                result.lines.push_back(*iterator2 * (-1.0));
                ++iterator2;
            }
            else
            {
                result.lines.push_back(*iterator1 - *iterator2);
                ++iterator1;
                ++iterator2;
            }
        }
        while (iterator1 != lines.end())
        {
            result.lines.push_back(*iterator1);
            ++iterator1;
        }
        while (iterator2 != other.lines.end())
        {
            result.lines.push_back(*iterator2 * (-1.0));
            ++iterator2;
        }
        return result;
    }

    Spectrum operator*(double scalar) const
    {
        Spectrum result;
        for (const auto &line : lines)
        {
            result.lines.push_back(line * scalar);
        }
        return result;
    }

    Spectrum &operator+=(const Spectrum &other)
    {
        *this = *this + other;
        return *this;
    }
    Spectrum &operator-=(const Spectrum &other)
    {
        *this = *this - other;
        return *this;
    }
    Spectrum &operator*=(double scalar)
    {
        for (auto &line : lines)
        {
            line *= scalar;
        }
        return *this;
    }

    friend Spectrum operator*(double scalar, const Spectrum &spectrum) { return spectrum * scalar; }

    friend std::ostream &operator<<(std::ostream &output, const Spectrum &spectrum)
    {
        for (const auto &line : spectrum.lines)
        {
            output << line << "\n";
        }
        return output;
    }
};

template <typename T>
class Domain
{
    static_assert(is_arithmetic_v<T>, "Domain<T>: T must be double, Complex<N>, or Quaternion<N>");

  private:
    std::deque<Spectral<T>> lines;

    Spectral<T> interpolate(const Spectral<T> &lower, const Spectral<T> &upper, double ratio) const
    {
        double clampedRatio = ::mathematics::limit(ratio, 0.0, 1.0);
        Frequency interpolatedFrequency(
            ::mathematics::add(lower.sigma(), ::mathematics::multiply(::mathematics::subtract(upper.sigma(), lower.sigma()), clampedRatio)),
            ::mathematics::add(lower.omega(), ::mathematics::multiply(::mathematics::subtract(upper.omega(), lower.omega()), clampedRatio)));
        T interpolatedPhasor = lower.getPhasor() + (upper.getPhasor() - lower.getPhasor()) * static_cast<T>(clampedRatio);
        double interpolatedResolution = ::mathematics::limit(lower.getResolution(), 0.0, upper.getResolution());
        double interpolatedConfidence = ::mathematics::limit(
            ::mathematics::add(lower.getConfidence(), ::mathematics::multiply(::mathematics::subtract(upper.getConfidence(), lower.getConfidence()), clampedRatio)),
            0.0, 1.0);
        return Spectral<T>(interpolatedFrequency, interpolatedPhasor, interpolatedResolution, interpolatedConfidence);
    }

    double magnitude(const Spectral<T> &value) const { return value.amplitude(); }

  public:
    Domain() = default;

    void addLine(const Spectral<T> &spectral)
    {
        auto iterator = std::lower_bound(lines.begin(), lines.end(), spectral);
        if (iterator != lines.end() && iterator->getFrequency() == spectral.getFrequency())
        {
            *iterator = spectral;
        }
        else
        {
            lines.insert(iterator, spectral);
        }
    }

    void addLine(const Frequency &frequency, const T &phasor, double resolution = 0.0, double confidence = 1.0)
    {
        addLine(Spectral<T>(frequency, phasor, resolution, confidence));
    }

    void addLine(double hertz, const T &phasor, double resolution = 0.0, double confidence = 1.0)
    {
        addLine(Spectral<T>(frequencyFromHertz(hertz), phasor, resolution, confidence));
    }

    void removeLine(const Frequency &frequency)
    {
        Spectral<T> spectral(frequency, T());
        auto iterator = std::lower_bound(lines.begin(), lines.end(), spectral);
        if (iterator != lines.end() && iterator->getFrequency() == frequency)
        {
            lines.erase(iterator);
        }
    }

    void removeLine(double hertz) { removeLine(frequencyFromHertz(hertz)); }

    Spectral<T> at(const Frequency &frequency) const
    {
        Spectral<T> spectral(frequency, T());
        if (lines.empty())
        {
            return Spectral<T>();
        }
        if (lines.size() == 1)
        {
            return lines.front();
        }
        auto upper = std::lower_bound(lines.begin(), lines.end(), spectral);
        if (upper == lines.begin())
        {
            return *upper;
        }
        if (upper == lines.end())
        {
            return lines.back();
        }
        auto lower = std::prev(upper);
        double ratio = ::mathematics::divide(frequency.hertz() - lower->hertz(), upper->hertz() - lower->hertz());
        return interpolate(*lower, *upper, ratio);
    }

    Spectral<T> at(double hertz) const { return at(frequencyFromHertz(hertz)); }

    size_t size() const { return lines.size(); }
    bool empty() const { return lines.empty(); }
    const std::deque<Spectral<T>> &getLines() const { return lines; }

    Domain bandpass(const Frequency &low, const Frequency &high) const
    {
        Domain result;
        for (const auto &line : lines)
        {
            if (line.getFrequency() >= low && line.getFrequency() <= high)
            {
                result.addLine(line);
            }
        }
        if (!lines.empty() && result.lines.empty())
        {
            result.addLine(at(low));
            result.addLine(at(high));
        }
        return result;
    }

    Domain bandpass(double lowHertz, double highHertz) const
    {
        return bandpass(frequencyFromHertz(lowHertz), frequencyFromHertz(highHertz));
    }

    double centroid() const
    {
        double weightedSum = 0.0;
        double totalWeight = 0.0;
        for (const auto &line : lines)
        {
            double weight = magnitude(line);
            weightedSum = ::mathematics::add(weightedSum, line.hertz() * weight);
            totalWeight = ::mathematics::add(totalWeight, weight);
        }
        if (totalWeight < ::mathematics::tolerance)
        {
            return 0.0;
        }
        return ::mathematics::divide(weightedSum, totalWeight);
    }

    double spread() const
    {
        double center = centroid();
        double totalWeight = 0.0;
        double variance = 0.0;
        for (const auto &line : lines)
        {
            double weight = magnitude(line);
            double difference = ::mathematics::subtract(line.hertz(), center);
            variance = ::mathematics::add(variance, difference * difference * weight);
            totalWeight = ::mathematics::add(totalWeight, weight);
        }
        if (totalWeight < ::mathematics::tolerance)
        {
            return 0.0;
        }
        return ::mathematics::squareRoot(::mathematics::divide(variance, totalWeight));
    }

    std::vector<double> peaks(double threshold = 0.0) const
    {
        if (lines.size() < 3)
        {
            return {};
        }
        std::vector<double> result;
        for (size_t index = 1; index < lines.size() - 1; ++index)
        {
            double previous = magnitude(lines[index - 1]);
            double current = magnitude(lines[index]);
            double next = magnitude(lines[index + 1]);
            if (current > previous && current > next && current >= threshold)
            {
                result.push_back(lines[index].hertz());
            }
        }
        return result;
    }

    Domain operator+(const Domain &other) const
    {
        Domain result;
        auto iterator1 = lines.begin();
        auto iterator2 = other.lines.begin();
        while (iterator1 != lines.end() && iterator2 != other.lines.end())
        {
            if (iterator1->getFrequency() < iterator2->getFrequency())
            {
                result.lines.push_back(*iterator1);
                ++iterator1;
            }
            else if (iterator2->getFrequency() < iterator1->getFrequency())
            {
                result.lines.push_back(*iterator2);
                ++iterator2;
            }
            else
            {
                result.lines.push_back(*iterator1 + *iterator2);
                ++iterator1;
                ++iterator2;
            }
        }
        while (iterator1 != lines.end())
        {
            result.lines.push_back(*iterator1);
            ++iterator1;
        }
        while (iterator2 != other.lines.end())
        {
            result.lines.push_back(*iterator2);
            ++iterator2;
        }
        return result;
    }

    Domain operator-(const Domain &other) const
    {
        Domain result;
        auto iterator1 = lines.begin();
        auto iterator2 = other.lines.begin();
        while (iterator1 != lines.end() && iterator2 != other.lines.end())
        {
            if (iterator1->getFrequency() < iterator2->getFrequency())
            {
                result.lines.push_back(*iterator1);
                ++iterator1;
            }
            else if (iterator2->getFrequency() < iterator1->getFrequency())
            {
                result.lines.push_back(*iterator2 * (-1.0));
                ++iterator2;
            }
            else
            {
                result.lines.push_back(*iterator1 - *iterator2);
                ++iterator1;
                ++iterator2;
            }
        }
        while (iterator1 != lines.end())
        {
            result.lines.push_back(*iterator1);
            ++iterator1;
        }
        while (iterator2 != other.lines.end())
        {
            result.lines.push_back(*iterator2 * (-1.0));
            ++iterator2;
        }
        return result;
    }

    Domain operator*(double scalar) const
    {
        Domain result;
        for (const auto &line : lines)
        {
            result.lines.push_back(line * scalar);
        }
        return result;
    }

    Domain &operator+=(const Domain &other)
    {
        *this = *this + other;
        return *this;
    }
    Domain &operator-=(const Domain &other)
    {
        *this = *this - other;
        return *this;
    }
    Domain &operator*=(double scalar)
    {
        for (auto &line : lines)
        {
            line *= scalar;
        }
        return *this;
    }

    friend Domain operator*(double scalar, const Domain &domain) { return domain * scalar; }

    friend std::ostream &operator<<(std::ostream &output, const Domain &domain)
    {
        for (const auto &line : domain.lines)
        {
            output << line << "\n";
        }
        return output;
    }
};

template <typename N>
class Wave
{
    static_assert(is_physical_v<N>, "Wave<N>: N must be arithmetic, scalar, phasor, or vector type");

  private:
    std::deque<Sample<N>> samples;
    double sampleRate;

    ::mathematics::Complex<double> to(const N &value) const
    {
        if constexpr (std::is_arithmetic_v<N>)
        {
            return ::mathematics::Complex<double>(static_cast<double>(value), 0.0);
        }
        else if constexpr (std::is_same_v<N, ::mathematics::Complex<double>>)
        {
            return value;
        }
        else if constexpr (std::is_same_v<N, ::mathematics::Quaternion<double>>)
        {
            return ::mathematics::Complex<double>(static_cast<double>(value.getReal()), static_cast<double>(value.getImagI()));
        }
        else if constexpr (is_phasor_v<N>)
        {
            return ::mathematics::Complex<double>(static_cast<double>(value.getReal()), static_cast<double>(value.getImag()));
        }
        else if constexpr (is_vector_v<N>)
        {
            return ::mathematics::Complex<double>(static_cast<double>(value.getReal()), static_cast<double>(value.getImagI()));
        }
        else if constexpr (std::is_same_v<N, Time>)
        {
            return ::mathematics::Complex<double>(value.seconds(), 0.0);
        }
        else if constexpr (std::is_same_v<N, Period>)
        {
            return ::mathematics::Complex<double>(value.getDuration(), 0.0);
        }
        else if constexpr (std::is_same_v<N, Resistance>)
        {
            return ::mathematics::Complex<double>(value.ohms(), 0.0);
        }
        else if constexpr (std::is_same_v<N, Conductance>)
        {
            return ::mathematics::Complex<double>(value.siemens(), 0.0);
        }
        else if constexpr (std::is_same_v<N, Capacitance>)
        {
            return ::mathematics::Complex<double>(value.farads(), 0.0);
        }
        else if constexpr (std::is_same_v<N, Inductance>)
        {
            return ::mathematics::Complex<double>(value.henries(), 0.0);
        }
        else if constexpr (std::is_same_v<N, Charge>)
        {
            return ::mathematics::Complex<double>(value.coulombs(), 0.0);
        }
        else if constexpr (std::is_same_v<N, Flux>)
        {
            return ::mathematics::Complex<double>(value.webers(), 0.0);
        }
        else if constexpr (std::is_same_v<N, Temperature>)
        {
            return ::mathematics::Complex<double>(value.kelvin(), 0.0);
        }
        else if constexpr (std::is_same_v<N, HeatCapacity>)
        {
            return ::mathematics::Complex<double>(value.joulesPerKelvin(), 0.0);
        }
        else if constexpr (std::is_same_v<N, Noise>)
        {
            return ::mathematics::Complex<double>(value.density(), 0.0);
        }
        else
        {
            return ::mathematics::Complex<double>(0.0, 0.0);
        }
    }

    N from(const ::mathematics::Complex<double> &complex) const
    {
        if constexpr (std::is_arithmetic_v<N>)
        {
            return static_cast<N>(complex.module());
        }
        else if constexpr (std::is_same_v<N, ::mathematics::Complex<double>>)
        {
            return complex;
        }
        else if constexpr (std::is_same_v<N, ::mathematics::Quaternion<double>>)
        {
            return ::mathematics::Quaternion<double>(0.0, complex.getReal(), complex.getImag(), 0.0);
        }
        else if constexpr (is_phasor_v<N>)
        {
            return N(static_cast<typename N::value_type>(complex.getReal()), static_cast<typename N::value_type>(complex.getImag()));
        }
        else if constexpr (is_vector_v<N>)
        {
            return N(static_cast<typename N::value_type>(0), static_cast<typename N::value_type>(complex.getReal()), static_cast<typename N::value_type>(complex.getImag()), static_cast<typename N::value_type>(0));
        }
        else if constexpr (std::is_same_v<N, Time>)
        {
            return Time(complex.getReal());
        }
        else if constexpr (std::is_same_v<N, Period>)
        {
            return Period(complex.getReal(), 1.0);
        }
        else if constexpr (std::is_same_v<N, Resistance>)
        {
            return Resistance(complex.module());
        }
        else if constexpr (std::is_same_v<N, Conductance>)
        {
            return Conductance(complex.module());
        }
        else if constexpr (std::is_same_v<N, Capacitance>)
        {
            return Capacitance(complex.module());
        }
        else if constexpr (std::is_same_v<N, Inductance>)
        {
            return Inductance(complex.module());
        }
        else if constexpr (std::is_same_v<N, Charge>)
        {
            return Charge(complex.module());
        }
        else if constexpr (std::is_same_v<N, Flux>)
        {
            return Flux(complex.module());
        }
        else if constexpr (std::is_same_v<N, Temperature>)
        {
            return Temperature(complex.module());
        }
        else if constexpr (std::is_same_v<N, HeatCapacity>)
        {
            return HeatCapacity(complex.module());
        }
        else if constexpr (std::is_same_v<N, Noise>)
        {
            return Noise(complex.module());
        }
        else
        {
            return N();
        }
    }

    void discreteCosineTransform(std::vector<double> &data) const
    {
        size_t size = data.size();
        if (size <= 1)
        {
            return;
        }
        std::vector<::mathematics::Complex<double>> complexData(size * 2, ::mathematics::Complex<double>(0.0, 0.0));
        for (size_t i = 0; i < size; ++i)
        {
            complexData[i] = ::mathematics::Complex<double>(data[i], 0.0);
            complexData[size * 2 - 1 - i] = ::mathematics::Complex<double>(data[i], 0.0);
        }
        ::mathematics::fastFourierTransform1D(complexData, false);
        for (size_t i = 0; i < size; ++i)
        {
            double angle = -::mathematics::pi * static_cast<double>(i) / static_cast<double>(size * 2);
            complexData[i] = complexData[i] * ::mathematics::Complex<double>(::mathematics::cosine(angle), ::mathematics::sine(angle)) * 2.0;
            data[i] = complexData[i].getReal();
        }
    }

  public:
    Wave(double inputSampleRate = 1.0) : sampleRate(inputSampleRate)
    {
        if (::mathematics::abnormal(sampleRate))
        {
            throw std::domain_error("Wave: NaN or Inf sample rate");
        }
        if (sampleRate <= 0.0)
        {
            throw std::domain_error("Wave: non-positive sample rate");
        }
    }

    Wave(std::deque<Sample<N>> inputSamples, double inputSampleRate)
        : samples(std::move(inputSamples)), sampleRate(inputSampleRate)
    {
        if (::mathematics::abnormal(sampleRate))
        {
            throw std::domain_error("Wave: NaN or Inf sample rate");
        }
        if (sampleRate <= 0.0)
        {
            throw std::domain_error("Wave: non-positive sample rate");
        }
    }

    void append(const N &value) { samples.push_back(Sample<N>(Frequency(), value)); }
    void append(const Sample<N> &sample) { samples.push_back(sample); }

    size_t size() const { return samples.size(); }
    bool empty() const { return samples.empty(); }
    double getSampleRate() const { return sampleRate; }
    double duration() const { return ::mathematics::divide(static_cast<double>(samples.size()), sampleRate); }
    double timeAt(size_t index) const { return ::mathematics::divide(static_cast<double>(index), sampleRate); }

    const N &operator[](size_t index) const { return samples[index].value; }
    N &operator[](size_t index) { return samples[index].value; }
    const std::deque<Sample<N>> &getSamples() const { return samples; }

    template <typename U = N>
    std::enable_if_t<!is_arithmetic_v<U>, Spectrum<U>> fastFourierTransform() const
    {
        Spectrum<U> spectrum;
        if (samples.empty())
        {
            return spectrum;
        }
        size_t padded = ::mathematics::nextPowerOfTwo(samples.size());
        std::vector<::mathematics::Complex<double>> data(padded, ::mathematics::Complex<double>(0.0, 0.0));
        double normalizationFactor = 1.0;
        if (samples.size() > 1)
        {
            double windowSum = 0.0;
            for (size_t i = 0; i < samples.size(); ++i)
            {
                double window = ::mathematics::smooth(static_cast<double>(i), 0.0, static_cast<double>(samples.size() - 1));
                window = 0.5 * (1.0 - ::mathematics::cosine(2.0 * ::mathematics::pi * window / static_cast<double>(samples.size() - 1)));
                data[i] = to(samples[i].value) * window;
                windowSum = ::mathematics::add(windowSum, window * window);
            }
            normalizationFactor = ::mathematics::squareRoot(::mathematics::divide(static_cast<double>(samples.size()), windowSum));
        }
        else
        {
            for (size_t i = 0; i < samples.size(); ++i)
            {
                data[i] = to(samples[i].value);
            }
        }
        ::mathematics::fastFourierTransform1D(data, false);
        double frequencyStep = ::mathematics::divide(sampleRate, static_cast<double>(padded));
        double normalization = ::mathematics::divide(normalizationFactor, static_cast<double>(padded));
        for (size_t i = 0; i <= padded / 2; ++i)
        {
            double hertz = static_cast<double>(i) * frequencyStep;
            spectrum.addLine(frequencyFromHertz(hertz), from(data[i] * normalization));
        }
        return spectrum;
    }

    template <typename U = N>
    std::enable_if_t<is_arithmetic_v<U>, Domain<::mathematics::Complex<double>>> fastFourierTransform() const
    {
        Domain<::mathematics::Complex<double>> domain;
        if (samples.empty())
        {
            return domain;
        }
        size_t padded = ::mathematics::nextPowerOfTwo(samples.size());
        std::vector<::mathematics::Complex<double>> data(padded, ::mathematics::Complex<double>(0.0, 0.0));
        double normalizationFactor = 1.0;
        if (samples.size() > 1)
        {
            double windowSum = 0.0;
            for (size_t i = 0; i < samples.size(); ++i)
            {
                double window = ::mathematics::smooth(static_cast<double>(i), 0.0, static_cast<double>(samples.size() - 1));
                window = 0.5 * (1.0 - ::mathematics::cosine(2.0 * ::mathematics::pi * window / static_cast<double>(samples.size() - 1)));
                data[i] = to(samples[i].value) * window;
                windowSum = ::mathematics::add(windowSum, window * window);
            }
            normalizationFactor = ::mathematics::squareRoot(::mathematics::divide(static_cast<double>(samples.size()), windowSum));
        }
        else
        {
            for (size_t i = 0; i < samples.size(); ++i)
            {
                data[i] = to(samples[i].value);
            }
        }
        ::mathematics::fastFourierTransform1D(data, false);
        double frequencyStep = ::mathematics::divide(sampleRate, static_cast<double>(padded));
        double normalization = ::mathematics::divide(normalizationFactor, static_cast<double>(padded));
        for (size_t i = 0; i <= padded / 2; ++i)
        {
            ::mathematics::Complex<double> complex = data[i] * normalization;
            double hertz = static_cast<double>(i) * frequencyStep;
            domain.addLine(frequencyFromHertz(hertz), complex, frequencyStep, 1.0);
        }
        return domain;
    }

    template <typename U = N>
    std::enable_if_t<!is_arithmetic_v<U>, Spectrum<U>> discreteFourierTransform() const
    {
        Spectrum<U> spectrum;
        if (samples.empty())
        {
            return spectrum;
        }
        size_t length = samples.size();
        double frequencyStep = ::mathematics::divide(sampleRate, static_cast<double>(length));
        for (size_t k = 0; k <= length / 2; ++k)
        {
            double hertz = static_cast<double>(k) * frequencyStep;
            ::mathematics::Complex<double> sum(0.0, 0.0);
            for (size_t t = 0; t < length; ++t)
            {
                double angle = -2.0 * ::mathematics::pi * static_cast<double>(k * t) / static_cast<double>(length);
                ::mathematics::Complex<double> twiddle(::mathematics::cosine(angle), ::mathematics::sine(angle));
                sum = sum + to(samples[t].value) * twiddle;
            }
            spectrum.addLine(frequencyFromHertz(hertz), from(sum / static_cast<double>(length)));
        }
        return spectrum;
    }

    template <typename U = N>
    std::enable_if_t<is_arithmetic_v<U>, Domain<::mathematics::Complex<double>>> discreteFourierTransform() const
    {
        Domain<::mathematics::Complex<double>> domain;
        if (samples.empty())
        {
            return domain;
        }
        size_t length = samples.size();
        double frequencyStep = ::mathematics::divide(sampleRate, static_cast<double>(length));
        for (size_t k = 0; k <= length / 2; ++k)
        {
            double hertz = static_cast<double>(k) * frequencyStep;
            ::mathematics::Complex<double> sum(0.0, 0.0);
            for (size_t t = 0; t < length; ++t)
            {
                double angle = -2.0 * ::mathematics::pi * static_cast<double>(k * t) / static_cast<double>(length);
                ::mathematics::Complex<double> twiddle(::mathematics::cosine(angle), ::mathematics::sine(angle));
                sum = sum + to(samples[t].value) * twiddle;
            }
            ::mathematics::Complex<double> value = sum / static_cast<double>(length);
            domain.addLine(frequencyFromHertz(hertz), value, frequencyStep, 1.0);
        }
        return domain;
    }

    template <typename U = N>
    std::enable_if_t<!is_arithmetic_v<U>, Wave<U>> inverseFastFourierTransform(const Spectrum<U> &spectrum) const
    {
        size_t spectrumSize = spectrum.size();
        if (spectrumSize == 0)
        {
            return Wave<U>(sampleRate);
        }
        size_t outputSize = (spectrumSize - 1) * 2;
        std::vector<::mathematics::Complex<double>> data(outputSize, ::mathematics::Complex<double>(0.0, 0.0));
        const auto &lines = spectrum.getLines();
        for (size_t i = 0; i < spectrumSize; ++i)
        {
            data[i] = to(lines[i].value);
            if (i > 0 && i < spectrumSize - 1)
            {
                data[outputSize - i] = data[i].conjugate();
            }
        }
        ::mathematics::fastFourierTransform1D(data, true);
        std::deque<Sample<U>> resultSamples;
        size_t keep = ::mathematics::limit(static_cast<double>(samples.size()), 0.0, static_cast<double>(outputSize));
        for (size_t i = 0; i < keep; ++i)
        {
            resultSamples.push_back(Sample<U>(Frequency(), from(data[i])));
        }
        return Wave<U>(resultSamples, sampleRate);
    }

    template <typename U = N>
    std::enable_if_t<is_arithmetic_v<U>, Wave<U>> inverseFastFourierTransform(const Domain<::mathematics::Complex<double>> &domain) const
    {
        size_t domainSize = domain.size();
        if (domainSize == 0)
        {
            return Wave<U>(sampleRate);
        }
        size_t outputSize = (domainSize - 1) * 2;
        std::vector<::mathematics::Complex<double>> data(outputSize, ::mathematics::Complex<double>(0.0, 0.0));
        const auto &lines = domain.getLines();
        for (size_t i = 0; i < domainSize; ++i)
        {
            data[i] = lines[i].complexCoefficient();
            if (i > 0 && i < domainSize - 1)
            {
                data[outputSize - i] = data[i].conjugate();
            }
        }
        ::mathematics::fastFourierTransform1D(data, true);
        std::deque<Sample<U>> resultSamples;
        size_t keep = ::mathematics::limit(static_cast<double>(samples.size()), 0.0, static_cast<double>(outputSize));
        for (size_t i = 0; i < keep; ++i)
        {
            resultSamples.push_back(Sample<U>(Frequency(), from(data[i])));
        }
        return Wave<U>(resultSamples, sampleRate);
    }

    template <typename U = N>
    std::enable_if_t<!is_arithmetic_v<U>, std::vector<Spectrum<U>>> shortTimeFourierTransform(size_t windowSize, size_t hopSize) const
    {
        std::vector<Spectrum<U>> spectrogram;
        if (windowSize == 0 || hopSize == 0 || samples.empty())
        {
            return spectrogram;
        }
        for (size_t start = 0; start + windowSize <= samples.size(); start += hopSize)
        {
            Wave frame(sampleRate);
            for (size_t i = 0; i < windowSize; ++i)
            {
                frame.append(samples[start + i].value);
            }
            spectrogram.push_back(frame.fastFourierTransform());
        }
        return spectrogram;
    }

    template <typename U = N>
    std::enable_if_t<is_arithmetic_v<U>, std::vector<Domain<::mathematics::Complex<double>>>> shortTimeFourierTransform(size_t windowSize, size_t hopSize) const
    {
        std::vector<Domain<::mathematics::Complex<double>>> spectrogram;
        if (windowSize == 0 || hopSize == 0 || samples.empty())
        {
            return spectrogram;
        }
        for (size_t start = 0; start + windowSize <= samples.size(); start += hopSize)
        {
            Wave frame(sampleRate);
            for (size_t i = 0; i < windowSize; ++i)
            {
                frame.append(samples[start + i].value);
            }
            spectrogram.push_back(frame.fastFourierTransform());
        }
        return spectrogram;
    }

    template <typename U = N>
    std::enable_if_t<!is_arithmetic_v<U>, Spectrum<U>> chirpZTransform(double startFrequency, double endFrequency, size_t numPoints) const
    {
        Spectrum<U> spectrum;
        if (numPoints < 2 || samples.empty())
        {
            return spectrum;
        }
        double startOmega = ::mathematics::divide(2.0 * ::mathematics::pi * startFrequency, sampleRate);
        double endOmega = ::mathematics::divide(2.0 * ::mathematics::pi * endFrequency, sampleRate);
        double deltaOmega = ::mathematics::divide(endOmega - startOmega, static_cast<double>(numPoints - 1));
        ::mathematics::Complex<double> a(::mathematics::cosine(startOmega), ::mathematics::sine(startOmega));
        ::mathematics::Complex<double> w(::mathematics::cosine(-deltaOmega), ::mathematics::sine(-deltaOmega));
        size_t n = samples.size();
        size_t m = numPoints;
        size_t l = ::mathematics::nextPowerOfTwo(n + m - 1);
        std::vector<::mathematics::Complex<double>> x(l, ::mathematics::Complex<double>(0.0, 0.0));
        std::vector<::mathematics::Complex<double>> y(l, ::mathematics::Complex<double>(0.0, 0.0));
        for (size_t i = 0; i < n; ++i)
        {
            x[i] = to(samples[i].value) * a.power(-static_cast<double>(i)) * w.power(static_cast<double>(i * i) * 0.5);
        }
        for (size_t k = 0; k < m; ++k)
        {
            y[k] = w.power(-static_cast<double>(k * k) * 0.5);
        }
        for (size_t k = 1; k < n; ++k)
        {
            y[l - k] = w.power(-static_cast<double>((l - k) * (l - k)) * 0.5);
        }
        ::mathematics::fastFourierTransform1D(x, false);
        ::mathematics::fastFourierTransform1D(y, false);
        for (size_t i = 0; i < l; ++i)
        {
            x[i] = x[i] * y[i];
        }
        ::mathematics::fastFourierTransform1D(x, true);
        for (size_t k = 0; k < m; ++k)
        {
            double hertz = startFrequency + static_cast<double>(k) * (endFrequency - startFrequency) / static_cast<double>(m - 1);
            spectrum.addLine(frequencyFromHertz(hertz), from(x[k] * w.power(static_cast<double>(k * k) * 0.5)));
        }
        return spectrum;
    }

    template <typename U = N>
    std::enable_if_t<is_arithmetic_v<U>, Domain<::mathematics::Complex<double>>> chirpZTransform(double startFrequency, double endFrequency, size_t numPoints) const
    {
        Domain<::mathematics::Complex<double>> domain;
        if (numPoints < 2 || samples.empty())
        {
            return domain;
        }
        double startOmega = ::mathematics::divide(2.0 * ::mathematics::pi * startFrequency, sampleRate);
        double endOmega = ::mathematics::divide(2.0 * ::mathematics::pi * endFrequency, sampleRate);
        double deltaOmega = ::mathematics::divide(endOmega - startOmega, static_cast<double>(numPoints - 1));
        ::mathematics::Complex<double> a(::mathematics::cosine(startOmega), ::mathematics::sine(startOmega));
        ::mathematics::Complex<double> w(::mathematics::cosine(-deltaOmega), ::mathematics::sine(-deltaOmega));
        size_t n = samples.size();
        size_t m = numPoints;
        size_t l = ::mathematics::nextPowerOfTwo(n + m - 1);
        std::vector<::mathematics::Complex<double>> x(l, ::mathematics::Complex<double>(0.0, 0.0));
        std::vector<::mathematics::Complex<double>> y(l, ::mathematics::Complex<double>(0.0, 0.0));
        for (size_t i = 0; i < n; ++i)
        {
            x[i] = to(samples[i].value) * a.power(-static_cast<double>(i)) * w.power(static_cast<double>(i * i) * 0.5);
        }
        for (size_t k = 0; k < m; ++k)
        {
            y[k] = w.power(-static_cast<double>(k * k) * 0.5);
        }
        for (size_t k = 1; k < n; ++k)
        {
            y[l - k] = w.power(-static_cast<double>((l - k) * (l - k)) * 0.5);
        }
        ::mathematics::fastFourierTransform1D(x, false);
        ::mathematics::fastFourierTransform1D(y, false);
        for (size_t i = 0; i < l; ++i)
        {
            x[i] = x[i] * y[i];
        }
        ::mathematics::fastFourierTransform1D(x, true);
        double frequencyStep = ::mathematics::divide(endFrequency - startFrequency, static_cast<double>(m - 1));
        for (size_t k = 0; k < m; ++k)
        {
            double hertz = startFrequency + static_cast<double>(k) * frequencyStep;
            ::mathematics::Complex<double> value = x[k] * w.power(static_cast<double>(k * k) * 0.5);
            domain.addLine(frequencyFromHertz(hertz), value, frequencyStep, 1.0);
        }
        return domain;
    }

    Domain<::mathematics::Complex<double>> powerSpectralDensity() const
    {
        Domain<::mathematics::Complex<double>> domain;
        if (samples.empty())
        {
            return domain;
        }
        auto transform = fastFourierTransform();
        double frequencyStep = 0.0;
        if constexpr (!is_arithmetic_v<N>)
        {
            const auto &lines = transform.template getLines();
            if (lines.size() >= 2)
            {
                frequencyStep = lines[1].frequency.hertz() - lines[0].frequency.hertz();
            }
            for (const auto &line : lines)
            {
                ::mathematics::Complex<double> complex = to(line.value);
                double power = ::mathematics::divide(complex.getReal() * complex.getReal() + complex.getImag() * complex.getImag(), static_cast<double>(samples.size()));
                domain.addLine(line.frequency, ::mathematics::Complex<double>(power, 0.0), frequencyStep, 1.0);
            }
        }
        else
        {
            const auto &lines = transform.getLines();
            if (lines.size() >= 2)
            {
                frequencyStep = lines[1].getFrequency().hertz() - lines[0].getFrequency().hertz();
            }
            for (const auto &line : lines)
            {
                ::mathematics::Complex<double> complex = line.complexCoefficient();
                double power = ::mathematics::divide(complex.getReal() * complex.getReal() + complex.getImag() * complex.getImag(), static_cast<double>(samples.size()));
                domain.addLine(line.getFrequency(), ::mathematics::Complex<double>(power, 0.0), frequencyStep, 1.0);
            }
        }
        return domain;
    }

    Domain<::mathematics::Complex<double>> amplitudeSpectrum() const
    {
        Domain<::mathematics::Complex<double>> domain;
        if (samples.empty())
        {
            return domain;
        }
        auto transform = fastFourierTransform();
        double frequencyStep = 0.0;
        if constexpr (!is_arithmetic_v<N>)
        {
            const auto &lines = transform.template getLines();
            if (lines.size() >= 2)
            {
                frequencyStep = lines[1].frequency.hertz() - lines[0].frequency.hertz();
            }
            for (const auto &line : lines)
            {
                domain.addLine(line.frequency, ::mathematics::Complex<double>(to(line.value).module(), 0.0), frequencyStep, 1.0);
            }
        }
        else
        {
            const auto &lines = transform.getLines();
            if (lines.size() >= 2)
            {
                frequencyStep = lines[1].getFrequency().hertz() - lines[0].getFrequency().hertz();
            }
            for (const auto &line : lines)
            {
                domain.addLine(line.getFrequency(), ::mathematics::Complex<double>(line.amplitude(), 0.0), frequencyStep, 1.0);
            }
        }
        return domain;
    }

    Domain<::mathematics::Complex<double>> phaseSpectrum() const
    {
        Domain<::mathematics::Complex<double>> domain;
        if (samples.empty())
        {
            return domain;
        }
        auto transform = fastFourierTransform();
        double frequencyStep = 0.0;
        if constexpr (!is_arithmetic_v<N>)
        {
            const auto &lines = transform.template getLines();
            if (lines.size() >= 2)
            {
                frequencyStep = lines[1].frequency.hertz() - lines[0].frequency.hertz();
            }
            for (const auto &line : lines)
            {
                ::mathematics::Complex<double> complex = to(line.value);
                domain.addLine(line.frequency, ::mathematics::Complex<double>(::mathematics::arctangent2(complex.getImag(), complex.getReal()), 0.0), frequencyStep, 1.0);
            }
        }
        else
        {
            const auto &lines = transform.getLines();
            if (lines.size() >= 2)
            {
                frequencyStep = lines[1].getFrequency().hertz() - lines[0].getFrequency().hertz();
            }
            for (const auto &line : lines)
            {
                domain.addLine(line.getFrequency(), ::mathematics::Complex<double>(line.phase(), 0.0), frequencyStep, 1.0);
            }
        }
        return domain;
    }

    Domain<::mathematics::Complex<double>> autoPowerSpectrum() const { return powerSpectralDensity(); }

    Domain<::mathematics::Complex<double>> crossPowerSpectrum(const Wave &other) const
    {
        Domain<::mathematics::Complex<double>> domain;
        if (samples.empty() || other.samples.empty())
        {
            return domain;
        }
        if (::mathematics::absolute(sampleRate - other.sampleRate) > ::mathematics::tolerance)
        {
            throw std::runtime_error("crossPowerSpectrum: rate mismatch");
        }
        auto transformA = fastFourierTransform();
        auto transformB = other.fastFourierTransform();
        double frequencyStep = 0.0;
        if constexpr (!is_arithmetic_v<N>)
        {
            const auto &linesA = transformA.template getLines();
            const auto &linesB = transformB.template getLines();
            size_t minSize = std::min(linesA.size(), linesB.size());
            if (minSize >= 2)
            {
                frequencyStep = linesA[1].frequency.hertz() - linesA[0].frequency.hertz();
            }
            for (size_t i = 0; i < minSize; ++i)
            {
                ::mathematics::Complex<double> complexA = to(linesA[i].value);
                ::mathematics::Complex<double> complexB = to(linesB[i].value);
                ::mathematics::Complex<double> cross = complexA * complexB.conjugate();
                double power = ::mathematics::divide(cross.module(), static_cast<double>(std::min(samples.size(), other.samples.size())));
                domain.addLine(linesA[i].frequency, ::mathematics::Complex<double>(power, 0.0), frequencyStep, 1.0);
            }
        }
        else
        {
            const auto &linesA = transformA.getLines();
            const auto &linesB = transformB.getLines();
            size_t minSize = std::min(linesA.size(), linesB.size());
            if (minSize >= 2)
            {
                frequencyStep = linesA[1].getFrequency().hertz() - linesA[0].getFrequency().hertz();
            }
            for (size_t i = 0; i < minSize; ++i)
            {
                ::mathematics::Complex<double> complexA = linesA[i].complexCoefficient();
                ::mathematics::Complex<double> complexB = linesB[i].complexCoefficient();
                ::mathematics::Complex<double> cross = complexA * complexB.conjugate();
                double power = ::mathematics::divide(cross.module(), static_cast<double>(std::min(samples.size(), other.samples.size())));
                domain.addLine(linesA[i].getFrequency(), ::mathematics::Complex<double>(power, 0.0), frequencyStep, 1.0);
            }
        }
        return domain;
    }

    Domain<::mathematics::Complex<double>> cepstrum() const
    {
        Domain<::mathematics::Complex<double>> domain;
        if (samples.empty())
        {
            return domain;
        }
        Domain<::mathematics::Complex<double>> psd = powerSpectralDensity();
        std::vector<double> logSpectrum;
        const auto &lines = psd.getLines();
        for (const auto &line : lines)
        {
            logSpectrum.push_back(::mathematics::logarithm(line.amplitude() + ::mathematics::tolerance));
        }
        Wave<double> logWave(sampleRate);
        for (double value : logSpectrum)
        {
            logWave.append(value);
        }
        domain = logWave.fastFourierTransform();
        return domain;
    }

    Wave hilbertTransform() const
    {
        Wave result(sampleRate);
        if (samples.empty())
        {
            return result;
        }
        size_t padded = ::mathematics::nextPowerOfTwo(samples.size());
        std::vector<::mathematics::Complex<double>> data(padded, ::mathematics::Complex<double>(0.0, 0.0));
        for (size_t i = 0; i < samples.size(); ++i)
        {
            data[i] = to(samples[i].value);
        }
        ::mathematics::fastFourierTransform1D(data, false);
        for (size_t i = 0; i < padded; ++i)
        {
            if (i == 0 || i == padded / 2)
            {
                data[i] = ::mathematics::Complex<double>(0.0, 0.0);
            }
            else if (i < padded / 2)
            {
                data[i] = data[i] * ::mathematics::Complex<double>(0.0, 1.0);
            }
            else
            {
                data[i] = data[i] * ::mathematics::Complex<double>(0.0, -1.0);
            }
        }
        ::mathematics::fastFourierTransform1D(data, true);
        for (size_t i = 0; i < samples.size(); ++i)
        {
            result.append(from(data[i]));
        }
        return result;
    }

    Domain<::mathematics::Complex<double>> discreteCosineTransform() const
    {
        Domain<::mathematics::Complex<double>> domain;
        if (samples.empty())
        {
            return domain;
        }
        size_t padded = ::mathematics::nextPowerOfTwo(samples.size());
        std::vector<double> data(padded, 0.0);
        for (size_t i = 0; i < samples.size(); ++i)
        {
            data[i] = to(samples[i].value).module();
        }
        discreteCosineTransform(data);
        double frequencyStep = ::mathematics::divide(sampleRate, static_cast<double>(padded * 2));
        for (size_t i = 0; i <= padded / 2; ++i)
        {
            double hertz = static_cast<double>(i) * frequencyStep;
            domain.addLine(frequencyFromHertz(hertz), ::mathematics::Complex<double>(data[i], 0.0), frequencyStep, 1.0);
        }
        return domain;
    }

    Domain<::mathematics::Complex<double>> discreteSineTransform() const
    {
        Domain<::mathematics::Complex<double>> domain;
        if (samples.empty())
        {
            return domain;
        }
        size_t padded = ::mathematics::nextPowerOfTwo(samples.size());
        std::vector<double> data(padded, 0.0);
        for (size_t i = 0; i < samples.size(); ++i)
        {
            data[i] = to(samples[i].value).module() * ::mathematics::sine(::mathematics::pi * static_cast<double>(i + 1) / static_cast<double>(padded * 2));
        }
        discreteCosineTransform(data);
        double frequencyStep = ::mathematics::divide(sampleRate, static_cast<double>(padded * 2));
        for (size_t i = 0; i <= padded / 2; ++i)
        {
            double hertz = static_cast<double>(i) * frequencyStep;
            domain.addLine(frequencyFromHertz(hertz), ::mathematics::Complex<double>(data[i], 0.0), frequencyStep, 1.0);
        }
        return domain;
    }

    Wave operator+(const Wave &other) const
    {
        if (::mathematics::absolute(sampleRate - other.sampleRate) > ::mathematics::tolerance)
        {
            throw std::runtime_error("Wave+: rate mismatch");
        }
        size_t minSize = std::min(samples.size(), other.samples.size());
        std::deque<Sample<N>> result;
        for (size_t i = 0; i < minSize; ++i)
        {
            result.push_back(Sample<N>(Frequency(), samples[i].value + other.samples[i].value));
        }
        return Wave(result, sampleRate);
    }

    Wave operator-(const Wave &other) const
    {
        if (::mathematics::absolute(sampleRate - other.sampleRate) > ::mathematics::tolerance)
        {
            throw std::runtime_error("Wave-: rate mismatch");
        }
        size_t minSize = std::min(samples.size(), other.samples.size());
        std::deque<Sample<N>> result;
        for (size_t i = 0; i < minSize; ++i)
        {
            result.push_back(Sample<N>(Frequency(), samples[i].value - other.samples[i].value));
        }
        return Wave(result, sampleRate);
    }

    Wave operator*(double scalar) const
    {
        std::deque<Sample<N>> result;
        for (const auto &sample : samples)
        {
            result.push_back(Sample<N>(Frequency(), sample.value * scalar));
        }
        return Wave(result, sampleRate);
    }

    Wave operator/(double scalar) const
    {
        std::deque<Sample<N>> result;
        for (const auto &sample : samples)
        {
            result.push_back(Sample<N>(Frequency(), sample.value / scalar));
        }
        return Wave(result, sampleRate);
    }

    Wave &operator+=(const Wave &other)
    {
        if (::mathematics::absolute(sampleRate - other.sampleRate) > ::mathematics::tolerance)
        {
            throw std::runtime_error("Wave+=: rate mismatch");
        }
        size_t minSize = std::min(samples.size(), other.samples.size());
        for (size_t i = 0; i < minSize; ++i)
        {
            samples[i].value = samples[i].value + other.samples[i].value;
        }
        return *this;
    }

    Wave &operator-=(const Wave &other)
    {
        if (::mathematics::absolute(sampleRate - other.sampleRate) > ::mathematics::tolerance)
        {
            throw std::runtime_error("Wave-=: rate mismatch");
        }
        size_t minSize = std::min(samples.size(), other.samples.size());
        for (size_t i = 0; i < minSize; ++i)
        {
            samples[i].value = samples[i].value - other.samples[i].value;
        }
        return *this;
    }

    Wave &operator*=(double scalar)
    {
        for (auto &sample : samples)
        {
            sample.value = sample.value * scalar;
        }
        return *this;
    }
    Wave &operator/=(double scalar)
    {
        for (auto &sample : samples)
        {
            sample.value = sample.value / scalar;
        }
        return *this;
    }

    friend Wave operator*(double scalar, const Wave &wave) { return wave * scalar; }

    friend std::ostream &operator<<(std::ostream &output, const Wave &wave)
    {
        for (size_t i = 0; i < wave.samples.size(); ++i)
        {
            output << wave.timeAt(i) << " s: " << wave.samples[i].value << "\n";
        }
        return output;
    }
};

inline Spectral<::mathematics::Complex<double>> spectralFromSample(const Sample<::mathematics::Complex<double>> &sample, double resolution = 0.0, double confidence = 1.0)
{
    return Spectral<::mathematics::Complex<double>>(sample.frequency, sample.value, resolution, confidence);
}

inline Sample<::mathematics::Complex<double>> sampleFromSpectral(const Spectral<::mathematics::Complex<double>> &spectral)
{
    return Sample<::mathematics::Complex<double>>(spectral.getFrequency(), spectral.getPhasor());
}

template <typename N, std::enable_if_t<is_physical_v<N>, int> = 0>
Domain<::mathematics::Complex<double>> domainFromSpectrum(const Spectrum<N> &spectrum)
{
    Domain<::mathematics::Complex<double>> domain;
    double frequencyStep = 0.0;
    const auto &lines = spectrum.getLines();
    if (lines.size() >= 2)
    {
        frequencyStep = lines[1].frequency.hertz() - lines[0].frequency.hertz();
    }
    for (const auto &line : lines)
    {
        ::mathematics::Complex<double> coefficient;
        if constexpr (is_scalar_v<N>)
        {
            coefficient = ::mathematics::Complex<double>(static_cast<double>(line.value), 0.0);
        }
        else if constexpr (is_phasor_v<N>)
        {
            coefficient = ::mathematics::Complex<double>(static_cast<double>(line.value.getReal()), static_cast<double>(line.value.getImag()));
        }
        else
        {
            coefficient = ::mathematics::Complex<double>(static_cast<double>(line.value.getReal()), static_cast<double>(line.value.getImagI()));
        }
        domain.addLine(Spectral<::mathematics::Complex<double>>(line.frequency, coefficient, frequencyStep, 1.0));
    }
    return domain;
}

inline Spectrum<::mathematics::Complex<double>> spectrumFromDomain(const Domain<::mathematics::Complex<double>> &domain)
{
    Spectrum<::mathematics::Complex<double>> spectrum;
    for (const auto &line : domain.getLines())
    {
        spectrum.addLine(line.getFrequency(), line.getPhasor());
    }
    return spectrum;
}

} // namespace physics