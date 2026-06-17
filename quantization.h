#pragma once
#include <complex>
#include <stdexcept>
#include <limits>
#include <iostream>
#include "mathematics.h"
#include "quaternion.h"

namespace physics
{
class Time
{
  private:
    double value;

  public:
    Time() : value(0.0) {}

    explicit Time(double seconds) : value(seconds)
    {
        if (::mathematics::abnormal(value))
        {
            throw std::domain_error("Time: NaN or Inf input");
        }
        value = ::mathematics::exact(value);
    }

    Time(const Time &other) : value(other.value) {}
    Time(Time &&other) noexcept : value(other.value) { other.value = 0.0; }

    Time &operator=(const Time &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Time &operator=(Time &&other) noexcept
    {
        if (this != &other)
        {
            value = other.value;
            other.value = 0.0;
        }
        return *this;
    }

    double seconds() const { return value; }
    double milliseconds() const { return ::mathematics::multiply(value, 1e3); }
    double microseconds() const { return ::mathematics::multiply(value, 1e6); }

    size_t samples(double sampleRate) const
    {
        if (::mathematics::abnormal(sampleRate))
        {
            throw std::domain_error("Time::samples: NaN or Inf sample rate");
        }
        if (sampleRate <= 0.0)
        {
            throw std::domain_error("Time::samples: non-positive sample rate");
        }
        double result = ::mathematics::multiply(value, sampleRate);
        if (::mathematics::abnormal(result))
        {
            throw std::overflow_error("Time::samples: result overflow");
        }
        return static_cast<size_t>(result);
    }

    Time operator+(const Time &other) const { return Time(::mathematics::add(value, other.value)); }
    Time operator-(const Time &other) const { return Time(::mathematics::subtract(value, other.value)); }
    Time operator*(double scalar) const { return Time(::mathematics::multiply(value, scalar)); }
    Time operator/(double scalar) const { return Time(::mathematics::divide(value, scalar)); }
    double operator/(const Time &other) const { return ::mathematics::divide(value, other.value); }
    Time operator-() const { return Time(-value); }

    Time &operator+=(const Time &other)
    {
        value = ::mathematics::add(value, other.value);
        return *this;
    }
    Time &operator-=(const Time &other)
    {
        value = ::mathematics::subtract(value, other.value);
        return *this;
    }
    Time &operator*=(double scalar)
    {
        value = ::mathematics::multiply(value, scalar);
        return *this;
    }
    Time &operator/=(double scalar)
    {
        value = ::mathematics::divide(value, scalar);
        return *this;
    }

    bool operator==(const Time &other) const { return ::mathematics::absolute(::mathematics::subtract(value, other.value)) < ::mathematics::tolerance; }
    bool operator!=(const Time &other) const { return !(*this == other); }
    bool operator<(const Time &other) const { return value < other.value; }
    bool operator<=(const Time &other) const { return value <= other.value; }
    bool operator>(const Time &other) const { return value > other.value; }
    bool operator>=(const Time &other) const { return value >= other.value; }

    friend Time operator*(double scalar, const Time &time) { return time * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Time &time) { return output << time.value << " s"; }
};

class Frequency
{
  private:
    ::mathematics::Complex<double> value;

  public:
    Frequency() : value(0.0, 0.0) {}
    Frequency(double sigma, double omega) : value(sigma, omega) {}

    Frequency(const Frequency &other) : value(other.value) {}
    Frequency(Frequency &&other) noexcept : value(std::move(other.value)) {}

    Frequency &operator=(const Frequency &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Frequency &operator=(Frequency &&other) noexcept
    {
        if (this != &other)
        {
            value = std::move(other.value);
        }
        return *this;
    }

    double sigma() const { return value.getReal(); }
    double omega() const { return value.getImag(); }
    double hertz() const { return ::mathematics::divide(value.getImag(), 2.0 * ::mathematics::pi); }
    double neper() const { return value.getReal(); }

    double dampingRatio() const
    {
        double magnitude = value.module();
        if (magnitude < ::mathematics::tolerance)
        {
            return 0.0;
        }
        return ::mathematics::divide(-value.getReal(), magnitude);
    }

    double qualityFactor() const
    {
        double damping = dampingRatio();
        if (damping < ::mathematics::tolerance)
        {
            return std::numeric_limits<double>::infinity();
        }
        return ::mathematics::divide(0.5, damping);
    }

    ::mathematics::Complex<double> evolution(double timeSeconds) const
    {
        ::mathematics::Complex<double> exponent(::mathematics::multiply(value.getReal(), timeSeconds), ::mathematics::multiply(value.getImag(), timeSeconds));
        return exponent.exponential();
    }

    const ::mathematics::Complex<double> &laplace() const { return value; }

    Frequency operator+(const Frequency &other) const
    {
        ::mathematics::Complex<double> result = value + other.value;
        return Frequency(result.getReal(), result.getImag());
    }
    Frequency operator-(const Frequency &other) const
    {
        ::mathematics::Complex<double> result = value - other.value;
        return Frequency(result.getReal(), result.getImag());
    }
    Frequency operator*(double scalar) const
    {
        ::mathematics::Complex<double> result = value * scalar;
        return Frequency(result.getReal(), result.getImag());
    }
    Frequency operator/(double scalar) const
    {
        ::mathematics::Complex<double> result = value / scalar;
        return Frequency(result.getReal(), result.getImag());
    }
    Frequency operator-() const
    {
        ::mathematics::Complex<double> result = -value;
        return Frequency(result.getReal(), result.getImag());
    }

    Frequency &operator+=(const Frequency &other)
    {
        value = value + other.value;
        return *this;
    }
    Frequency &operator-=(const Frequency &other)
    {
        value = value - other.value;
        return *this;
    }
    Frequency &operator*=(double scalar)
    {
        value = value * scalar;
        return *this;
    }
    Frequency &operator/=(double scalar)
    {
        value = value / scalar;
        return *this;
    }

    bool operator==(const Frequency &other) const { return (value - other.value).module() < ::mathematics::tolerance; }
    bool operator!=(const Frequency &other) const { return !(*this == other); }
    bool operator<(const Frequency &other) const
    {
        if (::mathematics::absolute(::mathematics::subtract(omega(), other.omega())) > ::mathematics::tolerance)
        {
            return omega() < other.omega();
        }
        return sigma() < other.sigma();
    }
    bool operator<=(const Frequency &other) const { return (*this < other) || (*this == other); }
    bool operator>(const Frequency &other) const { return !(*this <= other); }
    bool operator>=(const Frequency &other) const { return !(*this < other); }

    friend Frequency operator*(double scalar, const Frequency &frequency) { return frequency * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Frequency &frequency) { return output << "(" << frequency.sigma() << " + j" << frequency.omega() << ") rad/s"; }
};

class Period
{
  private:
    double duration;
    double bandwidth;

  public:
    Period(double durationSeconds, double bandwidthHertz) : duration(durationSeconds), bandwidth(bandwidthHertz)
    {
        if (::mathematics::abnormal(duration) || ::mathematics::abnormal(bandwidth))
        {
            throw std::domain_error("Period: NaN or Inf input");
        }
        if (duration <= 0.0)
        {
            throw std::domain_error("Period: non-positive duration");
        }
        if (bandwidth <= 0.0)
        {
            throw std::domain_error("Period: non-positive bandwidth");
        }
        if (::mathematics::multiply(duration, bandwidth) < ::mathematics::subtract(::mathematics::divide(1.0, 4.0 * ::mathematics::pi), ::mathematics::tolerance))
        {
            throw std::domain_error("Period: violates Gabor uncertainty principle");
        }
        duration = ::mathematics::exact(duration);
        bandwidth = ::mathematics::exact(bandwidth);
    }

    Period(const Period &other) : duration(other.duration), bandwidth(other.bandwidth) {}
    Period(Period &&other) noexcept : duration(other.duration), bandwidth(other.bandwidth)
    {
        other.duration = 1.0;
        other.bandwidth = 1.0;
    }

    Period &operator=(const Period &other)
    {
        if (this != &other)
        {
            duration = other.duration;
            bandwidth = other.bandwidth;
        }
        return *this;
    }
    Period &operator=(Period &&other) noexcept
    {
        if (this != &other)
        {
            duration = other.duration;
            bandwidth = other.bandwidth;
            other.duration = 1.0;
            other.bandwidth = 1.0;
        }
        return *this;
    }

    double getDuration() const { return duration; }
    double getBandwidth() const { return bandwidth; }
    double timeBandwidthProduct() const { return ::mathematics::multiply(duration, bandwidth); }
    double frequencyResolution() const { return ::mathematics::divide(1.0, duration); }
    double timeResolution() const { return ::mathematics::divide(1.0, bandwidth); }
    size_t nyquistSamples() const { return static_cast<size_t>(::mathematics::multiply(::mathematics::multiply(2.0, duration), bandwidth)); }
    size_t samples(double sampleRate) const { return static_cast<size_t>(::mathematics::multiply(duration, sampleRate)); }

    Period operator+(const Period &other) const { return Period(::mathematics::add(duration, other.duration), std::max(bandwidth, other.bandwidth)); }
    Period operator*(double scalar) const { return Period(::mathematics::multiply(duration, scalar), ::mathematics::divide(bandwidth, scalar)); }
    Period operator/(double scalar) const { return Period(::mathematics::divide(duration, scalar), ::mathematics::multiply(bandwidth, scalar)); }
    double operator/(const Period &other) const { return ::mathematics::divide(::mathematics::multiply(duration, bandwidth), ::mathematics::multiply(other.duration, other.bandwidth)); }

    Period &operator*=(double scalar)
    {
        duration = ::mathematics::multiply(duration, scalar);
        bandwidth = ::mathematics::divide(bandwidth, scalar);
        return *this;
    }
    Period &operator/=(double scalar)
    {
        duration = ::mathematics::divide(duration, scalar);
        bandwidth = ::mathematics::multiply(bandwidth, scalar);
        return *this;
    }

    bool operator==(const Period &other) const { return ::mathematics::absolute(::mathematics::subtract(duration, other.duration)) < ::mathematics::tolerance && ::mathematics::absolute(::mathematics::subtract(bandwidth, other.bandwidth)) < ::mathematics::tolerance; }
    bool operator!=(const Period &other) const { return !(*this == other); }

    friend Period operator*(double scalar, const Period &period) { return period * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Period &period) { return output << "T=" << period.duration << "s B=" << period.bandwidth << "Hz"; }
};

class Phasor
{
  private:
    ::mathematics::Complex<double> value;

  public:
    Phasor() : value(0.0, 0.0) {}
    explicit Phasor(::mathematics::Complex<double> complexValue) : value(complexValue) {}

    Phasor(double amplitude, double phaseRadians)
    {
        if (::mathematics::abnormal(amplitude) || ::mathematics::abnormal(phaseRadians))
        {
            throw std::domain_error("Phasor: NaN or Inf input");
        }
        if (amplitude < 0.0)
        {
            throw std::domain_error("Phasor: negative amplitude");
        }
        value = ::mathematics::polarComplex(amplitude, phaseRadians);
        value = value.exact();
    }

    Phasor(const Phasor &other) : value(other.value) {}
    Phasor(Phasor &&other) noexcept : value(std::move(other.value)) {}

    Phasor &operator=(const Phasor &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Phasor &operator=(Phasor &&other) noexcept
    {
        if (this != &other)
        {
            value = std::move(other.value);
        }
        return *this;
    }

    double amplitude() const { return value.module(); }
    double phase() const { return value.argument(); }
    double real() const { return value.getReal(); }
    double imaginary() const { return value.getImag(); }
    double power() const { return value.dot(value); }
    const ::mathematics::Complex<double> &complex() const { return value; }

    double decibels(double reference) const
    {
        double amplitudeValue = amplitude();
        if (amplitudeValue < ::mathematics::tolerance)
        {
            return -std::numeric_limits<double>::infinity();
        }
        return ::mathematics::multiply(20.0, ::mathematics::logarithm(10.0, ::mathematics::divide(amplitudeValue, reference)));
    }

    Phasor evolve(const Frequency &frequency, double timeSeconds) const { return Phasor(value * frequency.evolution(timeSeconds)); }
    Phasor evolve(const Frequency &frequency, const Period &period) const { return evolve(frequency, period.getDuration()); }

    double resonate(const Frequency &frequency) const
    {
        Phasor response = Phasor(value * frequency.evolution(0.0));
        double signalEnergy = power();
        if (signalEnergy < ::mathematics::tolerance)
        {
            return 0.0;
        }
        return ::mathematics::limit(::mathematics::divide(response.power(), signalEnergy), 0.0, 1.0);
    }

    Phasor operator+(const Phasor &other) const { return Phasor(value + other.value); }
    Phasor operator-(const Phasor &other) const { return Phasor(value - other.value); }
    Phasor operator*(const Phasor &other) const { return Phasor(value * other.value); }
    Phasor operator/(const Phasor &other) const { return Phasor(value / other.value); }
    Phasor operator*(double scalar) const { return Phasor(value * scalar); }
    Phasor operator/(double scalar) const { return Phasor(value / scalar); }
    Phasor operator-() const { return Phasor(-value); }

    Phasor &operator+=(const Phasor &other)
    {
        value = value + other.value;
        return *this;
    }
    Phasor &operator-=(const Phasor &other)
    {
        value = value - other.value;
        return *this;
    }
    Phasor &operator*=(const Phasor &other)
    {
        value = value * other.value;
        return *this;
    }
    Phasor &operator/=(const Phasor &other)
    {
        value = value / other.value;
        return *this;
    }
    Phasor &operator*=(double scalar)
    {
        value = value * scalar;
        return *this;
    }
    Phasor &operator/=(double scalar)
    {
        value = value / scalar;
        return *this;
    }

    bool operator==(const Phasor &other) const { return (value - other.value).module() < ::mathematics::tolerance; }
    bool operator!=(const Phasor &other) const { return !(*this == other); }

    friend Phasor operator*(double scalar, const Phasor &phasor) { return phasor * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Phasor &phasor) { return output << phasor.amplitude() << " angle " << phasor.phase() << " rad"; }
};

class Voltage
{
  private:
    Phasor value;

  public:
    Voltage() : value() {}
    Voltage(double amplitude, double phaseRadians = 0.0) : value(amplitude, phaseRadians) {}
    explicit Voltage(const Phasor &phasor) : value(phasor) {}

    Voltage(const Voltage &other) : value(other.value) {}
    Voltage(Voltage &&other) noexcept : value(std::move(other.value)) {}

    Voltage &operator=(const Voltage &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Voltage &operator=(Voltage &&other) noexcept
    {
        if (this != &other)
        {
            value = std::move(other.value);
        }
        return *this;
    }

    double amplitude() const { return value.amplitude(); }
    double phase() const { return value.phase(); }
    double rootMeanSquare() const { return ::mathematics::divide(value.amplitude(), ::mathematics::squareRoot(2.0)); }
    double decibels(double reference) const { return value.decibels(reference); }
    double decibelsMicrovolt() const { return value.decibels(1e-6); }
    const Phasor &phasor() const { return value; }

    Voltage operator+(const Voltage &other) const { return Voltage(value + other.value); }
    Voltage operator-(const Voltage &other) const { return Voltage(value - other.value); }
    Voltage operator*(double scalar) const { return Voltage(value * scalar); }
    Voltage operator/(double scalar) const { return Voltage(value / scalar); }
    Voltage operator-() const { return Voltage(-value); }

    Voltage &operator+=(const Voltage &other)
    {
        value += other.value;
        return *this;
    }
    Voltage &operator-=(const Voltage &other)
    {
        value -= other.value;
        return *this;
    }
    Voltage &operator*=(double scalar)
    {
        value *= scalar;
        return *this;
    }
    Voltage &operator/=(double scalar)
    {
        value /= scalar;
        return *this;
    }

    bool operator==(const Voltage &other) const { return value == other.value; }
    bool operator!=(const Voltage &other) const { return !(*this == other); }

    friend Voltage operator*(double scalar, const Voltage &voltage) { return voltage * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Voltage &voltage) { return output << voltage.amplitude() << " V angle " << voltage.phase() << " rad"; }
};

class Current
{
  private:
    Phasor value;

  public:
    Current() : value() {}
    Current(double amplitude, double phaseRadians = 0.0) : value(amplitude, phaseRadians) {}
    explicit Current(const Phasor &phasor) : value(phasor) {}

    Current(const Current &other) : value(other.value) {}
    Current(Current &&other) noexcept : value(std::move(other.value)) {}

    Current &operator=(const Current &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Current &operator=(Current &&other) noexcept
    {
        if (this != &other)
        {
            value = std::move(other.value);
        }
        return *this;
    }

    double amplitude() const { return value.amplitude(); }
    double phase() const { return value.phase(); }
    double rootMeanSquare() const { return ::mathematics::divide(value.amplitude(), ::mathematics::squareRoot(2.0)); }
    double decibels(double reference) const { return value.decibels(reference); }
    const Phasor &phasor() const { return value; }

    Current operator+(const Current &other) const { return Current(value + other.value); }
    Current operator-(const Current &other) const { return Current(value - other.value); }
    Current operator*(double scalar) const { return Current(value * scalar); }
    Current operator/(double scalar) const { return Current(value / scalar); }
    Current operator-() const { return Current(-value); }

    Current &operator+=(const Current &other)
    {
        value += other.value;
        return *this;
    }
    Current &operator-=(const Current &other)
    {
        value -= other.value;
        return *this;
    }
    Current &operator*=(double scalar)
    {
        value *= scalar;
        return *this;
    }
    Current &operator/=(double scalar)
    {
        value /= scalar;
        return *this;
    }

    bool operator==(const Current &other) const { return value == other.value; }
    bool operator!=(const Current &other) const { return !(*this == other); }

    friend Current operator*(double scalar, const Current &current) { return current * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Current &current) { return output << current.amplitude() << " A angle " << current.phase() << " rad"; }
};

class Resistance
{
  private:
    double value;

  public:
    Resistance() : value(0.0) {}

    explicit Resistance(double ohms) : value(ohms)
    {
        if (::mathematics::abnormal(value))
        {
            throw std::domain_error("Resistance: NaN or Inf input");
        }
        if (value < 0.0)
        {
            throw std::domain_error("Resistance: negative resistance");
        }
        value = ::mathematics::exact(value);
    }

    Resistance(const Resistance &other) : value(other.value) {}
    Resistance(Resistance &&other) noexcept : value(other.value) { other.value = 0.0; }

    Resistance &operator=(const Resistance &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Resistance &operator=(Resistance &&other) noexcept
    {
        if (this != &other)
        {
            value = other.value;
            other.value = 0.0;
        }
        return *this;
    }

    double ohms() const { return value; }

    Resistance operator+(const Resistance &other) const { return Resistance(::mathematics::add(value, other.value)); }
    Resistance operator*(double scalar) const { return Resistance(::mathematics::multiply(value, scalar)); }
    Resistance operator/(double scalar) const { return Resistance(::mathematics::divide(value, scalar)); }

    Resistance &operator+=(const Resistance &other)
    {
        value = ::mathematics::add(value, other.value);
        return *this;
    }
    Resistance &operator*=(double scalar)
    {
        value = ::mathematics::multiply(value, scalar);
        return *this;
    }
    Resistance &operator/=(double scalar)
    {
        value = ::mathematics::divide(value, scalar);
        return *this;
    }

    bool operator==(const Resistance &other) const { return ::mathematics::absolute(::mathematics::subtract(value, other.value)) < ::mathematics::tolerance; }
    bool operator!=(const Resistance &other) const { return !(*this == other); }
    bool operator<(const Resistance &other) const { return value < other.value; }
    bool operator<=(const Resistance &other) const { return value <= other.value; }
    bool operator>(const Resistance &other) const { return value > other.value; }
    bool operator>=(const Resistance &other) const { return value >= other.value; }

    friend Resistance operator*(double scalar, const Resistance &resistance) { return resistance * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Resistance &resistance) { return output << resistance.value << " Ohm"; }
};

class Conductance
{
  private:
    double value;

  public:
    Conductance() : value(0.0) {}

    explicit Conductance(double siemens) : value(siemens)
    {
        if (::mathematics::abnormal(value))
        {
            throw std::domain_error("Conductance: NaN or Inf input");
        }
        if (value < 0.0)
        {
            throw std::domain_error("Conductance: negative conductance");
        }
        value = ::mathematics::exact(value);
    }

    Conductance(const Conductance &other) : value(other.value) {}
    Conductance(Conductance &&other) noexcept : value(other.value) { other.value = 0.0; }

    Conductance &operator=(const Conductance &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Conductance &operator=(Conductance &&other) noexcept
    {
        if (this != &other)
        {
            value = other.value;
            other.value = 0.0;
        }
        return *this;
    }

    double siemens() const { return value; }

    Conductance operator+(const Conductance &other) const { return Conductance(::mathematics::add(value, other.value)); }
    Conductance operator*(double scalar) const { return Conductance(::mathematics::multiply(value, scalar)); }
    Conductance operator/(double scalar) const { return Conductance(::mathematics::divide(value, scalar)); }

    Conductance &operator+=(const Conductance &other)
    {
        value = ::mathematics::add(value, other.value);
        return *this;
    }
    Conductance &operator*=(double scalar)
    {
        value = ::mathematics::multiply(value, scalar);
        return *this;
    }
    Conductance &operator/=(double scalar)
    {
        value = ::mathematics::divide(value, scalar);
        return *this;
    }

    bool operator==(const Conductance &other) const { return ::mathematics::absolute(::mathematics::subtract(value, other.value)) < ::mathematics::tolerance; }
    bool operator!=(const Conductance &other) const { return !(*this == other); }
    bool operator<(const Conductance &other) const { return value < other.value; }
    bool operator<=(const Conductance &other) const { return value <= other.value; }
    bool operator>(const Conductance &other) const { return value > other.value; }
    bool operator>=(const Conductance &other) const { return value >= other.value; }

    friend Conductance operator*(double scalar, const Conductance &conductance) { return conductance * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Conductance &conductance) { return output << conductance.value << " S"; }
};

class Capacitance
{
  private:
    double value;

  public:
    Capacitance() : value(0.0) {}

    explicit Capacitance(double farads) : value(farads)
    {
        if (::mathematics::abnormal(value))
        {
            throw std::domain_error("Capacitance: NaN or Inf input");
        }
        if (value < 0.0)
        {
            throw std::domain_error("Capacitance: negative capacitance");
        }
        value = ::mathematics::exact(value);
    }

    Capacitance(const Capacitance &other) : value(other.value) {}
    Capacitance(Capacitance &&other) noexcept : value(other.value) { other.value = 0.0; }

    Capacitance &operator=(const Capacitance &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Capacitance &operator=(Capacitance &&other) noexcept
    {
        if (this != &other)
        {
            value = other.value;
            other.value = 0.0;
        }
        return *this;
    }

    double farads() const { return value; }

    Capacitance operator+(const Capacitance &other) const { return Capacitance(::mathematics::add(value, other.value)); }
    Capacitance operator*(double scalar) const { return Capacitance(::mathematics::multiply(value, scalar)); }
    Capacitance operator/(double scalar) const { return Capacitance(::mathematics::divide(value, scalar)); }

    Capacitance &operator+=(const Capacitance &other)
    {
        value = ::mathematics::add(value, other.value);
        return *this;
    }
    Capacitance &operator*=(double scalar)
    {
        value = ::mathematics::multiply(value, scalar);
        return *this;
    }
    Capacitance &operator/=(double scalar)
    {
        value = ::mathematics::divide(value, scalar);
        return *this;
    }

    bool operator==(const Capacitance &other) const { return ::mathematics::absolute(::mathematics::subtract(value, other.value)) < ::mathematics::tolerance; }
    bool operator!=(const Capacitance &other) const { return !(*this == other); }
    bool operator<(const Capacitance &other) const { return value < other.value; }
    bool operator<=(const Capacitance &other) const { return value <= other.value; }
    bool operator>(const Capacitance &other) const { return value > other.value; }
    bool operator>=(const Capacitance &other) const { return value >= other.value; }

    friend Capacitance operator*(double scalar, const Capacitance &capacitance) { return capacitance * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Capacitance &capacitance) { return output << capacitance.value << " F"; }
};

class Inductance
{
  private:
    double value;

  public:
    Inductance() : value(0.0) {}

    explicit Inductance(double henries) : value(henries)
    {
        if (::mathematics::abnormal(value))
        {
            throw std::domain_error("Inductance: NaN or Inf input");
        }
        if (value < 0.0)
        {
            throw std::domain_error("Inductance: negative inductance");
        }
        value = ::mathematics::exact(value);
    }

    Inductance(const Inductance &other) : value(other.value) {}
    Inductance(Inductance &&other) noexcept : value(other.value) { other.value = 0.0; }

    Inductance &operator=(const Inductance &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Inductance &operator=(Inductance &&other) noexcept
    {
        if (this != &other)
        {
            value = other.value;
            other.value = 0.0;
        }
        return *this;
    }

    double henries() const { return value; }

    Inductance operator+(const Inductance &other) const { return Inductance(::mathematics::add(value, other.value)); }
    Inductance operator*(double scalar) const { return Inductance(::mathematics::multiply(value, scalar)); }
    Inductance operator/(double scalar) const { return Inductance(::mathematics::divide(value, scalar)); }

    Inductance &operator+=(const Inductance &other)
    {
        value = ::mathematics::add(value, other.value);
        return *this;
    }
    Inductance &operator*=(double scalar)
    {
        value = ::mathematics::multiply(value, scalar);
        return *this;
    }
    Inductance &operator/=(double scalar)
    {
        value = ::mathematics::divide(value, scalar);
        return *this;
    }

    bool operator==(const Inductance &other) const { return ::mathematics::absolute(::mathematics::subtract(value, other.value)) < ::mathematics::tolerance; }
    bool operator!=(const Inductance &other) const { return !(*this == other); }
    bool operator<(const Inductance &other) const { return value < other.value; }
    bool operator<=(const Inductance &other) const { return value <= other.value; }
    bool operator>(const Inductance &other) const { return value > other.value; }
    bool operator>=(const Inductance &other) const { return value >= other.value; }

    friend Inductance operator*(double scalar, const Inductance &inductance) { return inductance * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Inductance &inductance) { return output << inductance.value << " H"; }
};

class Impedance
{
  private:
    ::mathematics::Complex<double> value;

  public:
    Impedance(double resistanceValue = 0.0, double reactanceValue = 0.0) : value(resistanceValue, reactanceValue)
    {
        if (::mathematics::abnormal(resistanceValue) || ::mathematics::abnormal(reactanceValue))
        {
            throw std::domain_error("Impedance: NaN or Inf input");
        }
        value = value.exact();
    }

    Impedance(const Impedance &other) : value(other.value) {}
    Impedance(Impedance &&other) noexcept : value(std::move(other.value)) {}

    Impedance &operator=(const Impedance &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Impedance &operator=(Impedance &&other) noexcept
    {
        if (this != &other)
        {
            value = std::move(other.value);
        }
        return *this;
    }

    double getResistance() const { return value.getReal(); }
    double getReactance() const { return value.getImag(); }
    double magnitude() const { return value.module(); }
    double phase() const { return value.argument(); }
    bool isInductive() const { return value.getImag() > ::mathematics::tolerance; }
    bool isCapacitive() const { return value.getImag() < -::mathematics::tolerance; }
    bool isResistive() const { return ::mathematics::absolute(value.getImag()) < ::mathematics::tolerance; }
    const ::mathematics::Complex<double> &complex() const { return value; }

    Impedance operator+(const Impedance &other) const
    {
        ::mathematics::Complex<double> result = value + other.value;
        return Impedance(result.getReal(), result.getImag());
    }
    Impedance operator-(const Impedance &other) const
    {
        ::mathematics::Complex<double> result = value - other.value;
        return Impedance(result.getReal(), result.getImag());
    }
    Impedance operator*(double scalar) const
    {
        ::mathematics::Complex<double> result = value * scalar;
        return Impedance(result.getReal(), result.getImag());
    }
    Impedance operator/(double scalar) const
    {
        ::mathematics::Complex<double> result = value / scalar;
        return Impedance(result.getReal(), result.getImag());
    }
    Impedance operator-() const
    {
        ::mathematics::Complex<double> result = -value;
        return Impedance(result.getReal(), result.getImag());
    }

    Impedance &operator+=(const Impedance &other)
    {
        value = value + other.value;
        return *this;
    }
    Impedance &operator-=(const Impedance &other)
    {
        value = value - other.value;
        return *this;
    }
    Impedance &operator*=(double scalar)
    {
        value = value * scalar;
        return *this;
    }
    Impedance &operator/=(double scalar)
    {
        value = value / scalar;
        return *this;
    }

    bool operator==(const Impedance &other) const { return (value - other.value).module() < ::mathematics::tolerance; }
    bool operator!=(const Impedance &other) const { return !(*this == other); }

    friend Impedance operator*(double scalar, const Impedance &impedance) { return impedance * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Impedance &impedance) { return output << impedance.getResistance() << " + j" << impedance.getReactance() << " Ohm"; }
};

class Admittance
{
  private:
    ::mathematics::Complex<double> value;

  public:
    Admittance(double conductanceValue = 0.0, double susceptanceValue = 0.0) : value(conductanceValue, susceptanceValue)
    {
        if (::mathematics::abnormal(conductanceValue) || ::mathematics::abnormal(susceptanceValue))
        {
            throw std::domain_error("Admittance: NaN or Inf input");
        }
        value = value.exact();
    }

    Admittance(const Admittance &other) : value(other.value) {}
    Admittance(Admittance &&other) noexcept : value(std::move(other.value)) {}

    Admittance &operator=(const Admittance &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Admittance &operator=(Admittance &&other) noexcept
    {
        if (this != &other)
        {
            value = std::move(other.value);
        }
        return *this;
    }

    double getConductance() const { return value.getReal(); }
    double getSusceptance() const { return value.getImag(); }
    double magnitude() const { return value.module(); }
    double phase() const { return value.argument(); }
    bool isCapacitive() const { return value.getImag() > ::mathematics::tolerance; }
    bool isInductive() const { return value.getImag() < -::mathematics::tolerance; }
    const ::mathematics::Complex<double> &complex() const { return value; }

    Admittance operator+(const Admittance &other) const
    {
        ::mathematics::Complex<double> result = value + other.value;
        return Admittance(result.getReal(), result.getImag());
    }
    Admittance operator-(const Admittance &other) const
    {
        ::mathematics::Complex<double> result = value - other.value;
        return Admittance(result.getReal(), result.getImag());
    }
    Admittance operator*(double scalar) const
    {
        ::mathematics::Complex<double> result = value * scalar;
        return Admittance(result.getReal(), result.getImag());
    }
    Admittance operator/(double scalar) const
    {
        ::mathematics::Complex<double> result = value / scalar;
        return Admittance(result.getReal(), result.getImag());
    }
    Admittance operator-() const
    {
        ::mathematics::Complex<double> result = -value;
        return Admittance(result.getReal(), result.getImag());
    }

    Admittance &operator+=(const Admittance &other)
    {
        value = value + other.value;
        return *this;
    }
    Admittance &operator-=(const Admittance &other)
    {
        value = value - other.value;
        return *this;
    }
    Admittance &operator*=(double scalar)
    {
        value = value * scalar;
        return *this;
    }
    Admittance &operator/=(double scalar)
    {
        value = value / scalar;
        return *this;
    }

    bool operator==(const Admittance &other) const { return (value - other.value).module() < ::mathematics::tolerance; }
    bool operator!=(const Admittance &other) const { return !(*this == other); }

    friend Admittance operator*(double scalar, const Admittance &admittance) { return admittance * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Admittance &admittance) { return output << admittance.getConductance() << " + j" << admittance.getSusceptance() << " S"; }
};

class Power
{
  private:
    ::mathematics::Complex<double> value;

  public:
    Power(double activeWatts = 0.0, double reactiveVoltAmperesReactive = 0.0) : value(activeWatts, reactiveVoltAmperesReactive)
    {
        if (::mathematics::abnormal(activeWatts) || ::mathematics::abnormal(reactiveVoltAmperesReactive))
        {
            throw std::domain_error("Power: NaN or Inf input");
        }
        value = value.exact();
    }

    Power(const Power &other) : value(other.value) {}
    Power(Power &&other) noexcept : value(std::move(other.value)) {}

    Power &operator=(const Power &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Power &operator=(Power &&other) noexcept
    {
        if (this != &other)
        {
            value = std::move(other.value);
        }
        return *this;
    }

    double active() const { return value.getReal(); }
    double reactive() const { return value.getImag(); }
    double apparent() const { return value.module(); }

    double powerFactor() const
    {
        double apparentPower = apparent();
        if (apparentPower < ::mathematics::tolerance)
        {
            return 0.0;
        }
        return ::mathematics::limit(::mathematics::divide(active(), apparentPower), 0.0, 1.0);
    }

    const ::mathematics::Complex<double> &complex() const { return value; }

    Power operator+(const Power &other) const
    {
        ::mathematics::Complex<double> result = value + other.value;
        return Power(result.getReal(), result.getImag());
    }
    Power operator-(const Power &other) const
    {
        ::mathematics::Complex<double> result = value - other.value;
        return Power(result.getReal(), result.getImag());
    }
    Power operator*(double scalar) const
    {
        ::mathematics::Complex<double> result = value * scalar;
        return Power(result.getReal(), result.getImag());
    }
    Power operator/(double scalar) const
    {
        ::mathematics::Complex<double> result = value / scalar;
        return Power(result.getReal(), result.getImag());
    }
    Power operator-() const
    {
        ::mathematics::Complex<double> result = -value;
        return Power(result.getReal(), result.getImag());
    }

    Power &operator+=(const Power &other)
    {
        value = value + other.value;
        return *this;
    }
    Power &operator-=(const Power &other)
    {
        value = value - other.value;
        return *this;
    }
    Power &operator*=(double scalar)
    {
        value = value * scalar;
        return *this;
    }
    Power &operator/=(double scalar)
    {
        value = value / scalar;
        return *this;
    }

    bool operator==(const Power &other) const { return (value - other.value).module() < ::mathematics::tolerance; }
    bool operator!=(const Power &other) const { return !(*this == other); }

    friend Power operator*(double scalar, const Power &power) { return power * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Power &power) { return output << power.active() << " W + j" << power.reactive() << " VAR"; }
};

class Energy
{
  private:
    ::mathematics::Complex<double> value;

  public:
    Energy(double storedJoules = 0.0, double dissipatedWatts = 0.0) : value(storedJoules, dissipatedWatts)
    {
        if (::mathematics::abnormal(storedJoules) || ::mathematics::abnormal(dissipatedWatts))
        {
            throw std::domain_error("Energy: NaN or Inf input");
        }
        if (storedJoules < 0.0)
        {
            throw std::domain_error("Energy: negative stored energy");
        }
        if (dissipatedWatts < 0.0)
        {
            throw std::domain_error("Energy: negative dissipated power");
        }
        value = value.exact();
    }

    Energy(const Energy &other) : value(other.value) {}
    Energy(Energy &&other) noexcept : value(std::move(other.value)) {}

    Energy &operator=(const Energy &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Energy &operator=(Energy &&other) noexcept
    {
        if (this != &other)
        {
            value = std::move(other.value);
        }
        return *this;
    }

    double stored() const { return value.getReal(); }
    double dissipated() const { return value.getImag(); }
    double total() const { return value.module(); }
    const ::mathematics::Complex<double> &complex() const { return value; }

    Energy operator+(const Energy &other) const
    {
        ::mathematics::Complex<double> result = value + other.value;
        return Energy(result.getReal(), result.getImag());
    }
    Energy operator-(const Energy &other) const
    {
        ::mathematics::Complex<double> result = value - other.value;
        return Energy(result.getReal(), result.getImag());
    }
    Energy operator*(double scalar) const
    {
        ::mathematics::Complex<double> result = value * scalar;
        return Energy(result.getReal(), result.getImag());
    }
    Energy operator/(double scalar) const
    {
        ::mathematics::Complex<double> result = value / scalar;
        return Energy(result.getReal(), result.getImag());
    }
    Energy operator-() const
    {
        ::mathematics::Complex<double> result = -value;
        return Energy(result.getReal(), result.getImag());
    }

    Energy &operator+=(const Energy &other)
    {
        value = value + other.value;
        return *this;
    }
    Energy &operator-=(const Energy &other)
    {
        value = value - other.value;
        return *this;
    }
    Energy &operator*=(double scalar)
    {
        value = value * scalar;
        return *this;
    }
    Energy &operator/=(double scalar)
    {
        value = value / scalar;
        return *this;
    }

    bool operator==(const Energy &other) const { return (value - other.value).module() < ::mathematics::tolerance; }
    bool operator!=(const Energy &other) const { return !(*this == other); }

    friend Energy operator*(double scalar, const Energy &energy) { return energy * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Energy &energy) { return output << "E(" << energy.stored() << " J, " << energy.dissipated() << " W)"; }
};

class Entropy
{
  private:
    ::mathematics::Complex<double> value;

  public:
    Entropy(double informationNats = 0.0, double productionRateNatsPerSec = 0.0) : value(informationNats, productionRateNatsPerSec)
    {
        if (::mathematics::abnormal(informationNats) || ::mathematics::abnormal(productionRateNatsPerSec))
        {
            throw std::domain_error("Entropy: NaN or Inf input");
        }
        if (informationNats < 0.0)
        {
            throw std::domain_error("Entropy: negative information");
        }
        if (productionRateNatsPerSec < 0.0)
        {
            throw std::domain_error("Entropy: negative production rate");
        }
        value = value.exact();
    }

    Entropy(const Entropy &other) : value(other.value) {}
    Entropy(Entropy &&other) noexcept : value(std::move(other.value)) {}

    Entropy &operator=(const Entropy &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Entropy &operator=(Entropy &&other) noexcept
    {
        if (this != &other)
        {
            value = std::move(other.value);
        }
        return *this;
    }

    double information() const { return value.getReal(); }
    double productionRate() const { return value.getImag(); }
    const ::mathematics::Complex<double> &complex() const { return value; }

    Entropy operator+(const Entropy &other) const
    {
        ::mathematics::Complex<double> result = value + other.value;
        return Entropy(result.getReal(), result.getImag());
    }
    Entropy operator-(const Entropy &other) const
    {
        ::mathematics::Complex<double> result = value - other.value;
        return Entropy(result.getReal(), result.getImag());
    }
    Entropy operator*(double scalar) const
    {
        ::mathematics::Complex<double> result = value * scalar;
        return Entropy(result.getReal(), result.getImag());
    }
    Entropy operator/(double scalar) const
    {
        ::mathematics::Complex<double> result = value / scalar;
        return Entropy(result.getReal(), result.getImag());
    }
    Entropy operator-() const
    {
        ::mathematics::Complex<double> result = -value;
        return Entropy(result.getReal(), result.getImag());
    }

    Entropy &operator+=(const Entropy &other)
    {
        value = value + other.value;
        return *this;
    }
    Entropy &operator-=(const Entropy &other)
    {
        value = value - other.value;
        return *this;
    }
    Entropy &operator*=(double scalar)
    {
        value = value * scalar;
        return *this;
    }
    Entropy &operator/=(double scalar)
    {
        value = value / scalar;
        return *this;
    }

    bool operator==(const Entropy &other) const { return (value - other.value).module() < ::mathematics::tolerance; }
    bool operator!=(const Entropy &other) const { return !(*this == other); }

    friend Entropy operator*(double scalar, const Entropy &entropy) { return entropy * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Entropy &entropy) { return output << "S(" << entropy.information() << " nat, " << entropy.productionRate() << " nat/s)"; }
};

class Noise
{
  private:
    double powerSpectralDensity;
    double temperature;
    double frequency;

  public:
    Noise(double powerSpectralDensityValue = 0.0, double temperatureKelvin = 290.0, double frequencyHertz = 0.0)
        : powerSpectralDensity(powerSpectralDensityValue), temperature(temperatureKelvin), frequency(frequencyHertz)
    {
        if (::mathematics::abnormal(powerSpectralDensity) || ::mathematics::abnormal(temperature) || ::mathematics::abnormal(frequency))
        {
            throw std::domain_error("Noise: NaN or Inf input");
        }
        if (powerSpectralDensity < 0.0)
        {
            throw std::domain_error("Noise: negative PSD");
        }
        if (temperature < 0.0)
        {
            throw std::domain_error("Noise: negative temperature");
        }
        if (frequency < 0.0)
        {
            throw std::domain_error("Noise: negative frequency");
        }
    }

    Noise(const Noise &other) : powerSpectralDensity(other.powerSpectralDensity), temperature(other.temperature), frequency(other.frequency) {}
    Noise(Noise &&other) noexcept : powerSpectralDensity(other.powerSpectralDensity), temperature(other.temperature), frequency(other.frequency)
    {
        other.powerSpectralDensity = 0.0;
        other.temperature = 290.0;
        other.frequency = 0.0;
    }

    Noise &operator=(const Noise &other)
    {
        if (this != &other)
        {
            powerSpectralDensity = other.powerSpectralDensity;
            temperature = other.temperature;
            frequency = other.frequency;
        }
        return *this;
    }
    Noise &operator=(Noise &&other) noexcept
    {
        if (this != &other)
        {
            powerSpectralDensity = other.powerSpectralDensity;
            temperature = other.temperature;
            frequency = other.frequency;
            other.powerSpectralDensity = 0.0;
            other.temperature = 290.0;
            other.frequency = 0.0;
        }
        return *this;
    }

    double density() const { return powerSpectralDensity; }

    double totalPower(double bandwidthHertz) const { return ::mathematics::multiply(powerSpectralDensity, bandwidthHertz); }

    double rootMeanSquareVoltage(double resistanceOhms, double bandwidthHertz) const
    {
        if (::mathematics::abnormal(resistanceOhms) || ::mathematics::abnormal(bandwidthHertz))
        {
            throw std::domain_error("Noise::rootMeanSquareVoltage: NaN or Inf input");
        }
        if (resistanceOhms < 0.0)
        {
            throw std::domain_error("Noise::rootMeanSquareVoltage: negative resistance");
        }
        if (bandwidthHertz < 0.0)
        {
            throw std::domain_error("Noise::rootMeanSquareVoltage: negative bandwidth");
        }
        double meanSquare = ::mathematics::multiply(::mathematics::multiply(powerSpectralDensity, resistanceOhms), bandwidthHertz);
        if (::mathematics::abnormal(meanSquare))
        {
            throw std::overflow_error("Noise::rootMeanSquareVoltage: overflow");
        }
        return ::mathematics::squareRoot(meanSquare);
    }

    double noiseTemperature(double resistanceOhms) const
    {
        const double boltzmannConstant = 1.380649e-23;
        return ::mathematics::divide(powerSpectralDensity, ::mathematics::multiply(::mathematics::multiply(4.0, boltzmannConstant), resistanceOhms));
    }

    Noise operator+(const Noise &other) const { return Noise(::mathematics::add(powerSpectralDensity, other.powerSpectralDensity)); }
    Noise operator-(const Noise &other) const { return Noise(::mathematics::subtract(powerSpectralDensity, other.powerSpectralDensity)); }
    Noise operator*(double scalar) const { return Noise(::mathematics::multiply(powerSpectralDensity, scalar)); }
    Noise operator/(double scalar) const { return Noise(::mathematics::divide(powerSpectralDensity, scalar)); }
    Noise operator-() const { return Noise(-powerSpectralDensity, temperature, frequency); }

    Noise &operator+=(const Noise &other)
    {
        powerSpectralDensity = ::mathematics::add(powerSpectralDensity, other.powerSpectralDensity);
        return *this;
    }
    Noise &operator-=(const Noise &other)
    {
        powerSpectralDensity = ::mathematics::subtract(powerSpectralDensity, other.powerSpectralDensity);
        return *this;
    }
    Noise &operator*=(double scalar)
    {
        powerSpectralDensity = ::mathematics::multiply(powerSpectralDensity, scalar);
        return *this;
    }
    Noise &operator/=(double scalar)
    {
        powerSpectralDensity = ::mathematics::divide(powerSpectralDensity, scalar);
        return *this;
    }

    bool operator==(const Noise &other) const { return ::mathematics::absolute(::mathematics::subtract(powerSpectralDensity, other.powerSpectralDensity)) < ::mathematics::tolerance; }
    bool operator!=(const Noise &other) const { return !(*this == other); }

    friend Noise operator*(double scalar, const Noise &noise) { return noise * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Noise &noise) { return output << noise.powerSpectralDensity << " V^2/Hz at " << noise.frequency << " Hz " << noise.temperature << " K"; }
};

class Charge
{
  private:
    double value;

  public:
    Charge() : value(0.0) {}

    explicit Charge(double coulombs) : value(coulombs)
    {
        if (::mathematics::abnormal(value))
        {
            throw std::domain_error("Charge: NaN or Inf input");
        }
        value = ::mathematics::exact(value);
    }

    Charge(const Charge &other) : value(other.value) {}
    Charge(Charge &&other) noexcept : value(other.value) { other.value = 0.0; }

    Charge &operator=(const Charge &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Charge &operator=(Charge &&other) noexcept
    {
        if (this != &other)
        {
            value = other.value;
            other.value = 0.0;
        }
        return *this;
    }

    double coulombs() const { return value; }
    double elementaryCharges() const { return ::mathematics::divide(value, 1.602176634e-19); }

    Charge operator+(const Charge &other) const { return Charge(::mathematics::add(value, other.value)); }
    Charge operator-(const Charge &other) const { return Charge(::mathematics::subtract(value, other.value)); }
    Charge operator*(double scalar) const { return Charge(::mathematics::multiply(value, scalar)); }
    Charge operator/(double scalar) const { return Charge(::mathematics::divide(value, scalar)); }
    Charge operator-() const { return Charge(-value); }

    Charge &operator+=(const Charge &other)
    {
        value = ::mathematics::add(value, other.value);
        return *this;
    }
    Charge &operator-=(const Charge &other)
    {
        value = ::mathematics::subtract(value, other.value);
        return *this;
    }
    Charge &operator*=(double scalar)
    {
        value = ::mathematics::multiply(value, scalar);
        return *this;
    }
    Charge &operator/=(double scalar)
    {
        value = ::mathematics::divide(value, scalar);
        return *this;
    }

    bool operator==(const Charge &other) const { return ::mathematics::absolute(::mathematics::subtract(value, other.value)) < ::mathematics::tolerance; }
    bool operator!=(const Charge &other) const { return !(*this == other); }
    bool operator<(const Charge &other) const { return value < other.value; }
    bool operator<=(const Charge &other) const { return value <= other.value; }
    bool operator>(const Charge &other) const { return value > other.value; }
    bool operator>=(const Charge &other) const { return value >= other.value; }

    friend Charge operator*(double scalar, const Charge &charge) { return charge * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Charge &charge) { return output << charge.value << " C"; }
};

class Flux
{
  private:
    double value;

  public:
    Flux() : value(0.0) {}

    explicit Flux(double webers) : value(webers)
    {
        if (::mathematics::abnormal(value))
        {
            throw std::domain_error("Flux: NaN or Inf input");
        }
        value = ::mathematics::exact(value);
    }

    Flux(const Flux &other) : value(other.value) {}
    Flux(Flux &&other) noexcept : value(other.value) { other.value = 0.0; }

    Flux &operator=(const Flux &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Flux &operator=(Flux &&other) noexcept
    {
        if (this != &other)
        {
            value = other.value;
            other.value = 0.0;
        }
        return *this;
    }

    double webers() const { return value; }
    double voltSeconds() const { return value; }
    double teslaSquareMeters() const { return value; }

    Flux operator+(const Flux &other) const { return Flux(::mathematics::add(value, other.value)); }
    Flux operator-(const Flux &other) const { return Flux(::mathematics::subtract(value, other.value)); }
    Flux operator*(double scalar) const { return Flux(::mathematics::multiply(value, scalar)); }
    Flux operator/(double scalar) const { return Flux(::mathematics::divide(value, scalar)); }
    double operator/(const Flux &other) const { return ::mathematics::divide(value, other.value); }
    Flux operator-() const { return Flux(-value); }

    Flux &operator+=(const Flux &other)
    {
        value = ::mathematics::add(value, other.value);
        return *this;
    }
    Flux &operator-=(const Flux &other)
    {
        value = ::mathematics::subtract(value, other.value);
        return *this;
    }
    Flux &operator*=(double scalar)
    {
        value = ::mathematics::multiply(value, scalar);
        return *this;
    }
    Flux &operator/=(double scalar)
    {
        value = ::mathematics::divide(value, scalar);
        return *this;
    }

    bool operator==(const Flux &other) const { return ::mathematics::absolute(::mathematics::subtract(value, other.value)) < ::mathematics::tolerance; }
    bool operator!=(const Flux &other) const { return !(*this == other); }
    bool operator<(const Flux &other) const { return value < other.value; }
    bool operator<=(const Flux &other) const { return value <= other.value; }
    bool operator>(const Flux &other) const { return value > other.value; }
    bool operator>=(const Flux &other) const { return value >= other.value; }

    friend Flux operator*(double scalar, const Flux &flux) { return flux * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Flux &flux) { return output << flux.value << " Wb"; }
};

class Displacement
{
  private:
    ::mathematics::Quaternion<double> value;

  public:
    Displacement() : value(0.0, 0.0, 0.0, 0.0) {}
    Displacement(double x, double y, double z) : value(0.0, x, y, z)
    {
        if (::mathematics::abnormal(x) || ::mathematics::abnormal(y) || ::mathematics::abnormal(z))
        {
            throw std::domain_error("Displacement: NaN or Inf input");
        }
    }

    explicit Displacement(const ::mathematics::Quaternion<double> &quaternionValue) : value(quaternionValue) {}

    Displacement(const Displacement &other) : value(other.value) {}
    Displacement(Displacement &&other) noexcept : value(std::move(other.value)) {}

    Displacement &operator=(const Displacement &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Displacement &operator=(Displacement &&other) noexcept
    {
        if (this != &other)
        {
            value = std::move(other.value);
        }
        return *this;
    }

    double x() const { return value.getImagI(); }
    double y() const { return value.getImagJ(); }
    double z() const { return value.getImagK(); }
    double magnitude() const { return value.module(); }
    const ::mathematics::Quaternion<double> &asQuaternion() const { return value; }

    double dot(const Displacement &other) const { return value.dot(other.value); }
    Displacement cross(const Displacement &other) const { return Displacement(value.cross(other.value)); }

    Displacement operator+(const Displacement &other) const { return Displacement(value + other.value); }
    Displacement operator-(const Displacement &other) const { return Displacement(value - other.value); }
    Displacement operator*(double scalar) const { return Displacement(value * scalar); }
    Displacement operator/(double scalar) const { return Displacement(value / scalar); }
    Displacement operator-() const { return Displacement(-value); }

    Displacement &operator+=(const Displacement &other)
    {
        value = value + other.value;
        return *this;
    }
    Displacement &operator-=(const Displacement &other)
    {
        value = value - other.value;
        return *this;
    }
    Displacement &operator*=(double scalar)
    {
        value = value * scalar;
        return *this;
    }
    Displacement &operator/=(double scalar)
    {
        value = value / scalar;
        return *this;
    }

    bool operator==(const Displacement &other) const { return value == other.value; }
    bool operator!=(const Displacement &other) const { return !(*this == other); }

    friend Displacement operator*(double scalar, const Displacement &displacement) { return displacement * scalar; }
    friend Displacement operator^(const Displacement &left, const Displacement &right) { return left.cross(right); }
    friend std::ostream &operator<<(std::ostream &output, const Displacement &displacement) { return output << "r(" << displacement.x() << ", " << displacement.y() << ", " << displacement.z() << ") m"; }
};

class Velocity
{
  private:
    ::mathematics::Quaternion<double> value;

  public:
    Velocity() : value(0.0, 0.0, 0.0, 0.0) {}
    Velocity(double x, double y, double z) : value(0.0, x, y, z)
    {
        if (::mathematics::abnormal(x) || ::mathematics::abnormal(y) || ::mathematics::abnormal(z))
        {
            throw std::domain_error("Velocity: NaN or Inf input");
        }
    }

    explicit Velocity(const ::mathematics::Quaternion<double> &quaternionValue) : value(quaternionValue) {}

    Velocity(const Velocity &other) : value(other.value) {}
    Velocity(Velocity &&other) noexcept : value(std::move(other.value)) {}

    Velocity &operator=(const Velocity &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Velocity &operator=(Velocity &&other) noexcept
    {
        if (this != &other)
        {
            value = std::move(other.value);
        }
        return *this;
    }

    double x() const { return value.getImagI(); }
    double y() const { return value.getImagJ(); }
    double z() const { return value.getImagK(); }
    double speed() const { return value.module(); }
    const ::mathematics::Quaternion<double> &asQuaternion() const { return value; }

    double dot(const Velocity &other) const { return value.dot(other.value); }
    Velocity cross(const Velocity &other) const { return Velocity(value.cross(other.value)); }

    Velocity operator+(const Velocity &other) const { return Velocity(value + other.value); }
    Velocity operator-(const Velocity &other) const { return Velocity(value - other.value); }
    Velocity operator*(double scalar) const { return Velocity(value * scalar); }
    Velocity operator/(double scalar) const { return Velocity(value / scalar); }
    Velocity operator-() const { return Velocity(-value); }

    Velocity &operator+=(const Velocity &other)
    {
        value = value + other.value;
        return *this;
    }
    Velocity &operator-=(const Velocity &other)
    {
        value = value - other.value;
        return *this;
    }
    Velocity &operator*=(double scalar)
    {
        value = value * scalar;
        return *this;
    }
    Velocity &operator/=(double scalar)
    {
        value = value / scalar;
        return *this;
    }

    bool operator==(const Velocity &other) const { return value == other.value; }
    bool operator!=(const Velocity &other) const { return !(*this == other); }

    friend Velocity operator*(double scalar, const Velocity &velocity) { return velocity * scalar; }
    friend Velocity operator^(const Velocity &left, const Velocity &right) { return left.cross(right); }
    friend std::ostream &operator<<(std::ostream &output, const Velocity &velocity) { return output << "v(" << velocity.x() << ", " << velocity.y() << ", " << velocity.z() << ") m/s"; }
};

class Acceleration
{
  private:
    ::mathematics::Quaternion<double> value;

  public:
    Acceleration() : value(0.0, 0.0, 0.0, 0.0) {}
    Acceleration(double x, double y, double z) : value(0.0, x, y, z)
    {
        if (::mathematics::abnormal(x) || ::mathematics::abnormal(y) || ::mathematics::abnormal(z))
        {
            throw std::domain_error("Acceleration: NaN or Inf input");
        }
    }

    explicit Acceleration(const ::mathematics::Quaternion<double> &quaternionValue) : value(quaternionValue) {}

    Acceleration(const Acceleration &other) : value(other.value) {}
    Acceleration(Acceleration &&other) noexcept : value(std::move(other.value)) {}

    Acceleration &operator=(const Acceleration &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Acceleration &operator=(Acceleration &&other) noexcept
    {
        if (this != &other)
        {
            value = std::move(other.value);
        }
        return *this;
    }

    double x() const { return value.getImagI(); }
    double y() const { return value.getImagJ(); }
    double z() const { return value.getImagK(); }
    double magnitude() const { return value.module(); }
    const ::mathematics::Quaternion<double> &asQuaternion() const { return value; }

    double dot(const Acceleration &other) const { return value.dot(other.value); }

    Acceleration operator+(const Acceleration &other) const { return Acceleration(value + other.value); }
    Acceleration operator-(const Acceleration &other) const { return Acceleration(value - other.value); }
    Acceleration operator*(double scalar) const { return Acceleration(value * scalar); }
    Acceleration operator/(double scalar) const { return Acceleration(value / scalar); }
    Acceleration operator-() const { return Acceleration(-value); }

    Acceleration &operator+=(const Acceleration &other)
    {
        value = value + other.value;
        return *this;
    }
    Acceleration &operator-=(const Acceleration &other)
    {
        value = value - other.value;
        return *this;
    }
    Acceleration &operator*=(double scalar)
    {
        value = value * scalar;
        return *this;
    }
    Acceleration &operator/=(double scalar)
    {
        value = value / scalar;
        return *this;
    }

    bool operator==(const Acceleration &other) const { return value == other.value; }
    bool operator!=(const Acceleration &other) const { return !(*this == other); }

    friend Acceleration operator*(double scalar, const Acceleration &acceleration) { return acceleration * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Acceleration &acceleration) { return output << "a(" << acceleration.x() << ", " << acceleration.y() << ", " << acceleration.z() << ") m/s^2"; }
};

class Force
{
  private:
    ::mathematics::Quaternion<double> value;

  public:
    Force() : value(0.0, 0.0, 0.0, 0.0) {}
    Force(double x, double y, double z) : value(0.0, x, y, z)
    {
        if (::mathematics::abnormal(x) || ::mathematics::abnormal(y) || ::mathematics::abnormal(z))
        {
            throw std::domain_error("Force: NaN or Inf input");
        }
    }

    explicit Force(const ::mathematics::Quaternion<double> &quaternionValue) : value(quaternionValue) {}

    Force(const Force &other) : value(other.value) {}
    Force(Force &&other) noexcept : value(std::move(other.value)) {}

    Force &operator=(const Force &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Force &operator=(Force &&other) noexcept
    {
        if (this != &other)
        {
            value = std::move(other.value);
        }
        return *this;
    }

    double x() const { return value.getImagI(); }
    double y() const { return value.getImagJ(); }
    double z() const { return value.getImagK(); }
    double magnitude() const { return value.module(); }
    const ::mathematics::Quaternion<double> &asQuaternion() const { return value; }

    double dot(const Force &other) const { return value.dot(other.value); }
    Force cross(const Force &other) const { return Force(value.cross(other.value)); }

    Force operator+(const Force &other) const { return Force(value + other.value); }
    Force operator-(const Force &other) const { return Force(value - other.value); }
    Force operator*(double scalar) const { return Force(value * scalar); }
    Force operator/(double scalar) const { return Force(value / scalar); }
    Force operator-() const { return Force(-value); }

    Force &operator+=(const Force &other)
    {
        value = value + other.value;
        return *this;
    }
    Force &operator-=(const Force &other)
    {
        value = value - other.value;
        return *this;
    }
    Force &operator*=(double scalar)
    {
        value = value * scalar;
        return *this;
    }
    Force &operator/=(double scalar)
    {
        value = value / scalar;
        return *this;
    }

    bool operator==(const Force &other) const { return value == other.value; }
    bool operator!=(const Force &other) const { return !(*this == other); }

    friend Force operator*(double scalar, const Force &force) { return force * scalar; }
    friend Force operator^(const Force &left, const Force &right) { return left.cross(right); }
    friend std::ostream &operator<<(std::ostream &output, const Force &force) { return output << "F(" << force.x() << ", " << force.y() << ", " << force.z() << ") N"; }
};

class Momentum
{
  private:
    ::mathematics::Quaternion<double> value;

  public:
    Momentum() : value(0.0, 0.0, 0.0, 0.0) {}
    Momentum(double x, double y, double z) : value(0.0, x, y, z)
    {
        if (::mathematics::abnormal(x) || ::mathematics::abnormal(y) || ::mathematics::abnormal(z))
        {
            throw std::domain_error("Momentum: NaN or Inf input");
        }
    }

    explicit Momentum(const ::mathematics::Quaternion<double> &quaternionValue) : value(quaternionValue) {}

    Momentum(const Momentum &other) : value(other.value) {}
    Momentum(Momentum &&other) noexcept : value(std::move(other.value)) {}

    Momentum &operator=(const Momentum &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Momentum &operator=(Momentum &&other) noexcept
    {
        if (this != &other)
        {
            value = std::move(other.value);
        }
        return *this;
    }

    double x() const { return value.getImagI(); }
    double y() const { return value.getImagJ(); }
    double z() const { return value.getImagK(); }
    double magnitude() const { return value.module(); }
    const ::mathematics::Quaternion<double> &asQuaternion() const { return value; }

    double dot(const Momentum &other) const { return value.dot(other.value); }
    Momentum cross(const Momentum &other) const { return Momentum(value.cross(other.value)); }

    Momentum operator+(const Momentum &other) const { return Momentum(value + other.value); }
    Momentum operator-(const Momentum &other) const { return Momentum(value - other.value); }
    Momentum operator*(double scalar) const { return Momentum(value * scalar); }
    Momentum operator/(double scalar) const { return Momentum(value / scalar); }
    Momentum operator-() const { return Momentum(-value); }

    Momentum &operator+=(const Momentum &other)
    {
        value = value + other.value;
        return *this;
    }
    Momentum &operator-=(const Momentum &other)
    {
        value = value - other.value;
        return *this;
    }
    Momentum &operator*=(double scalar)
    {
        value = value * scalar;
        return *this;
    }
    Momentum &operator/=(double scalar)
    {
        value = value / scalar;
        return *this;
    }

    bool operator==(const Momentum &other) const { return value == other.value; }
    bool operator!=(const Momentum &other) const { return !(*this == other); }

    friend Momentum operator*(double scalar, const Momentum &momentum) { return momentum * scalar; }
    friend Momentum operator^(const Momentum &left, const Momentum &right) { return left.cross(right); }
    friend std::ostream &operator<<(std::ostream &output, const Momentum &momentum) { return output << "p(" << momentum.x() << ", " << momentum.y() << ", " << momentum.z() << ") kg*m/s"; }
};

class Torque
{
  private:
    ::mathematics::Quaternion<double> value;

  public:
    Torque() : value(0.0, 0.0, 0.0, 0.0) {}
    Torque(double x, double y, double z) : value(0.0, x, y, z)
    {
        if (::mathematics::abnormal(x) || ::mathematics::abnormal(y) || ::mathematics::abnormal(z))
        {
            throw std::domain_error("Torque: NaN or Inf input");
        }
    }

    explicit Torque(const ::mathematics::Quaternion<double> &quaternionValue) : value(quaternionValue) {}

    Torque(const Torque &other) : value(other.value) {}
    Torque(Torque &&other) noexcept : value(std::move(other.value)) {}

    Torque &operator=(const Torque &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Torque &operator=(Torque &&other) noexcept
    {
        if (this != &other)
        {
            value = std::move(other.value);
        }
        return *this;
    }

    double x() const { return value.getImagI(); }
    double y() const { return value.getImagJ(); }
    double z() const { return value.getImagK(); }
    double magnitude() const { return value.module(); }
    const ::mathematics::Quaternion<double> &asQuaternion() const { return value; }

    double dot(const Torque &other) const { return value.dot(other.value); }
    Torque cross(const Torque &other) const { return Torque(value.cross(other.value)); }

    Torque operator+(const Torque &other) const { return Torque(value + other.value); }
    Torque operator-(const Torque &other) const { return Torque(value - other.value); }
    Torque operator*(double scalar) const { return Torque(value * scalar); }
    Torque operator/(double scalar) const { return Torque(value / scalar); }
    Torque operator-() const { return Torque(-value); }

    Torque &operator+=(const Torque &other)
    {
        value = value + other.value;
        return *this;
    }
    Torque &operator-=(const Torque &other)
    {
        value = value - other.value;
        return *this;
    }
    Torque &operator*=(double scalar)
    {
        value = value * scalar;
        return *this;
    }
    Torque &operator/=(double scalar)
    {
        value = value / scalar;
        return *this;
    }

    bool operator==(const Torque &other) const { return value == other.value; }
    bool operator!=(const Torque &other) const { return !(*this == other); }

    friend Torque operator*(double scalar, const Torque &torque) { return torque * scalar; }
    friend Torque operator^(const Torque &left, const Torque &right) { return left.cross(right); }
    friend std::ostream &operator<<(std::ostream &output, const Torque &torque) { return output << "tau(" << torque.x() << ", " << torque.y() << ", " << torque.z() << ") N*m"; }
};

class Temperature
{
  private:
    double value;

  public:
    Temperature() : value(0.0) {}

    explicit Temperature(double kelvin) : value(kelvin)
    {
        if (::mathematics::abnormal(value))
        {
            throw std::domain_error("Temperature: NaN or Inf input");
        }
        if (value < 0.0)
        {
            throw std::domain_error("Temperature: negative temperature");
        }
        value = ::mathematics::exact(value);
    }

    Temperature(const Temperature &other) : value(other.value) {}
    Temperature(Temperature &&other) noexcept : value(other.value) { other.value = 0.0; }

    Temperature &operator=(const Temperature &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    Temperature &operator=(Temperature &&other) noexcept
    {
        if (this != &other)
        {
            value = other.value;
            other.value = 0.0;
        }
        return *this;
    }

    double kelvin() const { return value; }
    double celsius() const { return ::mathematics::subtract(value, 273.15); }
    double fahrenheit() const { return ::mathematics::add(::mathematics::multiply(::mathematics::subtract(value, 273.15), 1.8), 32.0); }

    Temperature operator+(const Temperature &other) const { return Temperature(::mathematics::add(value, other.value)); }
    Temperature operator-(const Temperature &other) const { return Temperature(::mathematics::subtract(value, other.value)); }
    Temperature operator*(double scalar) const { return Temperature(::mathematics::multiply(value, scalar)); }
    Temperature operator/(double scalar) const { return Temperature(::mathematics::divide(value, scalar)); }
    double operator/(const Temperature &other) const { return ::mathematics::divide(value, other.value); }

    Temperature &operator+=(const Temperature &other)
    {
        value = ::mathematics::add(value, other.value);
        return *this;
    }
    Temperature &operator-=(const Temperature &other)
    {
        value = ::mathematics::subtract(value, other.value);
        return *this;
    }
    Temperature &operator*=(double scalar)
    {
        value = ::mathematics::multiply(value, scalar);
        return *this;
    }
    Temperature &operator/=(double scalar)
    {
        value = ::mathematics::divide(value, scalar);
        return *this;
    }

    bool operator==(const Temperature &other) const { return ::mathematics::absolute(::mathematics::subtract(value, other.value)) < ::mathematics::tolerance; }
    bool operator!=(const Temperature &other) const { return !(*this == other); }
    bool operator<(const Temperature &other) const { return value < other.value; }
    bool operator<=(const Temperature &other) const { return value <= other.value; }
    bool operator>(const Temperature &other) const { return value > other.value; }
    bool operator>=(const Temperature &other) const { return value >= other.value; }

    friend Temperature operator*(double scalar, const Temperature &temperature) { return temperature * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const Temperature &temperature) { return output << temperature.value << " K"; }
};

class HeatCapacity
{
  private:
    double value;

  public:
    HeatCapacity() : value(0.0) {}

    explicit HeatCapacity(double joulesPerKelvin) : value(joulesPerKelvin)
    {
        if (::mathematics::abnormal(value))
        {
            throw std::domain_error("HeatCapacity: NaN or Inf input");
        }
        if (value < 0.0)
        {
            throw std::domain_error("HeatCapacity: negative heat capacity");
        }
        value = ::mathematics::exact(value);
    }

    HeatCapacity(const HeatCapacity &other) : value(other.value) {}
    HeatCapacity(HeatCapacity &&other) noexcept : value(other.value) { other.value = 0.0; }

    HeatCapacity &operator=(const HeatCapacity &other)
    {
        if (this != &other)
        {
            value = other.value;
        }
        return *this;
    }
    HeatCapacity &operator=(HeatCapacity &&other) noexcept
    {
        if (this != &other)
        {
            value = other.value;
            other.value = 0.0;
        }
        return *this;
    }

    double joulesPerKelvin() const { return value; }

    HeatCapacity operator+(const HeatCapacity &other) const { return HeatCapacity(::mathematics::add(value, other.value)); }
    HeatCapacity operator*(double scalar) const { return HeatCapacity(::mathematics::multiply(value, scalar)); }
    HeatCapacity operator/(double scalar) const { return HeatCapacity(::mathematics::divide(value, scalar)); }

    HeatCapacity &operator+=(const HeatCapacity &other)
    {
        value = ::mathematics::add(value, other.value);
        return *this;
    }
    HeatCapacity &operator*=(double scalar)
    {
        value = ::mathematics::multiply(value, scalar);
        return *this;
    }
    HeatCapacity &operator/=(double scalar)
    {
        value = ::mathematics::divide(value, scalar);
        return *this;
    }

    bool operator==(const HeatCapacity &other) const { return ::mathematics::absolute(::mathematics::subtract(value, other.value)) < ::mathematics::tolerance; }
    bool operator!=(const HeatCapacity &other) const { return !(*this == other); }
    bool operator<(const HeatCapacity &other) const { return value < other.value; }
    bool operator<=(const HeatCapacity &other) const { return value <= other.value; }
    bool operator>(const HeatCapacity &other) const { return value > other.value; }
    bool operator>=(const HeatCapacity &other) const { return value >= other.value; }

    friend HeatCapacity operator*(double scalar, const HeatCapacity &heatCapacity) { return heatCapacity * scalar; }
    friend std::ostream &operator<<(std::ostream &output, const HeatCapacity &heatCapacity) { return output << heatCapacity.value << " J/K"; }
};

inline Frequency frequencyFromHertz(double hertzValue)
{
    if (::mathematics::abnormal(hertzValue))
    {
        throw std::domain_error("frequencyFromHertz: NaN or Inf input");
    }
    if (hertzValue < 0.0)
    {
        throw std::domain_error("frequencyFromHertz: negative frequency");
    }
    return Frequency(0.0, ::mathematics::multiply(2.0 * ::mathematics::pi, hertzValue));
}

inline Frequency frequencyFromAngularVelocity(double angularFrequency)
{
    if (::mathematics::abnormal(angularFrequency))
    {
        throw std::domain_error("frequencyFromAngularVelocity: NaN or Inf input");
    }
    if (angularFrequency < 0.0)
    {
        throw std::domain_error("frequencyFromAngularVelocity: negative frequency");
    }
    return Frequency(0.0, angularFrequency);
}

inline Frequency frequencyFromLaplace(double sigma, double omega)
{
    if (::mathematics::abnormal(sigma) || ::mathematics::abnormal(omega))
    {
        throw std::domain_error("frequencyFromLaplace: NaN or Inf input");
    }
    return Frequency(sigma, omega);
}

inline Frequency frequencyFromDampedNatural(double naturalHertz, double dampingRatio)
{
    if (::mathematics::abnormal(naturalHertz) || ::mathematics::abnormal(dampingRatio))
    {
        throw std::domain_error("frequencyFromDampedNatural: NaN or Inf input");
    }
    if (naturalHertz < 0.0)
    {
        throw std::domain_error("frequencyFromDampedNatural: negative natural frequency");
    }
    if (dampingRatio < 0.0 || dampingRatio > 1.0)
    {
        throw std::domain_error("frequencyFromDampedNatural: damping ratio out of [0,1]");
    }
    double angularNatural = ::mathematics::multiply(2.0 * ::mathematics::pi, naturalHertz);
    double angularDamped = ::mathematics::multiply(angularNatural, ::mathematics::squareRoot(::mathematics::subtract(1.0, ::mathematics::multiply(dampingRatio, dampingRatio))));
    double decayRate = ::mathematics::multiply(-dampingRatio, angularNatural);
    return Frequency(decayRate, angularDamped);
}

inline Period periodFromDuration(double durationSeconds) { return Period(durationSeconds, ::mathematics::divide(1.0, durationSeconds)); }

inline Period periodFromHarmonic(double frequencyHertz, int cycles)
{
    double ds = ::mathematics::divide(static_cast<double>(cycles), frequencyHertz);
    double bh = ::mathematics::divide(frequencyHertz, static_cast<double>(cycles));
    return Period(ds, bh);
}

inline Phasor phasorFromCartesian(double real, double imaginary)
{
    if (::mathematics::abnormal(real) || ::mathematics::abnormal(imaginary))
    {
        throw std::domain_error("phasorFromCartesian: NaN or Inf input");
    }
    return Phasor(::mathematics::Complex<double>(real, imaginary));
}

inline Phasor phasorZero() { return Phasor(); }
inline Phasor phasorUnity() { return Phasor(1.0, 0.0); }

inline Voltage voltageFromCartesian(double real, double imaginary) { return Voltage(phasorFromCartesian(real, imaginary)); }
inline Voltage voltageZero() { return Voltage(phasorZero()); }

inline Current currentFromCartesian(double real, double imaginary) { return Current(phasorFromCartesian(real, imaginary)); }
inline Current currentZero() { return Current(phasorZero()); }

inline Noise noiseFromThermal(double temperatureKelvin, double resistanceOhms, double frequencyHertz = 0.0)
{
    const double boltzmannConstant = 1.380649e-23;
    const double planckConstant = 6.62607015e-34;
    double classical = ::mathematics::multiply(::mathematics::multiply(::mathematics::multiply(4.0, boltzmannConstant), temperatureKelvin), resistanceOhms);
    if (frequencyHertz < 1e3)
    {
        return Noise(classical, temperatureKelvin, frequencyHertz);
    }
    double ratio = ::mathematics::divide(::mathematics::multiply(planckConstant, frequencyHertz), ::mathematics::multiply(boltzmannConstant, temperatureKelvin));
    double quantumFactor = ratio > ::mathematics::tolerance ? ::mathematics::divide(ratio, ::mathematics::subtract(::mathematics::exponential(ratio), 1.0)) : 1.0;
    double psdValue = ::mathematics::multiply(classical, quantumFactor);
    if (::mathematics::abnormal(psdValue))
    {
        throw std::overflow_error("noiseFromThermal: PSD overflow");
    }
    return Noise(psdValue, temperatureKelvin, frequencyHertz);
}

inline Noise noiseFromShot(double currentAmperes)
{
    const double elementaryCharge = 1.602176634e-19;
    return Noise(::mathematics::multiply(::mathematics::multiply(2.0, elementaryCharge), ::mathematics::absolute(currentAmperes)));
}

inline Temperature temperatureFromCelsius(double celsiusValue) { return Temperature(::mathematics::add(celsiusValue, 273.15)); }
inline Temperature temperatureFromFahrenheit(double fahrenheitValue) { return Temperature(::mathematics::divide(::mathematics::add(fahrenheitValue, 459.67), 1.8)); }

inline Velocity velocityFromDisplacementTime(const Displacement &displacement, const Time &time)
{
    double seconds = time.seconds();
    if (seconds < ::mathematics::tolerance)
    {
        throw std::domain_error("velocityFromDisplacementTime: zero time");
    }
    double scale = ::mathematics::divide(1.0, seconds);
    return Velocity(::mathematics::multiply(displacement.x(), scale), ::mathematics::multiply(displacement.y(), scale), ::mathematics::multiply(displacement.z(), scale));
}

inline Acceleration accelerationFromVelocityTime(const Velocity &velocity, const Time &time)
{
    double seconds = time.seconds();
    if (seconds < ::mathematics::tolerance)
    {
        throw std::domain_error("accelerationFromVelocityTime: zero time");
    }
    double scale = ::mathematics::divide(1.0, seconds);
    return Acceleration(::mathematics::multiply(velocity.x(), scale), ::mathematics::multiply(velocity.y(), scale), ::mathematics::multiply(velocity.z(), scale));
}

inline Force forceFromMassAcceleration(double mass, const Acceleration &acceleration)
{
    return Force(::mathematics::multiply(acceleration.x(), mass), ::mathematics::multiply(acceleration.y(), mass), ::mathematics::multiply(acceleration.z(), mass));
}

inline Momentum momentumFromMassVelocity(double mass, const Velocity &velocity)
{
    return Momentum(::mathematics::multiply(velocity.x(), mass), ::mathematics::multiply(velocity.y(), mass), ::mathematics::multiply(velocity.z(), mass));
}

inline Torque torqueFromDisplacementForce(const Displacement &displacement, const Force &force)
{
    ::mathematics::Quaternion<double> result = displacement.asQuaternion().cross(force.asQuaternion());
    return Torque(result.getImagI(), result.getImagJ(), result.getImagK());
}

inline Capacitance capacitanceFromChargeVoltage(const Charge &charge, const Voltage &voltage)
{
    double voltageAmplitude = voltage.amplitude();
    if (voltageAmplitude < ::mathematics::tolerance)
    {
        throw std::domain_error("capacitanceFromChargeVoltage: zero voltage");
    }
    return Capacitance(::mathematics::divide(charge.coulombs(), voltageAmplitude));
}

inline Inductance inductanceFromFluxCurrent(const Flux &flux, const Current &current)
{
    double currentAmplitude = current.amplitude();
    if (currentAmplitude < ::mathematics::tolerance)
    {
        throw std::domain_error("inductanceFromFluxCurrent: zero current");
    }
    return Inductance(::mathematics::divide(flux.webers(), currentAmplitude));
}

inline Resistance resistanceFromImpedance(const Impedance &impedance) { return Resistance(impedance.getResistance()); }
inline Conductance conductanceFromAdmittance(const Admittance &admittance) { return Conductance(admittance.getConductance()); }
inline Resistance resistanceFromConductance(const Conductance &conductance) { return Resistance(::mathematics::divide(1.0, conductance.siemens())); }
inline Conductance conductanceFromResistance(const Resistance &resistance) { return Conductance(::mathematics::divide(1.0, resistance.ohms())); }

inline Impedance impedanceFromResistance(const Resistance &resistance) { return Impedance(resistance.ohms(), 0.0); }
inline Impedance impedanceFromReactance(double reactanceValue) { return Impedance(0.0, reactanceValue); }
inline Impedance impedanceFromCapacitance(const Capacitance &capacitance, const Frequency &frequency) { return Impedance(0.0, ::mathematics::divide(-1.0, ::mathematics::multiply(frequency.omega(), capacitance.farads()))); }
inline Impedance impedanceFromInductance(const Inductance &inductance, const Frequency &frequency) { return Impedance(0.0, ::mathematics::multiply(frequency.omega(), inductance.henries())); }

inline Admittance admittanceFromConductance(const Conductance &conductance) { return Admittance(conductance.siemens(), 0.0); }

inline Voltage operator*(const Resistance &resistance, const Current &current) { return Voltage(::mathematics::multiply(current.amplitude(), resistance.ohms()), current.phase()); }
inline Voltage operator*(const Current &current, const Resistance &resistance) { return resistance * current; }

inline Current operator*(const Conductance &conductance, const Voltage &voltage) { return Current(::mathematics::multiply(voltage.amplitude(), conductance.siemens()), voltage.phase()); }
inline Current operator*(const Voltage &voltage, const Conductance &conductance) { return conductance * voltage; }

inline Charge operator*(const Capacitance &capacitance, const Voltage &voltage) { return Charge(::mathematics::multiply(capacitance.farads(), voltage.amplitude())); }
inline Charge operator*(const Voltage &voltage, const Capacitance &capacitance) { return capacitance * voltage; }

inline Flux operator*(const Inductance &inductance, const Current &current) { return Flux(::mathematics::multiply(inductance.henries(), current.amplitude())); }
inline Flux operator*(const Current &current, const Inductance &inductance) { return inductance * current; }

inline Power powerFromResistanceCurrent(const Resistance &resistance, const Current &current)
{
    double rmsCurrent = current.rootMeanSquare();
    return Power(::mathematics::multiply(::mathematics::multiply(rmsCurrent, rmsCurrent), resistance.ohms()), 0.0);
}

inline Power powerFromConductanceVoltage(const Conductance &conductance, const Voltage &voltage)
{
    double rmsVoltage = voltage.rootMeanSquare();
    return Power(::mathematics::multiply(::mathematics::multiply(rmsVoltage, rmsVoltage), conductance.siemens()), 0.0);
}

inline Energy energyFromCapacitanceVoltage(const Capacitance &capacitance, const Voltage &voltage)
{
    return Energy(::mathematics::multiply(0.5, ::mathematics::multiply(capacitance.farads(), ::mathematics::multiply(voltage.amplitude(), voltage.amplitude()))), 0.0);
}

inline Energy energyFromInductanceCurrent(const Inductance &inductance, const Current &current)
{
    return Energy(::mathematics::multiply(0.5, ::mathematics::multiply(inductance.henries(), ::mathematics::multiply(current.amplitude(), current.amplitude()))), 0.0);
}

inline Current currentFromCapacitanceFrequencyVoltage(const Capacitance &capacitance, const Frequency &frequency, const Voltage &voltage)
{
    double currentAmplitude = ::mathematics::multiply(::mathematics::multiply(capacitance.farads(), frequency.omega()), voltage.amplitude());
    return Current(currentAmplitude, ::mathematics::add(voltage.phase(), ::mathematics::pi / 2.0));
}

inline Voltage voltageFromInductanceFrequencyCurrent(const Inductance &inductance, const Frequency &frequency, const Current &current)
{
    double voltageAmplitude = ::mathematics::multiply(::mathematics::multiply(inductance.henries(), frequency.omega()), current.amplitude());
    return Voltage(voltageAmplitude, ::mathematics::add(current.phase(), ::mathematics::pi / 2.0));
}

inline Voltage operator*(const Current &current, const Impedance &impedance)
{
    ::mathematics::Complex<double> result = current.phasor().complex() * impedance.complex();
    if (::mathematics::abnormal(result))
    {
        throw std::overflow_error("Current * Impedance overflow");
    }
    return Voltage(Phasor(result));
}

inline Voltage operator*(const Impedance &impedance, const Current &current) { return current * impedance; }

inline Current operator/(const Voltage &voltage, const Impedance &impedance)
{
    ::mathematics::Complex<double> result = voltage.phasor().complex() / impedance.complex();
    if (::mathematics::abnormal(result))
    {
        throw std::overflow_error("Voltage / Impedance overflow");
    }
    return Current(Phasor(result));
}

inline Voltage operator*(const Current &current, const Admittance &admittance)
{
    ::mathematics::Complex<double> result = current.phasor().complex() * admittance.complex();
    if (::mathematics::abnormal(result))
    {
        throw std::overflow_error("Current * Admittance overflow");
    }
    return Voltage(Phasor(result));
}

inline Current operator*(const Voltage &voltage, const Admittance &admittance)
{
    ::mathematics::Complex<double> result = voltage.phasor().complex() * admittance.complex();
    if (::mathematics::abnormal(result))
    {
        throw std::overflow_error("Voltage * Admittance overflow");
    }
    return Current(Phasor(result));
}

inline Impedance operator~(const Admittance &admittance)
{
    ::mathematics::Complex<double> inverseValue = admittance.complex().inverse();
    return Impedance(inverseValue.getReal(), inverseValue.getImag());
}

inline Admittance operator~(const Impedance &impedance)
{
    ::mathematics::Complex<double> inverseValue = impedance.complex().inverse();
    return Admittance(inverseValue.getReal(), inverseValue.getImag());
}

inline Power operator*(const Voltage &voltage, const Current &current)
{
    ::mathematics::Complex<double> complexPower = voltage.phasor().complex() * current.phasor().complex().conjugate();
    if (::mathematics::abnormal(complexPower))
    {
        throw std::overflow_error("Voltage * Current overflow");
    }
    return Power(complexPower.getReal(), complexPower.getImag());
}

inline Energy operator*(const Power &power, const Time &time) { return Energy(::mathematics::multiply(power.active(), time.seconds()), 0.0); }
inline Energy operator*(const Time &time, const Power &power) { return power * time; }

inline Phasor operator*(const Frequency &frequency, const Time &time) { return Phasor(frequency.evolution(time.seconds())); }
inline Phasor operator*(const Time &time, const Frequency &frequency) { return frequency * time; }

inline Voltage operator*(const Noise &noise, const Impedance &impedance)
{
    double rmsVoltage = noise.rootMeanSquareVoltage(impedance.getResistance(), 1.0);
    return Voltage(::mathematics::multiply(rmsVoltage, ::mathematics::squareRoot(2.0)), 0.0);
}

inline Voltage operator*(const Impedance &impedance, const Noise &noise) { return noise * impedance; }

inline Energy operator*(const Voltage &voltage, const Charge &charge) { return Energy(::mathematics::multiply(voltage.amplitude(), charge.coulombs()), 0.0); }
inline Energy operator*(const Charge &charge, const Voltage &voltage) { return voltage * charge; }

inline Voltage operator/(const Flux &flux, const Time &time)
{
    double seconds = time.seconds();
    if (seconds < ::mathematics::tolerance)
    {
        throw std::domain_error("Flux / Time: zero time step");
    }
    return Voltage(::mathematics::negative(::mathematics::divide(flux.webers(), seconds)), 0.0);
}

inline Velocity operator/(const Displacement &displacement, const Time &time) { return velocityFromDisplacementTime(displacement, time); }
inline Acceleration operator/(const Velocity &velocity, const Time &time) { return accelerationFromVelocityTime(velocity, time); }

inline Energy operator*(const Force &force, const Displacement &displacement)
{
    double dotProduct = ::mathematics::add(::mathematics::add(::mathematics::multiply(force.x(), displacement.x()), ::mathematics::multiply(force.y(), displacement.y())), ::mathematics::multiply(force.z(), displacement.z()));
    return Energy(dotProduct, 0.0);
}

inline Power operator*(const Force &force, const Velocity &velocity)
{
    double dotProduct = ::mathematics::add(::mathematics::add(::mathematics::multiply(force.x(), velocity.x()), ::mathematics::multiply(force.y(), velocity.y())), ::mathematics::multiply(force.z(), velocity.z()));
    return Power(dotProduct, 0.0);
}

inline Force operator/(const Momentum &momentum, const Time &time)
{
    double seconds = time.seconds();
    if (seconds < ::mathematics::tolerance)
    {
        throw std::domain_error("Momentum / Time: zero time step");
    }
    double scale = ::mathematics::divide(1.0, seconds);
    return Force(::mathematics::multiply(momentum.x(), scale), ::mathematics::multiply(momentum.y(), scale), ::mathematics::multiply(momentum.z(), scale));
}

inline Torque operator^(const Displacement &displacement, const Force &force) { return torqueFromDisplacementForce(displacement, force); }

inline Impedance operator*(const Capacitance &capacitance, const Frequency &frequency) { return impedanceFromCapacitance(capacitance, frequency); }
inline Impedance operator*(const Frequency &frequency, const Capacitance &capacitance) { return capacitance * frequency; }

inline Impedance operator*(const Inductance &inductance, const Frequency &frequency) { return impedanceFromInductance(inductance, frequency); }
inline Impedance operator*(const Frequency &frequency, const Inductance &inductance) { return inductance * frequency; }

inline Energy operator*(const Temperature &temperature, const HeatCapacity &heatCapacity) { return Energy(::mathematics::multiply(temperature.kelvin(), heatCapacity.joulesPerKelvin()), 0.0); }
inline Energy operator*(const HeatCapacity &heatCapacity, const Temperature &temperature) { return temperature * heatCapacity; }

} // namespace physics

namespace std
{
template <>
struct hash<physics::Time>
{
    size_t operator()(const physics::Time &time) const { return hash<double>{}(time.seconds()); }
};

template <>
struct hash<physics::Frequency>
{
    size_t operator()(const physics::Frequency &frequency) const
    {
        size_t seed = 0;
        seed ^= hash<double>{}(frequency.sigma()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash<double>{}(frequency.omega()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

template <>
struct hash<physics::Period>
{
    size_t operator()(const physics::Period &period) const
    {
        size_t seed = 0;
        seed ^= hash<double>{}(period.getDuration()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash<double>{}(period.getBandwidth()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

template <>
struct hash<physics::Phasor>
{
    size_t operator()(const physics::Phasor &phasor) const
    {
        size_t seed = 0;
        seed ^= hash<double>{}(phasor.real()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash<double>{}(phasor.imaginary()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

template <>
struct hash<physics::Voltage>
{
    size_t operator()(const physics::Voltage &voltage) const { return hash<physics::Phasor>{}(voltage.phasor()); }
};
template <>
struct hash<physics::Current>
{
    size_t operator()(const physics::Current &current) const { return hash<physics::Phasor>{}(current.phasor()); }
};
template <>
struct hash<physics::Resistance>
{
    size_t operator()(const physics::Resistance &resistance) const { return hash<double>{}(resistance.ohms()); }
};
template <>
struct hash<physics::Conductance>
{
    size_t operator()(const physics::Conductance &conductance) const { return hash<double>{}(conductance.siemens()); }
};
template <>
struct hash<physics::Capacitance>
{
    size_t operator()(const physics::Capacitance &capacitance) const { return hash<double>{}(capacitance.farads()); }
};
template <>
struct hash<physics::Inductance>
{
    size_t operator()(const physics::Inductance &inductance) const { return hash<double>{}(inductance.henries()); }
};

template <>
struct hash<physics::Impedance>
{
    size_t operator()(const physics::Impedance &impedance) const
    {
        size_t seed = 0;
        seed ^= hash<double>{}(impedance.getResistance()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash<double>{}(impedance.getReactance()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

template <>
struct hash<physics::Admittance>
{
    size_t operator()(const physics::Admittance &admittance) const
    {
        size_t seed = 0;
        seed ^= hash<double>{}(admittance.getConductance()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash<double>{}(admittance.getSusceptance()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

template <>
struct hash<physics::Power>
{
    size_t operator()(const physics::Power &power) const
    {
        size_t seed = 0;
        seed ^= hash<double>{}(power.active()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash<double>{}(power.reactive()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

template <>
struct hash<physics::Energy>
{
    size_t operator()(const physics::Energy &energy) const
    {
        size_t seed = 0;
        seed ^= hash<double>{}(energy.stored()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash<double>{}(energy.dissipated()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

template <>
struct hash<physics::Entropy>
{
    size_t operator()(const physics::Entropy &entropy) const
    {
        size_t seed = 0;
        seed ^= hash<double>{}(entropy.information()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash<double>{}(entropy.productionRate()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

template <>
struct hash<physics::Noise>
{
    size_t operator()(const physics::Noise &noise) const { return hash<double>{}(noise.density()); }
};
template <>
struct hash<physics::Charge>
{
    size_t operator()(const physics::Charge &charge) const { return hash<double>{}(charge.coulombs()); }
};
template <>
struct hash<physics::Flux>
{
    size_t operator()(const physics::Flux &flux) const { return hash<double>{}(flux.webers()); }
};
template <>
struct hash<physics::Displacement>
{
    size_t operator()(const physics::Displacement &displacement) const { return hash<::mathematics::Quaternion<double>>{}(displacement.asQuaternion()); }
};
template <>
struct hash<physics::Velocity>
{
    size_t operator()(const physics::Velocity &velocity) const { return hash<::mathematics::Quaternion<double>>{}(velocity.asQuaternion()); }
};
template <>
struct hash<physics::Acceleration>
{
    size_t operator()(const physics::Acceleration &acceleration) const { return hash<::mathematics::Quaternion<double>>{}(acceleration.asQuaternion()); }
};
template <>
struct hash<physics::Force>
{
    size_t operator()(const physics::Force &force) const { return hash<::mathematics::Quaternion<double>>{}(force.asQuaternion()); }
};
template <>
struct hash<physics::Momentum>
{
    size_t operator()(const physics::Momentum &momentum) const { return hash<::mathematics::Quaternion<double>>{}(momentum.asQuaternion()); }
};
template <>
struct hash<physics::Torque>
{
    size_t operator()(const physics::Torque &torque) const { return hash<::mathematics::Quaternion<double>>{}(torque.asQuaternion()); }
};
template <>
struct hash<physics::Temperature>
{
    size_t operator()(const physics::Temperature &temperature) const { return hash<double>{}(temperature.kelvin()); }
};
template <>
struct hash<physics::HeatCapacity>
{
    size_t operator()(const physics::HeatCapacity &heatCapacity) const { return hash<double>{}(heatCapacity.joulesPerKelvin()); }
};

template <>
struct less<physics::Time>
{
    bool operator()(const physics::Time &left, const physics::Time &right) const { return left.seconds() < right.seconds(); }
};

template <>
struct less<physics::Frequency>
{
    bool operator()(const physics::Frequency &left, const physics::Frequency &right) const
    {
        if (left.sigma() != right.sigma())
        {
            return left.sigma() < right.sigma();
        }
        return left.omega() < right.omega();
    }
};

template <>
struct less<physics::Period>
{
    bool operator()(const physics::Period &left, const physics::Period &right) const
    {
        if (left.getDuration() != right.getDuration())
        {
            return left.getDuration() < right.getDuration();
        }
        return left.getBandwidth() < right.getBandwidth();
    }
};

template <>
struct less<physics::Phasor>
{
    bool operator()(const physics::Phasor &left, const physics::Phasor &right) const
    {
        if (left.real() != right.real())
        {
            return left.real() < right.real();
        }
        return left.imaginary() < right.imaginary();
    }
};

template <>
struct less<physics::Voltage>
{
    bool operator()(const physics::Voltage &left, const physics::Voltage &right) const
    {
        if (left.amplitude() != right.amplitude())
        {
            return left.amplitude() < right.amplitude();
        }
        return left.phase() < right.phase();
    }
};

template <>
struct less<physics::Current>
{
    bool operator()(const physics::Current &left, const physics::Current &right) const
    {
        if (left.amplitude() != right.amplitude())
        {
            return left.amplitude() < right.amplitude();
        }
        return left.phase() < right.phase();
    }
};

template <>
struct less<physics::Resistance>
{
    bool operator()(const physics::Resistance &left, const physics::Resistance &right) const { return left.ohms() < right.ohms(); }
};
template <>
struct less<physics::Conductance>
{
    bool operator()(const physics::Conductance &left, const physics::Conductance &right) const { return left.siemens() < right.siemens(); }
};
template <>
struct less<physics::Capacitance>
{
    bool operator()(const physics::Capacitance &left, const physics::Capacitance &right) const { return left.farads() < right.farads(); }
};
template <>
struct less<physics::Inductance>
{
    bool operator()(const physics::Inductance &left, const physics::Inductance &right) const { return left.henries() < right.henries(); }
};

template <>
struct less<physics::Impedance>
{
    bool operator()(const physics::Impedance &left, const physics::Impedance &right) const
    {
        if (left.getResistance() != right.getResistance())
        {
            return left.getResistance() < right.getResistance();
        }
        return left.getReactance() < right.getReactance();
    }
};

template <>
struct less<physics::Admittance>
{
    bool operator()(const physics::Admittance &left, const physics::Admittance &right) const
    {
        if (left.getConductance() != right.getConductance())
        {
            return left.getConductance() < right.getConductance();
        }
        return left.getSusceptance() < right.getSusceptance();
    }
};

template <>
struct less<physics::Power>
{
    bool operator()(const physics::Power &left, const physics::Power &right) const
    {
        if (left.active() != right.active())
        {
            return left.active() < right.active();
        }
        return left.reactive() < right.reactive();
    }
};

template <>
struct less<physics::Energy>
{
    bool operator()(const physics::Energy &left, const physics::Energy &right) const
    {
        if (left.stored() != right.stored())
        {
            return left.stored() < right.stored();
        }
        return left.dissipated() < right.dissipated();
    }
};

template <>
struct less<physics::Entropy>
{
    bool operator()(const physics::Entropy &left, const physics::Entropy &right) const
    {
        if (left.information() != right.information())
        {
            return left.information() < right.information();
        }
        return left.productionRate() < right.productionRate();
    }
};

template <>
struct less<physics::Noise>
{
    bool operator()(const physics::Noise &left, const physics::Noise &right) const { return left.density() < right.density(); }
};
template <>
struct less<physics::Charge>
{
    bool operator()(const physics::Charge &left, const physics::Charge &right) const { return left.coulombs() < right.coulombs(); }
};
template <>
struct less<physics::Flux>
{
    bool operator()(const physics::Flux &left, const physics::Flux &right) const { return left.webers() < right.webers(); }
};
template <>
struct less<physics::Displacement>
{
    bool operator()(const physics::Displacement &left, const physics::Displacement &right) const { return left.asQuaternion() < right.asQuaternion(); }
};
template <>
struct less<physics::Velocity>
{
    bool operator()(const physics::Velocity &left, const physics::Velocity &right) const { return left.asQuaternion() < right.asQuaternion(); }
};
template <>
struct less<physics::Acceleration>
{
    bool operator()(const physics::Acceleration &left, const physics::Acceleration &right) const { return left.asQuaternion() < right.asQuaternion(); }
};
template <>
struct less<physics::Force>
{
    bool operator()(const physics::Force &left, const physics::Force &right) const { return left.asQuaternion() < right.asQuaternion(); }
};
template <>
struct less<physics::Momentum>
{
    bool operator()(const physics::Momentum &left, const physics::Momentum &right) const { return left.asQuaternion() < right.asQuaternion(); }
};
template <>
struct less<physics::Torque>
{
    bool operator()(const physics::Torque &left, const physics::Torque &right) const { return left.asQuaternion() < right.asQuaternion(); }
};
template <>
struct less<physics::Temperature>
{
    bool operator()(const physics::Temperature &left, const physics::Temperature &right) const { return left.kelvin() < right.kelvin(); }
};
template <>
struct less<physics::HeatCapacity>
{
    bool operator()(const physics::HeatCapacity &left, const physics::HeatCapacity &right) const { return left.joulesPerKelvin() < right.joulesPerKelvin(); }
};

} // namespace std