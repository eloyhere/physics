#pragma once
#include <vector>
#include <deque>
#include <complex>
#include <functional>
#include "mathematics.h"
#include "quaternion.h"
#include "quantization.h"
#include "wave.h"

namespace physics
{
template <typename T>
class Field
{
    static_assert(is_arithmetic_v<T>,
                  "Field<T>: T must be arithmetic, Complex<N>, std::complex<N>, or Quaternion<N>");

  private:
    Wave<::mathematics::Quaternion<double>> wave;

    ::mathematics::Quaternion<double> to(const T &value) const
    {
        if constexpr (std::is_arithmetic_v<T>)
        {
            return ::mathematics::Quaternion<double>(static_cast<double>(value), 0.0, 0.0, 0.0);
        }
        else if constexpr (std::is_same_v<T, ::mathematics::Complex<double>>)
        {
            return ::mathematics::Quaternion<double>(static_cast<double>(value.getReal()), static_cast<double>(value.getImag()), 0.0, 0.0);
        }
        else if constexpr (std::is_same_v<T, std::complex<double>>)
        {
            return ::mathematics::Quaternion<double>(static_cast<double>(value.real()), static_cast<double>(value.imag()), 0.0, 0.0);
        }
        else if constexpr (std::is_same_v<T, ::mathematics::Quaternion<double>>)
        {
            return value;
        }
        else
        {
            return ::mathematics::Quaternion<double>(0.0, 0.0, 0.0, 0.0);
        }
    }

    T from(const ::mathematics::Quaternion<double> &quaternion) const
    {
        if constexpr (std::is_arithmetic_v<T>)
        {
            return static_cast<T>(quaternion.module());
        }
        else if constexpr (std::is_same_v<T, ::mathematics::Complex<double>>)
        {
            return ::mathematics::Complex<double>(quaternion.getReal(), quaternion.getImagI());
        }
        else if constexpr (std::is_same_v<T, std::complex<double>>)
        {
            return std::complex<double>(quaternion.getReal(), quaternion.getImagI());
        }
        else if constexpr (std::is_same_v<T, ::mathematics::Quaternion<double>>)
        {
            return quaternion;
        }
        else
        {
            return T();
        }
    }

    size_t indexAt(const Time &time) const
    {
        double timeSeconds = time.seconds();
        double sampleRate = wave.getSampleRate();
        double index = ::mathematics::multiply(timeSeconds, sampleRate);
        return static_cast<size_t>(::mathematics::limit(index, 0.0, static_cast<double>(wave.size() - 1)));
    }

  public:
    Field(double sampleRate, size_t size)
        : wave(sampleRate)
    {
        for (size_t i = 0; i < size; ++i)
        {
            wave.append(::mathematics::Quaternion<double>(0.0, 0.0, 0.0, 0.0));
        }
    }

    Field(double sampleRate, size_t size, const T &uniformValue)
        : wave(sampleRate)
    {
        auto q = to(uniformValue);
        for (size_t i = 0; i < size; ++i)
        {
            wave.append(q);
        }
    }

    T perturb(const T &particle, Time timePosition, Time timeRange, Time timeDuration)
    {
        ::mathematics::Quaternion<double> particleQ = to(particle);
        ::mathematics::Quaternion<double> scatteredQ = wave.perturb(particleQ, timePosition, timeRange, timeDuration);
        return from(scatteredQ);
    }

    Wave<T> perturb(const Wave<T> &incoming, Time timePosition, Time timeRange, Time timeDuration)
    {
        Wave<::mathematics::Quaternion<double>> incomingQ(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            incomingQ.append(to(incoming[i]));
        }
        Wave<::mathematics::Quaternion<double>> scatteredQ = wave.perturb(incomingQ, timePosition, timeRange, timeDuration);
        Wave<T> scattered(incoming.getSampleRate());
        for (size_t i = 0; i < scatteredQ.size(); ++i)
        {
            scattered.append(from(scatteredQ[i]));
        }
        return scattered;
    }

    Field &spread(double intensity, Time position, Time range, Time duration)
    {
        wave.spread(intensity, position, range, duration);
        return *this;
    }

    template <typename U = T>
    std::enable_if_t<is_arithmetic_v<U> && !std::is_fundamental_v<U>, Field &>
    spread(const T &complexIntensity, Time position, Time range, Time duration)
    {
        ::mathematics::Quaternion<double> intensityQ = to(complexIntensity);
        wave.spread(intensityQ, position, range, duration);
        return *this;
    }

    Field &spread(const Wave<T> &incoming, Time duration)
    {
        Wave<::mathematics::Quaternion<double>> incomingQ(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            incomingQ.append(to(incoming[i]));
        }
        wave.spread(incomingQ, duration);
        return *this;
    }

    T displacement(const Time &time) const
    {
        size_t index = indexAt(time);
        return from(wave[index]);
    }

    T velocity(const Time &time) const
    {
        size_t index = indexAt(time);
        if (index < 1 || index >= wave.size() - 1)
        {
            return T();
        }
        double deltaTime = ::mathematics::divide(1.0, wave.getSampleRate());
        ::mathematics::Quaternion<double> forward = wave[index + 1];
        ::mathematics::Quaternion<double> backward = wave[index - 1];
        double denominator = ::mathematics::multiply(2.0, deltaTime);
        return from(::mathematics::Quaternion<double>(
            0.0,
            ::mathematics::divide(::mathematics::subtract(forward.getImagI(), backward.getImagI()), denominator),
            ::mathematics::divide(::mathematics::subtract(forward.getImagJ(), backward.getImagJ()), denominator),
            ::mathematics::divide(::mathematics::subtract(forward.getImagK(), backward.getImagK()), denominator)));
    }

    T accelerate(const Time &time) const
    {
        size_t index = indexAt(time);
        if (index < 2 || index >= wave.size() - 2)
        {
            return T();
        }
        double deltaTime = ::mathematics::divide(1.0, wave.getSampleRate());
        double deltaTimeSquared = ::mathematics::multiply(deltaTime, deltaTime);

        ::mathematics::Quaternion<double> current = wave[index];
        ::mathematics::Quaternion<double> forward1 = wave[index + 1];
        ::mathematics::Quaternion<double> backward1 = wave[index - 1];
        ::mathematics::Quaternion<double> forward2 = wave[index + 2];
        ::mathematics::Quaternion<double> backward2 = wave[index - 2];

        double accelerationScalar = ::mathematics::divide(
            ::mathematics::subtract(
                ::mathematics::add(
                    ::mathematics::add(::mathematics::multiply(forward1.getReal(), 16.0), ::mathematics::multiply(backward1.getReal(), 16.0)),
                    ::mathematics::multiply(current.getReal(), -30.0)),
                ::mathematics::add(forward2.getReal(), backward2.getReal())),
            ::mathematics::multiply(12.0, deltaTimeSquared));

        double accelerationI = ::mathematics::divide(
            ::mathematics::subtract(
                ::mathematics::add(
                    ::mathematics::add(::mathematics::multiply(forward1.getImagI(), 16.0), ::mathematics::multiply(backward1.getImagI(), 16.0)),
                    ::mathematics::multiply(current.getImagI(), -30.0)),
                ::mathematics::add(forward2.getImagI(), backward2.getImagI())),
            ::mathematics::multiply(12.0, deltaTimeSquared));

        double accelerationJ = ::mathematics::divide(
            ::mathematics::subtract(
                ::mathematics::add(
                    ::mathematics::add(::mathematics::multiply(forward1.getImagJ(), 16.0), ::mathematics::multiply(backward1.getImagJ(), 16.0)),
                    ::mathematics::multiply(current.getImagJ(), -30.0)),
                ::mathematics::add(forward2.getImagJ(), backward2.getImagJ())),
            ::mathematics::multiply(12.0, deltaTimeSquared));

        double accelerationK = ::mathematics::divide(
            ::mathematics::subtract(
                ::mathematics::add(
                    ::mathematics::add(::mathematics::multiply(forward1.getImagK(), 16.0), ::mathematics::multiply(backward1.getImagK(), 16.0)),
                    ::mathematics::multiply(current.getImagK(), -30.0)),
                ::mathematics::add(forward2.getImagK(), backward2.getImagK())),
            ::mathematics::multiply(12.0, deltaTimeSquared));

        return from(::mathematics::Quaternion<double>(accelerationScalar, accelerationI, accelerationJ, accelerationK));
    }

    Wave<T> extractWave(const Time &duration, double sampleRate) const
    {
        Wave<T> result(sampleRate);
        size_t samples = static_cast<size_t>(::mathematics::multiply(duration.seconds(), sampleRate));
        double deltaTime = ::mathematics::divide(1.0, sampleRate);
        for (size_t k = 0; k < samples; ++k)
        {
            Time time(::mathematics::multiply(static_cast<double>(k), deltaTime));
            result.append(displacement(time));
        }
        return result;
    }

    Momentum momentum(const Time &startTime, const Time &endTime) const
    {
        double t1 = startTime.seconds();
        double t2 = endTime.seconds();
        if (t1 > t2)
        {
            throw std::domain_error("momentum: start time must be <= end time");
        }
        double deltaTime = 0.001;
        int steps = static_cast<int>(::mathematics::divide(::mathematics::subtract(t2, t1), deltaTime));
        if (steps < 1)
        {
            steps = 1;
            deltaTime = ::mathematics::subtract(t2, t1);
        }

        ::mathematics::Quaternion<double> sum(0.0, 0.0, 0.0, 0.0);
        for (int k = 0; k < steps; ++k)
        {
            double t = ::mathematics::add(t1, ::mathematics::multiply(static_cast<double>(k), deltaTime));
            Time time(t);
            auto velocityQuaternion = to(velocity(time));
            sum = sum + velocityQuaternion * deltaTime;
        }
        return Momentum(sum.getImagI(), sum.getImagJ(), sum.getImagK());
    }

    Energy energy() const
    {
        double kineticEnergy = 0.0;
        double potentialEnergy = 0.0;
        double deltaTime = ::mathematics::divide(1.0, wave.getSampleRate());

        for (size_t i = 1; i < wave.size() - 1; ++i)
        {
            ::mathematics::Quaternion<double> current = wave[i];
            ::mathematics::Quaternion<double> forward = wave[i + 1];
            ::mathematics::Quaternion<double> backward = wave[i - 1];

            double velocityDifference = ::mathematics::subtract(forward.getReal(), backward.getReal());
            double velocitySquared = ::mathematics::divide(
                ::mathematics::multiply(velocityDifference, velocityDifference),
                ::mathematics::multiply(4.0, ::mathematics::multiply(deltaTime, deltaTime)));

            kineticEnergy = ::mathematics::add(kineticEnergy, ::mathematics::multiply(0.5, velocitySquared));
            potentialEnergy = ::mathematics::add(potentialEnergy, ::mathematics::multiply(0.5, current.dot(current)));
        }
        double total = ::mathematics::add(kineticEnergy, potentialEnergy);
        return Energy(total, 0.0);
    }

    Entropy entropy() const
    {
        double totalPower = 0.0;
        for (size_t i = 0; i < wave.size(); ++i)
        {
            totalPower = ::mathematics::add(totalPower, wave[i].dot(wave[i]));
        }
        if (totalPower < ::mathematics::tolerance)
        {
            return Entropy(0.0, 0.0);
        }

        double information = 0.0;
        for (size_t i = 0; i < wave.size(); ++i)
        {
            double power = wave[i].dot(wave[i]);
            double powerFraction = ::mathematics::divide(power, totalPower);
            if (powerFraction > ::mathematics::tolerance)
            {
                information = ::mathematics::subtract(information, ::mathematics::multiply(powerFraction, ::mathematics::logarithm(powerFraction)));
            }
        }
        return Entropy(information, 0.0);
    }

    Field operator+(const Field &other) const
    {
        if (wave.size() != other.wave.size() || ::mathematics::absolute(::mathematics::subtract(wave.getSampleRate(), other.wave.getSampleRate())) > ::mathematics::tolerance)
        {
            throw std::invalid_argument("Field+: size or sample rate mismatch");
        }
        Field result(wave.getSampleRate(), wave.size());
        for (size_t i = 0; i < wave.size(); ++i)
        {
            result.wave[i] = wave[i] + other.wave[i];
        }
        return result;
    }

    Field operator-(const Field &other) const
    {
        if (wave.size() != other.wave.size() || ::mathematics::absolute(::mathematics::subtract(wave.getSampleRate(), other.wave.getSampleRate())) > ::mathematics::tolerance)
        {
            throw std::invalid_argument("Field-: size or sample rate mismatch");
        }
        Field result(wave.getSampleRate(), wave.size());
        for (size_t i = 0; i < wave.size(); ++i)
        {
            result.wave[i] = wave[i] - other.wave[i];
        }
        return result;
    }

    Field operator*(double scalar) const
    {
        Field result(wave.getSampleRate(), wave.size());
        for (size_t i = 0; i < wave.size(); ++i)
        {
            result.wave[i] = wave[i] * scalar;
        }
        return result;
    }

    Field operator/(double scalar) const
    {
        if (::mathematics::absolute(scalar) < ::mathematics::tolerance)
        {
            throw std::domain_error("Field/: division by zero");
        }
        Field result(wave.getSampleRate(), wave.size());
        for (size_t i = 0; i < wave.size(); ++i)
        {
            result.wave[i] = wave[i] / scalar;
        }
        return result;
    }

    Field &operator+=(const Field &other)
    {
        if (wave.size() != other.wave.size() || ::mathematics::absolute(::mathematics::subtract(wave.getSampleRate(), other.wave.getSampleRate())) > ::mathematics::tolerance)
        {
            throw std::invalid_argument("Field+=: size or sample rate mismatch");
        }
        for (size_t i = 0; i < wave.size(); ++i)
        {
            wave[i] = wave[i] + other.wave[i];
        }
        return *this;
    }

    Field &operator-=(const Field &other)
    {
        if (wave.size() != other.wave.size() || ::mathematics::absolute(::mathematics::subtract(wave.getSampleRate(), other.wave.getSampleRate())) > ::mathematics::tolerance)
        {
            throw std::invalid_argument("Field-=: size or sample rate mismatch");
        }
        for (size_t i = 0; i < wave.size(); ++i)
        {
            wave[i] = wave[i] - other.wave[i];
        }
        return *this;
    }

    Field &operator*=(double scalar)
    {
        for (size_t i = 0; i < wave.size(); ++i)
        {
            wave[i] = wave[i] * scalar;
        }
        return *this;
    }

    Field &operator/=(double scalar)
    {
        if (::mathematics::absolute(scalar) < ::mathematics::tolerance)
        {
            throw std::domain_error("Field/=: division by zero");
        }
        for (size_t i = 0; i < wave.size(); ++i)
        {
            wave[i] = wave[i] / scalar;
        }
        return *this;
    }

    bool operator==(const Field &other) const
    {
        if (wave.size() != other.wave.size())
        {
            return false;
        }
        for (size_t i = 0; i < wave.size(); ++i)
        {
            if (!(wave[i] == other.wave[i]))
            {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const Field &other) const { return !(*this == other); }

    bool operator<(const Field &other) const
    {
        double thisEnergy = energy().total();
        double otherEnergy = other.energy().total();
        return thisEnergy < otherEnergy;
    }

    bool operator<=(const Field &other) const { return (*this < other) || (*this == other); }
    bool operator>(const Field &other) const { return !(*this <= other); }
    bool operator>=(const Field &other) const { return !(*this < other); }

    friend Field operator*(double scalar, const Field &field) { return field * scalar; }

    size_t size() const { return wave.size(); }
    double getSampleRate() const { return wave.getSampleRate(); }
    double duration() const { return wave.duration(); }

    const Wave<::mathematics::Quaternion<double>> &getWave() const { return wave; }
    Wave<::mathematics::Quaternion<double>> &getWave() { return wave; }

    const ::mathematics::Quaternion<double> &operator[](size_t index) const { return wave[index]; }
    ::mathematics::Quaternion<double> &operator[](size_t index) { return wave[index]; }
};

template <>
class Field<Voltage>
{
  private:
    Field<::mathematics::Complex<double>> field;

  public:
    Field(double sampleRate, size_t size) : field(sampleRate, size) {}
    Field(double sampleRate, size_t size, const Voltage &value) : field(sampleRate, size, value.phasor().complex()) {}

    Voltage perturb(const Voltage &particle, Time timePosition, Time timeRange, Time timeDuration)
    {
        ::mathematics::Complex<double> scattered = field.perturb(particle.phasor().complex(), timePosition, timeRange, timeDuration);
        return Voltage(Phasor(scattered));
    }

    Wave<Voltage> perturb(const Wave<Voltage> &incoming, Time timePosition, Time timeRange, Time timeDuration)
    {
        Wave<::mathematics::Complex<double>> incomingComplex(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            incomingComplex.append(incoming[i].phasor().complex());
        }
        Wave<::mathematics::Complex<double>> scatteredComplex = field.perturb(incomingComplex, timePosition, timeRange, timeDuration);
        Wave<Voltage> scattered(scatteredComplex.getSampleRate());
        for (size_t i = 0; i < scatteredComplex.size(); ++i)
        {
            scattered.append(Voltage(Phasor(scatteredComplex[i])));
        }
        return scattered;
    }

    Field<Voltage> &spread(double intensity, Time position, Time range, Time duration)
    {
        field.spread(intensity, position, range, duration);
        return *this;
    }

    Field<Voltage> &spread(const Voltage &complexIntensity, Time position, Time range, Time duration)
    {
        field.spread(complexIntensity.phasor().complex(), position, range, duration);
        return *this;
    }

    Field<Voltage> &spread(const Wave<Voltage> &incoming, Time duration)
    {
        Wave<::mathematics::Complex<double>> incomingComplex(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            incomingComplex.append(incoming[i].phasor().complex());
        }
        field.spread(incomingComplex, duration);
        return *this;
    }

    Voltage displacement(const Time &time) const { return Voltage(Phasor(field.displacement(time))); }
    Voltage velocity(const Time &time) const { return Voltage(Phasor(field.velocity(time))); }
    Voltage accelerate(const Time &time) const { return Voltage(Phasor(field.accelerate(time))); }

    Wave<Voltage> extractWave(const Time &duration, double sampleRate) const
    {
        Wave<Voltage> result(sampleRate);
        size_t samples = static_cast<size_t>(::mathematics::multiply(duration.seconds(), sampleRate));
        double deltaTime = ::mathematics::divide(1.0, sampleRate);
        for (size_t k = 0; k < samples; ++k)
        {
            Time time(::mathematics::multiply(static_cast<double>(k), deltaTime));
            result.append(displacement(time));
        }
        return result;
    }

    Momentum momentum(const Time &startTime, const Time &endTime) const { return field.momentum(startTime, endTime); }
    Energy energy() const { return field.energy(); }
    Entropy entropy() const { return field.entropy(); }

    Field<Voltage> operator+(const Field<Voltage> &other) const
    {
        Field<Voltage> result(field.getSampleRate(), field.size());
        result.field = field + other.field;
        return result;
    }
    Field<Voltage> operator-(const Field<Voltage> &other) const
    {
        Field<Voltage> result(field.getSampleRate(), field.size());
        result.field = field - other.field;
        return result;
    }
    Field<Voltage> operator*(double scalar) const
    {
        Field<Voltage> result(field.getSampleRate(), field.size());
        result.field = field * scalar;
        return result;
    }
    Field<Voltage> operator/(double scalar) const
    {
        Field<Voltage> result(field.getSampleRate(), field.size());
        result.field = field / scalar;
        return result;
    }
    Field<Voltage> &operator+=(const Field<Voltage> &other)
    {
        field += other.field;
        return *this;
    }
    Field<Voltage> &operator-=(const Field<Voltage> &other)
    {
        field -= other.field;
        return *this;
    }
    Field<Voltage> &operator*=(double scalar)
    {
        field *= scalar;
        return *this;
    }
    Field<Voltage> &operator/=(double scalar)
    {
        field /= scalar;
        return *this;
    }
    bool operator==(const Field<Voltage> &other) const { return field == other.field; }
    bool operator!=(const Field<Voltage> &other) const { return !(*this == other); }
    bool operator<(const Field<Voltage> &other) const { return field < other.field; }
    bool operator<=(const Field<Voltage> &other) const { return field <= other.field; }
    bool operator>(const Field<Voltage> &other) const { return field > other.field; }
    bool operator>=(const Field<Voltage> &other) const { return field >= other.field; }
    friend Field<Voltage> operator*(double scalar, const Field<Voltage> &f) { return f * scalar; }

    size_t size() const { return field.size(); }
    double getSampleRate() const { return field.getSampleRate(); }
    double duration() const { return field.duration(); }
    const Wave<::mathematics::Quaternion<double>> &getWave() const { return field.getWave(); }
    Wave<::mathematics::Quaternion<double>> &getWave() { return field.getWave(); }
};

template <>
class Field<Current>
{
  private:
    Field<::mathematics::Complex<double>> field;

  public:
    Field(double sampleRate, size_t size) : field(sampleRate, size) {}
    Field(double sampleRate, size_t size, const Current &value) : field(sampleRate, size, value.phasor().complex()) {}

    Current perturb(const Current &particle, Time timePosition, Time timeRange, Time timeDuration)
    {
        ::mathematics::Complex<double> scattered = field.perturb(particle.phasor().complex(), timePosition, timeRange, timeDuration);
        return Current(Phasor(scattered));
    }

    Wave<Current> perturb(const Wave<Current> &incoming, Time timePosition, Time timeRange, Time timeDuration)
    {
        Wave<::mathematics::Complex<double>> incomingComplex(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            incomingComplex.append(incoming[i].phasor().complex());
        }
        Wave<::mathematics::Complex<double>> scatteredComplex = field.perturb(incomingComplex, timePosition, timeRange, timeDuration);
        Wave<Current> scattered(scatteredComplex.getSampleRate());
        for (size_t i = 0; i < scatteredComplex.size(); ++i)
        {
            scattered.append(Current(Phasor(scatteredComplex[i])));
        }
        return scattered;
    }

    Field<Current> &spread(double intensity, Time position, Time range, Time duration)
    {
        field.spread(intensity, position, range, duration);
        return *this;
    }

    Field<Current> &spread(const Current &complexIntensity, Time position, Time range, Time duration)
    {
        field.spread(complexIntensity.phasor().complex(), position, range, duration);
        return *this;
    }

    Field<Current> &spread(const Wave<Current> &incoming, Time duration)
    {
        Wave<::mathematics::Complex<double>> incomingComplex(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            incomingComplex.append(incoming[i].phasor().complex());
        }
        field.spread(incomingComplex, duration);
        return *this;
    }

    Current displacement(const Time &time) const { return Current(Phasor(field.displacement(time))); }
    Current velocity(const Time &time) const { return Current(Phasor(field.velocity(time))); }
    Current accelerate(const Time &time) const { return Current(Phasor(field.accelerate(time))); }

    Wave<Current> extractWave(const Time &duration, double sampleRate) const
    {
        Wave<Current> result(sampleRate);
        size_t samples = static_cast<size_t>(::mathematics::multiply(duration.seconds(), sampleRate));
        double deltaTime = ::mathematics::divide(1.0, sampleRate);
        for (size_t k = 0; k < samples; ++k)
        {
            Time time(::mathematics::multiply(static_cast<double>(k), deltaTime));
            result.append(displacement(time));
        }
        return result;
    }

    Momentum momentum(const Time &startTime, const Time &endTime) const { return field.momentum(startTime, endTime); }
    Energy energy() const { return field.energy(); }
    Entropy entropy() const { return field.entropy(); }

    Field<Current> operator+(const Field<Current> &other) const
    {
        Field<Current> result(field.getSampleRate(), field.size());
        result.field = field + other.field;
        return result;
    }
    Field<Current> operator-(const Field<Current> &other) const
    {
        Field<Current> result(field.getSampleRate(), field.size());
        result.field = field - other.field;
        return result;
    }
    Field<Current> operator*(double scalar) const
    {
        Field<Current> result(field.getSampleRate(), field.size());
        result.field = field * scalar;
        return result;
    }
    Field<Current> operator/(double scalar) const
    {
        Field<Current> result(field.getSampleRate(), field.size());
        result.field = field / scalar;
        return result;
    }
    Field<Current> &operator+=(const Field<Current> &other)
    {
        field += other.field;
        return *this;
    }
    Field<Current> &operator-=(const Field<Current> &other)
    {
        field -= other.field;
        return *this;
    }
    Field<Current> &operator*=(double scalar)
    {
        field *= scalar;
        return *this;
    }
    Field<Current> &operator/=(double scalar)
    {
        field /= scalar;
        return *this;
    }
    bool operator==(const Field<Current> &other) const { return field == other.field; }
    bool operator!=(const Field<Current> &other) const { return !(*this == other); }
    bool operator<(const Field<Current> &other) const { return field < other.field; }
    bool operator<=(const Field<Current> &other) const { return field <= other.field; }
    bool operator>(const Field<Current> &other) const { return field > other.field; }
    bool operator>=(const Field<Current> &other) const { return field >= other.field; }
    friend Field<Current> operator*(double scalar, const Field<Current> &f) { return f * scalar; }

    size_t size() const { return field.size(); }
    double getSampleRate() const { return field.getSampleRate(); }
    double duration() const { return field.duration(); }
    const Wave<::mathematics::Quaternion<double>> &getWave() const { return field.getWave(); }
    Wave<::mathematics::Quaternion<double>> &getWave() { return field.getWave(); }
};

template <>
class Field<Force>
{
  private:
    Field<::mathematics::Quaternion<double>> field;

  public:
    Field(double sampleRate, size_t size) : field(sampleRate, size) {}
    Field(double sampleRate, size_t size, const Force &value) : field(sampleRate, size, ::mathematics::Quaternion<double>(0.0, static_cast<double>(value.x()), static_cast<double>(value.y()), static_cast<double>(value.z()))) {}

    Force perturb(const Force &particle, Time timePosition, Time timeRange, Time timeDuration)
    {
        ::mathematics::Quaternion<double> particleQ(0.0, particle.x(), particle.y(), particle.z());
        ::mathematics::Quaternion<double> scatteredQ = field.perturb(particleQ, timePosition, timeRange, timeDuration);
        return Force(scatteredQ.getImagI(), scatteredQ.getImagJ(), scatteredQ.getImagK());
    }

    Wave<Force> perturb(const Wave<Force> &incoming, Time timePosition, Time timeRange, Time timeDuration)
    {
        Wave<::mathematics::Quaternion<double>> incomingQ(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            Force f = incoming[i];
            incomingQ.append(::mathematics::Quaternion<double>(0.0, f.x(), f.y(), f.z()));
        }
        Wave<::mathematics::Quaternion<double>> scatteredQ = field.perturb(incomingQ, timePosition, timeRange, timeDuration);
        Wave<Force> scattered(scatteredQ.getSampleRate());
        for (size_t i = 0; i < scatteredQ.size(); ++i)
        {
            scattered.append(Force(scatteredQ[i].getImagI(), scatteredQ[i].getImagJ(), scatteredQ[i].getImagK()));
        }
        return scattered;
    }

    Field<Force> &spread(double intensity, Time position, Time range, Time duration)
    {
        field.spread(intensity, position, range, duration);
        return *this;
    }

    Field<Force> &spread(const Force &complexIntensity, Time position, Time range, Time duration)
    {
        ::mathematics::Quaternion<double> intensityQ(0.0, complexIntensity.x(), complexIntensity.y(), complexIntensity.z());
        field.spread(intensityQ, position, range, duration);
        return *this;
    }

    Field<Force> &spread(const Wave<Force> &incoming, Time duration)
    {
        Wave<::mathematics::Quaternion<double>> incomingQ(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            Force f = incoming[i];
            incomingQ.append(::mathematics::Quaternion<double>(0.0, f.x(), f.y(), f.z()));
        }
        field.spread(incomingQ, duration);
        return *this;
    }

    Force displacement(const Time &time) const
    {
        auto q = field.displacement(time);
        return Force(q.getImagI(), q.getImagJ(), q.getImagK());
    }
    Force velocity(const Time &time) const
    {
        auto q = field.velocity(time);
        return Force(q.getImagI(), q.getImagJ(), q.getImagK());
    }
    Force accelerate(const Time &time) const
    {
        auto q = field.accelerate(time);
        return Force(q.getImagI(), q.getImagJ(), q.getImagK());
    }

    Wave<Force> extractWave(const Time &duration, double sampleRate) const
    {
        Wave<Force> result(sampleRate);
        size_t samples = static_cast<size_t>(::mathematics::multiply(duration.seconds(), sampleRate));
        double deltaTime = ::mathematics::divide(1.0, sampleRate);
        for (size_t k = 0; k < samples; ++k)
        {
            Time time(::mathematics::multiply(static_cast<double>(k), deltaTime));
            result.append(displacement(time));
        }
        return result;
    }

    Momentum momentum(const Time &startTime, const Time &endTime) const { return field.momentum(startTime, endTime); }
    Energy energy() const { return field.energy(); }
    Entropy entropy() const { return field.entropy(); }

    Field<Force> operator+(const Field<Force> &other) const
    {
        Field<Force> result(field.getSampleRate(), field.size());
        result.field = field + other.field;
        return result;
    }
    Field<Force> operator-(const Field<Force> &other) const
    {
        Field<Force> result(field.getSampleRate(), field.size());
        result.field = field - other.field;
        return result;
    }
    Field<Force> operator*(double scalar) const
    {
        Field<Force> result(field.getSampleRate(), field.size());
        result.field = field * scalar;
        return result;
    }
    Field<Force> operator/(double scalar) const
    {
        Field<Force> result(field.getSampleRate(), field.size());
        result.field = field / scalar;
        return result;
    }
    Field<Force> &operator+=(const Field<Force> &other)
    {
        field += other.field;
        return *this;
    }
    Field<Force> &operator-=(const Field<Force> &other)
    {
        field -= other.field;
        return *this;
    }
    Field<Force> &operator*=(double scalar)
    {
        field *= scalar;
        return *this;
    }
    Field<Force> &operator/=(double scalar)
    {
        field /= scalar;
        return *this;
    }
    bool operator==(const Field<Force> &other) const { return field == other.field; }
    bool operator!=(const Field<Force> &other) const { return !(*this == other); }
    bool operator<(const Field<Force> &other) const { return field < other.field; }
    bool operator<=(const Field<Force> &other) const { return field <= other.field; }
    bool operator>(const Field<Force> &other) const { return field > other.field; }
    bool operator>=(const Field<Force> &other) const { return field >= other.field; }
    friend Field<Force> operator*(double scalar, const Field<Force> &f) { return f * scalar; }

    size_t size() const { return field.size(); }
    double getSampleRate() const { return field.getSampleRate(); }
    double duration() const { return field.duration(); }
    const Wave<::mathematics::Quaternion<double>> &getWave() const { return field.getWave(); }
    Wave<::mathematics::Quaternion<double>> &getWave() { return field.getWave(); }
};

template <>
class Field<Displacement>
{
  private:
    Field<::mathematics::Quaternion<double>> field;

  public:
    Field(double sampleRate, size_t size) : field(sampleRate, size) {}
    Field(double sampleRate, size_t size, const Displacement &value) : field(sampleRate, size, ::mathematics::Quaternion<double>(0.0, static_cast<double>(value.x()), static_cast<double>(value.y()), static_cast<double>(value.z()))) {}

    Displacement perturb(const Displacement &particle, Time timePosition, Time timeRange, Time timeDuration)
    {
        ::mathematics::Quaternion<double> particleQ(0.0, particle.x(), particle.y(), particle.z());
        ::mathematics::Quaternion<double> scatteredQ = field.perturb(particleQ, timePosition, timeRange, timeDuration);
        return Displacement(scatteredQ.getImagI(), scatteredQ.getImagJ(), scatteredQ.getImagK());
    }

    Wave<Displacement> perturb(const Wave<Displacement> &incoming, Time timePosition, Time timeRange, Time timeDuration)
    {
        Wave<::mathematics::Quaternion<double>> incomingQ(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            Displacement d = incoming[i];
            incomingQ.append(::mathematics::Quaternion<double>(0.0, d.x(), d.y(), d.z()));
        }
        Wave<::mathematics::Quaternion<double>> scatteredQ = field.perturb(incomingQ, timePosition, timeRange, timeDuration);
        Wave<Displacement> scattered(scatteredQ.getSampleRate());
        for (size_t i = 0; i < scatteredQ.size(); ++i)
        {
            scattered.append(Displacement(scatteredQ[i].getImagI(), scatteredQ[i].getImagJ(), scatteredQ[i].getImagK()));
        }
        return scattered;
    }

    Field<Displacement> &spread(double intensity, Time position, Time range, Time duration)
    {
        field.spread(intensity, position, range, duration);
        return *this;
    }

    Field<Displacement> &spread(const Displacement &complexIntensity, Time position, Time range, Time duration)
    {
        ::mathematics::Quaternion<double> intensityQ(0.0, complexIntensity.x(), complexIntensity.y(), complexIntensity.z());
        field.spread(intensityQ, position, range, duration);
        return *this;
    }

    Field<Displacement> &spread(const Wave<Displacement> &incoming, Time duration)
    {
        Wave<::mathematics::Quaternion<double>> incomingQ(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            Displacement d = incoming[i];
            incomingQ.append(::mathematics::Quaternion<double>(0.0, d.x(), d.y(), d.z()));
        }
        field.spread(incomingQ, duration);
        return *this;
    }

    Displacement displacement(const Time &time) const
    {
        auto q = field.displacement(time);
        return Displacement(q.getImagI(), q.getImagJ(), q.getImagK());
    }
    Displacement velocity(const Time &time) const
    {
        auto q = field.velocity(time);
        return Displacement(q.getImagI(), q.getImagJ(), q.getImagK());
    }
    Displacement accelerate(const Time &time) const
    {
        auto q = field.accelerate(time);
        return Displacement(q.getImagI(), q.getImagJ(), q.getImagK());
    }

    Wave<Displacement> extractWave(const Time &duration, double sampleRate) const
    {
        Wave<Displacement> result(sampleRate);
        size_t samples = static_cast<size_t>(::mathematics::multiply(duration.seconds(), sampleRate));
        double deltaTime = ::mathematics::divide(1.0, sampleRate);
        for (size_t k = 0; k < samples; ++k)
        {
            Time time(::mathematics::multiply(static_cast<double>(k), deltaTime));
            result.append(displacement(time));
        }
        return result;
    }

    Momentum momentum(const Time &startTime, const Time &endTime) const { return field.momentum(startTime, endTime); }
    Energy energy() const { return field.energy(); }
    Entropy entropy() const { return field.entropy(); }

    Field<Displacement> operator+(const Field<Displacement> &other) const
    {
        Field<Displacement> result(field.getSampleRate(), field.size());
        result.field = field + other.field;
        return result;
    }
    Field<Displacement> operator-(const Field<Displacement> &other) const
    {
        Field<Displacement> result(field.getSampleRate(), field.size());
        result.field = field - other.field;
        return result;
    }
    Field<Displacement> operator*(double scalar) const
    {
        Field<Displacement> result(field.getSampleRate(), field.size());
        result.field = field * scalar;
        return result;
    }
    Field<Displacement> operator/(double scalar) const
    {
        Field<Displacement> result(field.getSampleRate(), field.size());
        result.field = field / scalar;
        return result;
    }
    Field<Displacement> &operator+=(const Field<Displacement> &other)
    {
        field += other.field;
        return *this;
    }
    Field<Displacement> &operator-=(const Field<Displacement> &other)
    {
        field -= other.field;
        return *this;
    }
    Field<Displacement> &operator*=(double scalar)
    {
        field *= scalar;
        return *this;
    }
    Field<Displacement> &operator/=(double scalar)
    {
        field /= scalar;
        return *this;
    }
    bool operator==(const Field<Displacement> &other) const { return field == other.field; }
    bool operator!=(const Field<Displacement> &other) const { return !(*this == other); }
    bool operator<(const Field<Displacement> &other) const { return field < other.field; }
    bool operator<=(const Field<Displacement> &other) const { return field <= other.field; }
    bool operator>(const Field<Displacement> &other) const { return field > other.field; }
    bool operator>=(const Field<Displacement> &other) const { return field >= other.field; }
    friend Field<Displacement> operator*(double scalar, const Field<Displacement> &f) { return f * scalar; }

    size_t size() const { return field.size(); }
    double getSampleRate() const { return field.getSampleRate(); }
    double duration() const { return field.duration(); }
    const Wave<::mathematics::Quaternion<double>> &getWave() const { return field.getWave(); }
    Wave<::mathematics::Quaternion<double>> &getWave() { return field.getWave(); }
};

template <>
class Field<Velocity>
{
  private:
    Field<::mathematics::Quaternion<double>> field;

  public:
    Field(double sampleRate, size_t size) : field(sampleRate, size) {}
    Field(double sampleRate, size_t size, const Velocity &value) : field(sampleRate, size, ::mathematics::Quaternion<double>(0.0, static_cast<double>(value.x()), static_cast<double>(value.y()), static_cast<double>(value.z()))) {}

    Velocity perturb(const Velocity &particle, Time timePosition, Time timeRange, Time timeDuration)
    {
        ::mathematics::Quaternion<double> particleQ(0.0, particle.x(), particle.y(), particle.z());
        ::mathematics::Quaternion<double> scatteredQ = field.perturb(particleQ, timePosition, timeRange, timeDuration);
        return Velocity(scatteredQ.getImagI(), scatteredQ.getImagJ(), scatteredQ.getImagK());
    }

    Wave<Velocity> perturb(const Wave<Velocity> &incoming, Time timePosition, Time timeRange, Time timeDuration)
    {
        Wave<::mathematics::Quaternion<double>> incomingQ(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            Velocity v = incoming[i];
            incomingQ.append(::mathematics::Quaternion<double>(0.0, v.x(), v.y(), v.z()));
        }
        Wave<::mathematics::Quaternion<double>> scatteredQ = field.perturb(incomingQ, timePosition, timeRange, timeDuration);
        Wave<Velocity> scattered(scatteredQ.getSampleRate());
        for (size_t i = 0; i < scatteredQ.size(); ++i)
        {
            scattered.append(Velocity(scatteredQ[i].getImagI(), scatteredQ[i].getImagJ(), scatteredQ[i].getImagK()));
        }
        return scattered;
    }

    Field<Velocity> &spread(double intensity, Time position, Time range, Time duration)
    {
        field.spread(intensity, position, range, duration);
        return *this;
    }

    Field<Velocity> &spread(const Velocity &complexIntensity, Time position, Time range, Time duration)
    {
        ::mathematics::Quaternion<double> intensityQ(0.0, complexIntensity.x(), complexIntensity.y(), complexIntensity.z());
        field.spread(intensityQ, position, range, duration);
        return *this;
    }

    Field<Velocity> &spread(const Wave<Velocity> &incoming, Time duration)
    {
        Wave<::mathematics::Quaternion<double>> incomingQ(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            Velocity v = incoming[i];
            incomingQ.append(::mathematics::Quaternion<double>(0.0, v.x(), v.y(), v.z()));
        }
        field.spread(incomingQ, duration);
        return *this;
    }

    Velocity displacement(const Time &time) const
    {
        auto q = field.displacement(time);
        return Velocity(q.getImagI(), q.getImagJ(), q.getImagK());
    }
    Velocity velocity(const Time &time) const
    {
        auto q = field.velocity(time);
        return Velocity(q.getImagI(), q.getImagJ(), q.getImagK());
    }
    Velocity accelerate(const Time &time) const
    {
        auto q = field.accelerate(time);
        return Velocity(q.getImagI(), q.getImagJ(), q.getImagK());
    }

    Wave<Velocity> extractWave(const Time &duration, double sampleRate) const
    {
        Wave<Velocity> result(sampleRate);
        size_t samples = static_cast<size_t>(::mathematics::multiply(duration.seconds(), sampleRate));
        double deltaTime = ::mathematics::divide(1.0, sampleRate);
        for (size_t k = 0; k < samples; ++k)
        {
            Time time(::mathematics::multiply(static_cast<double>(k), deltaTime));
            result.append(displacement(time));
        }
        return result;
    }

    Momentum momentum(const Time &startTime, const Time &endTime) const { return field.momentum(startTime, endTime); }
    Energy energy() const { return field.energy(); }
    Entropy entropy() const { return field.entropy(); }

    Field<Velocity> operator+(const Field<Velocity> &other) const
    {
        Field<Velocity> result(field.getSampleRate(), field.size());
        result.field = field + other.field;
        return result;
    }
    Field<Velocity> operator-(const Field<Velocity> &other) const
    {
        Field<Velocity> result(field.getSampleRate(), field.size());
        result.field = field - other.field;
        return result;
    }
    Field<Velocity> operator*(double scalar) const
    {
        Field<Velocity> result(field.getSampleRate(), field.size());
        result.field = field * scalar;
        return result;
    }
    Field<Velocity> operator/(double scalar) const
    {
        Field<Velocity> result(field.getSampleRate(), field.size());
        result.field = field / scalar;
        return result;
    }
    Field<Velocity> &operator+=(const Field<Velocity> &other)
    {
        field += other.field;
        return *this;
    }
    Field<Velocity> &operator-=(const Field<Velocity> &other)
    {
        field -= other.field;
        return *this;
    }
    Field<Velocity> &operator*=(double scalar)
    {
        field *= scalar;
        return *this;
    }
    Field<Velocity> &operator/=(double scalar)
    {
        field /= scalar;
        return *this;
    }
    bool operator==(const Field<Velocity> &other) const { return field == other.field; }
    bool operator!=(const Field<Velocity> &other) const { return !(*this == other); }
    bool operator<(const Field<Velocity> &other) const { return field < other.field; }
    bool operator<=(const Field<Velocity> &other) const { return field <= other.field; }
    bool operator>(const Field<Velocity> &other) const { return field > other.field; }
    bool operator>=(const Field<Velocity> &other) const { return field >= other.field; }
    friend Field<Velocity> operator*(double scalar, const Field<Velocity> &f) { return f * scalar; }

    size_t size() const { return field.size(); }
    double getSampleRate() const { return field.getSampleRate(); }
    double duration() const { return field.duration(); }
    const Wave<::mathematics::Quaternion<double>> &getWave() const { return field.getWave(); }
    Wave<::mathematics::Quaternion<double>> &getWave() { return field.getWave(); }
};

template <>
class Field<Acceleration>
{
  private:
    Field<::mathematics::Quaternion<double>> field;

  public:
    Field(double sampleRate, size_t size) : field(sampleRate, size) {}
    Field(double sampleRate, size_t size, const Acceleration &value) : field(sampleRate, size, ::mathematics::Quaternion<double>(0.0, static_cast<double>(value.x()), static_cast<double>(value.y()), static_cast<double>(value.z()))) {}

    Acceleration perturb(const Acceleration &particle, Time timePosition, Time timeRange, Time timeDuration)
    {
        ::mathematics::Quaternion<double> particleQ(0.0, particle.x(), particle.y(), particle.z());
        ::mathematics::Quaternion<double> scatteredQ = field.perturb(particleQ, timePosition, timeRange, timeDuration);
        return Acceleration(scatteredQ.getImagI(), scatteredQ.getImagJ(), scatteredQ.getImagK());
    }

    Wave<Acceleration> perturb(const Wave<Acceleration> &incoming, Time timePosition, Time timeRange, Time timeDuration)
    {
        Wave<::mathematics::Quaternion<double>> incomingQ(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            Acceleration a = incoming[i];
            incomingQ.append(::mathematics::Quaternion<double>(0.0, a.x(), a.y(), a.z()));
        }
        Wave<::mathematics::Quaternion<double>> scatteredQ = field.perturb(incomingQ, timePosition, timeRange, timeDuration);
        Wave<Acceleration> scattered(scatteredQ.getSampleRate());
        for (size_t i = 0; i < scatteredQ.size(); ++i)
        {
            scattered.append(Acceleration(scatteredQ[i].getImagI(), scatteredQ[i].getImagJ(), scatteredQ[i].getImagK()));
        }
        return scattered;
    }

    Field<Acceleration> &spread(double intensity, Time position, Time range, Time duration)
    {
        field.spread(intensity, position, range, duration);
        return *this;
    }

    Field<Acceleration> &spread(const Acceleration &complexIntensity, Time position, Time range, Time duration)
    {
        ::mathematics::Quaternion<double> intensityQ(0.0, complexIntensity.x(), complexIntensity.y(), complexIntensity.z());
        field.spread(intensityQ, position, range, duration);
        return *this;
    }

    Field<Acceleration> &spread(const Wave<Acceleration> &incoming, Time duration)
    {
        Wave<::mathematics::Quaternion<double>> incomingQ(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            Acceleration a = incoming[i];
            incomingQ.append(::mathematics::Quaternion<double>(0.0, a.x(), a.y(), a.z()));
        }
        field.spread(incomingQ, duration);
        return *this;
    }

    Acceleration displacement(const Time &time) const
    {
        auto q = field.displacement(time);
        return Acceleration(q.getImagI(), q.getImagJ(), q.getImagK());
    }
    Acceleration velocity(const Time &time) const
    {
        auto q = field.velocity(time);
        return Acceleration(q.getImagI(), q.getImagJ(), q.getImagK());
    }
    Acceleration accelerate(const Time &time) const
    {
        auto q = field.accelerate(time);
        return Acceleration(q.getImagI(), q.getImagJ(), q.getImagK());
    }

    Wave<Acceleration> extractWave(const Time &duration, double sampleRate) const
    {
        Wave<Acceleration> result(sampleRate);
        size_t samples = static_cast<size_t>(::mathematics::multiply(duration.seconds(), sampleRate));
        double deltaTime = ::mathematics::divide(1.0, sampleRate);
        for (size_t k = 0; k < samples; ++k)
        {
            Time time(::mathematics::multiply(static_cast<double>(k), deltaTime));
            result.append(displacement(time));
        }
        return result;
    }

    Momentum momentum(const Time &startTime, const Time &endTime) const { return field.momentum(startTime, endTime); }
    Energy energy() const { return field.energy(); }
    Entropy entropy() const { return field.entropy(); }

    Field<Acceleration> operator+(const Field<Acceleration> &other) const
    {
        Field<Acceleration> result(field.getSampleRate(), field.size());
        result.field = field + other.field;
        return result;
    }
    Field<Acceleration> operator-(const Field<Acceleration> &other) const
    {
        Field<Acceleration> result(field.getSampleRate(), field.size());
        result.field = field - other.field;
        return result;
    }
    Field<Acceleration> operator*(double scalar) const
    {
        Field<Acceleration> result(field.getSampleRate(), field.size());
        result.field = field * scalar;
        return result;
    }
    Field<Acceleration> operator/(double scalar) const
    {
        Field<Acceleration> result(field.getSampleRate(), field.size());
        result.field = field / scalar;
        return result;
    }
    Field<Acceleration> &operator+=(const Field<Acceleration> &other)
    {
        field += other.field;
        return *this;
    }
    Field<Acceleration> &operator-=(const Field<Acceleration> &other)
    {
        field -= other.field;
        return *this;
    }
    Field<Acceleration> &operator*=(double scalar)
    {
        field *= scalar;
        return *this;
    }
    Field<Acceleration> &operator/=(double scalar)
    {
        field /= scalar;
        return *this;
    }
    bool operator==(const Field<Acceleration> &other) const { return field == other.field; }
    bool operator!=(const Field<Acceleration> &other) const { return !(*this == other); }
    bool operator<(const Field<Acceleration> &other) const { return field < other.field; }
    bool operator<=(const Field<Acceleration> &other) const { return field <= other.field; }
    bool operator>(const Field<Acceleration> &other) const { return field > other.field; }
    bool operator>=(const Field<Acceleration> &other) const { return field >= other.field; }
    friend Field<Acceleration> operator*(double scalar, const Field<Acceleration> &f) { return f * scalar; }

    size_t size() const { return field.size(); }
    double getSampleRate() const { return field.getSampleRate(); }
    double duration() const { return field.duration(); }
    const Wave<::mathematics::Quaternion<double>> &getWave() const { return field.getWave(); }
    Wave<::mathematics::Quaternion<double>> &getWave() { return field.getWave(); }
};

template <>
class Field<Momentum>
{
  private:
    Field<::mathematics::Quaternion<double>> field;

  public:
    Field(double sampleRate, size_t size) : field(sampleRate, size) {}
    Field(double sampleRate, size_t size, const Momentum &value) : field(sampleRate, size, ::mathematics::Quaternion<double>(0.0, static_cast<double>(value.x()), static_cast<double>(value.y()), static_cast<double>(value.z()))) {}

    Momentum perturb(const Momentum &particle, Time timePosition, Time timeRange, Time timeDuration)
    {
        ::mathematics::Quaternion<double> particleQ(0.0, particle.x(), particle.y(), particle.z());
        ::mathematics::Quaternion<double> scatteredQ = field.perturb(particleQ, timePosition, timeRange, timeDuration);
        return Momentum(scatteredQ.getImagI(), scatteredQ.getImagJ(), scatteredQ.getImagK());
    }

    Wave<Momentum> perturb(const Wave<Momentum> &incoming, Time timePosition, Time timeRange, Time timeDuration)
    {
        Wave<::mathematics::Quaternion<double>> incomingQ(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            Momentum p = incoming[i];
            incomingQ.append(::mathematics::Quaternion<double>(0.0, p.x(), p.y(), p.z()));
        }
        Wave<::mathematics::Quaternion<double>> scatteredQ = field.perturb(incomingQ, timePosition, timeRange, timeDuration);
        Wave<Momentum> scattered(scatteredQ.getSampleRate());
        for (size_t i = 0; i < scatteredQ.size(); ++i)
        {
            scattered.append(Momentum(scatteredQ[i].getImagI(), scatteredQ[i].getImagJ(), scatteredQ[i].getImagK()));
        }
        return scattered;
    }

    Field<Momentum> &spread(double intensity, Time position, Time range, Time duration)
    {
        field.spread(intensity, position, range, duration);
        return *this;
    }

    Field<Momentum> &spread(const Momentum &complexIntensity, Time position, Time range, Time duration)
    {
        ::mathematics::Quaternion<double> intensityQ(0.0, complexIntensity.x(), complexIntensity.y(), complexIntensity.z());
        field.spread(intensityQ, position, range, duration);
        return *this;
    }

    Field<Momentum> &spread(const Wave<Momentum> &incoming, Time duration)
    {
        Wave<::mathematics::Quaternion<double>> incomingQ(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            Momentum p = incoming[i];
            incomingQ.append(::mathematics::Quaternion<double>(0.0, p.x(), p.y(), p.z()));
        }
        field.spread(incomingQ, duration);
        return *this;
    }

    Momentum displacement(const Time &time) const
    {
        auto q = field.displacement(time);
        return Momentum(q.getImagI(), q.getImagJ(), q.getImagK());
    }
    Momentum velocity(const Time &time) const
    {
        auto q = field.velocity(time);
        return Momentum(q.getImagI(), q.getImagJ(), q.getImagK());
    }
    Momentum accelerate(const Time &time) const
    {
        auto q = field.accelerate(time);
        return Momentum(q.getImagI(), q.getImagJ(), q.getImagK());
    }

    Wave<Momentum> extractWave(const Time &duration, double sampleRate) const
    {
        Wave<Momentum> result(sampleRate);
        size_t samples = static_cast<size_t>(::mathematics::multiply(duration.seconds(), sampleRate));
        double deltaTime = ::mathematics::divide(1.0, sampleRate);
        for (size_t k = 0; k < samples; ++k)
        {
            Time time(::mathematics::multiply(static_cast<double>(k), deltaTime));
            result.append(displacement(time));
        }
        return result;
    }

    Momentum momentum(const Time &startTime, const Time &endTime) const { return field.momentum(startTime, endTime); }
    Energy energy() const { return field.energy(); }
    Entropy entropy() const { return field.entropy(); }

    Field<Momentum> operator+(const Field<Momentum> &other) const
    {
        Field<Momentum> result(field.getSampleRate(), field.size());
        result.field = field + other.field;
        return result;
    }
    Field<Momentum> operator-(const Field<Momentum> &other) const
    {
        Field<Momentum> result(field.getSampleRate(), field.size());
        result.field = field - other.field;
        return result;
    }
    Field<Momentum> operator*(double scalar) const
    {
        Field<Momentum> result(field.getSampleRate(), field.size());
        result.field = field * scalar;
        return result;
    }
    Field<Momentum> operator/(double scalar) const
    {
        Field<Momentum> result(field.getSampleRate(), field.size());
        result.field = field / scalar;
        return result;
    }
    Field<Momentum> &operator+=(const Field<Momentum> &other)
    {
        field += other.field;
        return *this;
    }
    Field<Momentum> &operator-=(const Field<Momentum> &other)
    {
        field -= other.field;
        return *this;
    }
    Field<Momentum> &operator*=(double scalar)
    {
        field *= scalar;
        return *this;
    }
    Field<Momentum> &operator/=(double scalar)
    {
        field /= scalar;
        return *this;
    }
    bool operator==(const Field<Momentum> &other) const { return field == other.field; }
    bool operator!=(const Field<Momentum> &other) const { return !(*this == other); }
    bool operator<(const Field<Momentum> &other) const { return field < other.field; }
    bool operator<=(const Field<Momentum> &other) const { return field <= other.field; }
    bool operator>(const Field<Momentum> &other) const { return field > other.field; }
    bool operator>=(const Field<Momentum> &other) const { return field >= other.field; }
    friend Field<Momentum> operator*(double scalar, const Field<Momentum> &f) { return f * scalar; }

    size_t size() const { return field.size(); }
    double getSampleRate() const { return field.getSampleRate(); }
    double duration() const { return field.duration(); }
    const Wave<::mathematics::Quaternion<double>> &getWave() const { return field.getWave(); }
    Wave<::mathematics::Quaternion<double>> &getWave() { return field.getWave(); }
};

template <>
class Field<Torque>
{
  private:
    Field<::mathematics::Quaternion<double>> field;

  public:
    Field(double sampleRate, size_t size) : field(sampleRate, size) {}
    Field(double sampleRate, size_t size, const Torque &value) : field(sampleRate, size, ::mathematics::Quaternion<double>(0.0, static_cast<double>(value.x()), static_cast<double>(value.y()), static_cast<double>(value.z()))) {}

    Torque perturb(const Torque &particle, Time timePosition, Time timeRange, Time timeDuration)
    {
        ::mathematics::Quaternion<double> particleQ(0.0, particle.x(), particle.y(), particle.z());
        ::mathematics::Quaternion<double> scatteredQ = field.perturb(particleQ, timePosition, timeRange, timeDuration);
        return Torque(scatteredQ.getImagI(), scatteredQ.getImagJ(), scatteredQ.getImagK());
    }

    Wave<Torque> perturb(const Wave<Torque> &incoming, Time timePosition, Time timeRange, Time timeDuration)
    {
        Wave<::mathematics::Quaternion<double>> incomingQ(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            Torque t = incoming[i];
            incomingQ.append(::mathematics::Quaternion<double>(0.0, t.x(), t.y(), t.z()));
        }
        Wave<::mathematics::Quaternion<double>> scatteredQ = field.perturb(incomingQ, timePosition, timeRange, timeDuration);
        Wave<Torque> scattered(scatteredQ.getSampleRate());
        for (size_t i = 0; i < scatteredQ.size(); ++i)
        {
            scattered.append(Torque(scatteredQ[i].getImagI(), scatteredQ[i].getImagJ(), scatteredQ[i].getImagK()));
        }
        return scattered;
    }

    Field<Torque> &spread(double intensity, Time position, Time range, Time duration)
    {
        field.spread(intensity, position, range, duration);
        return *this;
    }

    Field<Torque> &spread(const Torque &complexIntensity, Time position, Time range, Time duration)
    {
        ::mathematics::Quaternion<double> intensityQ(0.0, complexIntensity.x(), complexIntensity.y(), complexIntensity.z());
        field.spread(intensityQ, position, range, duration);
        return *this;
    }

    Field<Torque> &spread(const Wave<Torque> &incoming, Time duration)
    {
        Wave<::mathematics::Quaternion<double>> incomingQ(incoming.getSampleRate());
        for (size_t i = 0; i < incoming.size(); ++i)
        {
            Torque t = incoming[i];
            incomingQ.append(::mathematics::Quaternion<double>(0.0, t.x(), t.y(), t.z()));
        }
        field.spread(incomingQ, duration);
        return *this;
    }

    Torque displacement(const Time &time) const
    {
        auto q = field.displacement(time);
        return Torque(q.getImagI(), q.getImagJ(), q.getImagK());
    }
    Torque velocity(const Time &time) const
    {
        auto q = field.velocity(time);
        return Torque(q.getImagI(), q.getImagJ(), q.getImagK());
    }
    Torque accelerate(const Time &time) const
    {
        auto q = field.accelerate(time);
        return Torque(q.getImagI(), q.getImagJ(), q.getImagK());
    }

    Wave<Torque> extractWave(const Time &duration, double sampleRate) const
    {
        Wave<Torque> result(sampleRate);
        size_t samples = static_cast<size_t>(::mathematics::multiply(duration.seconds(), sampleRate));
        double deltaTime = ::mathematics::divide(1.0, sampleRate);
        for (size_t k = 0; k < samples; ++k)
        {
            Time time(::mathematics::multiply(static_cast<double>(k), deltaTime));
            result.append(displacement(time));
        }
        return result;
    }

    Momentum momentum(const Time &startTime, const Time &endTime) const { return field.momentum(startTime, endTime); }
    Energy energy() const { return field.energy(); }
    Entropy entropy() const { return field.entropy(); }

    Field<Torque> operator+(const Field<Torque> &other) const
    {
        Field<Torque> result(field.getSampleRate(), field.size());
        result.field = field + other.field;
        return result;
    }
    Field<Torque> operator-(const Field<Torque> &other) const
    {
        Field<Torque> result(field.getSampleRate(), field.size());
        result.field = field - other.field;
        return result;
    }
    Field<Torque> operator*(double scalar) const
    {
        Field<Torque> result(field.getSampleRate(), field.size());
        result.field = field * scalar;
        return result;
    }
    Field<Torque> operator/(double scalar) const
    {
        Field<Torque> result(field.getSampleRate(), field.size());
        result.field = field / scalar;
        return result;
    }
    Field<Torque> &operator+=(const Field<Torque> &other)
    {
        field += other.field;
        return *this;
    }
    Field<Torque> &operator-=(const Field<Torque> &other)
    {
        field -= other.field;
        return *this;
    }
    Field<Torque> &operator*=(double scalar)
    {
        field *= scalar;
        return *this;
    }
    Field<Torque> &operator/=(double scalar)
    {
        field /= scalar;
        return *this;
    }
    bool operator==(const Field<Torque> &other) const { return field == other.field; }
    bool operator!=(const Field<Torque> &other) const { return !(*this == other); }
    bool operator<(const Field<Torque> &other) const { return field < other.field; }
    bool operator<=(const Field<Torque> &other) const { return field <= other.field; }
    bool operator>(const Field<Torque> &other) const { return field > other.field; }
    bool operator>=(const Field<Torque> &other) const { return field >= other.field; }
    friend Field<Torque> operator*(double scalar, const Field<Torque> &f) { return f * scalar; }

    size_t size() const { return field.size(); }
    double getSampleRate() const { return field.getSampleRate(); }
    double duration() const { return field.duration(); }
    const Wave<::mathematics::Quaternion<double>> &getWave() const { return field.getWave(); }
    Wave<::mathematics::Quaternion<double>> &getWave() { return field.getWave(); }
};

} // namespace physics

namespace std
{
template <typename T>
struct hash<physics::Field<T>>
{
    size_t operator()(const physics::Field<T> &field) const
    {
        const physics::Wave<::mathematics::Quaternion<double>> &wave = field.getWave();
        size_t seed = 0;
        for (size_t i = 0; i < wave.size(); ++i)
        {
            const ::mathematics::Quaternion<double> &q = wave[i];
            seed ^= hash<double>{}(q.getReal()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= hash<double>{}(q.getImagI()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= hash<double>{}(q.getImagJ()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= hash<double>{}(q.getImagK()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

template <typename T>
struct less<physics::Field<T>>
{
    bool operator()(const physics::Field<T> &left, const physics::Field<T> &right) const
    {
        return left < right;
    }
};
} // namespace std