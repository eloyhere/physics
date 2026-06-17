#pragma once
#include <tuple>
#include <cmath>
#include <vector>
#include <complex>
#include <initializer_list>
#include <stdexcept>
#include <limits>
#include <iostream>
#include "mathematics.h"

namespace mathematics
{
template <typename N>
class Complex
{
  private:
    N real;
    N imag;

  public:
    Complex() : real(N(0)), imag(N(0)) {}

    explicit Complex(const N &valueReal) : real(valueReal), imag(N(0)) {}

    explicit Complex(const N &valueReal, const N &valueImag) : real(valueReal), imag(valueImag) {}

    explicit Complex(const std::initializer_list<N> &list)
    {
        switch (list.size())
        {
        case 0:
            real = N(0);
            imag = N(0);
            break;
        case 1:
            real = *(list.begin());
            imag = N(0);
            break;
        default:
            real = *(list.begin());
            imag = *(list.begin() + 1);
            break;
        }
    }

    Complex(const std::complex<N> &complexValue) : real(complexValue.real()), imag(complexValue.imag()) {}

    operator std::complex<N>() const
    {
        return std::complex<N>(real, imag);
    }

    const N &getReal() const { return real; }
    N &getReal() { return real; }
    const N &getImag() const { return imag; }
    N &getImag() { return imag; }

    N realPart() const { return real; }
    void realPart(const N &value) { real = value; }
    N imagPart() const { return imag; }
    void imagPart(const N &value) { imag = value; }

    Complex<N> exact(N epsilon = N(::mathematics::tolerance)) const
    {
        return Complex<N>(::mathematics::exact(real, epsilon), ::mathematics::exact(imag, epsilon));
    }

    Complex<N> limit(N start, N end) const
    {
        return Complex<N>(::mathematics::limit(real, start, end), ::mathematics::limit(imag, start, end));
    }

    Complex<N> smooth(N threshold, N steepness) const
    {
        return Complex<N>(::mathematics::smooth(real, threshold, steepness), ::mathematics::smooth(imag, threshold, steepness));
    }

    Complex<N> conjugate() const
    {
        return Complex<N>(real, ::mathematics::negative(imag));
    }

    Complex<N> inverse() const
    {
        N denominator = ::mathematics::add(::mathematics::multiply(real, real), ::mathematics::multiply(imag, imag));
        if (denominator < N(1e-11))
        {
            return Complex<N>(N(1), N(0));
        }
        N scale = ::mathematics::divide(N(1), denominator);
        return Complex<N>(::mathematics::multiply(real, scale), ::mathematics::multiply(::mathematics::negative(imag), scale));
    }

    N module() const
    {
        N maximumValue = std::max(::mathematics::absolute(real), ::mathematics::absolute(imag));
        if (maximumValue < N(1e-11))
        {
            return N(0);
        }
        N scale = ::mathematics::divide(N(1), maximumValue);
        N scaledReal = ::mathematics::multiply(real, scale);
        N scaledImag = ::mathematics::multiply(imag, scale);
        N sumSquares = ::mathematics::add(::mathematics::multiply(scaledReal, scaledReal), ::mathematics::multiply(scaledImag, scaledImag));
        return ::mathematics::multiply(::mathematics::squareRoot(sumSquares), maximumValue);
    }

    N argument() const
    {
        return ::mathematics::arctangent2(imag, real);
    }

    Complex<N> normalize() const
    {
        N magnitude = module();
        if (magnitude < N(1e-11))
        {
            return Complex<N>(N(1), N(0));
        }
        N scale = ::mathematics::divide(N(1), magnitude);
        return Complex<N>(::mathematics::multiply(real, scale), ::mathematics::multiply(imag, scale));
    }

    Complex<N> exponential() const
    {
        N expReal = ::mathematics::exponential(real);
        if (::mathematics::absolute(imag) < N(1e-11))
        {
            return Complex<N>(expReal, N(0));
        }
        N cosValue = ::mathematics::cosine(imag);
        N sinValue = ::mathematics::sine(imag);
        return Complex<N>(::mathematics::multiply(expReal, cosValue), ::mathematics::multiply(expReal, sinValue));
    }

    Complex<N> logarithm() const
    {
        N magnitude = module();
        if (magnitude < N(1e-11))
        {
            return Complex<N>(N(0), N(0));
        }
        N logMagnitude = ::mathematics::logarithm(magnitude);
        N arg = argument();
        return Complex<N>(logMagnitude, arg);
    }

    Complex<N> logarithm(const Complex<N> &base) const
    {
        return logarithm() * base.logarithm().inverse();
    }

    Complex<N> power(const Complex<N> &exponent) const
    {
        N magnitude = module();
        if (magnitude < N(1e-11))
        {
            if (::mathematics::absolute(exponent.real) < N(1e-11) && ::mathematics::absolute(exponent.imag) < N(1e-11))
            {
                return Complex<N>(N(1), N(0));
            }
            return Complex<N>(N(0), N(0));
        }
        return (logarithm() * exponent).exponential();
    }

    Complex<N> power(const N &exponent) const
    {
        return power(Complex<N>(exponent, N(0)));
    }

    Complex<N> squareRoot() const
    {
        N magnitude = module();
        if (magnitude < N(1e-11))
        {
            return Complex<N>(N(0), N(0));
        }
        N arg = argument();
        N sqrtMagnitude = ::mathematics::squareRoot(magnitude);
        N halfArg = ::mathematics::divide(arg, N(2));
        return Complex<N>(::mathematics::multiply(sqrtMagnitude, ::mathematics::cosine(halfArg)), ::mathematics::multiply(sqrtMagnitude, ::mathematics::sine(halfArg)));
    }

    Complex<N> sine() const
    {
        N sinReal = ::mathematics::sine(real);
        N cosReal = ::mathematics::cosine(real);
        N sinhImag = ::mathematics::hyperbolicSine(imag);
        N coshImag = ::mathematics::hyperbolicCosine(imag);
        return Complex<N>(::mathematics::multiply(sinReal, coshImag), ::mathematics::multiply(cosReal, sinhImag));
    }

    Complex<N> cosine() const
    {
        N cosReal = ::mathematics::cosine(real);
        N sinReal = ::mathematics::sine(real);
        N sinhImag = ::mathematics::hyperbolicSine(imag);
        N coshImag = ::mathematics::hyperbolicCosine(imag);
        return Complex<N>(::mathematics::multiply(cosReal, coshImag), ::mathematics::multiply(::mathematics::negative(sinReal), sinhImag));
    }

    Complex<N> tangent() const
    {
        return sine() / cosine();
    }

    Complex<N> arcsine() const
    {
        Complex<N> i(N(0), N(1));
        return ::mathematics::negative(i) * ((i * (*this) + (Complex<N>(N(1), N(0)) - (*this) * (*this)).squareRoot()).logarithm());
    }

    Complex<N> arccosine() const
    {
        Complex<N> i(N(0), N(1));
        return ::mathematics::negative(i) * ((*this + i * (Complex<N>(N(1), N(0)) - (*this) * (*this)).squareRoot()).logarithm());
    }

    Complex<N> arctangent() const
    {
        Complex<N> i(N(0), N(1));
        return ::mathematics::multiply(i, N(0.5)) * ((Complex<N>(N(1), N(0)) - i * (*this)).logarithm() - (Complex<N>(N(1), N(0)) + i * (*this)).logarithm());
    }

    Complex<N> rotate(const N &angle) const
    {
        Complex<N> rotation(::mathematics::cosine(angle), ::mathematics::sine(angle));
        return (*this) * rotation;
    }

    Complex<N> rotate(const Complex<N> &rotation) const
    {
        return (*this) * rotation.normalize();
    }

    Complex<N> lerp(const Complex<N> &other, N t) const
    {
        N oneMinusT = ::mathematics::subtract(N(1), t);
        return Complex<N>(::mathematics::add(::mathematics::multiply(real, oneMinusT), ::mathematics::multiply(other.real, t)),
                          ::mathematics::add(::mathematics::multiply(imag, oneMinusT), ::mathematics::multiply(other.imag, t)));
    }

    Complex<N> nlerp(const Complex<N> &other, N t) const
    {
        Complex<N> interpolated = lerp(other, t);
        N magnitude = interpolated.module();
        if (magnitude < N(1e-11))
        {
            return *this;
        }
        N scale = ::mathematics::divide(N(1), magnitude);
        return Complex<N>(::mathematics::multiply(interpolated.real, scale), ::mathematics::multiply(interpolated.imag, scale));
    }

    Complex<N> slerp(const Complex<N> &other, N t) const
    {
        Complex<N> normalizedThis = normalize();
        Complex<N> normalizedOther = other.normalize();
        N dotProduct = ::mathematics::add(::mathematics::multiply(normalizedThis.real, normalizedOther.real), ::mathematics::multiply(normalizedThis.imag, normalizedOther.imag));
        Complex<N> target = normalizedOther;
        if (dotProduct < N(0))
        {
            dotProduct = ::mathematics::negative(dotProduct);
            target = target.negative();
        }
        if (dotProduct > N(1))
        {
            dotProduct = N(1);
        }
        N theta = ::mathematics::arccosine(dotProduct);
        N sinTheta = ::mathematics::sine(theta);
        if (::mathematics::absolute(sinTheta) < N(1e-11))
        {
            return lerp(normalizedOther, t);
        }
        N factorA = ::mathematics::divide(::mathematics::sine(::mathematics::multiply(::mathematics::subtract(N(1), t), theta)), sinTheta);
        N factorB = ::mathematics::divide(::mathematics::sine(::mathematics::multiply(t, theta)), sinTheta);
        return Complex<N>(::mathematics::add(::mathematics::multiply(normalizedThis.real, factorA), ::mathematics::multiply(target.real, factorB)),
                          ::mathematics::add(::mathematics::multiply(normalizedThis.imag, factorA), ::mathematics::multiply(target.imag, factorB)));
    }

    N dot(const Complex<N> &other) const
    {
        return ::mathematics::add(::mathematics::multiply(real, other.real), ::mathematics::multiply(imag, other.imag));
    }

    N distance(const Complex<N> &other) const
    {
        Complex<N> diff = (*this) - other;
        return diff.module();
    }

    N angle(const Complex<N> &other) const
    {
        Complex<N> normalizedThis = normalize();
        Complex<N> normalizedOther = other.normalize();
        N dotProduct = ::mathematics::add(::mathematics::multiply(normalizedThis.real, normalizedOther.real), ::mathematics::multiply(normalizedThis.imag, normalizedOther.imag));
        if (dotProduct > N(1))
        {
            dotProduct = N(1);
        }
        if (dotProduct < N(-1))
        {
            dotProduct = N(-1);
        }
        return ::mathematics::arccosine(dotProduct);
    }

    std::tuple<std::complex<N>, std::complex<N>> toComplex() const
    {
        return {std::complex<N>(real, imag), std::complex<N>(real, ::mathematics::negative(imag))};
    }

    Complex<N> &operator=(const Complex<N> &) = default;

    Complex<N> operator+(const Complex<N> &other) const
    {
        return Complex<N>(::mathematics::add(real, other.real), ::mathematics::add(imag, other.imag));
    }

    Complex<N> operator-(const Complex<N> &other) const
    {
        return Complex<N>(::mathematics::subtract(real, other.real), ::mathematics::subtract(imag, other.imag));
    }

    Complex<N> operator*(const Complex<N> &other) const
    {
        return Complex<N>(::mathematics::subtract(::mathematics::multiply(real, other.real), ::mathematics::multiply(imag, other.imag)),
                          ::mathematics::add(::mathematics::multiply(real, other.imag), ::mathematics::multiply(imag, other.real)));
    }

    Complex<N> operator/(const Complex<N> &other) const
    {
        return (*this) * other.inverse();
    }

    Complex<N> operator*(const N &scalar) const
    {
        return Complex<N>(::mathematics::multiply(real, scalar), ::mathematics::multiply(imag, scalar));
    }

    Complex<N> operator/(const N &scalar) const
    {
        return Complex<N>(::mathematics::divide(real, scalar), ::mathematics::divide(imag, scalar));
    }

    Complex<N> &operator+=(const Complex<N> &other)
    {
        real = ::mathematics::add(real, other.real);
        imag = ::mathematics::add(imag, other.imag);
        return *this;
    }

    Complex<N> &operator-=(const Complex<N> &other)
    {
        real = ::mathematics::subtract(real, other.real);
        imag = ::mathematics::subtract(imag, other.imag);
        return *this;
    }

    Complex<N> &operator*=(const Complex<N> &other)
    {
        *this = *this * other;
        return *this;
    }

    Complex<N> &operator/=(const Complex<N> &other)
    {
        *this = *this / other;
        return *this;
    }

    Complex<N> &operator*=(const N &scalar)
    {
        real = ::mathematics::multiply(real, scalar);
        imag = ::mathematics::multiply(imag, scalar);
        return *this;
    }

    Complex<N> &operator/=(const N &scalar)
    {
        real = ::mathematics::divide(real, scalar);
        imag = ::mathematics::divide(imag, scalar);
        return *this;
    }

    Complex<N> operator-() const
    {
        return Complex<N>(::mathematics::negative(real), ::mathematics::negative(imag));
    }

    bool operator==(const Complex<N> &other) const
    {
        return ::mathematics::absolute(::mathematics::subtract(real, other.real)) < N(1e-9) && ::mathematics::absolute(::mathematics::subtract(imag, other.imag)) < N(1e-9);
    }

    bool operator!=(const Complex<N> &other) const { return !(*this == other); }

    bool operator<(const Complex<N> &other) const
    {
        N magnitude1 = module();
        N magnitude2 = other.module();
        if (::mathematics::absolute(::mathematics::subtract(magnitude1, magnitude2)) >= N(1e-9))
        {
            return magnitude1 < magnitude2;
        }
        if (::mathematics::absolute(::mathematics::subtract(real, other.real)) >= N(1e-9))
        {
            return real < other.real;
        }
        return imag < other.imag;
    }

    bool operator<=(const Complex<N> &other) const { return (*this < other) || (*this == other); }
    bool operator>(const Complex<N> &other) const { return !(*this <= other); }
    bool operator>=(const Complex<N> &other) const { return !(*this < other); }

    N &operator[](std::size_t index)
    {
        switch (index)
        {
        case 0:
            return real;
        case 1:
            return imag;
        default:
            throw std::out_of_range("Index out of range");
        }
    }

    const N &operator[](std::size_t index) const
    {
        switch (index)
        {
        case 0:
            return real;
        case 1:
            return imag;
        default:
            throw std::out_of_range("Index out of range");
        }
    }

    friend Complex<N> operator*(const N &scalar, const Complex<N> &complexValue) { return complexValue * scalar; }
    friend Complex<N> operator/(const N &scalar, const Complex<N> &complexValue) { return Complex<N>(scalar, N(0)) / complexValue; }
    friend Complex<N> operator+(const N &scalar, const Complex<N> &complexValue) { return Complex<N>(::mathematics::add(scalar, complexValue.real), complexValue.imag); }
    friend Complex<N> operator-(const N &scalar, const Complex<N> &complexValue) { return Complex<N>(::mathematics::subtract(scalar, complexValue.real), ::mathematics::negative(complexValue.imag)); }

    friend std::ostream &operator<<(std::ostream &output, const Complex<N> &complexValue)
    {
        output << "(" << complexValue.real;
        if (complexValue.imag >= N(0))
        {
            output << " + " << complexValue.imag << "i)";
        }
        else
        {
            output << " - " << ::mathematics::absolute(complexValue.imag) << "i)";
        }
        return output;
    }
};

template <typename N>
Complex<N> polarComplex(const N &magnitude, const N &angle)
{
    return Complex<N>(::mathematics::multiply(magnitude, ::mathematics::cosine(angle)), ::mathematics::multiply(magnitude, ::mathematics::sine(angle)));
}

template <typename N>
class Quaternion
{
  private:
    N real;
    N imagI;
    N imagJ;
    N imagK;

  public:
    Quaternion() : real(N(1)), imagI(N(0)), imagJ(N(0)), imagK(N(0)) {}

    explicit Quaternion(const N &valueReal) : real(valueReal), imagI(N(0)), imagJ(N(0)), imagK(N(0)) {}

    explicit Quaternion(const N &valueReal, const N &valueImagI) : real(valueReal), imagI(valueImagI), imagJ(N(0)), imagK(N(0)) {}

    explicit Quaternion(const N &valueReal, const N &valueImagI, const N &valueImagJ) : real(valueReal), imagI(valueImagI), imagJ(valueImagJ), imagK(N(0)) {}

    explicit Quaternion(const N &valueReal, const N &valueImagI, const N &valueImagJ, const N &valueImagK) : real(valueReal), imagI(valueImagI), imagJ(valueImagJ), imagK(valueImagK) {}

    explicit Quaternion(const std::initializer_list<N> &list)
    {
        switch (list.size())
        {
        case 0:
            real = N(1);
            imagI = N(0);
            imagJ = N(0);
            imagK = N(0);
            break;
        case 1:
            real = *(list.begin());
            imagI = N(0);
            imagJ = N(0);
            imagK = N(0);
            break;
        case 2:
            real = *(list.begin());
            imagI = *(list.begin() + 1);
            imagJ = N(0);
            imagK = N(0);
            break;
        case 3:
            real = *(list.begin());
            imagI = *(list.begin() + 1);
            imagJ = *(list.begin() + 2);
            imagK = N(0);
            break;
        default:
            real = *(list.begin());
            imagI = *(list.begin() + 1);
            imagJ = *(list.begin() + 2);
            imagK = *(list.begin() + 3);
            break;
        }
    }

    const N &getReal() const { return real; }
    N &getReal() { return real; }
    const N &getImagI() const { return imagI; }
    N &getImagI() { return imagI; }
    const N &getImagJ() const { return imagJ; }
    N &getImagJ() { return imagJ; }
    const N &getImagK() const { return imagK; }
    N &getImagK() { return imagK; }

    N scalar() const { return real; }
    Quaternion<N> vector() const { return Quaternion<N>(N(0), imagI, imagJ, imagK); }

    Quaternion<N> exact(N epsilon = N(::mathematics::tolerance)) const
    {
        return Quaternion<N>(::mathematics::exact(real, epsilon), ::mathematics::exact(imagI, epsilon), ::mathematics::exact(imagJ, epsilon), ::mathematics::exact(imagK, epsilon));
    }

    Quaternion<N> cross(const Quaternion<N> &other) const
    {
        return Quaternion<N>(N(0),
                             ::mathematics::subtract(::mathematics::multiply(imagJ, other.imagK), ::mathematics::multiply(imagK, other.imagJ)),
                             ::mathematics::subtract(::mathematics::multiply(imagK, other.imagI), ::mathematics::multiply(imagI, other.imagK)),
                             ::mathematics::subtract(::mathematics::multiply(imagI, other.imagJ), ::mathematics::multiply(imagJ, other.imagI)));
    }

    N module() const
    {
        N maximumValue = std::max({::mathematics::absolute(real), ::mathematics::absolute(imagI), ::mathematics::absolute(imagJ), ::mathematics::absolute(imagK)});
        if (maximumValue < N(1e-11))
        {
            return N(0);
        }
        N scale = ::mathematics::divide(N(1), maximumValue);
        N scaledReal = ::mathematics::multiply(real, scale);
        N scaledImagI = ::mathematics::multiply(imagI, scale);
        N scaledImagJ = ::mathematics::multiply(imagJ, scale);
        N scaledImagK = ::mathematics::multiply(imagK, scale);
        N sumSquares = ::mathematics::add(::mathematics::add(::mathematics::add(::mathematics::multiply(scaledReal, scaledReal), ::mathematics::multiply(scaledImagI, scaledImagI)), ::mathematics::multiply(scaledImagJ, scaledImagJ)), ::mathematics::multiply(scaledImagK, scaledImagK));
        return ::mathematics::multiply(::mathematics::squareRoot(sumSquares), maximumValue);
    }

    Quaternion<N> pure() const
    {
        N imaginaryNorm = ::mathematics::squareRoot(::mathematics::add(::mathematics::add(::mathematics::multiply(imagI, imagI), ::mathematics::multiply(imagJ, imagJ)), ::mathematics::multiply(imagK, imagK)));
        if (imaginaryNorm < N(1e-11))
        {
            return Quaternion<N>(N(0), N(0), N(0), N(1));
        }
        N scale = ::mathematics::divide(N(1), imaginaryNorm);
        return Quaternion<N>(N(0), ::mathematics::multiply(imagI, scale), ::mathematics::multiply(imagJ, scale), ::mathematics::multiply(imagK, scale));
    }

    Quaternion<N> conjugate() const
    {
        return Quaternion<N>(real, ::mathematics::negative(imagI), ::mathematics::negative(imagJ), ::mathematics::negative(imagK));
    }

    Quaternion<N> inverse() const
    {
        N magnitudeSquared = ::mathematics::add(::mathematics::add(::mathematics::add(::mathematics::multiply(real, real), ::mathematics::multiply(imagI, imagI)), ::mathematics::multiply(imagJ, imagJ)), ::mathematics::multiply(imagK, imagK));
        if (magnitudeSquared < N(1e-11))
        {
            return Quaternion<N>(N(1), N(0), N(0), N(0));
        }
        N scale = ::mathematics::divide(N(1), magnitudeSquared);
        return Quaternion<N>(::mathematics::multiply(real, scale), ::mathematics::multiply(::mathematics::negative(imagI), scale), ::mathematics::multiply(::mathematics::negative(imagJ), scale), ::mathematics::multiply(::mathematics::negative(imagK), scale));
    }

    Quaternion<N> normalize() const
    {
        N magnitude = this->module();
        if (magnitude < N(1e-11))
        {
            return Quaternion<N>(N(1), N(0), N(0), N(0));
        }
        N scale = ::mathematics::divide(N(1), magnitude);
        return Quaternion<N>(::mathematics::multiply(real, scale), ::mathematics::multiply(imagI, scale), ::mathematics::multiply(imagJ, scale), ::mathematics::multiply(imagK, scale));
    }

    Quaternion<N> squareRoot() const
    {
        N magnitude = this->module();
        if (magnitude < N(1e-11))
        {
            return Quaternion<N>(N(0), N(0), N(0), N(0));
        }
        N scalarPart = real;
        N denominator = ::mathematics::add(magnitude, scalarPart);
        N factor = ::mathematics::divide(N(1), denominator);
        N half = N(0.5);
        N sqrtMagnitude = ::mathematics::squareRoot(::mathematics::multiply(magnitude, half));
        return Quaternion<N>(::mathematics::squareRoot(::mathematics::multiply(::mathematics::add(magnitude, scalarPart), half)),
                             ::mathematics::multiply(::mathematics::multiply(imagI, factor), sqrtMagnitude),
                             ::mathematics::multiply(::mathematics::multiply(imagJ, factor), sqrtMagnitude),
                             ::mathematics::multiply(::mathematics::multiply(imagK, factor), sqrtMagnitude));
    }

    Quaternion<N> lerp(const Quaternion<N> &other, N t) const
    {
        N oneMinusT = ::mathematics::subtract(N(1), t);
        return Quaternion<N>(::mathematics::add(::mathematics::multiply(real, oneMinusT), ::mathematics::multiply(other.real, t)),
                             ::mathematics::add(::mathematics::multiply(imagI, oneMinusT), ::mathematics::multiply(other.imagI, t)),
                             ::mathematics::add(::mathematics::multiply(imagJ, oneMinusT), ::mathematics::multiply(other.imagJ, t)),
                             ::mathematics::add(::mathematics::multiply(imagK, oneMinusT), ::mathematics::multiply(other.imagK, t)));
    }

    Quaternion<N> nlerp(const Quaternion<N> &other, N t) const
    {
        Quaternion<N> interpolated = lerp(other, t);
        N magnitude = interpolated.module();
        if (magnitude < N(1e-11))
        {
            return *this;
        }
        N scale = ::mathematics::divide(N(1), magnitude);
        return Quaternion<N>(::mathematics::multiply(interpolated.real, scale), ::mathematics::multiply(interpolated.imagI, scale), ::mathematics::multiply(interpolated.imagJ, scale), ::mathematics::multiply(interpolated.imagK, scale));
    }

    N dot(const Quaternion<N> &other) const
    {
        return ::mathematics::add(::mathematics::add(::mathematics::add(::mathematics::multiply(real, other.real), ::mathematics::multiply(imagI, other.imagI)), ::mathematics::multiply(imagJ, other.imagJ)), ::mathematics::multiply(imagK, other.imagK));
    }

    N distance(const Quaternion<N> &other) const
    {
        Quaternion<N> diff = (*this) - other;
        return diff.module();
    }

    N angle(const Quaternion<N> &other) const
    {
        N moduleA = this->module();
        N moduleB = other.module();
        if (moduleA < N(1e-11) || moduleB < N(1e-11))
        {
            return N(0);
        }
        Quaternion<N> normalizedA = *this / moduleA;
        Quaternion<N> normalizedB = other / moduleB;
        N dotProduct = ::mathematics::absolute(normalizedA.dot(normalizedB));
        if (dotProduct > N(1))
        {
            dotProduct = N(1);
        }
        N oneMinusDotSquared = ::mathematics::subtract(N(1), ::mathematics::multiply(dotProduct, dotProduct));
        return ::mathematics::multiply(N(2), ::mathematics::arctangent2(::mathematics::squareRoot(oneMinusDotSquared), dotProduct));
    }

    std::tuple<std::complex<N>, std::complex<N>, std::complex<N>, std::complex<N>, std::complex<N>, std::complex<N>> toComplex() const
    {
        return {std::complex<N>(real, imagI), std::complex<N>(real, imagJ), std::complex<N>(real, imagK),
                std::complex<N>(imagI, imagJ), std::complex<N>(imagI, imagK), std::complex<N>(imagJ, imagK)};
    }

    Quaternion<N> exponential() const
    {
        N imaginaryNorm = ::mathematics::squareRoot(::mathematics::add(::mathematics::add(::mathematics::multiply(imagI, imagI), ::mathematics::multiply(imagJ, imagJ)), ::mathematics::multiply(imagK, imagK)));
        N expReal = ::mathematics::exponential(real);
        if (imaginaryNorm < N(1e-11))
        {
            return Quaternion<N>(expReal, N(0), N(0), N(0));
        }
        N sinValue = ::mathematics::sine(imaginaryNorm);
        N cosValue = ::mathematics::cosine(imaginaryNorm);
        N factor = ::mathematics::divide(sinValue, imaginaryNorm);
        return Quaternion<N>(::mathematics::multiply(expReal, cosValue), ::mathematics::multiply(::mathematics::multiply(expReal, factor), imagI), ::mathematics::multiply(::mathematics::multiply(expReal, factor), imagJ), ::mathematics::multiply(::mathematics::multiply(expReal, factor), imagK));
    }

    Quaternion<N> logarithm() const
    {
        N magnitude = this->module();
        if (magnitude < N(1e-11))
        {
            return Quaternion<N>(N(0), N(0), N(0), N(0));
        }
        N imaginaryNorm = ::mathematics::squareRoot(::mathematics::add(::mathematics::add(::mathematics::multiply(imagI, imagI), ::mathematics::multiply(imagJ, imagJ)), ::mathematics::multiply(imagK, imagK)));
        N logMagnitude = ::mathematics::logarithm(magnitude);
        if (imaginaryNorm < N(1e-11))
        {
            if (real < N(0))
            {
                return Quaternion<N>(logMagnitude, N(::mathematics::pi), N(0), N(0));
            }
            return Quaternion<N>(logMagnitude, N(0), N(0), N(0));
        }
        N theta = ::mathematics::arctangent2(imaginaryNorm, real);
        N factor = ::mathematics::divide(theta, imaginaryNorm);
        return Quaternion<N>(logMagnitude, ::mathematics::multiply(factor, imagI), ::mathematics::multiply(factor, imagJ), ::mathematics::multiply(factor, imagK));
    }

    Quaternion<N> logarithm(const Quaternion<N> &base) const
    {
        return this->logarithm() * base.logarithm().inverse();
    }

    Quaternion<N> power(N t) const
    {
        N magnitude = this->module();
        if (magnitude < N(1e-11))
        {
            throw std::domain_error("power: zero quaternion cannot be raised to a power");
        }
        N imaginaryNorm = ::mathematics::squareRoot(::mathematics::add(::mathematics::add(::mathematics::multiply(imagI, imagI), ::mathematics::multiply(imagJ, imagJ)), ::mathematics::multiply(imagK, imagK)));
        N logMagnitude = ::mathematics::logarithm(magnitude);
        N newReal = ::mathematics::multiply(logMagnitude, t);
        if (imaginaryNorm < N(1e-11))
        {
            if (real < N(0))
            {
                N adjustedAngle = ::mathematics::multiply(N(::mathematics::pi), t);
                return Quaternion<N>(::mathematics::multiply(::mathematics::exponential(newReal), ::mathematics::cosine(adjustedAngle)), ::mathematics::multiply(::mathematics::exponential(newReal), ::mathematics::sine(adjustedAngle)), N(0), N(0));
            }
            return Quaternion<N>(::mathematics::exponential(newReal), N(0), N(0), N(0));
        }
        N theta = ::mathematics::arctangent2(imaginaryNorm, real);
        N newAngle = ::mathematics::multiply(theta, t);
        N factor = ::mathematics::divide(newAngle, imaginaryNorm);
        N expReal = ::mathematics::exponential(newReal);
        return Quaternion<N>(::mathematics::multiply(expReal, ::mathematics::cosine(newAngle)), ::mathematics::multiply(::mathematics::multiply(expReal, factor), imagI), ::mathematics::multiply(::mathematics::multiply(expReal, factor), imagJ), ::mathematics::multiply(::mathematics::multiply(expReal, factor), imagK));
    }

    Quaternion<N> slerp(const Quaternion<N> &other, N t) const
    {
        N dotProduct = this->dot(other);
        Quaternion<N> target = other;
        if (dotProduct < N(0))
        {
            dotProduct = ::mathematics::negative(dotProduct);
            target = -target;
        }
        if (dotProduct > N(1))
        {
            dotProduct = N(1);
        }
        N oneMinusDotSquared = ::mathematics::subtract(N(1), ::mathematics::multiply(dotProduct, dotProduct));
        N theta = ::mathematics::arctangent2(::mathematics::squareRoot(oneMinusDotSquared), dotProduct);
        N sinTheta = ::mathematics::sine(theta);
        if (::mathematics::absolute(sinTheta) < N(1e-11))
        {
            return lerp(target, t);
        }
        N factorA = ::mathematics::divide(::mathematics::sine(::mathematics::multiply(::mathematics::subtract(N(1), t), theta)), sinTheta);
        N factorB = ::mathematics::divide(::mathematics::sine(::mathematics::multiply(t, theta)), sinTheta);
        return Quaternion<N>(::mathematics::add(::mathematics::multiply(real, factorA), ::mathematics::multiply(target.real, factorB)),
                             ::mathematics::add(::mathematics::multiply(imagI, factorA), ::mathematics::multiply(target.imagI, factorB)),
                             ::mathematics::add(::mathematics::multiply(imagJ, factorA), ::mathematics::multiply(target.imagJ, factorB)),
                             ::mathematics::add(::mathematics::multiply(imagK, factorA), ::mathematics::multiply(target.imagK, factorB)));
    }

    Quaternion<N> around(const Quaternion<N> &rotation) const
    {
        return rotation * (*this) * rotation.inverse();
    }

    Quaternion<N> rotate(const Quaternion<N> &v) const { return (*this) * v * inverse(); }

    Quaternion<N> rotate(const std::initializer_list<Quaternion<N>> &rotations) const
    {
        Quaternion<N> result = *this;
        for (const auto &r : rotations)
        {
            result = r * result * r.inverse();
        }
        return result;
    }

    Quaternion<N> rotate(const std::vector<Quaternion<N>> &rotations) const
    {
        Quaternion<N> result = *this;
        for (const auto &r : rotations)
        {
            result = r * result * r.inverse();
        }
        return result;
    }

    Quaternion<N> derivative(const Quaternion<N> &angularVelocity) const
    {
        return Quaternion<N>(N(0.5)) * angularVelocity * (*this);
    }

    Quaternion<N> integrate(const Quaternion<N> &angularVelocity, N timeStep) const
    {
        Quaternion<N> halfAngularVelocity = angularVelocity * N(0.5) * timeStep;
        return (*this) * halfAngularVelocity.exponential();
    }

    Quaternion<N> &operator=(const Quaternion<N> &) = default;

    Quaternion<N> operator+(const Quaternion<N> &q) const
    {
        return Quaternion<N>(::mathematics::add(real, q.real), ::mathematics::add(imagI, q.imagI), ::mathematics::add(imagJ, q.imagJ), ::mathematics::add(imagK, q.imagK));
    }

    Quaternion<N> operator-(const Quaternion<N> &q) const
    {
        return Quaternion<N>(::mathematics::subtract(real, q.real), ::mathematics::subtract(imagI, q.imagI), ::mathematics::subtract(imagJ, q.imagJ), ::mathematics::subtract(imagK, q.imagK));
    }

    Quaternion<N> operator*(const N &scalar) const
    {
        return Quaternion<N>(::mathematics::multiply(real, scalar), ::mathematics::multiply(imagI, scalar), ::mathematics::multiply(imagJ, scalar), ::mathematics::multiply(imagK, scalar));
    }

    Quaternion<N> operator/(const N &scalar) const
    {
        return Quaternion<N>(::mathematics::divide(real, scalar), ::mathematics::divide(imagI, scalar), ::mathematics::divide(imagJ, scalar), ::mathematics::divide(imagK, scalar));
    }

    Quaternion<N> operator*(const Quaternion<N> &q) const
    {
        return Quaternion<N>(
            ::mathematics::subtract(::mathematics::subtract(::mathematics::subtract(::mathematics::multiply(real, q.real), ::mathematics::multiply(imagI, q.imagI)), ::mathematics::multiply(imagJ, q.imagJ)), ::mathematics::multiply(imagK, q.imagK)),
            ::mathematics::add(::mathematics::add(::mathematics::subtract(::mathematics::multiply(real, q.imagI), ::mathematics::multiply(imagK, q.imagJ)), ::mathematics::multiply(imagI, q.real)), ::mathematics::multiply(imagJ, q.imagK)),
            ::mathematics::add(::mathematics::subtract(::mathematics::add(::mathematics::multiply(real, q.imagJ), ::mathematics::multiply(imagK, q.imagI)), ::mathematics::multiply(imagI, q.imagK)), ::mathematics::multiply(imagJ, q.real)),
            ::mathematics::add(::mathematics::subtract(::mathematics::add(::mathematics::multiply(real, q.imagK), ::mathematics::multiply(imagI, q.imagJ)), ::mathematics::multiply(imagJ, q.imagI)), ::mathematics::multiply(imagK, q.real)));
    }

    Quaternion<N> operator/(const Quaternion<N> &q) const
    {
        return (*this) * q.inverse();
    }

    Quaternion<N> &operator+=(const Quaternion<N> &q)
    {
        real = ::mathematics::add(real, q.real);
        imagI = ::mathematics::add(imagI, q.imagI);
        imagJ = ::mathematics::add(imagJ, q.imagJ);
        imagK = ::mathematics::add(imagK, q.imagK);
        return *this;
    }

    Quaternion<N> &operator-=(const Quaternion<N> &q)
    {
        real = ::mathematics::subtract(real, q.real);
        imagI = ::mathematics::subtract(imagI, q.imagI);
        imagJ = ::mathematics::subtract(imagJ, q.imagJ);
        imagK = ::mathematics::subtract(imagK, q.imagK);
        return *this;
    }

    Quaternion<N> &operator*=(const Quaternion<N> &q)
    {
        *this = *this * q;
        return *this;
    }

    Quaternion<N> &operator*=(const N &scalar)
    {
        real = ::mathematics::multiply(real, scalar);
        imagI = ::mathematics::multiply(imagI, scalar);
        imagJ = ::mathematics::multiply(imagJ, scalar);
        imagK = ::mathematics::multiply(imagK, scalar);
        return *this;
    }

    Quaternion<N> &operator/=(const N &scalar)
    {
        real = ::mathematics::divide(real, scalar);
        imagI = ::mathematics::divide(imagI, scalar);
        imagJ = ::mathematics::divide(imagJ, scalar);
        imagK = ::mathematics::divide(imagK, scalar);
        return *this;
    }

    Quaternion<N> operator-() const
    {
        return Quaternion<N>(::mathematics::negative(real), ::mathematics::negative(imagI), ::mathematics::negative(imagJ), ::mathematics::negative(imagK));
    }

    bool operator==(const Quaternion<N> &q) const
    {
        return ::mathematics::absolute(::mathematics::subtract(real, q.real)) < N(1e-9) && ::mathematics::absolute(::mathematics::subtract(imagI, q.imagI)) < N(1e-9) && ::mathematics::absolute(::mathematics::subtract(imagJ, q.imagJ)) < N(1e-9) && ::mathematics::absolute(::mathematics::subtract(imagK, q.imagK)) < N(1e-9);
    }

    bool operator!=(const Quaternion<N> &q) const { return !(*this == q); }

    bool operator<(const Quaternion<N> &q) const
    {
        N magnitude1 = module();
        N magnitude2 = q.module();
        if (::mathematics::absolute(::mathematics::subtract(magnitude1, magnitude2)) >= N(1e-9))
            return magnitude1 < magnitude2;
        if (::mathematics::absolute(::mathematics::subtract(real, q.real)) >= N(1e-9))
            return real < q.real;
        if (::mathematics::absolute(::mathematics::subtract(imagI, q.imagI)) >= N(1e-9))
            return imagI < q.imagI;
        if (::mathematics::absolute(::mathematics::subtract(imagJ, q.imagJ)) >= N(1e-9))
            return imagJ < q.imagJ;
        return imagK < q.imagK;
    }

    bool operator<=(const Quaternion<N> &q) const { return (*this < q) || (*this == q); }
    bool operator>(const Quaternion<N> &q) const { return !(*this <= q); }
    bool operator>=(const Quaternion<N> &q) const { return !(*this < q); }

    N &operator[](std::size_t index)
    {
        switch (index)
        {
        case 0:
            return real;
        case 1:
            return imagI;
        case 2:
            return imagJ;
        case 3:
            return imagK;
        default:
            throw std::out_of_range("Index out of range");
        }
    }

    const N &operator[](std::size_t index) const
    {
        switch (index)
        {
        case 0:
            return real;
        case 1:
            return imagI;
        case 2:
            return imagJ;
        case 3:
            return imagK;
        default:
            throw std::out_of_range("Index out of range");
        }
    }

    friend Quaternion<N> operator*(const N &scalar, const Quaternion<N> &q) { return q * scalar; }
    friend Quaternion<N> operator/(const N &scalar, const Quaternion<N> &q) { return scalar * q.inverse(); }
    friend Quaternion<N> operator+(const N &scalar, const Quaternion<N> &q) { return Quaternion<N>(::mathematics::add(scalar, q.real), q.imagI, q.imagJ, q.imagK); }
    friend Quaternion<N> operator-(const N &scalar, const Quaternion<N> &q) { return Quaternion<N>(::mathematics::subtract(scalar, q.real), ::mathematics::negative(q.imagI), ::mathematics::negative(q.imagJ), ::mathematics::negative(q.imagK)); }

    friend std::ostream &operator<<(std::ostream &output, const Quaternion &q)
    {
        output << "(" << q.real << ", " << q.imagI << ", " << q.imagJ << ", " << q.imagK << ")";
        return output;
    }
};

template <typename N>
inline bool abnormal(const Complex<N> &value)
{
    return ::mathematics::abnormal(value.getReal()) || ::mathematics::abnormal(value.getImag());
}

template <typename N>
inline bool abnormal(const Quaternion<N> &value)
{
    return ::mathematics::abnormal(value.getReal()) || ::mathematics::abnormal(value.getImagI()) || ::mathematics::abnormal(value.getImagJ()) || ::mathematics::abnormal(value.getImagK());
}

template <typename N>
inline bool abnormal(const std::complex<N> &value)
{
    return ::mathematics::abnormal(value.real()) || ::mathematics::abnormal(value.imag());
}

template <typename N>
N module(const Complex<N> &value) { return value.module(); }

template <typename N>
N module(const Quaternion<N> &value) { return value.module(); }

template <typename N>
N module(const std::complex<N> &value)
{
    N maximumValue = std::max(::mathematics::absolute(value.real()), ::mathematics::absolute(value.imag()));
    if (maximumValue < N(1e-11))
    {
        return N(0);
    }
    N scale = ::mathematics::divide(N(1), maximumValue);
    N scaledReal = ::mathematics::multiply(value.real(), scale);
    N scaledImag = ::mathematics::multiply(value.imag(), scale);
    N sumSquares = ::mathematics::add(::mathematics::multiply(scaledReal, scaledReal), ::mathematics::multiply(scaledImag, scaledImag));
    return ::mathematics::multiply(::mathematics::squareRoot(sumSquares), maximumValue);
}

template <typename N>
Complex<N> conjugate(const Complex<N> &value) { return value.conjugate(); }

template <typename N>
Quaternion<N> conjugate(const Quaternion<N> &value) { return value.conjugate(); }

template <typename N>
std::complex<N> conjugate(const std::complex<N> &value)
{
    return std::complex<N>(value.real(), ::mathematics::negative(value.imag()));
}

template <typename N>
Complex<N> inverse(const Complex<N> &value) { return value.inverse(); }

template <typename N>
Quaternion<N> inverse(const Quaternion<N> &value) { return value.inverse(); }

template <typename N>
std::complex<N> inverse(const std::complex<N> &value)
{
    N denominator = ::mathematics::add(::mathematics::multiply(value.real(), value.real()), ::mathematics::multiply(value.imag(), value.imag()));
    if (denominator < N(1e-11))
    {
        return std::complex<N>(N(1), N(0));
    }
    N scale = ::mathematics::divide(N(1), denominator);
    return std::complex<N>(::mathematics::multiply(value.real(), scale), ::mathematics::multiply(::mathematics::negative(value.imag()), scale));
}

template <typename N>
Complex<N> exact(const Complex<N> &value, N epsilon = N(::mathematics::tolerance)) { return value.exact(epsilon); }

template <typename N>
Quaternion<N> exact(const Quaternion<N> &value, N epsilon = N(::mathematics::tolerance)) { return value.exact(epsilon); }

template <typename N>
std::complex<N> exact(const std::complex<N> &value, N epsilon = N(::mathematics::tolerance))
{
    return std::complex<N>(::mathematics::exact(value.real(), epsilon), ::mathematics::exact(value.imag(), epsilon));
}

template <typename N>
Complex<N> limit(const Complex<N> &value, N start, N end) { return value.limit(start, end); }

template <typename N>
Quaternion<N> limit(const Quaternion<N> &value, N start, N end)
{
    return Quaternion<N>(::mathematics::limit(value.getReal(), start, end), ::mathematics::limit(value.getImagI(), start, end), ::mathematics::limit(value.getImagJ(), start, end), ::mathematics::limit(value.getImagK(), start, end));
}

template <typename N>
std::complex<N> limit(const std::complex<N> &value, N start, N end)
{
    return std::complex<N>(::mathematics::limit(value.real(), start, end), ::mathematics::limit(value.imag(), start, end));
}

template <typename N>
Complex<N> smooth(const Complex<N> &value, N threshold, N steepness) { return value.smooth(threshold, steepness); }

template <typename N>
Quaternion<N> smooth(const Quaternion<N> &value, N threshold, N steepness)
{
    return Quaternion<N>(::mathematics::smooth(value.getReal(), threshold, steepness), ::mathematics::smooth(value.getImagI(), threshold, steepness), ::mathematics::smooth(value.getImagJ(), threshold, steepness), ::mathematics::smooth(value.getImagK(), threshold, steepness));
}

template <typename N>
std::complex<N> smooth(const std::complex<N> &value, N threshold, N steepness)
{
    return std::complex<N>(::mathematics::smooth(value.real(), threshold, steepness), ::mathematics::smooth(value.imag(), threshold, steepness));
}

template <typename N>
Complex<N> negative(const Complex<N> &value) { return -value; }

template <typename N>
Quaternion<N> negative(const Quaternion<N> &value) { return -value; }

template <typename N>
std::complex<N> negative(const std::complex<N> &value)
{
    return std::complex<N>(::mathematics::negative(value.real()), ::mathematics::negative(value.imag()));
}

template <typename N>
Complex<N> exponential(const Complex<N> &value) { return value.exponential(); }

template <typename N>
Quaternion<N> exponential(const Quaternion<N> &value) { return value.exponential(); }

template <typename N>
std::complex<N> exponential(const std::complex<N> &value)
{
    N expReal = ::mathematics::exponential(value.real());
    return std::complex<N>(::mathematics::multiply(expReal, ::mathematics::cosine(value.imag())), ::mathematics::multiply(expReal, ::mathematics::sine(value.imag())));
}

template <typename N>
Complex<N> logarithm(const Complex<N> &value) { return value.logarithm(); }

template <typename N>
Quaternion<N> logarithm(const Quaternion<N> &value) { return value.logarithm(); }

template <typename N>
std::complex<N> logarithm(const std::complex<N> &value)
{
    return std::complex<N>(::mathematics::logarithm(module(value)), ::mathematics::arctangent2(value.imag(), value.real()));
}

template <typename N>
Complex<N> power(const Complex<N> &base, const Complex<N> &exponent) { return base.power(exponent); }

template <typename N>
Complex<N> power(const Complex<N> &base, const N &exponent) { return base.power(exponent); }

template <typename N>
Quaternion<N> power(const Quaternion<N> &base, N exponent) { return base.power(exponent); }

template <typename N>
std::complex<N> power(const std::complex<N> &base, const N &exponent)
{
    N magnitude = module(base);
    N arg = ::mathematics::arctangent2(base.imag(), base.real());
    N newMagnitude = ::mathematics::power(magnitude, exponent);
    N newArg = ::mathematics::multiply(arg, exponent);
    return std::complex<N>(::mathematics::multiply(newMagnitude, ::mathematics::cosine(newArg)), ::mathematics::multiply(newMagnitude, ::mathematics::sine(newArg)));
}

template <typename N>
std::complex<N> power(const std::complex<N> &base, const std::complex<N> &exponent)
{
    return exponential(logarithm(base) * exponent);
}

template <typename N>
Complex<N> squareRoot(const Complex<N> &value) { return value.squareRoot(); }

template <typename N>
Quaternion<N> squareRoot(const Quaternion<N> &value) { return value.squareRoot(); }

template <typename N>
std::complex<N> squareRoot(const std::complex<N> &value)
{
    N magnitude = module(value);
    N arg = ::mathematics::arctangent2(value.imag(), value.real());
    N sqrtMagnitude = ::mathematics::squareRoot(magnitude);
    N halfArg = ::mathematics::divide(arg, N(2));
    return std::complex<N>(::mathematics::multiply(sqrtMagnitude, ::mathematics::cosine(halfArg)), ::mathematics::multiply(sqrtMagnitude, ::mathematics::sine(halfArg)));
}

template <typename N>
Complex<N> sine(const Complex<N> &value) { return value.sine(); }

template <typename N>
std::complex<N> sine(const std::complex<N> &value)
{
    N sinReal = ::mathematics::sine(value.real());
    N cosReal = ::mathematics::cosine(value.real());
    N sinhImag = ::mathematics::hyperbolicSine(value.imag());
    N coshImag = ::mathematics::hyperbolicCosine(value.imag());
    return std::complex<N>(::mathematics::multiply(sinReal, coshImag), ::mathematics::multiply(cosReal, sinhImag));
}

template <typename N>
Complex<N> cosine(const Complex<N> &value) { return value.cosine(); }

template <typename N>
std::complex<N> cosine(const std::complex<N> &value)
{
    N cosReal = ::mathematics::cosine(value.real());
    N sinReal = ::mathematics::sine(value.real());
    N sinhImag = ::mathematics::hyperbolicSine(value.imag());
    N coshImag = ::mathematics::hyperbolicCosine(value.imag());
    return std::complex<N>(::mathematics::multiply(cosReal, coshImag), ::mathematics::multiply(::mathematics::negative(sinReal), sinhImag));
}

template <typename N>
Complex<N> tangent(const Complex<N> &value) { return value.tangent(); }

template <typename N>
std::complex<N> tangent(const std::complex<N> &value)
{
    return sine(value) / cosine(value);
}

template <typename N>
Complex<N> arcsine(const Complex<N> &value) { return value.arcsine(); }

template <typename N>
std::complex<N> arcsine(const std::complex<N> &value)
{
    std::complex<N> i(N(0), N(1));
    return negative(i) * logarithm(i * value + squareRoot(std::complex<N>(N(1), N(0)) - value * value));
}

template <typename N>
Complex<N> arccosine(const Complex<N> &value) { return value.arccosine(); }

template <typename N>
std::complex<N> arccosine(const std::complex<N> &value)
{
    std::complex<N> i(N(0), N(1));
    return negative(i) * logarithm(value + i * squareRoot(std::complex<N>(N(1), N(0)) - value * value));
}

template <typename N>
Complex<N> arctangent(const Complex<N> &value) { return value.arctangent(); }

template <typename N>
std::complex<N> arctangent(const std::complex<N> &value)
{
    std::complex<N> i(N(0), N(1));
    return ::mathematics::multiply(i, N(0.5)) * (logarithm(std::complex<N>(N(1), N(0)) - i * value) - logarithm(std::complex<N>(N(1), N(0)) + i * value));
}

template <typename N>
N dot(const Complex<N> &left, const Complex<N> &right) { return left.dot(right); }

template <typename N>
N dot(const Quaternion<N> &left, const Quaternion<N> &right) { return left.dot(right); }

template <typename N>
N dot(const std::complex<N> &left, const std::complex<N> &right)
{
    return ::mathematics::add(::mathematics::multiply(left.real(), right.real()), ::mathematics::multiply(left.imag(), right.imag()));
}

template <typename N>
N distance(const Complex<N> &left, const Complex<N> &right) { return left.distance(right); }

template <typename N>
N distance(const Quaternion<N> &left, const Quaternion<N> &right) { return left.distance(right); }

template <typename N>
N distance(const std::complex<N> &left, const std::complex<N> &right)
{
    return module(left - right);
}

template <typename N>
N angle(const Complex<N> &left, const Complex<N> &right) { return left.angle(right); }

template <typename N>
N angle(const Quaternion<N> &left, const Quaternion<N> &right) { return left.angle(right); }

template <typename N>
N angle(const std::complex<N> &left, const std::complex<N> &right)
{
    N moduleA = module(left);
    N moduleB = module(right);
    if (moduleA < N(1e-11) || moduleB < N(1e-11))
    {
        return N(0);
    }
    std::complex<N> normalizedA = left / moduleA;
    std::complex<N> normalizedB = right / moduleB;
    N dotProduct = ::mathematics::absolute(dot(normalizedA, normalizedB));
    if (dotProduct > N(1))
    {
        dotProduct = N(1);
    }
    N oneMinusDotSquared = ::mathematics::subtract(N(1), ::mathematics::multiply(dotProduct, dotProduct));
    return ::mathematics::multiply(N(2), ::mathematics::arctangent2(::mathematics::squareRoot(oneMinusDotSquared), dotProduct));
}

} // namespace mathematics

namespace std
{
template <typename N>
struct hash<mathematics::Complex<N>>
{
    size_t operator()(const mathematics::Complex<N> &complexValue) const
    {
        size_t seed = 0;
        seed ^= hash<N>{}(complexValue.getReal()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash<N>{}(complexValue.getImag()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

template <typename N>
struct less<mathematics::Complex<N>>
{
    bool operator()(const mathematics::Complex<N> &left, const mathematics::Complex<N> &right) const
    {
        if (left.getReal() != right.getReal())
        {
            return left.getReal() < right.getReal();
        }
        return left.getImag() < right.getImag();
    }
};

template <typename N>
struct hash<mathematics::Quaternion<N>>
{
    size_t operator()(const mathematics::Quaternion<N> &q) const
    {
        size_t seed = 0;
        seed ^= hash<N>{}(q.getReal()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash<N>{}(q.getImagI()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash<N>{}(q.getImagJ()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash<N>{}(q.getImagK()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

template <typename N>
struct less<mathematics::Quaternion<N>>
{
    bool operator()(const mathematics::Quaternion<N> &left, const mathematics::Quaternion<N> &right) const
    {
        if (left.getReal() != right.getReal())
        {
            return left.getReal() < right.getReal();
        }
        if (left.getImagI() != right.getImagI())
        {
            return left.getImagI() < right.getImagI();
        }
        if (left.getImagJ() != right.getImagJ())
        {
            return left.getImagJ() < right.getImagJ();
        }
        return left.getImagK() < right.getImagK();
    }
};

} // namespace std