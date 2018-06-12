/*
Copyright (c) 2017 Vincent Ducharme

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef _UNIT_H_
#define _UNIT_H_

#include <cmath>
#include <cstdint>
#include <type_traits>

using int64  = std::int64_t;
using int32  = std::int32_t;
using uint64 = std::uint64_t;
using uint32 = std::uint32_t;

////////////////////////////////////////////
// Internal utilities
////////////////////////////////////////////
namespace internals
{
    ////////////////////////////////////////////
    // PGCD (Plus grand commun diviseur)
    ////////////////////////////////////////////
    template <int32 M, int32 N>
    struct PGCD
    {
        constexpr static int32 Value = PGCD<N, M%N>::Value;
    };

    template <int32 N>
    struct PGCD<N, 0>
    {
        constexpr static int32 Value = N;
    };

    template<>
    struct PGCD<0, 0>
    {
        constexpr static int32 Value = 1;
    };


    ////////////////////////////////////////////
    // Same Types
    ////////////////////////////////////////////
    template<class T, class T1>
    struct IsSameTypes
    {
        constexpr static bool Result = false;
    };

    template<class T>
    struct IsSameTypes<T, T>
    {
        constexpr static bool Result = true;
    };


    ////////////////////////////////////////////
    // Enable methods if a condition
    ////////////////////////////////////////////
    template<bool COND, typename T = void>
    struct EnableIf
    {
        using EnableType = T;
    };

    template<typename T>
    struct EnableIf<false, T> {};


    ////////////////////////////////////////////
    // Static If Else with Int
    ////////////////////////////////////////////
    template<bool COND, int64 T, int64 F>
    struct IfElseValue;

    template<int64 T, int64 F>
    struct IfElseValue<true, T, F>
    {
        constexpr static int64 Value = T;
    };

    template<int64 T, int64 F>
    struct IfElseValue<false, T, F>
    {
        constexpr static int64 Value = F;
    };


    ////////////////////////////////////////////
    // Static If Else with Types
    ////////////////////////////////////////////
    template<bool COND, typename T, typename F>
    struct IfElseType;

    template<typename T, typename F>
    struct IfElseType<true, T, F>
    {
        using Type = T;
    };

    template<typename T, typename F>
    struct IfElseType<false, T, F>
    {
        using Type = F;
    };

    ////////////////////////////////////////////
    // Test is a type is an arithmetic type
    ////////////////////////////////////////////
    template<typename V>
    struct IsArithmetic
    {
        constexpr static bool Result = std::is_arithmetic<V>::value;
    };
}


////////////////////////////////////////////
// Forward declaration
////////////////////////////////////////////

// Create unit in the form of U^(N/D)
template<typename U, int32 N, int32 D, typename UBase>
class PoweredUnit;

// Create unit in the form of U*(N/D)
template<typename U, int64 N, int64 D, typename UBase>
class ScaledUnit;

// Create unit in the form of UV
template<typename U, typename V>
class ComposedUnit;

// Create a vector with component of type U
template<typename U>
class Vector;

// Base class for Unit
template<typename U>
class Unit;

// Wrapper of double for unitless data
class Scalar;


////////////////////////////////////////////
// Scalar - Wrapper to a double
////////////////////////////////////////////
class Scalar final
{
public:
    constexpr Scalar(double val = 0.0)
        : m_value(val)
    {
    }

    // Special case when Power of 0
    template<typename U, int32 D, typename UBase>
    constexpr Scalar(const PoweredUnit<U, 0, D, UBase>& u)
        : m_value(u.Value())
    {
    }

    constexpr Scalar& operator+=(const Scalar& u)
    {
        m_value += u.m_value;
        return *this;
    }

    constexpr Scalar& operator-=(const Scalar& u)
    {
        m_value -= u.m_value;
        return *this;
    }

    constexpr Scalar& operator*=(Scalar u)
    {
        m_value *= u.m_value;
        return *this;
    }

    template<typename F>
    constexpr typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, Scalar&>::EnableType operator*=(const F& u)
    {
        m_value *= u;
        return *this;
    }

    constexpr Scalar& operator/=(Scalar u)
    {
        m_value /= u.m_value;
        return *this;
    }

    template<typename F>
    constexpr typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, Scalar&>::EnableType operator/=(const F& u)
    {
        m_value /= u;
        return *this;
    }

    constexpr operator double() const
    {
        return m_value;
    }

    constexpr double Value() const
    {
        return m_value;
    }

    template<typename U>
    constexpr U ToUnit() const
    {
        return U(m_value);
    }

    template<typename S>
    friend S& operator >> (S& in, Scalar& s)
    {
        return (in >> s.m_value);
    }

protected:
    double m_value;
};


////////////////////////////////////////////
// Unit
////////////////////////////////////////////
template<typename U>
class Unit
{
public:
    constexpr Unit()
        : m_value()
    {
    }

    constexpr explicit Unit(double val)
        : m_value(val)
    {
    }

    constexpr Unit(const Unit&) = default;

    ///////////////////////////
    // Assignation operators
    ///////////////////////////
    constexpr U& operator=(const Unit<U>& u)
    {
        m_value = u.Value();

        return static_cast<U&>(*this);
    }

    ///////////////////////////
    // Negation operators
    ///////////////////////////
    constexpr U operator-() const
    {
        return U(-m_value);
    }

    ///////////////////////////
    // Auto in(de)crement operators
    ///////////////////////////
    constexpr U operator++(int)
    {
        U temp(m_value);
        ++m_value;
        return temp;
    }

    constexpr U& operator++()
    {
        ++m_value;
        return static_cast<U&>(*this);
    }

    constexpr U operator--(int)
    {
        U temp(m_value);
        --m_value;
        return temp;
    }

    constexpr U& operator--()
    {
        --m_value;
        return static_cast<U&>(*this);
    }

    ///////////////////////////
    // Add/Sub operators
    ///////////////////////////
    constexpr U operator-(const U& u) const
    {
        return U(m_value - u.m_value);
    }

    constexpr U operator+(const U& u) const
    {
        return U(m_value + u.m_value);
    }

    ///////////////////////////
    // Multiplication operators
    ///////////////////////////
    // Multiply by Scalar -> Unit
    constexpr U operator*(Scalar u) const
    {
        return U(m_value * u.Value());
    }

    // Multiply by an arithmetic value -> Unit
    template<typename F>
    constexpr typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, U>::EnableType operator*(const F& u) const
    {
        return U(m_value * u);
    }

    ///////////////////////////
    // Division operators
    ///////////////////////////
    // Divide by same type -> Scalar
    constexpr Scalar operator/(const U& u) const
    {
        return Scalar(m_value / u.m_value);
    }

    // Divide by Scalar -> Unit
    constexpr U operator/(Scalar u) const
    {
        return U(m_value / u.Value());
    }

    // Divide by an arithmetic value -> Unit
    template<typename F>
    constexpr typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, U>::EnableType operator/(const F& u) const
    {
        return U(m_value / u);
    }

    ///////////////////////////
    // Add/Sub equal operator
    ///////////////////////////
    constexpr U& operator+=(const U& u)
    {
        m_value += u.m_value;
        return static_cast<U&>(*this);
    }

    constexpr U& operator-=(const U& u)
    {
        m_value -= u.m_value;
        return static_cast<U&>(*this);
    }

    ///////////////////////////
    // Multiply equal operator
    ///////////////////////////
    constexpr U& operator*=(Scalar u)
    {
        m_value *= u.Value();
        return static_cast<U&>(*this);
    }

    template<typename F>
    constexpr typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, U&>::EnableType operator*=(const F& u)
    {
        m_value *= u;
        return static_cast<U&>(*this);
    }

    ///////////////////////////
    // Divide equal operator
    ///////////////////////////
    constexpr U& operator/=(Scalar u)
    {
        m_value /= u.Value();
        return static_cast<U&>(*this);
    }

    template<typename F>
    constexpr typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, U&>::EnableType operator/=(const F& u)
    {
        m_value /= u;
        return static_cast<U&>(*this);
    }

    ///////////////////////////
    // Logical operators
    ///////////////////////////
    constexpr bool operator==(const U& u) const
    {
        return m_value == u.m_value;
    }

    constexpr bool operator!=(const U& u) const
    {
        return !(*this == u);
    }

    constexpr bool operator<(const U& u) const
    {
        return m_value < u.m_value;
    }

    constexpr bool operator>(const U& u) const
    {
        return m_value > u.m_value;
    }

    constexpr bool operator<=(const U& u) const
    {
        return !(*this > u);
    }

    constexpr bool operator>=(const U& u) const
    {
        return !(*this < u);
    }

    ///////////////////////////
    // Other functions
    ///////////////////////////
    constexpr double Value() const { return m_value; }

    constexpr explicit operator double() const
    {
        return m_value;
    }

    template<typename V>
    constexpr V ToUnit() const { return V(m_value); }

    template<typename S>
    friend S& operator >> (S& is, Unit<U>& u)
    {
        return (is >> u.m_value);
    }

protected:
    double m_value;
};


// Scalar multiply unit
template<typename U>
constexpr U operator*(const Scalar d, const Unit<U>& u)
{
    return u * d;
}

// Arithmetic value multiply unit
template<typename F, typename U>
constexpr typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, U>::EnableType operator*(const F& f, const Unit<U>& u)
{
    return u * f;
}

////////////////////////////////////////////
// Internal utilities
////////////////////////////////////////////
namespace internals
{
    ////////////////////////////////////////////
    // Trait for Unit type
    ////////////////////////////////////////////
    template<typename U>
    struct IsUnit
    {
        constexpr static bool Result = false;
    };

    template<typename U, int64 SN, int64 SD, typename UBase>
    struct IsUnit< ::ScaledUnit<U, SN, SD, UBase> >
    {
        constexpr static bool Result = true;
    };

    template<typename U, int32 PN, int32 PD, typename UBase>
    struct IsUnit< ::PoweredUnit<U, PN, PD, UBase> >
    {
        constexpr static bool Result = true;
    };

    template<typename U, typename V>
    struct IsUnit< ::ComposedUnit<U, V> >
    {
        constexpr static bool Result = true;
    };

    ////////////////////////////////////////////
    // Trait for Scalar type
    ////////////////////////////////////////////
    template<typename U>
    struct IsScalar
    {
        constexpr static bool Result = false;
    };

    template<>
    struct IsScalar< ::Scalar >
    {
        constexpr static bool Result = true;
    };

    ////////////////////////////////////////////
    // Get the Scale of the Unit
    ////////////////////////////////////////////
    template<class T>
    struct GetScale
    {
        constexpr static uint64 Num = 1;
        constexpr static uint64 Den = 1;
    };

    template<class U, int64 N, int64 D, typename UBase>
    struct GetScale< ::ScaledUnit<U, N, D, UBase> >
    {
        constexpr static uint64 Num = N * GetScale<U>::Num;
        constexpr static uint64 Den = D * GetScale<U>::Den;
    };

    template<class U, int32 N, int32 D, typename UBase>
    struct GetScale< ::PoweredUnit<U, N, D, UBase> >
    {
        constexpr static uint64 Num = GetScale<U>::Num;
        constexpr static uint64 Den = GetScale<U>::Den;
    };


    ////////////////////////////////////////////
    // Get the Power of the Unit
    ////////////////////////////////////////////
    template<class T>
    struct GetUnitPower
    {
        constexpr static int32 Num = 0;
        constexpr static int32 Den = 0;
    };

    template<class U, int32 N, int32 D, typename UBase>
    struct GetUnitPower< ::PoweredUnit<U, N, D, UBase> >
    {
        constexpr static int32 Num = N + GetUnitPower<U>::Num;
        constexpr static int32 Den = D + GetUnitPower<U>::Den;
    };


    template<class T>
    struct GetPower
    {
        constexpr static int32 Num = 1;
        constexpr static int32 Den = 1;
    };

    template<class U, int32 N, int32 D, typename UBase>
    struct GetPower< ::PoweredUnit<U, N, D, UBase> >
    {
        constexpr static int32 Num = N * GetPower<U>::Num;
        constexpr static int32 Den = D * GetPower<U>::Den;
    };


    ////////////////////////////////////////////
    // Get Base Unit of a Unit
    ////////////////////////////////////////////
    template<typename U>
    struct GetBaseUnit
    {
        using BaseUnit = U;
    };

    template<typename SV, int64 SN, int64 SD, typename SUBase>
    struct GetBaseUnit< ::ScaledUnit<SV, SN, SD, SUBase> >
    {
        using BaseUnit = SUBase;
    };

    template<typename PV, int32 PN, int32 PD, typename PUBase>
    struct GetBaseUnit< ::PoweredUnit<PV, PN, PD, PUBase> >
    {
        using BaseUnit = PUBase;
    };


    ////////////////////////////////////////////
    // Trait for same base unit
    ////////////////////////////////////////////
    template<typename U, typename V>
    struct IsSameBaseUnit
    {
        constexpr static bool Result = IsSameTypes<typename GetBaseUnit<U>::BaseUnit, typename GetBaseUnit<V>::BaseUnit>::Result;
    };


    ////////////////////////////////////////////
    // Get the change in scale between two
    // scales
    ////////////////////////////////////////////
    template<int64 SNFrom, int64 SDFrom, int64 SNTarget, int64 SDTarget>
    struct GetScaleChange
    {
        constexpr static double Get() { return (static_cast<double>(SDTarget)* static_cast<double>(SNFrom)) / (static_cast<double>(SNTarget)* static_cast<double>(SDFrom)); }
    };

    // Particular case when the same scale, do not make the division for nothing
    template<int64 N, int64 D>
    struct GetScaleChange<N, D, N, D>
    {
        constexpr static double Get() { return 1.0; }
    };


    ////////////////////////////////////////////
    // Get the change factor between two units
    // taking into account the power and the
    // scale
    ////////////////////////////////////////////
    template<int64 SNTarget, int64 SDTarget, int64 SNFrom, int64 SDFrom, int32 PN, int32 PD>
    struct GetChangeFactorBase
    {
        constexpr static double GetFactor() { return std::pow(GetScaleChange<SNFrom, SDFrom, SNTarget, SDTarget>::Get(), (PN / static_cast<double>(PD))); }
    };

    // Particular case when Power of 1, do not make the pow function for nothing
    template<int64 SNTarget, int64 SDTarget, int64 SNFrom, int64 SDFrom, int32 F>
    struct GetChangeFactorBase<SNTarget, SDTarget, SNFrom, SDFrom, F, F>
    {
        constexpr static double GetFactor() { return GetScaleChange<SNFrom, SDFrom, SNTarget, SDTarget>::Get(); }
    };


    ////////////////////////////////////////////
    // Change the value of a unit V to
    // another unit U
    ////////////////////////////////////////////
    template < typename U, typename V >
    struct ChangeFactor
    {
        constexpr static double GetFactor()
        {
            return GetChangeFactorBase<GetScale<U>::Num, GetScale<U>::Den, GetScale<V>::Num, GetScale<V>::Den, GetPower<V>::Num, GetPower<V>::Den>::GetFactor();
        }

        constexpr static double GetValue(double value)
        {
            return value * GetFactor();
        }
    };

    template<typename U>
    struct ChangeFactor<U, U>
    {
        constexpr static double GetFactor()
        {
            return 1.0;
        }

        constexpr static double GetValue(double value)
        {
            return value;
        }
    };


    ////////////////////////////////////////////
    // Inverse Power result type
    ////////////////////////////////////////////
    template<typename U>
    struct InversePower
    {
        using Type = ::PoweredUnit<U, -1, 1, typename GetBaseUnit<U>::BaseUnit>;
    };

    template<>
    struct InversePower< ::Scalar>
    {
        using Type = ::Scalar;
    };

    template<typename V, int64 SN, int64 SD, typename UBase>
    struct InversePower< ::ScaledUnit<V, SN, SD, UBase> >
    {
        using Type = ::PoweredUnit< ::ScaledUnit<V, SN, SD, UBase>, -1, 1, UBase>;
    };

    template<typename V, int32 PN, int32 PD, typename UBase>
    struct InversePower< ::PoweredUnit<V, PN, PD, UBase> >
    {
        using Type = ::PoweredUnit<V, -PN, PD, UBase>;
    };

    template<typename V, int32 PD, typename UBase>
    struct InversePower< ::PoweredUnit<V, -1, PD, UBase> >
    {
        using Type = V;
    };

    template<typename U1, typename U2>
    struct InversePower< ::ComposedUnit<U1, U2> >
    {
        using Type = ::ComposedUnit<typename InversePower<U1>::Type, typename InversePower<U2>::Type>;
    };


    ////////////////////////////////////////////
    // Square root result type
    ////////////////////////////////////////////
    template<typename U>
    struct SqrtTransform
    {
        using Type = ::PoweredUnit<U, 1, 2, typename GetBaseUnit<U>::BaseUnit>;
    };

    template<>
    struct SqrtTransform< ::Scalar>
    {
        using Type = ::Scalar;
    };

    template<typename U, int32 N, int32 D, typename UBase>
    struct SqrtTransform< ::PoweredUnit<U, N, D, UBase> >
    {
        constexpr static int32 RN = IfElseValue<N % (D * 2) == 0, N / (D * 2), N>::Value;
        constexpr static int32 RD = IfElseValue<N % (D * 2) == 0, 1, D * 2>::Value;

        using Type = typename IfElseType<RN == RD, U, ::PoweredUnit<U, RN, RD, UBase> >::Type;
    };

    template<typename U1, typename U2>
    struct SqrtTransform< ::ComposedUnit<U1, U2> >
    {
        using Type = ::ComposedUnit<typename SqrtTransform<U1>::Type, typename SqrtTransform<U2>::Type>;
    };


    ////////////////////////////////////////////
    // Power result type
    ////////////////////////////////////////////
    template<typename U, int32 NMult, int32 DMult>
    struct PowerTransform
    {
        constexpr static int32 RN = GetPower<U>::Num * NMult;
        constexpr static int32 RD = GetPower<U>::Den * DMult;

        using Type = ::PoweredUnit<U, RN / PGCD<RN, RD>::Value, RD / PGCD<RN, RD>::Value, typename GetBaseUnit<U>::BaseUnit>;
    };

    template<int32 NMult, int32 DMult>
    struct PowerTransform< ::Scalar, NMult, DMult>
    {
        using Type = ::Scalar;
    };

    template<typename U, int32 N, int32 D, typename UBase, int NMult, int DMult>
    struct PowerTransform< ::PoweredUnit<U, N, D, UBase>, NMult, DMult >
    {
        constexpr static int32 RN = GetPower< ::PoweredUnit<U, N, D, UBase> >::Num * NMult;
        constexpr static int32 RD = GetPower< ::PoweredUnit<U, N, D, UBase> >::Den * DMult;

        using Type = ::PoweredUnit<U, RN / PGCD<RN, RD>::Value, RD / PGCD<RN, RD>::Value, UBase>;
    };

    template<typename U1, typename U2, int32 NMult, int32 DMult>
    struct PowerTransform< ::ComposedUnit<U1, U2>, NMult, DMult>
    {
        using Type = ::ComposedUnit<typename PowerTransform<U1, NMult, DMult>::Type, typename PowerTransform<U2, NMult, DMult>::Type >;
    };


    ////////////////////////////////////////////
    // Transform two units together
    ////////////////////////////////////////////
    template<typename U, typename V>
    struct Transform;

    template<typename V>
    struct Transform< ::Scalar, V>
    {
        using ReturnTypeMultiply = V;
        using ReturnTypeDivide = typename InversePower<V>::Type;

        constexpr static bool Find = true;

        constexpr static double GetChangeFactor()
        {
            return 1.0;
        }
    };

    template<typename U, int64 N, int64 D, typename UBase, typename V>
    struct Transform< ::ScaledUnit<U, N, D, UBase>, V>
    {
        // If the operation is valid, return the operator ReturnType,
        // Otherwise, return BaseType itself
        // Used for the construction of the main resulting type
        using ReturnTypeMultiply = typename ::ScaledUnit<U, N, D, UBase>::template OperatorResultType<V>::MultiplyType;
        using ReturnTypeDivide = typename ::ScaledUnit<U, N, D, UBase>::template OperatorResultType<V>::DivideType;

        // If the baseUnit of V is UBase, we find a match
        constexpr static bool VIsScalar = IsScalar<V>::Result;
        constexpr static bool Find = VIsScalar || IsSameBaseUnit<V, UBase>::Result;

        constexpr static double GetChangeFactor()
        {
            return ((Find && !VIsScalar) ? ChangeFactor< ::ScaledUnit<U, N, D, UBase>, V>::GetFactor() : 1.0);
        }
    };

    template<typename U, int32 N, int32 D, typename UBase, typename V>
    struct Transform< ::PoweredUnit<U, N, D, UBase>, V>
    {
        // If the operation is valid, return the operator ReturnType,
        // Otherwise, return BaseType itself
        // Used for the construction of the main resulting type
        using ReturnTypeMultiply = typename ::PoweredUnit<U, N, D, UBase>::template OperatorResultType<V>::MultiplyType;
        using ReturnTypeDivide = typename ::PoweredUnit<U, N, D, UBase>::template OperatorResultType<V>::DivideType;

        // We find a match if the two Unit are of the same baseUnit type
        constexpr static bool VIsScalar = IsScalar<V>::Result;
        constexpr static bool Find = VIsScalar || IsSameBaseUnit<V, UBase>::Result;

        constexpr static double GetChangeFactor()
        {
            return ((Find && !VIsScalar) ? ChangeFactor< ::PoweredUnit<U, N, D, UBase>, V>::GetFactor() : 1.0);
        }
    };

    // Will return ComposedUnit<U1, U2> minus V, or ComposedUnit<U1, U2>
    template<typename U1, typename U2, typename V>
    struct Transform< ::ComposedUnit<U1, U2>, V>
    {
        // Try to transform each unit of the ComposedUnit by the unit V
        using U1MultiplyType = typename Transform<U1, V>::ReturnTypeMultiply;
        using U2MultiplyType = typename Transform<U2, V>::ReturnTypeMultiply;

        using U1DivideType = typename Transform<U1, V>::ReturnTypeDivide;
        using U2DivideType = typename Transform<U2, V>::ReturnTypeDivide;

        constexpr static bool Find = Transform<U1, V>::Find || Transform<U2, V>::Find || IsSameTypes< ::ComposedUnit<U1, U2>, V>::Result;

        // ResultType = ComposedUnit<TransformU1, TransformU2>;
        // If we find a match
        //    If (U1*V==Scalar) ResultType = U2
        //    elseif (U2*V==Scalar) ResultType = U1
        //    else ResultType = ComposedUnit<TransformU1, TransformU2>
        // Otherwise, ResultType = ComposedUnit<TransformU1, TransformU2>
        //If there is a match, the result should be the composition of each transform
        using MultiplyTypeNoScalar = typename ::ComposedUnit<U1MultiplyType, U2MultiplyType>;
        using MultiplyResultStep1 = typename IfElseType<IsScalar<U1MultiplyType>::Result, U2, MultiplyTypeNoScalar>::Type;
        using MultiplyResultStep2 = typename IfElseType<IsScalar<U2MultiplyType>::Result, U1, MultiplyResultStep1>::Type;

        using ReturnTypeMultiply = typename IfElseType<Find, MultiplyResultStep2, ::ComposedUnit<U1, U2> >::Type;


        using DivideTypeNoScalar = typename ::ComposedUnit<U1DivideType, U2DivideType>;
        using DivideResultStep1 = typename IfElseType<IsScalar<U1DivideType>::Result, U2, DivideTypeNoScalar>::Type;
        using DivideResultStep2 = typename IfElseType<IsScalar<U2DivideType>::Result, U1, DivideResultStep1>::Type;

        using ReturnTypeDivide = typename IfElseType<Find, DivideResultStep2, ::ComposedUnit<U1, U2> >::Type;

        constexpr static double GetChangeFactor()
        {
            // GetChangeFactor return 1 if no match found
            return Transform<U1, V>::GetChangeFactor() * Transform<U2, V>::GetChangeFactor();
        }
    };


    template<typename U, typename V>
    struct TransformBase
    {
        constexpr static bool FindMultiply = Transform<U, V>::Find;
        constexpr static bool FindDivide = Transform<U, V>::Find;

        // If we find a match, return the Transform result
        // Otherwise, composed the two units
        using ReturnTypeMultiply = typename IfElseType<FindMultiply, typename Transform<U, V>::ReturnTypeMultiply, ::ComposedUnit<U, V> >::Type;
        using ReturnTypeDivide = typename IfElseType<FindDivide, typename Transform<U, V>::ReturnTypeDivide, ::ComposedUnit<U, typename InversePower<V>::Type> >::Type;

        constexpr static double GetChangeFactorMultiply()
        {
            return Transform<U, V>::GetChangeFactor();
        }

        constexpr static double GetChangeFactorDivide()
        {
            return Transform<U, V>::GetChangeFactor();
        }
    };

    template<typename U, typename U1, typename U2>
    struct TransformBase<U, ::ComposedUnit<U1, U2> >
    {
        using MultiplyType = TransformBase<typename TransformBase<U, U2>::ReturnTypeMultiply, U1>;
        using ReturnTypeMultiply = typename MultiplyType::ReturnTypeMultiply;

        using DivideType = TransformBase<typename TransformBase<U, U2>::ReturnTypeDivide, U1>;
        using ReturnTypeDivide = typename DivideType::ReturnTypeDivide;

        constexpr static bool FindMultiply = MultiplyType::FindMultiply;
        constexpr static bool FindDivide = DivideType::FindDivide;

        constexpr static double GetChangeFactorMultiply()
        {
            return MultiplyType::GetChangeFactorMultiply() * TransformBase<U, U2>::GetChangeFactorMultiply();
        }

        constexpr static double GetChangeFactorDivide()
        {
            return DivideType::GetChangeFactorDivide() * TransformBase<U, U2>::GetChangeFactorDivide();
        }
    };


    ////////////////////////////////////////////
    // Main utility struct to Transform two 
    // units together
    ////////////////////////////////////////////
    template<typename U, typename V>
    struct TransformUnit
    {
        using MultiplyType = typename TransformBase<U, V>::ReturnTypeMultiply;
        using DivideType = typename TransformBase<U, V>::ReturnTypeDivide;

        constexpr static MultiplyType Multiply(const U& u, const V& v)
        {
            return MultiplyType(u.Value() * (v.Value() * TransformBase<U, V>::GetChangeFactorMultiply()));
        }

        constexpr static DivideType Divide(const U& u, const V& v)
        {
            return DivideType(u.Value() / (v.Value() * TransformBase<U, V>::GetChangeFactorDivide()));
        }
    };

    template<typename U>
    struct TransformUnit<U, U>
    {
        using MultiplyType = typename TransformBase<U, U>::ReturnTypeMultiply;
        using DivideType = typename TransformBase<U, U>::ReturnTypeDivide;

        constexpr static MultiplyType Multiply(const U& u, const U& v)
        {
            return MultiplyType(u.Value() * v.Value());
        }

        constexpr static DivideType Divide(const U& u, const U& v)
        {
            return DivideType(u.Value() / v.Value());
        }
    };


    ////////////////////////////////////////////
    // To test if two Units are equivalent
    ////////////////////////////////////////////
    template<typename U, typename V>
    struct IsTransformableTo
    {
        constexpr static bool Result = IsSameBaseUnit<U, V>::Result && GetPower<U>::Num == GetPower<V>::Num && GetPower<U>::Den == GetPower<V>::Den;

        constexpr static double GetChangeFactor()
        {
            return Result ? ChangeFactor<U, V>::GetFactor() : 1.0;
        }
    };

    template<typename U, typename V>
    struct IsInSet
    {
        constexpr static bool Result = IsTransformableTo<U, V>::Result;

        constexpr static double GetChangeFactor()
        {
            return Result ? IsTransformableTo<U, V>::GetChangeFactor() : 1.0;
        }
    };

    template<typename U, typename V1, typename V2>
    struct IsInSet< U, ::ComposedUnit<V1, V2> >
    {
        constexpr static bool Result = IsInSet<U, V1>::Result || IsInSet<U, V2>::Result;

        constexpr static double GetChangeFactor()
        {
            return IsInSet<U, V1>::GetChangeFactor() * IsInSet<U, V2>::GetChangeFactor();
        }
    };

    template<typename U1, typename U2, typename V1, typename V2>
    struct IsInSet< ::ComposedUnit<U1, U2>, ::ComposedUnit<V1, V2> >
    {
        constexpr static bool Result = IsInSet<U1, ::ComposedUnit<V1, V2> >::Result && IsInSet<U2, ::ComposedUnit<V1, V2> >::Result;

        constexpr static double GetChangeFactor()
        {
            return IsInSet<U1, ::ComposedUnit<V1, V2> >::GetChangeFactor() * IsInSet<U2, ::ComposedUnit<V1, V2> >::GetChangeFactor();
        }
    };

    ////////////////////////////////////////////
    // To test if two Units are equivalent
    ////////////////////////////////////////////
    template<typename U, typename V>
    struct IsEquivalent
    {
        constexpr static bool Result = IsTransformableTo<U, V>::Result;

        constexpr static double GetFactor()
        {
            return Result ? IsTransformableTo<U, V>::GetChangeFactor() : 1.0;
        }
    };

    template<typename U1, typename U2, typename V1, typename V2>
    struct IsEquivalent < ::ComposedUnit<U1, U2>, ::ComposedUnit<V1, V2> >
    {
        constexpr static bool Result = IsInSet< ::ComposedUnit<U1, U2>, ::ComposedUnit<V1, V2> >::Result && IsInSet< ::ComposedUnit<V1, V2>, ::ComposedUnit<U1, U2> >::Result;

        constexpr static double GetFactor()
        {
            return Result ? IsInSet< ::ComposedUnit<U1, U2>, ::ComposedUnit<V1, V2> >::GetChangeFactor() : 1.0;
        }
    };


    ////////////////////////////////////////////
    // Write the unit output
    ////////////////////////////////////////////
    // General case
    template<typename Unit>
    struct outputUnit
    {
        template<typename S>
        static void Output(S& os)
        {
        }
    };

    template<typename U, int64 Num, int64 Den>
    struct outputUnitScaledBase;

    // Particular cases for ScaledUnit
    template<typename U, int64 Num, int64 Den>
    struct outputUnitScaledBase
    {
        template<typename S>
        static void Output(S& os)
        {
            os << Num;
            if (Den != 1)
                os << '/' << Den;
            os << '.';
            outputUnit<U>::Output(os);
        }
    };

    // TODO : Check to transform these traits for autogenerated with macro when defining ScaleFactor later?
    template<typename U>
    struct outputUnitScaledBase<U, 1, 1>
    {
        template<typename S>
        static void Output(S& os)
        {
            outputUnit<U>::Output(os);
        }
    };

    template<typename U>
    struct outputUnitScaledBase<U, 1, 10>
    {
        template<typename S>
        static void Output(S& os)
        {
            os << 'd';
            outputUnit<U>::Output(os);
        }
    };

    template<typename U>
    struct outputUnitScaledBase<U, 1, 100>
    {
        template<typename S>
        static void Output(S& os)
        {
            os << 'c';
            outputUnit<U>::Output(os);
        }
    };

    template<typename U>
    struct outputUnitScaledBase<U, 1, 1000>
    {
        template<typename S>
        static void Output(S& os)
        {
            os << 'm';
            outputUnit<U>::Output(os);
        }
    };

    template<typename U>
    struct outputUnitScaledBase<U, 1, 1000000000>
    {
        template<typename S>
        static void Output(S& os)
        {
            os << 'n';
            outputUnit<U>::Output(os);
        }
    };

    template<typename U>
    struct outputUnitScaledBase<U, 1, 1000000000000>
    {
        template<typename S>
        static void Output(S& os)
        {
            os << 'p';
            outputUnit<U>::Output(os);
        }
    };

    template<typename U>
    struct outputUnitScaledBase<U, 1, 1000000000000000>
    {
        template<typename S>
        static void Output(S& os)
        {
            os << 'f';
            outputUnit<U>::Output(os);
        }
    };

    template<typename U>
    struct outputUnitScaledBase<U, 1, 1000000000000000000>
    {
        template<typename S>
        static void Output(S& os)
        {
            os << 'a';
            outputUnit<U>::Output(os);
        }
    };

    template<typename U>
    struct outputUnitScaledBase<U, 1000, 1>
    {
        template<typename S>
        static void Output(S& os)
        {
            os << 'k';
            outputUnit<U>::Output(os);
        }
    };

    // General case for ScaledUnit
    template<typename U, int64 SN, int64 SD, typename UBase>
    struct outputUnit< ::ScaledUnit<U, SN, SD, UBase> >
    {
        template<typename S>
        static void Output(S& os)
        {
            outputUnitScaledBase<typename GetBaseUnit<U>::BaseUnit, GetScale< ::ScaledUnit<U, SN, SD, UBase> >::Num, GetScale< ::ScaledUnit<U, SN, SD, UBase> >::Den>::Output(os);
        }
    };

    // General case for PoweredUnit
    template<typename U, int32 PN, int32 PD, typename UBase>
    struct outputUnit< ::PoweredUnit<U, PN, PD, UBase> >
    {
        template<typename S>
        static void Output(S& os)
        {
            outputUnit<U>::Output(os);
            if (PN != PD)
            {
                os << '^';
                if (PN%PD)
                {
                    os << '(';
                    os << PN << '/' << PD;
                    os << ')';
                }
                else
                {
                    os << PN;
                }
            }
        }
    };

    // General case for ComposedUnit
    template<typename U1, typename U2>
    struct outputUnit< ::ComposedUnit<U1, U2> >
    {
        template<typename S>
        static void Output(S& os)
        {
            outputUnit<U1>::Output(os);
            os << '.';
            outputUnit<U2>::Output(os);
        }
    };
}


///////////////////////////
// Scaled Unit
///////////////////////////
template<typename U, int64 N = 1, int64 D = 1, typename UBase = typename ::internals::GetBaseUnit<U>::BaseUnit>
class ScaledUnit final : public Unit<ScaledUnit<U, N, D, UBase> >
{
public:

    ///////////////////////////
    // Traits for Operators Result Type
    ///////////////////////////
    template<typename T>
    struct OperatorResultType // Otherwise, handle by ComposedUnit
    {
        using MultiplyType = ScaledUnit<U, N, D, UBase>;
        using DivideType = ScaledUnit<U, N, D, UBase>;
    };

    template<typename V, int64 SN, int64 SD>
    struct OperatorResultType<ScaledUnit<V, SN, SD, UBase> > // Operator by a scale of the same base
    {
        using MultiplyType = PoweredUnit<ScaledUnit<U, N, D, UBase>, 2, 1, UBase>;
        using DivideType = Scalar;
    };

    template<typename V, int32 PN, int32 PD>
    struct OperatorResultType<PoweredUnit<V, PN, PD, UBase> > // Operator by Powered of Self
    {
        using MultiplyType = typename ::internals::IfElseType<PN + PD == 0, Scalar, PoweredUnit<ScaledUnit<U, N, D, UBase>, PN + PD, PD, UBase> >::Type;
        using DivideType = typename ::internals::IfElseType<PD - PN == 0, Scalar, PoweredUnit<ScaledUnit<U, N, D, UBase>, PD - PN, PD, UBase> >::Type;
    };


    ///////////////////////////
    // Constructors
    ///////////////////////////

    // From a double value
    constexpr explicit ScaledUnit(const double& val = 0.0)
        : Unit<ScaledUnit<U, N, D, UBase> >(val)
    {
    }

    // From a ScaledUnit with the same baseUnit, different scale -> Scale it!
    template<typename V, int64 SN, int64 SD>
    constexpr ScaledUnit(const ScaledUnit<V, SN, SD, UBase>& u)
        : Unit<ScaledUnit<U, N, D, UBase> >(::internals::ChangeFactor<ScaledUnit<U, N, D, UBase>, ScaledUnit<V, SN, SD, UBase> >::GetValue(u.Value()))
    {
    }

    // From a PoweredUnit (1/1) of a ScaledUnit same as this one -> Just the scale unit
    constexpr ScaledUnit(const PoweredUnit<ScaledUnit<U, N, D, UBase>, 1, 1, UBase>& u)
        : Unit<ScaledUnit<U, N, D, UBase> >(u.Value())
    {
    }

    // From a PoweredUnit (1/1) of a ScaledUnit with the same baseUnit, but not the same scale -> Scale it!
    template<typename V, int64 SN, int64 SD>
    constexpr ScaledUnit(const PoweredUnit<ScaledUnit<V, SN, SD, UBase>, 1, 1, UBase>& u)
        : Unit<ScaledUnit<U, N, D, UBase> >(::internals::ChangeFactor<ScaledUnit<U, N, D, UBase>, ScaledUnit<V, SN, SD, UBase> >::GetValue(u.Value()))
    {
    }

    template<typename V1, typename V2>
    constexpr ScaledUnit(const ComposedUnit<V1, V2>& u, typename ::internals::EnableIf< ::internals::IsSameTypes<UBase, ComposedUnit<V1, V2> >::Result>::EnableType* dummy = 0)
        : Unit<ScaledUnit<U, N, D, UBase> >(::internals::ChangeFactor<ScaledUnit<U, N, D, UBase>, ComposedUnit<V1, V2> >::GetValue(u.Value()))
    {
    }


    ///////////////////////////
    // Assignation operators
    ///////////////////////////
    // Assignation from the same Scaled Unit
    constexpr ScaledUnit<U, N, D, UBase>& operator=(const ScaledUnit<U, N, D, UBase>& u)
    {
        if (this == &u)
            return *this;

        Unit<ScaledUnit<U, N, D, UBase> >::m_value = u.Value();

        return *this;
    }

    // Assignation from a ScaledUnit with the same baseUnit, but not the same scale -> Scale it!
    template<typename V, int64 SN, int64 SD>
    constexpr ScaledUnit<U, N, D, UBase>& operator=(const ScaledUnit<V, SN, SD, UBase>& u)
    {
        using ChangeFactorType = ::internals::ChangeFactor<ScaledUnit<U, N, D, UBase>, ScaledUnit<V, SN, SD, UBase> >;

        Unit<ScaledUnit<U, N, D, UBase> >::m_value = ChangeFactorType::GetValue(u.Value());
        return *this;
    }

    template<typename U1, typename U2>
    constexpr typename ::internals::EnableIf< ::internals::IsSameTypes<UBase, ComposedUnit<U1, U2> >::Result, ScaledUnit<U, N, D, UBase> >::EnableType& operator=(const ComposedUnit<U1, U2>& u)
    {
        using ChangeFactorType = ::internals::ChangeFactor<ScaledUnit<U, N, D, UBase>, ComposedUnit<U1, U2> >;

        Unit<ScaledUnit<U, N, D, UBase> >::m_value = ChangeFactorType::GetValue(u.Value());
        return *this;
    }

    ///////////////////////////
    // Multiplication operators
    ///////////////////////////
    // Make accessible operator* from base class (Unit<U>)
    using Unit<ScaledUnit<U, N, D, UBase> >::operator*;

    // Multiply by self -> PoweredUnit
    constexpr typename OperatorResultType<ScaledUnit<U, N, D, UBase> >::MultiplyType operator*(const ScaledUnit<U, N, D, UBase>& u) const
    {
        using ReturnType = typename OperatorResultType<ScaledUnit<U, N, D, UBase> >::MultiplyType;
        return ReturnType(Unit<ScaledUnit<U, N, D, UBase> >::Value() * u.Value());
    }

    // Multiply by same baseUnit, Scaled -> PoweredUnit
    template<typename V, int64 SN, int64 SD>
    constexpr typename OperatorResultType<ScaledUnit<V, SN, SD, UBase> >::MultiplyType operator*(const ScaledUnit<V, SN, SD, UBase>& u) const
    {
        using ChangeFactorType = ::internals::ChangeFactor<ScaledUnit<U, N, D, UBase>, ScaledUnit<V, SN, SD, UBase> >;
        using ReturnType = typename OperatorResultType<ScaledUnit<V, SN, SD, UBase> >::MultiplyType;

        return ReturnType(Unit<ScaledUnit<U, N, D, UBase> >::Value() * ChangeFactorType::GetValue(u.Value()));
    }

    // Multiply by PoweredUnit, same baseUnit, same Scaled -> PoweredUnit'
    template<int32 PN, int32 PD>
    constexpr typename OperatorResultType<PoweredUnit<ScaledUnit<U, N, D, UBase>, PN, PD, UBase> >::MultiplyType operator*(const PoweredUnit<ScaledUnit<U, N, D, UBase>, PN, PD, UBase>& u) const
    {
        using ReturnType = typename OperatorResultType<PoweredUnit<ScaledUnit<U, N, D, UBase>, PN, PD, UBase> >::MultiplyType;

        return ReturnType(Unit<ScaledUnit<U, N, D, UBase> >::Value() * u.Value());
    }

    // Multiply by PoweredUnit, same baseUnit, different Scaled -> PoweredUnit'
    template<typename V, int64 SN, int64 SD, int32 PN, int32 PD>
    constexpr typename OperatorResultType<PoweredUnit<ScaledUnit<V, SN, SD, UBase>, PN, PD, UBase> >::MultiplyType operator*(const PoweredUnit<ScaledUnit<V, SN, SD, UBase>, PN, PD, UBase>& u) const
    {
        using ChangeFactorType = ::internals::ChangeFactor<ScaledUnit<U, N, D, UBase>, PoweredUnit<ScaledUnit<V, SN, SD, UBase>, PN, PD, UBase> >;
        using ReturnType = typename OperatorResultType<PoweredUnit<ScaledUnit<V, SN, SD, UBase>, PN, PD, UBase> >::MultiplyType;

        return ReturnType(Unit<ScaledUnit<U, N, D, UBase> >::Value() * ChangeFactorType::GetValue(u.Value()));
    }

    // Multiply by ScaledUnit, different baseUnit -> ComposedUnit
    template<typename V, int64 N2, int64 D2, typename UBase2>
    constexpr typename ::internals::TransformUnit<ScaledUnit<U, N, D, UBase>, ScaledUnit<V, N2, D2, UBase2> >::MultiplyType operator*(const ScaledUnit<V, N2, D2, UBase2>& u) const
    {
        using ReturnType = ::internals::TransformUnit<ScaledUnit<U, N, D, UBase>, ScaledUnit<V, N2, D2, UBase2> >;

        return ReturnType::Multiply(*this, u);
    }

    // Multiply by PoweredUnit, different unit -> ComposedUnit
    template<typename V, int32 N2, int32 D2, typename UBase2>
    constexpr typename ::internals::TransformUnit<ScaledUnit<U, N, D, UBase>, PoweredUnit<V, N2, D2, UBase2> >::MultiplyType operator*(const PoweredUnit<V, N2, D2, UBase2>& u) const
    {
        using ReturnType = ::internals::TransformUnit<ScaledUnit<U, N, D, UBase>, PoweredUnit<V, N2, D2, UBase2> >;

        return ReturnType::Multiply(*this, u);
    }

    // Multiply by ComposedUnit
    template<typename U1, typename U2>
    constexpr typename ::internals::TransformUnit<ComposedUnit<U1, U2>, ScaledUnit<U, N, D, UBase> >::MultiplyType operator*(const ComposedUnit<U1, U2>& u) const
    {
        using ReturnType = ::internals::TransformUnit<ComposedUnit<U1, U2>, ScaledUnit<U, N, D, UBase> >;

        return ReturnType::Multiply(u, *this);
    }


    ///////////////////////////
    // Division operators
    ///////////////////////////
    // Make accessible operator/ from base class (Unit<U>)
    using Unit<ScaledUnit<U, N, D, UBase> >::operator/;

    // Divide by ScaledUnit, same baseUnit, different scale -> double
    template<typename V, int64 SN, int64 SD>
    constexpr typename OperatorResultType<ScaledUnit<V, SN, SD, UBase> >::DivideType operator/(const ScaledUnit<V, SN, SD, UBase>& u) const
    {
        using ChangeFactorType = ::internals::ChangeFactor<ScaledUnit<U, N, D, UBase>, ScaledUnit<V, SN, SD, UBase> >;
        using ReturnType = typename OperatorResultType<ScaledUnit<V, SN, SD, UBase> >::DivideType;

        return ReturnType(Unit<ScaledUnit<U, N, D, UBase> >::Value() / ChangeFactorType::GetValue(u.Value()));
    }

    // Divide by PoweredUnit, same baseUnit, different scale -> Powered ScaledUnit
    template<typename V, int32 PN, int32 PD>
    constexpr typename OperatorResultType<PoweredUnit<V, PN, PD, UBase> >::DivideType operator/(const PoweredUnit<V, PN, PD, UBase>& u) const
    {
        using ChangeFactorType = ::internals::ChangeFactor<ScaledUnit<U, N, D, UBase>, PoweredUnit<V, PN, PD, UBase> >;
        using ReturnType = typename OperatorResultType<PoweredUnit<V, PN, PD, UBase> >::DivideType;

        return ReturnType(Unit<ScaledUnit<U, N, D, UBase> >::Value() / ChangeFactorType::GetValue(u.Value()));
    }

    // Divide by PoweredUnit, same baseUnit, same scale -> Powered ScaledUnit
    template<int32 PN, int32 PD>
    constexpr typename OperatorResultType<PoweredUnit<ScaledUnit<U, N, D, UBase>, PN, PD, UBase> >::DivideType operator/(const PoweredUnit<ScaledUnit<U, N, D, UBase>, PN, PD, UBase>& u) const
    {
        using ReturnType = typename OperatorResultType<PoweredUnit<ScaledUnit<U, N, D, UBase>, PN, PD, UBase> >::DivideType;

        return ReturnType(Unit<ScaledUnit<U, N, D, UBase> >::Value() / u.Value());
    }

    // Divide by ScaledUnit, different baseUnit -> ComposedUnit and Powered -1 for other ScaledUnit
    template<typename V, int64 SN, int64 SD, typename UBase2>
    constexpr typename ::internals::TransformUnit<ScaledUnit<U, N, D, UBase>, ScaledUnit<V, SN, SD, UBase2> >::DivideType operator/(const ScaledUnit<V, SN, SD, UBase2>& u) const
    {
        using ReturnType = ::internals::TransformUnit<ScaledUnit<U, N, D, UBase>, ScaledUnit<V, SN, SD, UBase2> >;

        return ReturnType::Divide(*this, u);
    }

    // Divide by ScaledUnit, different baseUnit -> ComposedUnit and Powered -1 for other ScaledUnit
    template<typename V, int32 PN, int32 PD, typename UBase2>
    constexpr typename ::internals::TransformUnit<ScaledUnit<U, N, D, UBase>, PoweredUnit<V, PN, PD, UBase2> >::DivideType operator/(const PoweredUnit<V, PN, PD, UBase2>& u) const
    {
        using ReturnType = ::internals::TransformUnit<ScaledUnit<U, N, D, UBase>, PoweredUnit<V, PN, PD, UBase2> >;

        return ReturnType::Divide(*this, u);
    }

    // Divide by a ComposedUnit -> ComposedUnit
    template<typename U1, typename U2>
    constexpr typename ::internals::TransformUnit<ScaledUnit<U, N, D, UBase>, ComposedUnit<U1, U2> >::DivideType operator/(const ComposedUnit<U1, U2>& u) const
    {
        using ReturnType = ::internals::TransformUnit<ScaledUnit<U, N, D, UBase>, ComposedUnit<U1, U2> >;

        return ReturnType::Divide(*this, u);
    }
};


//////////////////////////////////
// Powered Unit
//////////////////////////////////
template<typename U, int32 N, int32 D = 1, typename UBase = typename ::internals::GetBaseUnit<U>::BaseUnit>
class PoweredUnit final : public Unit<PoweredUnit<U, N, D, UBase> >
{
public:

    ///////////////////////////
    // Traits for Operators Result Type
    ///////////////////////////
    template<typename T>
    struct OperatorResultType // Otherwise, handle by the ComposedUnit
    {
        // Because we can't explicitly specialize the template into another template class, we must test if the T is the same as the U typename
        using MultiplyType = typename ::internals::IfElseType< ::internals::IsSameTypes<T, U>::Result, typename ::internals::IfElseType<N + D == 0, Scalar, PoweredUnit<U, N + D, D, UBase> >::Type, PoweredUnit<U, N, D, UBase> >::Type;
        using DivideType = typename ::internals::IfElseType< ::internals::IsSameTypes<T, U>::Result, typename ::internals::IfElseType<N - D == 1, U, typename ::internals::IfElseType<N - D == 0, Scalar, PoweredUnit<U, N - D, D, UBase> >::Type>::Type, PoweredUnit<U, N, D, UBase> >::Type;
    };

    template<typename V, int64 SN, int64 SD>
    struct OperatorResultType<ScaledUnit<V, SN, SD, UBase> > // Operator by ScaledUnit
    {
        using MultiplyType = typename ::internals::IfElseType<N + D == 0, Scalar, PoweredUnit<U, N + D, D, UBase> >::Type;
        using DivideType = typename ::internals::IfElseType<N - D == 1, U, typename ::internals::IfElseType<N - D == 0, Scalar, PoweredUnit<U, N - D, D, UBase> >::Type>::Type;
    };

    template<typename V, int32 PN, int32 PD>
    struct OperatorResultType<PoweredUnit<V, PN, PD, UBase> > // Operator by Powered of the same unit
    {
        constexpr static int32 PGCDMultiply = ::internals::PGCD<(N*PD + D*PN), D*PD>::Value;
        constexpr static int32 PGCDDivide = ::internals::PGCD<(N*PD - D*PN), D*PD>::Value;

        using MultiplyType = typename ::internals::IfElseType<(N*PD + D*PN) % (D*PD) == 0, typename ::internals::IfElseType<(N*PD) + (D*PN) == 0, Scalar, typename ::internals::IfElseType<(N*PD) + (PN*D) == (D*PD), U, PoweredUnit<U, ((N*PD) + (D*PN)) / (D*PD), 1, UBase> >::Type>::Type, PoweredUnit<U, ((N*PD) + (PN*D)) / PGCDMultiply, (D*PD) / PGCDMultiply, UBase> >::Type;
        using DivideType = typename ::internals::IfElseType<(N*PD - D*PN) % (D*PD) == 0, typename ::internals::IfElseType<(N*PD) - (D*PN) == 0, Scalar, typename ::internals::IfElseType<(N*PD) - (PN*D) == (D*PD), U, PoweredUnit<U, ((N*PD) - (D*PN)) / (D*PD), 1, UBase> >::Type>::Type, PoweredUnit<U, ((N*PD) - (PN*D)) / PGCDMultiply, (D*PD) / PGCDMultiply, UBase> >::Type;
    };


    ///////////////////////////
    // Constructors
    ///////////////////////////
    constexpr explicit PoweredUnit(const double& val = 0.0)
        : Unit<PoweredUnit<U, N, D, UBase> >(val)
    {
    }

    constexpr PoweredUnit(const PoweredUnit<U, N, D, UBase>& u)
        : Unit<PoweredUnit<U, N, D, UBase> >(u.Value())
    {
    }

    template<typename V>
    constexpr PoweredUnit(const PoweredUnit<V, N, D, UBase>& u)
        : Unit<PoweredUnit<U, N, D, UBase> >(::internals::ChangeFactor< PoweredUnit<U, N, D, UBase>, PoweredUnit<V, N, D, UBase> >::GetValue(u.Value()))
    {
    }


    ///////////////////////////
    // Multiplication operators
    ///////////////////////////
    // Make accessible operator* from base class (Unit<U>)
    using Unit<PoweredUnit<U, N, D, UBase> >::operator*;

    // Multiply by the inverse power -> double
    constexpr typename OperatorResultType<PoweredUnit<U, -N, D, UBase> >::MultiplyType operator*(const PoweredUnit<U, -N, D, UBase>& u) const
    {
        using ReturnType = typename OperatorResultType<PoweredUnit<U, -N, D, UBase> >::MultiplyType;

        return ReturnType(Unit<PoweredUnit<U, N, D, UBase> >::Value() * u.Value());
    }

    // Multiply by the unit itself (same scale, power = 1)
    constexpr typename OperatorResultType<U>::MultiplyType operator*(const U& u) const
    {
        using ReturnType = typename OperatorResultType<U>::MultiplyType;

        return ReturnType(Unit<PoweredUnit<U, N, D, UBase> >::Value() * u.Value());
    }

    // Multiply by a ScaledUnit, with same baseUnit, different scale
    template<typename V, int64 SN, int64 SD>
    constexpr typename OperatorResultType<ScaledUnit<V, SN, SD, UBase> >::MultiplyType operator*(const ScaledUnit<V, SN, SD, UBase>& u) const
    {
        using ChangeFactorType = ::internals::ChangeFactor<U, ScaledUnit<V, SN, SD, UBase> >;
        using ReturnType = typename OperatorResultType<ScaledUnit<V, SN, SD, UBase> >::MultiplyType;

        return ReturnType(Unit<PoweredUnit<U, N, D, UBase> >::Value() * ChangeFactorType::GetValue(u.Value()));
    }

    // Multiply by another PoweredUnit, with the same baseUnit, but different scale and different power
    template<typename V, int PN2, int PD2>
    constexpr typename OperatorResultType<PoweredUnit<V, PN2, PD2, UBase> >::MultiplyType operator*(const PoweredUnit<V, PN2, PD2, UBase>& u) const
    {
        using ChangeFactorType = ::internals::ChangeFactor<U, PoweredUnit<V, PN2, PD2, UBase> >;
        using ReturnType = typename OperatorResultType<PoweredUnit<V, PN2, PD2, UBase> >::MultiplyType;

        return ReturnType(Unit<PoweredUnit<U, N, D, UBase> >::Value() * ChangeFactorType::GetValue(u.Value()));
    }

    // Multiply by Another PoweredUnit -> ComposedUnit
    template<typename V, int32 PN2, int32 PD2, typename UBase2>
    constexpr typename ::internals::TransformUnit<PoweredUnit<U, N, D, UBase>, PoweredUnit<V, PN2, PD2, UBase2> >::MultiplyType operator*(const PoweredUnit<V, PN2, PD2, UBase2>& u) const
    {
        using ReturnType = ::internals::TransformUnit<PoweredUnit<U, N, D, UBase>, PoweredUnit<V, PN2, PD2, UBase2> >;

        return ReturnType::Multiply(*this, u);
    }

    // Multiply by Another ScaledUnit -> ComposedUnit
    template<typename V, int64 SN, int64 SD, typename UBase2>
    constexpr typename ::internals::TransformUnit<PoweredUnit<U, N, D, UBase>, ScaledUnit<V, SN, SD, UBase2> >::MultiplyType operator*(const ScaledUnit<V, SN, SD, UBase2>& u) const
    {
        using ReturnType = ::internals::TransformUnit<PoweredUnit<U, N, D, UBase>, ScaledUnit<V, SN, SD, UBase2> >;

        return ReturnType::Multiply(*this, u);
    }


    // Multiply by Another PoweredUnit -> ComposedUnit
    template<typename U1, typename U2>
    constexpr typename ::internals::TransformUnit<ComposedUnit<U1, U2>, PoweredUnit<U, N, D, UBase> >::MultiplyType operator*(const ComposedUnit<U1, U2>& u) const
    {
        using ReturnType = ::internals::TransformUnit<ComposedUnit<U1, U2>, PoweredUnit<U, N, D, UBase> >;

        return ReturnType::Multiply(u, *this);
    }


    ///////////////////////////
    // Division Operators
    ///////////////////////////
    // Make accessible operator/ from base class (Unit<U>)
    using Unit<PoweredUnit<U, N, D, UBase> >::operator/;

    // Divide by a transformed unit of the same power
    template<typename V>
    constexpr typename OperatorResultType<PoweredUnit<V, N, D, UBase> >::DivideType operator/(const PoweredUnit<V, N, D, UBase>& u) const
    {
        using ChangeFactorType = ::internals::ChangeFactor<U, PoweredUnit<V, N, D, UBase> >;
        using ReturnType = typename OperatorResultType<PoweredUnit<V, N, D, UBase> >::DivideType;

        return ReturnType(Unit<PoweredUnit<U, N, D, UBase> >::Value() / ChangeFactorType::GetValue(u.Value()));
    }

    // Divide by a different power of the same unit
    template<int32 N2, int32 D2>
    constexpr typename OperatorResultType<PoweredUnit<U, N2, D2, UBase> >::DivideType operator/(const PoweredUnit<U, N2, D2, UBase>& u) const
    {
        using ReturnType = typename OperatorResultType<PoweredUnit<U, N2, D2, UBase> >::DivideType;

        return ReturnType(Unit<PoweredUnit<U, N, D, UBase> >::Value() / u.Value());
    }

    // Divide by a transformed unit of different power
    template<typename V, int32 N2, int32 D2>
    constexpr typename OperatorResultType<PoweredUnit<V, N2, D2, UBase> >::DivideType operator/(const PoweredUnit<V, N2, D2, UBase>& u) const
    {
        using ChangeFactorType = ::internals::ChangeFactor<U, PoweredUnit<V, N2, D2, UBase> >;
        using ReturnType = typename OperatorResultType<PoweredUnit<V, N2, D2, UBase> >::DivideType;

        return ReturnType(Unit<PoweredUnit<U, N, D, UBase> >::Value() / ChangeFactorType::GetValue(u.Value()));
    }

    // Divide by a ScaledUnit of the same baseUnit, different scale -> Scale it, and divide
    template<typename V, int64 N2, int64 D2>
    constexpr typename OperatorResultType<ScaledUnit<V, N2, D2, UBase> >::DivideType operator/(const ScaledUnit<V, N2, D2, UBase>& u) const
    {
        using ChangeFactorType = ::internals::ChangeFactor<U, ScaledUnit<V, N2, D2, UBase> >;
        using ReturnType = typename OperatorResultType<ScaledUnit<V, N2, D2, UBase> >::DivideType;

        return ReturnType(Unit<PoweredUnit<U, N, D, UBase> >::Value() / ChangeFactorType::GetValue(u.Value()));
    }

    // Divide by a ScaledUnit of the same baseUnit -> Scale it, and divide
    constexpr typename OperatorResultType<U>::DivideType operator/(const U& u) const
    {
        using ReturnType = typename OperatorResultType<U>::DivideType;

        return ReturnType(Unit<PoweredUnit<U, N, D, UBase> >::Value() / u.Value());
    }

    // Divide by a PoweredUnit of another baseUnit -> Compose it
    template<typename V, int32 N2, int32 D2, typename UBase2>
    constexpr typename ::internals::TransformUnit<PoweredUnit<U, N, D, UBase>, PoweredUnit<V, N2, D2, UBase2> >::DivideType operator/(const PoweredUnit<V, N2, D2, UBase2>& u) const
    {
        using ReturnType = ::internals::TransformUnit<PoweredUnit<U, N, D, UBase>, PoweredUnit<V, N2, D2, UBase2> >;

        return ReturnType::Divide(*this, u);
    }

    // Divide by a ScaledUnit of another baseUnit -> Compose it
    template<typename V, int64 N2, int64 D2, typename UBase2>
    constexpr typename ::internals::TransformUnit<PoweredUnit<U, N, D, UBase>, ScaledUnit<V, N2, D2, UBase2> >::DivideType operator/(const ScaledUnit<V, N2, D2, UBase2>& u) const
    {
        using ReturnType = ::internals::TransformUnit<PoweredUnit<U, N, D, UBase>, ScaledUnit<V, N2, D2, UBase2> >;

        return ReturnType::Divide(*this, u);
    }

    // Divide by a ComposedUnit -> ComposedUnit
    template<typename U1, typename U2>
    constexpr typename ::internals::TransformUnit<PoweredUnit<U, N, D, UBase>, ComposedUnit<U1, U2> >::DivideType operator/(const ComposedUnit<U1, U2>& u) const
    {
        using ReturnType = ::internals::TransformUnit<PoweredUnit<U, N, D, UBase>, ComposedUnit<U1, U2> >;

        return ReturnType::Divide(*this, u);
    }
};


////////////////////////////////////////////
// Composed Unit
////////////////////////////////////////////
template<typename U1, typename U2>
class ComposedUnit final : public Unit<ComposedUnit<U1, U2> >
{
public:

    ///////////////////////////
    // Constructors
    ///////////////////////////
    constexpr explicit ComposedUnit(const double& val = 0.0)
        : Unit<ComposedUnit<U1, U2> >(val)
    {
    }

    constexpr ComposedUnit(const ComposedUnit<U1, U2>& u)
        : Unit<ComposedUnit<U1, U2> >(u.Value())
    {
    }

    template<typename V1, typename V2>
    constexpr ComposedUnit(const ComposedUnit<V1, V2>& u, typename ::internals::EnableIf< ::internals::IsEquivalent<ComposedUnit<U1, U2>, ComposedUnit<V1, V2> >::Result>::EnableType* dummy = 0)
        : Unit<ComposedUnit<U1, U2> >(u.Value() * ::internals::IsEquivalent<ComposedUnit<U1, U2>, ComposedUnit<V1, V2> >::GetFactor())
    {
    }

    // From a PoweredUnit (1/1) of a ComposedUnit same as this one -> Just the ComposedUnit
    constexpr ComposedUnit(const PoweredUnit<ComposedUnit<U1, U2>, 1, 1, ComposedUnit<U1, U2> >& u)
        : Unit<ComposedUnit<U1, U2> >(u.Value())
    {
    }

    template<typename U, int64 N, int64 D>
    constexpr ComposedUnit(const ScaledUnit<U, N, D, ComposedUnit<U1, U2> >& u)
        : Unit<ComposedUnit<U1, U2> >(::internals::ChangeFactor< ComposedUnit<U1, U2>, ScaledUnit<U, N, D, ComposedUnit<U1, U2> > >::GetValue(u.Value()))
    {
    }


    ///////////////////////////
    // Assignation operators
    ///////////////////////////
    using Unit<ComposedUnit<U1, U2> >::operator=;

    template<typename V1, typename V2>
    constexpr typename ::internals::EnableIf< ::internals::IsEquivalent<ComposedUnit<U1, U2>, ComposedUnit<V1, V2> >::Result, ComposedUnit<U1, U2> >::EnableType& operator=(const ComposedUnit<V1, V2>& v)
    {
        double changeFactor = ::internals::IsEquivalent<ComposedUnit<U1, U2>, ComposedUnit<V1, V2> >::GetFactor();
        Unit<ComposedUnit<U1, U2> >::m_value = v.Value() * changeFactor;
        return *this;
    }


    ///////////////////////////
    // Multiplication operators
    ///////////////////////////
    // Make accessible operator* from base class (Unit<U>)
    using Unit<ComposedUnit<U1, U2> >::operator*;

    // Multiply by Another Unit -> ComposedUnit
    template<typename V>
    constexpr typename ::internals::EnableIf< ::internals::IsUnit<V>::Result, typename ::internals::TransformUnit<ComposedUnit<U1, U2>, V>::MultiplyType>::EnableType operator*(const V& u) const
    {
        using ReturnType = ::internals::TransformUnit<ComposedUnit<U1, U2>, V>;

        return ReturnType::Multiply(*this, u);
    }

    // Multiply by self -> PoweredUnit
    constexpr typename ::internals::TransformUnit<ComposedUnit<U1, U2>, ComposedUnit<U1, U2> >::MultiplyType operator*(const ComposedUnit<U1, U2>& u) const
    {
        using ReturnType = typename ::internals::TransformUnit<ComposedUnit<U1, U2>, ComposedUnit<U1, U2> >::MultiplyType;

        return ReturnType(Unit<ComposedUnit<U1, U2> >::Value() * u.Value());
    }


    ///////////////////////////
    // Division operators
    ///////////////////////////
    // Make accessible operator/ from base class (Unit<U>)
    using Unit<ComposedUnit<U1, U2> >::operator/;

    // Divide by Another Unit -> ComposedUnit
    template<typename V>
    constexpr typename ::internals::EnableIf< ::internals::IsUnit<V>::Result, typename ::internals::TransformUnit<ComposedUnit<U1, U2>, V>::DivideType>::EnableType operator/(const V& u) const
    {
        using ReturnType = ::internals::TransformUnit<ComposedUnit<U1, U2>, V>;

        return ReturnType::Divide(*this, u);
    }
};


////////////////////////////////////////////
// Global Unit operators functions
////////////////////////////////////////////
template<typename U>
constexpr typename ::internals::InversePower<U>::Type operator/(const Scalar& d, const Unit<U>& u)
{
    return typename ::internals::InversePower<U>::Type(d.Value() / u.Value());
}

template<typename F, typename U>
constexpr typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, typename ::internals::InversePower<U>::Type>::EnableType operator/(const F& f, const Unit<U>& u)
{
    return typename ::internals::InversePower<U>::Type(f / u.Value());
}

template<typename U>
constexpr typename ::internals::EnableIf< ::internals::IsUnit<U>::Result, typename ::internals::SqrtTransform<U>::Type>::EnableType sqrt(const U& u)
{
    return typename ::internals::SqrtTransform<U>::Type(std::sqrt(u.Value()));
}

template<int32 PowNum, int32 PowDen = 1, typename U>
constexpr typename ::internals::EnableIf< ::internals::IsUnit<U>::Result, typename ::internals::PowerTransform<U, PowNum, PowDen>::Type>::EnableType pow(const Unit<U>& u)
{
    return typename ::internals::PowerTransform<U, PowNum, PowDen>::Type(std::pow(u.Value(), static_cast<double>(PowNum) / static_cast<double>(PowDen)));
}


////////////////////////////////////////////
// Vector with component Unit
////////////////////////////////////////////
template<typename U>
class Vector
{
public:
    constexpr static uint32 DIMENSION = 3;

    constexpr Vector()
    {
        v[0] = U(0.0);
        v[1] = U(0.0);
        v[2] = U(0.0);
    }

    constexpr Vector(U xpos, U ypos, U zpos = U())
    {
        v[0] = xpos;
        v[1] = ypos;
        v[2] = zpos;
    }

    template<typename V>
    constexpr Vector(const Vector<V>& u, typename ::internals::EnableIf< ::internals::IsEquivalent<U, V>::Result>::EnableType* dummy = 0)
    {
        double changeFactor = ::internals::IsEquivalent<U, V>::GetFactor();
        v[0] = U(u.x().Value() * changeFactor);
        v[1] = U(u.y().Value() * changeFactor);
        v[2] = U(u.z().Value() * changeFactor);
    }

    constexpr const U* constData() const { return v; }

    constexpr bool isNull() const
    {
        U tresshold = U(0.0000000001);
        return abs(v[0]) < tresshold && abs(v[1]) < tresshold && abs(v[2]) < tresshold;
    }

    constexpr inline U x() const { return v[0]; }
    constexpr inline U y() const { return v[1]; }
    constexpr inline U z() const { return v[2]; }

    constexpr inline void setX(U aX) { v[0] = aX; }
    constexpr inline void setY(U aY) { v[1] = aY; }
    constexpr inline void setZ(U aZ) { v[2] = aZ; }

    constexpr U& operator[](unsigned int idx)
    {
        return v[idx];
    }

    constexpr U operator[](unsigned int idx) const
    {
        return v[idx];
    }

    constexpr U length() const
    {
        return sqrt(lengthSquared());
    }

    constexpr typename ::internals::TransformUnit<U, U>::MultiplyType lengthSquared() const
    {
        return dotProduct(*this, *this);
    }

    constexpr Vector<Scalar> normalized() const
    {
        double len = dLenSquare();

        if (len - 1.0 == 0.0)
            return Vector<Scalar>(Scalar((double)v[0]), Scalar((double)v[1]), Scalar((double)v[2]));
        else if (len > 0.0)
            return *this / length();
        else
            return Vector<Scalar>();
    }

    constexpr void correct()
    {
        if (isNull())
        {
            v[0] = U(0.0);
            v[1] = U(0.0);
            v[2] = U(0.0);
        }
    }

    template<typename U2>
    constexpr static typename ::internals::TransformUnit<U, U2>::MultiplyType dotProduct(const Vector<U>& v1, const Vector<U2>& v2)
    {
        return v1.dotProduct(v2);
    }

    template<typename U2>
    constexpr typename ::internals::TransformUnit<U, U2>::MultiplyType dotProduct(const Vector<U2>& v2) const
    {
        using UnitTransformType = ::internals::TransformUnit<U, U2>;
        using ReturnType = typename UnitTransformType::MultiplyType;

        return ReturnType(UnitTransformType::Multiply(v[0], v2.x()) + UnitTransformType::Multiply(v[1], v2.y()) + UnitTransformType::Multiply(v[2], v2.z()));
    }

    template<typename U2>
    constexpr static Vector<typename ::internals::TransformUnit<U, U2>::MultiplyType> crossProduct(const Vector<U>& v1, const Vector<U2>& v2)
    {
        return v1.crossProduct(v2);
    }

    template<typename U2>
    constexpr Vector<typename ::internals::TransformUnit<U, U2>::MultiplyType> crossProduct(const Vector<U2>& v2)
    {
        return Vector<typename ::internals::TransformUnit<U, U2>::MultiplyType>(v[1] * v2.z() - v[2] * v2.y(),
            v[2] * v2.x() - v[0] * v2.z(),
            v[0] * v2.y() - v[1] * v2.x());
    }

    template<typename U2>
    constexpr static Vector<Scalar> normal(const Vector<U>& v1, const Vector<U2>& v2)
    {
        return crossProduct(v1, v2).normalized();
    }

    //////////////////////////////////////
    // Assignation operators
    //////////////////////////////////////
    template<typename U2>
    constexpr typename ::internals::EnableIf< ::internals::IsEquivalent<U, U2>::Result, Vector<U> >::EnableType& operator=(const Vector<U2>& vec)
    {
        v[0] = vec.x();
        v[1] = vec.y();
        v[2] = vec.z();
        return *this;
    }

    //////////////////////////////////////
    // Coumpound assigment operators
    //////////////////////////////////////
    template<typename U2>
    constexpr typename ::internals::EnableIf< ::internals::IsEquivalent<U, U2>::Result, Vector<U> >::EnableType &operator+=(const Vector<U2> &vector)
    {
        v[0] += vector.x();
        v[1] += vector.y();
        v[2] += vector.z();
        return *this;
    }

    template<typename U2>
    constexpr typename ::internals::EnableIf< ::internals::IsEquivalent<U, U2>::Result, Vector<U> >::EnableType &operator-=(const Vector<U2> &vector)
    {
        v[0] -= vector.x();
        v[1] -= vector.y();
        v[2] -= vector.z();
        return *this;
    }

    template<typename F>
    constexpr typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, Vector<U>&>::EnableType operator*=(const F& factor)
    {
        v[0] *= factor;
        v[1] *= factor;
        v[2] *= factor;
        return *this;
    }

    template<typename F>
    constexpr typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, Vector<U>&>::EnableType operator/=(const F& divisor)
    {
        v[0] /= divisor;
        v[1] /= divisor;
        v[2] /= divisor;
        return *this;
    }

    //////////////////////////////////////
    // Equality operators
    //////////////////////////////////////
    template<typename U2>
    constexpr typename ::internals::EnableIf< ::internals::IsEquivalent<U, U2>::Result, bool>::EnableType operator==(const Vector<U2> &v2) const
    {
        return v[0] == v2.x() && v[1] == v2.y() && v[2] == v2.z();
    }

    template<typename U2>
    constexpr typename ::internals::EnableIf< ::internals::IsEquivalent<U, U2>::Result, bool>::EnableType operator!=(const Vector<U2> &v2) const
    {
        return !(*this == v2);
    }


    //////////////////////////////////////
    // Addition/Substraction operators
    //////////////////////////////////////
    template<typename U2>
    constexpr typename ::internals::EnableIf< ::internals::IsEquivalent<U, U2>::Result, Vector<U> >::EnableType operator+(const Vector<U2> &v2) const
    {
        return Vector<U>(v[0] + v2.x(), v[1] + v2.y(), v[2] + v2.z());
    }

    template<typename U2>
    constexpr typename ::internals::EnableIf< ::internals::IsEquivalent<U, U2>::Result, Vector<U> >::EnableType operator-(const Vector<U2> &v2) const
    {
        return Vector<U>(v[0] - v2.x(), v[1] - v2.y(), v[2] - v2.z());
    }


    //////////////////////////////////////
    // Multiplication operators
    //////////////////////////////////////
    // Vector<U> *  double
    template<typename F>
    constexpr typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, Vector<U> >::EnableType operator*(const F& factor) const
    {
        return Vector<U>(v[0] * factor, v[1] * factor, v[2] * factor);
    }

    // Vector<U> * V
    template<typename V>
    constexpr typename ::internals::EnableIf< ::internals::IsUnit<V>::Result || ::internals::IsScalar<V>::Result, Vector<typename ::internals::TransformUnit<U, V>::MultiplyType> >::EnableType operator*(const V& u) const
    {
        using UnitTransformType = ::internals::TransformUnit<U, V>;
        using ReturnType = typename ::internals::TransformUnit<U, V>::MultiplyType;

        return Vector<ReturnType>(UnitTransformType::Multiply(v[0], u), UnitTransformType::Multiply(v[1], u), UnitTransformType::Multiply(v[2], u));
    }

    // Vector<U> * Vector<double> Component multiplication
    constexpr Vector<U> operator*(const Vector<Scalar> &v2) const
    {
        return Vector<U>(v[0] * v2.x(), v[1] * v2.y(), v[2] * v2.z());
    }

    //////////////////////////////////////
    // Division operators
    //////////////////////////////////////
    // Vector<U> / V
    template<typename V>
    constexpr typename ::internals::EnableIf< ::internals::IsUnit<V>::Result || ::internals::IsScalar<V>::Result, Vector<typename ::internals::TransformUnit<U, V>::DivideType> >::EnableType operator/(const V& u) const
    {
        using UnitTransformType = ::internals::TransformUnit<U, V>;
        using ReturnType = typename ::internals::TransformUnit<U, V>::DivideType;

        return Vector<ReturnType>(UnitTransformType::Divide(v[0], u), UnitTransformType::Divide(v[1], u), UnitTransformType::Divide(v[2], u));
    }

    // Vector<U> / double
    template<typename F>
    constexpr typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, Vector<U> >::EnableType operator/(const F& divisor)
    {
        return Vector<U>(v[0] / divisor, v[1] / divisor, v[2] / divisor);
    }

    //////////////////////////////////////
    // Other operator functions
    //////////////////////////////////////
    constexpr Vector<U> operator-() const
    {
        return Vector<U>(-v[0], -v[1], -v[2]);
    }

    //////////////////////////////////////
    // Components functions
    //////////////////////////////////////
    template <typename U2>
    constexpr Vector<U> ParallelComponent(const Vector<U2>& v1) const
    {
        return (dotProduct(v1) / v1.lengthSquared()) * v1;
    }

    template <typename U2>
    constexpr Vector<U> PerpendicularComponent(const Vector<U2>& v1) const
    {
        return *this - ParallelComponent(v1);
    }

    // Transform to scalar vector (remove units)
    constexpr Vector<Scalar> ToScalar() const
    {
        return Vector<Scalar>(Scalar((double)v[0]), Scalar((double)v[1]), Scalar((double)v[2]));
    }

private:
    constexpr Vector(double x, double y, double z, int dummy)
    {
        v[0] = U(x);
        v[1] = U(y);
        v[2] = U(z);
    }

    constexpr double dLenSquare() const
    {
        return (double)v[0] * (double)v[0] + (double)v[1] * (double)v[1] + (double)v[2] * (double)v[2];
    }

    U v[DIMENSION];
};


////////////////////////////////////////////
// Vector global operators
////////////////////////////////////////////
// double * Vector<U>
template<typename F, typename U>
constexpr typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, Vector<U> >::EnableType operator*(const F& factor, const Vector<U> &vector)
{
    return vector * factor;
}

// Unit * Vector<U>
template<typename U, typename V>
constexpr typename ::internals::EnableIf< ::internals::IsUnit<U>::Result || ::internals::IsScalar<U>::Result, Vector<typename ::internals::TransformUnit<V, U>::MultiplyType> >::EnableType operator*(const U& u, const Vector<V>& v)
{
    return v * u;
}


////////////////////////////////////////////
// Transform to absolute value (abs functions)
//    For Unit and Vector
////////////////////////////////////////////
template<typename U>
constexpr typename ::internals::EnableIf< ::internals::IsUnit<U>::Result, U>::EnableType abs(const U& u)
{
    return U(std::abs((double)u));
}

template<typename U>
constexpr Vector<U> abs(const Vector<U>& v)
{
    return Vector<U>(abs(v.x()), abs(v.y()), abs(v.z()));
}


/*template<typename U>
class Matrix
{
public:
    Matrix()
    {
        m[0][0] = U(0.0);
        m[1][0] = U(0.0);
        m[2][0] = U(0.0);

        m[0][1] = U(0.0);
        m[1][1] = U(0.0);
        m[2][1] = U(0.0);

        m[0][2] = U(0.0);
        m[1][2] = U(0.0);
        m[2][2] = U(0.0);
    }

    Matrix(U m11, U m12, U m13, U m21, U m22, U m23, U m31, U m32, U m33)
    {
        m[0][0] = m11;
        m[1][0] = m21;
        m[2][0] = m31;

        m[0][1] = m12;
        m[1][1] = m22;
        m[2][1] = m32;

        m[0][2] = m13;
        m[1][2] = m23;
        m[2][2] = m33;
    }

    template<typename V>
    Matrix(const Matrix<V>& u, typename ::internals::EnableIf< ::internals::IsEquivalent<U, V>::Result>::EnableType* dummy = 0)
    {
        m[0][0] = U(u[0][0]);
        m[1][0] = U(u[1][0]);
        m[2][0] = U(u[2][0]);

        m[0][1] = U(u[0][1]);
        m[1][1] = U(u[1][1]);
        m[2][1] = U(u[2][1]);

        m[0][2] = U(u[0][2]);
        m[1][2] = U(u[1][2]);
        m[2][2] = U(u[2][2]);
    }

    const U* constData() const { return m; }

    bool isIdentity() const
    {
        U tresshold = U(0.0000000001);
        return abs(U(1) - m[0][0]) < tresshold && abs(U(1) - m[1][1]) < tresshold && abs(U(1) - m[2][2]) < tresshold;
    }

    U[]& operator[](unsigned int idx)
    {
        return m[idx];
    }

    const U[]& operator[](unsigned int idx) const
    {
        return m[idx];
    }

    //////////////////////////////////////
    // Assignation operators
    //////////////////////////////////////
    template<typename U2>
    typename ::internals::EnableIf< ::internals::IsEquivalent<U, U2>::Result, Matrix<U> >::EnableType& operator=(const Matrix<U2>& mat)
    {
        m[0][0] = mat[0][0];
        m[1][0] = mat[1][0];
        m[2][0] = mat[2][0];

        m[0][1] = mat[0][1];
        m[1][1] = mat[1][1];
        m[2][1] = mat[2][1];

        m[0][2] = mat[0][2];
        m[1][2] = mat[1][2];
        m[2][2] = mat[2][2];

        return *this;
    }

    //////////////////////////////////////
    // Coumpound assigment operators
    //////////////////////////////////////
    template<typename U2>
    typename ::internals::EnableIf< ::internals::IsEquivalent<U, U2>::Result, Matrix<U> >::EnableType &operator+=(const Matrix<U2> &mat)
    {
        m[0][0] += mat[0][0];
        m[1][0] += mat[1][0];
        m[2][0] += mat[2][0];

        m[0][1] += mat[0][1];
        m[1][1] += mat[1][1];
        m[2][1] += mat[2][1];

        m[0][2] += mat[0][2];
        m[1][2] += mat[1][2];
        m[2][2] += mat[2][2];
        return *this;
    }

    template<typename U2>
    typename ::internals::EnableIf< ::internals::IsEquivalent<U, U2>::Result, Matrix<U> >::EnableType &operator-=(const Matrix<U2> &mat)
    {
        m[0][0] -= mat[0][0];
        m[1][0] -= mat[1][0];
        m[2][0] -= mat[2][0];

        m[0][1] -= mat[0][1];
        m[1][1] -= mat[1][1];
        m[2][1] -= mat[2][1];

        m[0][2] -= mat[0][2];
        m[1][2] -= mat[1][2];
        m[2][2] -= mat[2][2];
        return *this;
    }

    template<typename F>
    typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, Matrix<U>&>::EnableType operator*=(const F& factor)
    {
        m[0][0] *= factor;
        m[1][0] *= factor;
        m[2][0] *= factor;

        m[0][1] *= factor;
        m[1][1] *= factor;
        m[2][1] *= factor;

        m[0][2] *= factor;
        m[1][2] *= factor;
        m[2][2] *= factor;
        return *this;
    }

    template<typename F>
    typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, Vector<U>&>::EnableType operator/=(const F& divisor)
    {
        m[0][0] /= factor;
        m[1][0] /= factor;
        m[2][0] /= factor;

        m[0][1] /= factor;
        m[1][1] /= factor;
        m[2][1] /= factor;

        m[0][2] /= factor;
        m[1][2] /= factor;
        m[2][2] /= factor;
        return *this;
    }

    //////////////////////////////////////
    // Equality operators
    //////////////////////////////////////
    template<typename U2>
    typename ::internals::EnableIf< ::internals::IsEquivalent<U, U2>::Result, bool>::EnableType operator==(const Matrix<U2> &mat) const
    {
        return m[0][0] == mat[0][0] &&
            m[1][0] == mat[1][0] &&
            m[2][0] == mat[2][0] &&

            m[0][1] == mat[0][1] &&
            m[1][1] == mat[1][1] &&
            m[2][1] == mat[2][1] &&

            m[0][2] == mat[0][2] &&
            m[1][2] == mat[1][2] &&
            m[2][2] == mat[2][2];
    }

    template<typename U2>
    typename ::internals::EnableIf< ::internals::IsEquivalent<U, U2>::Result, bool>::EnableType operator!=(const Matrix<U2> &mat) const
    {
        return !(*this == mat);
    }


    //////////////////////////////////////
    // Addition/Substraction operators
    //////////////////////////////////////
    template<typename U2>
    typename ::internals::EnableIf< ::internals::IsEquivalent<U, U2>::Result, Matrix<U> >::EnableType operator+(const Matrix<U2> &mat) const
    {
        return Matrix<U>(m[0][0] + mat[0][0],
            m[1][0] + mat[1][0],
            m[2][0] + mat[2][0],
            m[0][1] + mat[0][1],
            m[1][1] + mat[1][1],
            m[2][1] + mat[2][1],
            m[0][2] + mat[0][2],
            m[1][2] + mat[1][2],
            m[2][2] + mat[2][2]);
    }

    template<typename U2>
    typename ::internals::EnableIf< ::internals::IsEquivalent<U, U2>::Result, Matrix<U> >::EnableType operator-(const Matrix<U2> &v2) const
    {
        return Matrix<U>(m[0][0] - mat[0][0],
            m[1][0] - mat[1][0],
            m[2][0] - mat[2][0],
            m[0][1] - mat[0][1],
            m[1][1] - mat[1][1],
            m[2][1] - mat[2][1],
            m[0][2] - mat[0][2],
            m[1][2] - mat[1][2],
            m[2][2] - mat[2][2]);
    }


    //////////////////////////////////////
    // Multiplication operators
    //////////////////////////////////////
    // Matrix<U> *  double
    template<typename F>
    typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, Vector<U> >::EnableType operator*(const F& factor) const
    {
        return Matrix<U>(m[0][0] * factor,
        m[1][0] * factor,
        m[2][0] * factor,
        m[0][1] * factor,
        m[1][1] * factor,
        m[2][1] * factor,
        m[0][2] * factor,
        m[1][2] * factor,
        m[2][2] * factor);
    }

    // Matrix<U> * V
    template<typename V>
    typename ::internals::EnableIf< ::internals::IsUnit<V>::Result || ::internals::IsScalar<V>::Result, Matrix<typename ::internals::TransformUnit<U, V>::MultiplyType> >::EnableType operator*(const V& u) const
    {
        using UnitTransformType = ::internals::TransformUnit<U, V>;
        using ReturnType = typename ::internals::TransformUnit<U, V>::MultiplyType;

        return Matrix<ReturnType>(UnitTransformType::Multiply(m[0][0], u), UnitTransformType::Multiply(m[0][1], u), UnitTransformType::Multiply(m[0][2], u),
            UnitTransformType::Multiply(m[1][0], u), UnitTransformType::Multiply(m[1][1], u), UnitTransformType::Multiply(m[1][2], u),
            UnitTransformType::Multiply(m[2][0], u), UnitTransformType::Multiply(m[2][1], u), UnitTransformType::Multiply(m[2][2], u));
    }

    // Matrix<U> * Vector<U> 
    Vector<typename ::internals::TransformUnit<U, U>::MultiplyType> operator*(const Vector<U> &v2) const
    {
        using UnitTransformType = ::internals::TransformUnit<U, V>;
        using ReturnType = typename ::internals::TransformUnit<U, V>::MultiplyType;
        return Vector<ReturnType>(UnitTransformType::Multiply(v2.x(), m[0][0]) + UnitTransformType::Multiply(v2.y(), m[0][1]) + UnitTransformType::Multiply(v2.z(), m[0][2]), 
                                  UnitTransformType::Multiply(v2.x(), m[1][0]) + UnitTransformType::Multiply(v2.y(), m[1][1]) + UnitTransformType::Multiply(v2.z(), m[1][2]),
                                  UnitTransformType::Multiply(v2.x(), m[2][0]) + UnitTransformType::Multiply(v2.y(), m[2][1]) + UnitTransformType::Multiply(v2.z(), m[2][2]));
    }

    //////////////////////////////////////
    // Division operators
    //////////////////////////////////////
    // Vector<U> / V
    template<typename V>
    typename ::internals::EnableIf< ::internals::IsUnit<V>::Result || ::internals::IsScalar<V>::Result, Matrix<typename ::internals::TransformUnit<U, V>::DivideType> >::EnableType operator/(const V& u) const
    {
        using UnitTransformType = ::internals::TransformUnit<U, V>;
        using ReturnType = typename ::internals::TransformUnit<U, V>::DivideType;

        return Matrix<ReturnType>(UnitTransformType::Divide(m[0][0], u), UnitTransformType::Divide(m[0][1], u), UnitTransformType::Divide(m[0][2], u),
            UnitTransformType::Divide(m[1][0], u), UnitTransformType::Divide(m[1][1], u), UnitTransformType::Divide(m[1][2], u),
            UnitTransformType::Divide(m[2][0], u), UnitTransformType::Divide(m[2][1], u), UnitTransformType::Divide(m[2][2], u));
    }

    // Matrix<U> / double
    template<typename F>
    typename ::internals::EnableIf< ::internals::IsArithmetic<F>::Result, Matrix<U> >::EnableType operator/(const F& divisor)
    {
        return Matrix<U>(m[0][0] / divisor,
            m[1][0] / divisor,
            m[2][0] / divisor,
            m[0][1] / divisor,
            m[1][1] / divisor,
            m[2][1] / divisor,
            m[0][2] / divisor,
            m[1][2] / divisor,
            m[2][2] / divisor);
    }

    //////////////////////////////////////
    // Other operator functions
    //////////////////////////////////////

    
    // Transform to scalar vector (remove units)
    Matrix<Scalar> ToScalar() const
    {
        return Vector<Scalar>(Scalar((double)m[0][0]), Scalar((double)m[0][1]), Scalar((double)m[0][2]),
            Scalar((double)m[1][0]), Scalar((double)m[1][1]), Scalar((double)m[1][2]),
            Scalar((double)m[2][0]), Scalar((double)m[2][1]), Scalar((double)m[2][2]));
    }

private:
    U[3][3] m;
};*/


////////////////////////////////////////////
// Write to stream
////////////////////////////////////////////
template<typename S, typename U>
S& operator<<(S& os, const Unit<U>& u)
{
    os << u.Value() << ' ';
    ::internals::outputUnit<U>::Output(os);
    return os;
}


template<typename S, typename V>
S& operator<<(S& os, const Vector<V>& v)
{
    os << '(' << (double)v.x();
    os << ", " << (double)v.y();
    os << ", " << (double)v.z();
    os << ')';
    return os;
}

#define UNIT_DISPLAY_NAME( unit, s ) \
    namespace internals { \
        template<> struct outputUnit< ::unit> { \
            template<typename Stream> \
            static void Output(Stream& os) { os << s; } \
        }; \
    }



// Used to create a type multiplying an arbitrary number of types
template<typename U, typename ... Values>
struct Multiply
{
    using Type = typename ::internals::TransformUnit<U, typename Multiply<Values...>::Type >::MultiplyType;
};

template<typename U>
struct Multiply<U>
{
    using Type = U;
};

// Used to create a type that divide two resulting types (can be used with Multiply)
template<typename U, typename V>
struct Divide
{
    using Type = typename ::internals::TransformUnit<U, V>::DivideType;
};

// Used to create a power of a type
template<typename U, int32 Num, int32 Den = 1>
struct Pow
{
    using Type = typename ::internals::PowerTransform<U, Num, Den>::Type;
};

template<typename U, int64 Num, int64 Den = 1>
struct Scale
{
    using Type = ScaledUnit<U, Num, Den>;
};



////////////////////////////////////////////
// SI Prefix definitions
////////////////////////////////////////////
template<typename Unit> struct standard { typedef ScaledUnit<Unit> Type; };

template<typename Unit> struct deca { typedef ScaledUnit<Unit, 10> Type; };
template<typename Unit> struct hecto { typedef ScaledUnit<Unit, 100> Type; };
template<typename Unit> struct kilo { typedef ScaledUnit<Unit, 1000> Type; };
template<typename Unit> struct mega { typedef ScaledUnit<typename kilo<Unit>::Type, 1000> Type; };
template<typename Unit> struct giga { typedef ScaledUnit<typename mega<Unit>::Type, 1000> Type; };
template<typename Unit> struct tera { typedef ScaledUnit<typename giga<Unit>::Type, 1000> Type; };
template<typename Unit> struct peta { typedef ScaledUnit<typename tera<Unit>::Type, 1000> Type; };
template<typename Unit> struct exa { typedef ScaledUnit<typename peta<Unit>::Type, 1000> Type; };
template<typename Unit> struct zetta { typedef ScaledUnit<typename exa<Unit>::Type, 1000> Type; };
template<typename Unit> struct yotta { typedef ScaledUnit<typename zetta<Unit>::Type, 1000> Type; };

template<typename Unit> struct deci { typedef ScaledUnit<Unit, 1, 10> Type; };
template<typename Unit> struct centi { typedef ScaledUnit<Unit, 1, 100> Type; };
template<typename Unit> struct milli { typedef ScaledUnit<Unit, 1, 1000> Type; };
template<typename Unit> struct micro { typedef ScaledUnit<typename milli<Unit>::Type, 1, 1000> Type; };
template<typename Unit> struct nano { typedef ScaledUnit<typename micro<Unit>::Type, 1, 1000> Type; };
template<typename Unit> struct pico { typedef ScaledUnit<typename nano<Unit>::Type, 1, 1000> Type; };
template<typename Unit> struct femto { typedef ScaledUnit<typename pico<Unit>::Type, 1, 1000> Type; };
template<typename Unit> struct atto { typedef ScaledUnit<typename femto<Unit>::Type, 1, 1000> Type; };
template<typename Unit> struct zepto { typedef ScaledUnit<typename atto<Unit>::Type, 1, 1000> Type; };
template<typename Unit> struct yocto { typedef ScaledUnit<typename zepto<Unit>::Type, 1, 1000> Type; };


////////////////////////////////////////////
// SI Base Unit - DO NOT USE DIRECTLY
////////////////////////////////////////////
struct metreBase;
struct secondBase;
struct gBase;
struct moleBase;
struct kelvinBase;
struct ampereBase;
struct candelaBase;

// Generate display
UNIT_DISPLAY_NAME(metreBase, "m")
UNIT_DISPLAY_NAME(secondBase, "s")
UNIT_DISPLAY_NAME(moleBase, "mol")
UNIT_DISPLAY_NAME(kelvinBase, "K")
UNIT_DISPLAY_NAME(gBase, "g")
UNIT_DISPLAY_NAME(ampereBase, "A")
UNIT_DISPLAY_NAME(candelaBase, "cd")

////////////////////////////////////////////
// SI Unit - The one to use
////////////////////////////////////////////


// Distance
using Metre = standard<metreBase>::Type;
// Time
using Second = standard<secondBase>::Type;
// Amount of substance
using Mole = standard<moleBase>::Type;
// Temperature
using Kelvin = standard<kelvinBase>::Type;
// Mass
using Gram = standard<gBase>::Type;
// Electric current
using Ampere = standard<ampereBase>::Type;
// Luminous intensity
using Candela = standard<candelaBase>::Type;

////////////////////////////////////////////
// Scaled SI Unit
////////////////////////////////////////////
// Length units
using Kilometre = kilo<Metre>::Type;
using Decimetre = deci<Metre>::Type;
using Centimetre = centi<Metre>::Type;
using Millimetre = milli<Metre>::Type;
using Micrometre = micro<Metre>::Type;
using Nanometre = nano<Metre>::Type;
using Picometre = pico<Metre>::Type;

// Mass units
using Kilogram = kilo<Gram>::Type;
using Centigram = centi<Gram>::Type;
using Milligram = milli<Gram>::Type;

// Time units
using Attosecond = atto<Second>::Type;
using Femtosecond = femto<Second>::Type;
using Picosecond = pico<Second>::Type;
using Nanosecond = nano<Second>::Type;
using Microsecond = micro<Second>::Type;
using Millisecond = milli<Second>::Type;
using Minute = ScaledUnit<Second, 60, 1>;
using Hour = ScaledUnit<Minute, 60, 1>;
using Day = ScaledUnit<Hour, 24, 1>;

// Electric current
using Kiloampere = kilo<Ampere>::Type;
using Milliampere = milli<Ampere>::Type;

// Temperature
using Kilokelvin = kilo<Kelvin>::Type;
using Millikelvin = milli<Kelvin>::Type;

// Amount of substance
using Kilomole = kilo<Mole>::Type;
using Millimole = milli<Mole>::Type;

// Luminous intensity
using Millicandela = milli<Candela>::Type;


////////////////////////////////////////////
// Powered SI Unit
////////////////////////////////////////////
using Centimetre2 = Pow<Centimetre, 2>::Type;
using Centimetre3 = Pow<Centimetre, 3>::Type;
using Centimetre4 = Pow<Centimetre, 4>::Type;

using CentimetreM1 = Pow<Centimetre, -1>::Type;
using CentimetreM2 = Pow<Centimetre, -2>::Type;
using CentimetreM3 = Pow<Centimetre, -3>::Type;
using CentimetreM4 = Pow<Centimetre, -4>::Type;

using Decimetre2 = Pow<Decimetre, 2>::Type;
using Decimetre3 = Pow<Decimetre, 3>::Type;
using Decimetre4 = Pow<Decimetre, 4>::Type;

using DecimetreM1 = Pow<Decimetre, -1>::Type;
using DecimetreM2 = Pow<Decimetre, -2>::Type;
using DecimetreM3 = Pow<Decimetre, -3>::Type;
using DecimetreM4 = Pow<Decimetre, -4>::Type;

using Metre2 = Pow<Metre, 2>::Type;
using Metre3 = Pow<Metre, 3>::Type;
using Metre4 = Pow<Metre, 4>::Type;

using MetreM1 = Pow<Metre, -1>::Type;
using MetreM2 = Pow<Metre, -2>::Type;
using MetreM3 = Pow<Metre, -3>::Type;
using MetreM4 = Pow<Metre, -4>::Type;

using Kilometre2 = Pow<Kilometre, 2>::Type;
using Kilometre3 = Pow<Kilometre, 3>::Type;
using Kilometre4 = Pow<Kilometre, 4>::Type;

using KilometreM1 = Pow<Kilometre, -1>::Type;
using KilometreM2 = Pow<Kilometre, -2>::Type;
using KilometreM3 = Pow<Kilometre, -3>::Type;
using KilometreM4 = Pow<Kilometre, -4>::Type;

using Millimetre2 = Pow<Millimetre, 2>::Type;
using Millimetre3 = Pow<Millimetre, 3>::Type;
using Millimetre4 = Pow<Millimetre, 4>::Type;

using MillimetreM1 = Pow<Millimetre, -1>::Type;
using MillimetreM2 = Pow<Millimetre, -2>::Type;
using MillimetreM3 = Pow<Millimetre, -3>::Type;
using MillimetreM4 = Pow<Millimetre, -4>::Type;



using Second2 = Pow<Second, 2>::Type;
using Second3 = Pow<Second, 3>::Type;
using Second4 = Pow<Second, 4>::Type;

using SecondM1 = Pow<Second, -1>::Type;
using SecondM2 = Pow<Second, -2>::Type;
using SecondM3 = Pow<Second, -3>::Type;
using SecondM4 = Pow<Second, -4>::Type;

using Millisecond2 = Pow<Millisecond, 2>::Type;
using Millisecond3 = Pow<Millisecond, 3>::Type;

using MillisecondM1 = Pow<Millisecond, -1>::Type;
using MillisecondM2 = Pow<Millisecond, -2>::Type;



using Ampere2 = Pow<Ampere, 2>::Type;
using Ampere3 = Pow<Ampere, 3>::Type;
using Ampere4 = Pow<Ampere, 4>::Type;

using AmpereM1 = Pow<Ampere, -1>::Type;
using AmpereM2 = Pow<Ampere, -2>::Type;
using AmpereM3 = Pow<Ampere, -3>::Type;

////////////////////////////////////////////
// Named Derived SI Unit
////////////////////////////////////////////
using Newton = Multiply<Kilogram, Metre, SecondM2>::Type;   // Force
using Pascal = Multiply<Kilogram, MetreM1, SecondM2>::Type; // Pressure, stress

using Joule = Multiply<Kilogram, Metre2, SecondM2>::Type;   // Energy, work, heat
using Watt = Multiply<Kilogram, Metre2, SecondM3>::Type;   // Power, radiant flux

using Coulomb = Multiply<Ampere, Second>::Type;                                                    // Electric charge or quantity of electricity
using Volt = Divide<Multiply<Kilogram, Metre2>::Type, Multiply<Second3, Ampere>::Type>::Type;   // Voltage (Electrical potential difference), Electromotive force (Joule/Coulomb)
using Farad = Divide<Multiply<Ampere2, Second4>::Type, Multiply<Kilogram, Metre2>::Type>::Type;  // Electric capacitance (Coulomb/Volt)
using Ohm = Divide<Multiply<Kilogram, Metre2>::Type, Multiply<Second3, Ampere2>::Type>::Type;  // Electric resistance, impedance, reactance (Volt/Ampere)
using Siemens = Divide<Multiply<Ampere2, Second3>::Type, Multiply<Kilogram, Metre2>::Type>::Type;  // Electrical conductance (Ampere/Volt)

using Weber = Divide<Kilogram, Multiply<Metre2, Second2, Ampere>::Type>::Type;                    // Magnetic flux
using Tesla = Divide<Kilogram, Multiply<Second2, Ampere>::Type>::Type;                            // Magnetic field strength
using Henry = Divide<Multiply<Kilogram, Metre2>::Type, Multiply<Second2, Ampere2>::Type>::Type;   // Inductance

using Lumen = Candela;                         // Luminous flux
using Lux = Divide<Candela, Metre2>::Type;   // Illuminance

using Becquerel = SecondM1;                      // Radioactivity
using Gray = Divide<Metre2, Second2>::Type; // Absorbed dose (of ionizing radiation)
using Sievert = Divide<Metre2, Second2>::Type; // Equivalent dose (of ionizing radiation)

using Katal = Divide<Mole, Second>::Type;      // Catalytic Activity
using Hertz = SecondM1;                        // Frequency



UNIT_DISPLAY_NAME(Newton, "N")
UNIT_DISPLAY_NAME(Pascal, "Pa")
UNIT_DISPLAY_NAME(Joule, "J")
UNIT_DISPLAY_NAME(Watt, "W")
UNIT_DISPLAY_NAME(Coulomb, "C")
UNIT_DISPLAY_NAME(Volt, "V")
UNIT_DISPLAY_NAME(Farad, "F")
UNIT_DISPLAY_NAME(Ohm, "Ohm")
UNIT_DISPLAY_NAME(Siemens, "S")
UNIT_DISPLAY_NAME(Weber, "Wb")
UNIT_DISPLAY_NAME(Tesla, "T")
UNIT_DISPLAY_NAME(Henry, "H")
UNIT_DISPLAY_NAME(Lux, "lx")
UNIT_DISPLAY_NAME(Gray, "Gy")
UNIT_DISPLAY_NAME(Katal, "kat")

//UNIT_DISPLAY_NAME(Hertz, "Hz")
//UNIT_DISPLAY_NAME(Becquerel, "Bq")
//UNIT_DISPLAY_NAME(Sievert, "Sv")

////////////////////////////////////////////
// Other Composed Unit
////////////////////////////////////////////
using Speed = Divide<Metre, Second>::Type;
using Acceleration = Divide<Metre, Second2>::Type;
using Momentum = Divide<Multiply<Kilogram, Metre>::Type, Second>::Type;
using Density = Divide<Kilogram, Metre3>::Type;

using Angstrom = ScaledUnit<Metre, 1, 10000000000>;


////////////////////////////////////////////
// Angle unit
////////////////////////////////////////////
struct angleBase;

using Radian = standard<angleBase>::Type;
using Degree = ScaledUnit<Radian, 31415926535897932, 1800000000000000000>;

UNIT_DISPLAY_NAME(Radian, "Radian")
UNIT_DISPLAY_NAME(Degree, "Degree")

#endif //_UNIT_H_
