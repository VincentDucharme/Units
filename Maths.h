#ifndef _MATHHELPER_H_
#define _MATHHELPER_H_

#include "Units.h"

#include <cmath>
#include <algorithm>

class Maths
{
public:
    template<typename T>
    static T Min(const T& a, const T& b)
    {
        return std::min(a, b);
    }

    template<typename T1, typename T2, typename = StaticUtilities::EnableIf<units::IsEquivalent<T1, T2> > >
    static T1 Min(const T1& a, const T2& b)
    {
        return Min(a, T1(b));
    }

    template<typename T>
    static T Max(const T& a, const T& b)
    {
        return std::max(a, b);
    }

    template<typename T1, typename T2, typename = StaticUtilities::EnableIf<units::IsEquivalent<T1, T2> > >
    static T1 Max(const T1& a, const T2& b)
    {
        return Max(a, T1(b));
    }

    template<typename T>
    static T Pow(const T& v, double p)
    {
        return std::pow(v, p);
    }

    template<typename T, typename = StaticUtilities::EnableIf<StaticUtilities::IsArithmetic<T> > >
    static T Sqrt(const T& v)
    {
        return std::sqrt(v);
    }

    template<int PowNum, int PowDen = 1, typename T>
    static units::Pow<T, PowNum, PowDen> Pow(const Unit<T>& v)
    {
        return pow<PowNum, PowDen>(v);
    }

    template<typename T>
    static units::Sqrt<T> Sqrt(const Unit<T>& v)
    {
        return sqrt(v);
    }

    template<typename T>
    static units::Pow<T, 1, 3> Cbrt(const Unit<T>& v)
    {
        return units::Pow<T, 1, 3>(std::cbrt(v.Value()));
    }

    static units::GetValueType<Radian> Sin(const Radian& r)
    {
        return std::sin(r.Value());
    }

    static Radian Arcsin(const Metre& op, const Metre& hypo)
    {
        return Radian(std::asin(op / hypo));
    }

    static Radian Arcsin(units::GetValueType<Radian> ratio)
    {
        return Radian(std::asin(ratio));
    }

    static units::GetValueType<Radian> Cos(const Radian& r)
    {
        return std::cos(r.Value());
    }

    static Radian Arccos(const Metre& adj, const Metre& hypo)
    {
        return Radian(std::acos(adj / hypo));
    }

    static Radian Arccos(units::GetValueType<Radian> ratio)
    {
        return Radian(std::acos(ratio));
    }

    static units::GetValueType<Radian> Tan(const Radian& r)
    {
        return std::tan(r.Value());
    }

    static Radian Arctan(units::GetValueType<Radian> ratio)
    {
        return Radian(std::atan(ratio));
    }

    static Radian Arctan(const Metre& op, const Metre& adj)
    {
        return Radian(std::atan2(op.Value(), adj.Value()));
    }

    template<typename T>
    static T Ceil(const T& v)
    {
        return std::ceil(v);
    }

    template<typename T>
    static T Ceil(const Unit<T>& v)
    {
        return T(std::ceil(v.Value()));
    }

    template<typename T>
    static T Floor(const T& v)
    {
        return std::floor(v);
    }

    template<typename T>
    static T Floor(const Unit<T>& v)
    {
        return T(std::floor(v.Value()));
    }
};

#endif //_MATHHELPER_H_