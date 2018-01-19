#ifndef VEC_H
#define VEC_H

#include <cassert>
#include <cmath>
#include <utility>
#include <algorithm>
#include <numeric>
#include <functional>
#include "hashtable.h"
#include <iostream>
#include "util.h"

// Defines a thin wrapper around fixed size C-style arrays, using template parameters,
// which is useful for dealing with vectors of different dimensions.
// For example, float[3] is equivalent to Vec<3,float>.
// Entries in the vector are accessed with the overloaded [] operator, so
// for example if x is a Vec<3,float>, then the middle entry is x[1].
// For convenience, there are a number of typedefs for abbreviation:
//   Vec<3,float> -> Vec3f
//   Vec<2,int>   -> Vec2i
// and so on.
// Arithmetic operators are appropriately overloaded, and functions are defined
// for additional operations (such as dot-products, norms, cross-products, etc.)

template<unsigned int N, class T>
struct Vec {
    T v[N];
    
    template <typename... Args, std::enable_if_t<sizeof...(Args) == N,int> = 0, std::enable_if_t<(std::is_convertible_v<Args,T> && ... ),int> = 0>
    Vec(Args&&... args): v{std::forward<Args>(args)...} {}
    Vec(void) {}
    template<class S>
    explicit Vec(const S *source) { _assign(_is(),source); }
    template <class S>
    explicit Vec(const Vec<N,S>& source) { _assign(_is(),source.v); }
    explicit Vec(const T& v) { nullaryOpAssign([&v]() -> T {return v;});}

    template <int... M>
        using IS = typename std::integer_sequence<int,M...>;
    constexpr static auto _is() { return std::make_integer_sequence<int,N>(); }


    template <int... M, typename... Args, std::enable_if_t<sizeof...(Args) == N,int> = 0>
    Vec& _assign(IS<M...>, Args... args) { (std::exchange(v[M],args),...); return *this; }
    template <int... M, typename S>
    Vec& _assign(IS<M...>, const S* s) { (std::exchange(v[M],s[M]),...); return *this; }

    template <typename... Args, std::enable_if_t<sizeof...(Args) == N,int> = 0, std::enable_if_t<(std::is_convertible_v<Args,T> && ... ),int> = 0>
    Vec& assign(Args... args) { return _assign(_is(),std::forward<Args>(args)...); }
    template <typename S, std::enable_if_t<std::is_convertible_v<S,T>,int> = 0>
    Vec<N,T>& assign(const Vec<N,S>& o) { return _assign(_is(),o.v);}

    //template <int... M, typename UnaryOp>
    //Vec& nullaryOpAssign(IS<M...>, UnaryOp&& op) { return _assign(op(),...); }
    template <typename UnaryOp>
    //Vec& nullaryOpAssign(UnaryOp&& op) { return nullaryOpAssign(_is(),std::forward<UnaryOp>(op)); }
    Vec& nullaryOpAssign(UnaryOp&& op) { for(int i = 0; i < N; ++i) { v[i] = op(); }; return *this; }



    template <int... M, typename UnaryOp>
    Vec unaryOp(IS<M...>, UnaryOp&& op) const { return Vec(op(v[M])...); }
    template <typename UnaryOp>
    Vec unaryOp(UnaryOp&& op) const { return unaryOp(_is(),std::forward<UnaryOp>(op)); }
    template <int... M, typename UnaryOp>
    Vec& unaryOpAssign(IS<M...>, UnaryOp&& op) { return assign(op(v[M])...); }
    template <typename UnaryOp>
    Vec& unaryOpAssign(UnaryOp&& op) { return unaryOpAssign(_is(),std::forward<UnaryOp>(op)); }


    template <int... M, typename BinaryOp>
    Vec binaryOp(IS<M...>, BinaryOp&& op, const Vec& o) const { return Vec(op(v[M],o[M])...); }
    template <typename BinaryOp>
    Vec binaryOp(BinaryOp&& op, const Vec& o) const { return binaryOp(_is(),std::forward<BinaryOp>(op),o); }

    template <int... M, typename BinaryOp>
    Vec& binaryOpAssign(IS<M...>, BinaryOp&& op, const Vec& o) { return assign(op(v[M],o[M])...); }
    template <typename BinaryOp>
    Vec& binaryOpAssign(BinaryOp&& op, const Vec& o) { return binaryOpAssign(_is(),std::forward<BinaryOp>(op),o); }



    
    
    
    
    
    T &operator[](int index)
    {
        return v[index];
    }
    
    const T &operator[](int index) const
    {
        return v[index];
    }
    
    bool nonzero(void) const
    {
        for(unsigned int i=0; i<N; ++i) if(v[i]) return true;
        return false;
    }
    
    Vec operator+(const Vec& o) const { return binaryOp(std::plus<T>(), o); }
    Vec& operator+=(const Vec& o) { return binaryOpAssign(std::plus<T>(), o); }

    Vec operator-(const Vec& o) const { return binaryOp(std::minus<T>(), o); }
    Vec& operator-=(const Vec& o) { return binaryOpAssign(std::minus<T>(), o); }

    
    Vec operator-(void) const { return unaryOp(std::negate<T>()); }
    
    Vec operator*(const T& s) const { return unaryOp([&s](const T& a) -> T { return s*a; }); }
    Vec& operator*=(const T& s) { return unaryOpAssign([&s](const T& a)-> T { return s*a; }); }

    
    Vec operator*(const Vec& o) const { return binaryOp(std::multiplies<T>(), o); }

    Vec& operator*=(const Vec& o) { return binaryOpAssign(std::multiplies<T>(), o); }
    
    
    Vec operator/(const T& s) const { return unaryOp([&s](const T& a)-> T { return a/s; }); }
    Vec& operator/=(const T& s) { return unaryOpAssign([&s](const T& a)-> T { return a/s; }); }
};

typedef Vec<2,double>         Vec2d;
typedef Vec<2,float>          Vec2f;
typedef Vec<2,int>            Vec2i;
typedef Vec<2,unsigned int>   Vec2ui;
typedef Vec<2,short>          Vec2s;
typedef Vec<2,unsigned short> Vec2us;
typedef Vec<2,char>           Vec2c;
typedef Vec<2,unsigned char>  Vec2uc;
typedef Vec<2,size_t>         Vec2st;

typedef Vec<3,double>         Vec3d;
typedef Vec<3,float>          Vec3f;
typedef Vec<3,int>            Vec3i;
typedef Vec<3,unsigned int>   Vec3ui;
typedef Vec<3,short>          Vec3s;
typedef Vec<3,unsigned short> Vec3us;
typedef Vec<3,char>           Vec3c;
typedef Vec<3,unsigned char>  Vec3uc;
typedef Vec<3,size_t>         Vec3st;

typedef Vec<4,double>         Vec4d;
typedef Vec<4,float>          Vec4f;
typedef Vec<4,int>            Vec4i;
typedef Vec<4,unsigned int>   Vec4ui;
typedef Vec<4,short>          Vec4s;
typedef Vec<4,unsigned short> Vec4us;
typedef Vec<4,char>           Vec4c;
typedef Vec<4,unsigned char>  Vec4uc;
typedef Vec<4,size_t>         Vec4st;

typedef Vec<6,double>         Vec6d;
typedef Vec<6,float>          Vec6f;
typedef Vec<6,unsigned int>   Vec6ui;
typedef Vec<6,int>            Vec6i;
typedef Vec<6,short>          Vec6s;
typedef Vec<6,unsigned short> Vec6us;
typedef Vec<6,char>           Vec6c;
typedef Vec<6,unsigned char>  Vec6uc;


template<unsigned int N, class T>
T mag2(const Vec<N,T> &a)
{
    T l=sqr(a.v[0]);
    for(unsigned int i=1; i<N; ++i) l+=sqr(a.v[i]);
    return l;
}

template<unsigned int N, class T>
T mag(const Vec<N,T> &a)
{ return (T)sqrt(mag2(a)); }

template<unsigned int N, class T> 
inline T dist2(const Vec<N,T> &a, const Vec<N,T> &b)
{ 
    T d=sqr(a.v[0]-b.v[0]);
    for(unsigned int i=1; i<N; ++i) d+=sqr(a.v[i]-b.v[i]);
    return d;
}

template<unsigned int N, class T> 
inline T dist(const Vec<N,T> &a, const Vec<N,T> &b)
{ return std::sqrt(dist2(a,b)); }

template<unsigned int N, class T> 
inline void normalize(Vec<N,T> &a)
{ a/=mag(a); }

template<unsigned int N, class T> 
inline Vec<N,T> normalized(const Vec<N,T> &a)
{ return a/mag(a); }

template<unsigned int N, class T> 
inline T infnorm(const Vec<N,T> &a)
{
    T d=std::fabs(a.v[0]);
    for(unsigned int i=1; i<N; ++i) d=max(std::fabs(a.v[i]),d);
    return d;
}

template<unsigned int N, class T>
void zero(Vec<N,T> &a)
{ 
    for(unsigned int i=0; i<N; ++i)
        a.v[i] = 0;
}

template<unsigned int N, class T>
std::ostream &operator<<(std::ostream &out, const Vec<N,T> &v)
{
    out<<v.v[0];
    for(unsigned int i=1; i<N; ++i)
        out<<' '<<v.v[i];
    return out;
}

template<unsigned int N, class T>
std::istream &operator>>(std::istream &in, Vec<N,T> &v)
{
    in>>v.v[0];
    for(unsigned int i=1; i<N; ++i)
        in>>v.v[i];
    return in;
}

template<unsigned int N, class T> 
inline bool operator==(const Vec<N,T> &a, const Vec<N,T> &b)
{ 
    bool t = (a.v[0] == b.v[0]);
    unsigned int i=1;
    while(i<N && t) {
        t = t && (a.v[i]==b.v[i]); 
        ++i;
    }
    return t;
}

template<unsigned int N, class T> 
inline bool operator!=(const Vec<N,T> &a, const Vec<N,T> &b)
{ 
    bool t = (a.v[0] != b.v[0]);
    unsigned int i=1;
    while(i<N && !t) {
        t = t || (a.v[i]!=b.v[i]); 
        ++i;
    }
    return t;
}

template<unsigned int N, class T>
inline Vec<N,T> operator*(T a, const Vec<N,T> &v)
{
    return v.unaryOp([&a](const T& v)-> T { return a*v; });
}

template<unsigned int N, class T>
inline T min(const Vec<N,T> &a)
{
    return std::min_element(a.v,a.v+N);
}

template<unsigned int N, class T>
inline Vec<N,T> min_union(const Vec<N,T> &a, const Vec<N,T> &b)
{
    return b.binaryOp([](const T& a, const T& b) -> T{ return std::min(a,b); },a);
}

template<unsigned int N, class T>
inline Vec<N,T> max_union(const Vec<N,T> &a, const Vec<N,T> &b)
{
    return b.binaryOp([](const T& a, const T& b) -> T { return std::max(a,b); },a);
}

template<unsigned int N, class T>
inline T max(const Vec<N,T> &a)
{
    return std::max_element(a.v,a.v+N);
}

template<unsigned int N, class T>
inline T dot(const Vec<N,T> &a, const Vec<N,T> &b)
{
    return std::inner_product(a.v,a.v+N,b.v,T{0});
}

template<class T> 
inline Vec<2,T> rotate(const Vec<2,T>& a, const T& angle) 
{
    T c = cos(angle);
    T s = sin(angle);
    return Vec<2,T>(c*a[0] - s*a[1],s*a[0] + c*a[1]); // counter-clockwise rotation
}

// Rotate the point (x,y,z) around the vector (u,v,w)
template<class T> 
inline Vec<3,T> rotate( const Vec<3,T>& x, const T& angle, const Vec<3,T>& u )
{
    T ux = u[0]*x[0];
    T uy = u[0]*x[1];
    T uz = u[0]*x[2];
    T vx = u[1]*x[0];
    T vy = u[1]*x[1];
    T vz = u[1]*x[2];
    T wx = u[2]*x[0];
    T wy = u[2]*x[1];
    T wz = u[2]*x[2];
    
    T sa = (T) sin(angle);
    T ca = (T) cos(angle);
    
    return Vec<3,T> ( u[0] * (ux+vy+wz) + (x[0]*(u[1]*u[1]+u[2]*u[2])-u[0]*(vy+wz))*ca+(-wy+vz)*sa,
                     u[1] * (ux+vy+wz) + (x[1]*(u[0]*u[0]+u[2]*u[2])-u[1]*(ux+wz))*ca+(wx-uz)*sa,
                     u[2] * (ux+vy+wz) + (x[2]*(u[0]*u[0]+u[1]*u[1])-u[2]*(ux+vy))*ca+(-vx+uy)*sa );
    
}

template<class T>
inline Vec<2,T> perp(const Vec<2,T> &a)
{ return Vec<2,T>(-a.v[1], a.v[0]); } // counter-clockwise rotation by 90 degrees

template<class T>
inline T cross(const Vec<2,T> &a, const Vec<2,T> &b)
{ return a.v[0]*b.v[1]-a.v[1]*b.v[0]; }

template<class T>
inline Vec<3,T> cross(const Vec<3,T> &a, const Vec<3,T> &b)
{ return Vec<3,T>(a.v[1]*b.v[2]-a.v[2]*b.v[1], a.v[2]*b.v[0]-a.v[0]*b.v[2], a.v[0]*b.v[1]-a.v[1]*b.v[0]); }

template<class T>
inline T triple(const Vec<3,T> &a, const Vec<3,T> &b, const Vec<3,T> &c)
{ return a.v[0]*(b.v[1]*c.v[2]-b.v[2]*c.v[1])
    +a.v[1]*(b.v[2]*c.v[0]-b.v[0]*c.v[2])
    +a.v[2]*(b.v[0]*c.v[1]-b.v[1]*c.v[0]); }

template<unsigned int N, class T>
inline unsigned int hash(const Vec<N,T> &a)
{
    unsigned int h=a.v[0];
    for(unsigned int i=1; i<N; ++i)
        h=hash(h ^ a.v[i]);
    return h;
}

template<unsigned int N, class T>
inline void assign(const Vec<N,T> &a, T &a0, T &a1)
{ 
    assert(N==2);
    a0=a.v[0]; a1=a.v[1];
}

template<unsigned int N, class T>
inline void assign(const Vec<N,T> &a, T &a0, T &a1, T &a2)
{ 
    assert(N==3);
    a0=a.v[0]; a1=a.v[1]; a2=a.v[2];
}

template<unsigned int N, class T>
inline void assign(const Vec<N,T> &a, T &a0, T &a1, T &a2, T &a3)
{ 
    assert(N==4);
    a0=a.v[0]; a1=a.v[1]; a2=a.v[2]; a3=a.v[3];
}

template<unsigned int N, class T>
inline void assign(const Vec<N,T> &a, T &a0, T &a1, T &a2, T &a3, T &a4, T &a5)
{ 
    assert(N==6);
    a0=a.v[0]; a1=a.v[1]; a2=a.v[2]; a3=a.v[3]; a4=a.v[4]; a5=a.v[5];
}

template<unsigned int N, class T>
inline Vec<N,int> round(const Vec<N,T> &a)
{ 
    Vec<N,int> rounded;
    for(unsigned int i=0; i<N; ++i)
        rounded.v[i]=lround(a.v[i]);
    return rounded; 
}

template<unsigned int N, class T>
inline Vec<N,int> floor(const Vec<N,T> &a)
{ 
    Vec<N,int> rounded;
    for(unsigned int i=0; i<N; ++i)
        rounded.v[i]=(int)floor(a.v[i]);
    return rounded; 
}

template<unsigned int N, class T>
inline Vec<N,int> ceil(const Vec<N,T> &a)
{ 
    Vec<N,int> rounded;
    for(unsigned int i=0; i<N; ++i)
        rounded.v[i]=(int)ceil(a.v[i]);
    return rounded; 
}

template<unsigned int N, class T>
inline Vec<N,T> fabs(const Vec<N,T> &a)
{ 
    Vec<N,T> result;
    for(unsigned int i=0; i<N; ++i)
        result.v[i]=(T)fabs(a.v[i]);
    return result; 
}

template<unsigned int N, class T>
inline void minmax(const Vec<N,T> &x0, const Vec<N,T> &x1, Vec<N,T> &xmin, Vec<N,T> &xmax)
{
    for(unsigned int i=0; i<N; ++i)
        minmax(x0.v[i], x1.v[i], xmin.v[i], xmax.v[i]);
}

template<unsigned int N, class T>
inline void minmax(const Vec<N,T> &x0, const Vec<N,T> &x1, const Vec<N,T> &x2, Vec<N,T> &xmin, Vec<N,T> &xmax)
{
    for(unsigned int i=0; i<N; ++i)
        minmax(x0.v[i], x1.v[i], x2.v[i], xmin.v[i], xmax.v[i]);
}

template<unsigned int N, class T>
inline void minmax(const Vec<N,T> &x0, const Vec<N,T> &x1, const Vec<N,T> &x2, const Vec<N,T> &x3,
                   Vec<N,T> &xmin, Vec<N,T> &xmax)
{
    for(unsigned int i=0; i<N; ++i)
        minmax(x0.v[i], x1.v[i], x2.v[i], x3.v[i], xmin.v[i], xmax.v[i]);
}

template<unsigned int N, class T>
inline void minmax(const Vec<N,T> &x0, const Vec<N,T> &x1, const Vec<N,T> &x2, const Vec<N,T> &x3, const Vec<N,T> &x4,
                   Vec<N,T> &xmin, Vec<N,T> &xmax)
{
    for(unsigned int i=0; i<N; ++i)
        minmax(x0.v[i], x1.v[i], x2.v[i], x3.v[i], x4.v[i], xmin.v[i], xmax.v[i]);
}

template<unsigned int N, class T>
inline void minmax(const Vec<N,T> &x0, const Vec<N,T> &x1, const Vec<N,T> &x2, const Vec<N,T> &x3, const Vec<N,T> &x4,
                   const Vec<N,T> &x5, Vec<N,T> &xmin, Vec<N,T> &xmax)
{
    for(unsigned int i=0; i<N; ++i)
        minmax(x0.v[i], x1.v[i], x2.v[i], x3.v[i], x4.v[i], x5.v[i], xmin.v[i], xmax.v[i]);
}

template<unsigned int N, class T>
inline void update_minmax(const Vec<N,T> &x, Vec<N,T> &xmin, Vec<N,T> &xmax)
{
    for(unsigned int i=0; i<N; ++i) update_minmax(x[i], xmin[i], xmax[i]);
}

#endif
