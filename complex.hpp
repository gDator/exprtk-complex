#ifndef COMPLEX_HPP
#define COMPLEX_HPP


#include <math.h>
#include <complex>
#include <string> 

const long double pi = std::atan(1) * 4;
const long double e = std::exp(1.0);

template<typename T>
class ComplexWrapper;
template <typename T> inline ComplexWrapper<T> operator+= (ComplexWrapper<T>  &lhs, ComplexWrapper<T>  const &rhs); 
template <typename T> inline ComplexWrapper<T> operator-= (ComplexWrapper<T>  &lhs, ComplexWrapper<T>  const &rhs);
template <typename T> inline ComplexWrapper<T> operator~ (ComplexWrapper<T>  const &lhs);  //konjugiert komplex
template <typename T> inline ComplexWrapper<T> operator+ (ComplexWrapper<T>  const &lhs, ComplexWrapper<T>  const &rhs);
template <typename T> inline ComplexWrapper<T> operator+ (ComplexWrapper<T>  const &lhs, T const &rhs);  
template <typename T> inline ComplexWrapper<T> operator+ (ComplexWrapper<T>  const &lhs);  
template <typename T> inline ComplexWrapper<T> operator- (ComplexWrapper<T>  const &lhs, ComplexWrapper<T>  const &rhs);
template <typename T> inline ComplexWrapper<T> operator- (ComplexWrapper<T>  const &lhs);
template <typename T> inline ComplexWrapper<T> operator- (ComplexWrapper<T>  const &lhs, T const &rhs);
template <typename T> inline ComplexWrapper<T> operator* (ComplexWrapper<T>  const &lhs, ComplexWrapper<T>  const &rhs);
template <typename T> inline ComplexWrapper<T> operator* (ComplexWrapper<T>  const &lhs, T const &rhs);
template <typename T> inline ComplexWrapper<T> operator* (T const &lhs, ComplexWrapper<T>  const &rhs);  
template <typename T> inline ComplexWrapper<T> operator/ (ComplexWrapper<T>  const &lhs, ComplexWrapper<T>  const &rhs);
template <typename T> inline ComplexWrapper<T> operator/ (ComplexWrapper<T>  const &lhs, T const &rhs); 
template <typename T> inline ComplexWrapper<T> operator/ (T const &lhs, ComplexWrapper<T>  const &rhs);

template <typename T> inline bool operator== (ComplexWrapper<T>  const &lhs, ComplexWrapper<T>  const &rhs);
template <typename T> inline bool operator> (ComplexWrapper<T>  const &lhs, ComplexWrapper<T>  const &rhs);
template <typename T> inline bool operator< (ComplexWrapper<T>  const &lhs, ComplexWrapper<T>  const &rhs);
template <typename T> bool operator>= (ComplexWrapper<T>  const &lhs, ComplexWrapper<T>  const &rhs);
template <typename T> bool operator<= (ComplexWrapper<T>  const &lhs, ComplexWrapper<T>  const &rhs);



template <typename T> void operator<< (std::ostream &s, ComplexWrapper<T> &lhs);

template<typename T>
class ComplexWrapper
{
    private:


    public:
        std::complex<T> std;
    ComplexWrapper(T re = 0, T im = 0) {std.real(re); std.imag(im);}
    // ComplexWrapper(const T re = 0, const T im = 0) {std.real(re); std.imag(im);}
    ComplexWrapper(std::complex<T> c) : std(c) {}
    ComplexWrapper(std::string s)
    {
        std::complex<T> t;
        std::istringstream is('(' + s + ')'); 
        is >> t;    //convert string to complex number
    }
    T real() const {return std.real();}
    T imag() const {return std.imag();}
    void real(T re) {std.real(re);}
    void imag(T im) {std.imag(im);}
    //ComplexWrapper(const ComplexWrapper &c) {return *this} //Kopoerkonstruktor
    //ComplexWrapper(const ComplexWrapper<T> &&c) {this->real = std::move(c.real); this->imag = std::move(c.imag); } //Movekonstruktor
    ~ComplexWrapper() {};

friend inline ComplexWrapper<T>& operator ~(ComplexWrapper<T>  &lhs)  {lhs.imag(-lhs.imag()); return lhs;} //konjugiert komplex
friend inline ComplexWrapper operator+= (ComplexWrapper  &lhs, ComplexWrapper  const &rhs)  {lhs.std += rhs.std; return lhs;}
friend inline ComplexWrapper operator-= (ComplexWrapper  &lhs, ComplexWrapper  const &rhs)  {lhs.std -= rhs.std; return lhs;}
friend inline ComplexWrapper operator+ (ComplexWrapper  const &lhs, ComplexWrapper  const &rhs)  {return ComplexWrapper{lhs.std + rhs.std};}
friend inline ComplexWrapper operator+ (ComplexWrapper  const &lhs, T const &rhs)  {return ComplexWrapper{lhs.real() + rhs, lhs.imag()};}
friend inline ComplexWrapper operator+ (T const &lhs, ComplexWrapper  const &rhs)  {return ComplexWrapper{rhs.real() + lhs, rhs.imag()};}
friend inline ComplexWrapper operator+ (ComplexWrapper  const &lhs)  {return lhs;}
friend inline ComplexWrapper operator- (ComplexWrapper  const &lhs, ComplexWrapper  const &rhs)  {return ComplexWrapper{lhs.std - rhs.std};}
friend inline ComplexWrapper operator- (ComplexWrapper  const &lhs)  {return ComplexWrapper{-lhs.std};}
friend inline ComplexWrapper operator- (ComplexWrapper  const &lhs, T const &rhs)  {return ComplexWrapper{lhs.std - rhs};}
friend inline ComplexWrapper operator* (ComplexWrapper  const &lhs, ComplexWrapper  const &rhs)  {return ComplexWrapper{lhs.std*rhs.std};}
friend inline ComplexWrapper operator* (ComplexWrapper  const &lhs, T const &rhs)  {return ComplexWrapper{lhs.std*rhs};}
friend inline ComplexWrapper operator* (T const &lhs, ComplexWrapper  const &rhs)  {return ComplexWrapper{rhs.std*lhs};}
friend inline ComplexWrapper operator/ (ComplexWrapper const &lhs, ComplexWrapper const &rhs)  {return ComplexWrapper{lhs.std/rhs.std};}
friend inline ComplexWrapper operator/ (ComplexWrapper const &lhs, T const &rhs)  {return ComplexWrapper{lhs.std/rhs};}
friend inline ComplexWrapper operator/ (T const &lhs, ComplexWrapper  const &rhs)  {return ComplexWrapper{lhs/rhs.std};}

friend inline ComplexWrapper& operator*= (ComplexWrapper  &lhs, ComplexWrapper  const &rhs)  {lhs.std *= rhs.std; return lhs;}
friend inline ComplexWrapper& operator*= (ComplexWrapper  &lhs, T const &rhs)  {lhs.std *=rhs; return lhs;}
friend inline ComplexWrapper& operator*= (T &lhs, ComplexWrapper  const &rhs)  {rhs.std *= lhs; return rhs;}

friend inline ComplexWrapper& operator/= (ComplexWrapper  &lhs, ComplexWrapper  const &rhs)  {lhs.std /= rhs.std; return lhs;}
friend inline ComplexWrapper& operator/= (ComplexWrapper  &lhs, T const &rhs)  {lhs.std /= rhs; return lhs;}
friend inline ComplexWrapper& operator/= (T &lhs, ComplexWrapper  const &rhs)  {rhs.std /= lhs; return rhs;}

friend inline bool operator== (ComplexWrapper  const &lhs, ComplexWrapper  const &rhs)  {return lhs.std == rhs.std;}
friend inline bool operator!= (ComplexWrapper  const &lhs, ComplexWrapper  const &rhs)  {return lhs.std != rhs.std;}
friend inline bool operator> (ComplexWrapper  const &lhs, ComplexWrapper const &rhs)  {return abs(lhs.std)>abs(rhs.std);}
friend inline bool operator< (ComplexWrapper  const &lhs, ComplexWrapper const &rhs) {return abs(lhs.std)<abs(rhs.std);}
friend bool operator>= (ComplexWrapper  const &lhs, ComplexWrapper  const &rhs) {return abs(lhs.std)>=abs(rhs.std);}
friend bool operator<= (ComplexWrapper const &lhs, ComplexWrapper  const &rhs) {return abs(lhs.std)<=abs(rhs.std); }
inline ComplexWrapper<T>& operator=(const ComplexWrapper<T> &c) = default;
//inline ComplexWrapper<T>& operator=(const ComplexWrapper<T> &&c) {*this = std::move(c); return *this;}

friend void operator<< (std::ostream &s, ComplexWrapper<T> &lhs) {s << lhs;}
friend void operator>> (std::istringstream &s, ComplexWrapper<T> &rhs) {s >> rhs.std;}
};

namespace complex
{
    namespace details
    {
        namespace constant
        {
            template <typename T> static const T e       =  T(2.718281828459045235360);
            template <typename T> static const T pi      =  T(3.141592653589793238462);
            template <typename T> static const T pi_2    =  T(1.570796326794896619231);
            template <typename T> static const T pi_4    =  T(0.785398163397448309616);
            template <typename T> static const T pi_180  =  T(0.017453292519943295769);
            template <typename T> static const T _1_pi   =  T(0.318309886183790671538);
            template <typename T> static const T _2_pi   =  T(0.636619772367581343076);
            template <typename T> static const T _180_pi = T(57.295779513082320876798);
            template <typename T> static const T log2    =  T(0.693147180559945309417);
            template <typename T> static const T sqrt2   =  T(1.414213562373095048801);
        }
    }
}


template <typename T> static inline T abs(ComplexWrapper<T> v) {return std::abs(v.std);}             //!
template <typename T> static inline T arg(ComplexWrapper<T> v) {return std::arg(v.std);}                        //! cos(phi) = x/r   
template <typename T> static inline ComplexWrapper<T> pow(const ComplexWrapper<T> v, const ComplexWrapper<T> w) {return ComplexWrapper<T>(std::pow(v.std, w.std));} //! (c+id)^(a+ib) = exp(a - ln(r) -b*phi)* exp(j(a*phi + b* ln r))
template <typename T> static inline ComplexWrapper<T> pow(const ComplexWrapper<T> v, const T w) {return ComplexWrapper<T>(std::pow(v.std, w));} //! (c+id)^(a+ib) = exp(a - ln(r) -b*phi)* exp(j(a*phi + b* ln r))
template <typename T> static inline ComplexWrapper<T> pow(const ComplexWrapper<T> v, const int w) {return ComplexWrapper<T>(std::pow(v.std, w));} //! (c+id)^(a+ib) = exp(a - ln(r) -b*phi)* exp(j(a*phi + b* ln r))
template <typename T> static inline ComplexWrapper<T> pow2(const ComplexWrapper<T> v) {return ComplexWrapper<T>(std::pow(v.std, 2)); } //
template <typename T> static inline ComplexWrapper<T>   acos(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::acos(v.std)); }          //!
template <typename T> static inline ComplexWrapper<T>  acosh(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::acosh(v.std)); }  //!
template <typename T> static inline ComplexWrapper<T>   asin(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::asin(v.std)); } //!
template <typename T> static inline ComplexWrapper<T>  asinh(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::asinh(v.std)); } //!
template <typename T> static inline ComplexWrapper<T>   atan(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::atan(v.std));} //!
template <typename T> static inline ComplexWrapper<T>  atanh(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::atanh(v.std));} //!
template <typename T> static inline ComplexWrapper<T>   ceil(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::ceil (v.real()), std::ceil(v.imag()));} //!
template <typename T> static inline ComplexWrapper<T>    cos(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::cos(v.std));} //!
template <typename T> static inline ComplexWrapper<T>   cosh(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::cosh(v.std)); }  //!cosh z = cos(iz)
template <typename T> static inline ComplexWrapper<T>    exp(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::exp(v.std));}//pow(ComplexWrapper<T>{complex::details::constant::e<T>, 0} , v); } //!
template <typename T> static inline ComplexWrapper<T>  floor(const ComplexWrapper<T> v) { return ComplexWrapper{std::floor(v.real()), std::floor(v.imag())};} //!
template <typename T> static inline ComplexWrapper<T>    ln(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::log(v.std)); } //!
template <typename T> static inline ComplexWrapper<T>  log_10(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::log10(v.std)); } //!
template <typename T> static inline ComplexWrapper<T>   log_2(const ComplexWrapper<T> v) { return ln(v)/std::log(T(2)); } //!
template <typename T> static inline ComplexWrapper<T>   modulus(const ComplexWrapper<T> v) { return ComplexWrapper<T>{std::abs(v.real()), std::abs(v.imag())}; } //!
template <typename T> static inline ComplexWrapper<T>   neg(const ComplexWrapper<T> v) { return -v;            }        //!
template <typename T> static inline ComplexWrapper<T>    pos(const ComplexWrapper<T> v) { return  v;            }       //!
template <typename T> static inline ComplexWrapper<T>    sin(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::sin(v.std)); }       //!
template <typename T> static inline ComplexWrapper<T>   sinh(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::sinh(v.std)); } //!
template <typename T> static inline ComplexWrapper<T>   sqrt(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::sqrt(v.std));} //sqrt(r) * (cos(phi/2)+i*sin(phi/2)) //!
template <typename T> static inline ComplexWrapper<T>    tan(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::tan(v.std)); } //!
template <typename T> static inline ComplexWrapper<T>   tanh(const ComplexWrapper<T> v) { return ComplexWrapper<T>(std::tanh(v.std));} //!
template <typename T> static inline ComplexWrapper<T>    cot(const ComplexWrapper<T> v) { return ComplexWrapper<T>{1, 0} / tan(v.std); } //!
template <typename T> static inline ComplexWrapper<T>    sec(const ComplexWrapper<T> v) { return ComplexWrapper<T>{1, 0} / cos(v.std); } //!
template <typename T> static inline ComplexWrapper<T>    csc(const ComplexWrapper<T> v) { return ComplexWrapper<T>{1, 0} / sin(v.std); } //!
template <typename T> static inline ComplexWrapper<T>    r2d(const ComplexWrapper<T> v) { return arg(v.std)*pi/180;} //!
template <typename T> static inline ComplexWrapper<T>    d2r(const ComplexWrapper<T> v) { return arg(v.std);} //?
template <typename T> static inline ComplexWrapper<T>    d2g(const ComplexWrapper<T> v) { return arg(v.std); } //?
template <typename T> static inline ComplexWrapper<T>    g2d(const ComplexWrapper<T> v) { return arg(v.std); } //?
template <typename T> static inline ComplexWrapper<T>   notl(const ComplexWrapper<T> v) { return (v != ComplexWrapper<T>{0, 0} ? ComplexWrapper<T>{0, 0} : ComplexWrapper<T>{1, 0}); } //!
template <typename T> static inline ComplexWrapper<T>   frac(const ComplexWrapper<T> v) { return (v - ComplexWrapper<T>{std::trunc(v.std.real()), std::trunc(v.std.imag())});}  //!marked as const because std::complex real and imag are const
template <typename T> static inline ComplexWrapper<T>  trunc(const ComplexWrapper<T> v) { return ComplexWrapper<T>{std::trunc(v.std.real()), std::trunc(v.std.imag())};} //!
template <typename T> static inline std::string      asString(const ComplexWrapper<T> v) {return std::string (std::to_string(v.real()) + " + " +  std::to_string(v.imag()) + "i");} //!
template <typename T> static inline bool             is_true(const ComplexWrapper<T> v) {return v.abs() != 0;} //!
template <typename T> static inline bool             is_false(const ComplexWrapper<T> v) {return v.abs() == 0;} //!


namespace std
{
   template <typename T>
   class numeric_limits<ComplexWrapper<T>>
   {
    
   public:

      static const bool is_specialized = true;
      static ComplexWrapper<T> (min) ()        { return ComplexWrapper<T>(std::numeric_limits<double>::min(),std::numeric_limits<double>::min()) ;}
      static ComplexWrapper<T> (max) ()        { return ComplexWrapper<T>(std::numeric_limits<double>::max(), std::numeric_limits<double>::max()) ;}
      static ComplexWrapper<T> lowest()        { return -(max)();}
      static ComplexWrapper<T> epsilon()       { return  ComplexWrapper<T>(std::numeric_limits<double>::epsilon(), std::numeric_limits<double>::epsilon());}
      static ComplexWrapper<T> round_error()   { return ComplexWrapper<T>(std::numeric_limits<double>::round_error(), std::numeric_limits<double>::round_error());}
      static ComplexWrapper<T> infinity()      { return ComplexWrapper<T>(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());}
      static ComplexWrapper<T> quiet_NaN()     { return ComplexWrapper<T>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());}
      static ComplexWrapper<T> signaling_NaN() { return ComplexWrapper<T>(std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN()); }
      static ComplexWrapper<T> denorm_min()    { return ComplexWrapper<T>(std::numeric_limits<double>::denorm_min(), std::numeric_limits<double>::denorm_min());}
      static const int digits             = std::numeric_limits<double>::digits;
      static const int digits10           = std::numeric_limits<double>::digits10;
      static const int radix              = std::numeric_limits<double>::radix;
      static const int min_exponent       = std::numeric_limits<double>::min_exponent;
      static const int min_exponent10     = std::numeric_limits<double>::min_exponent10;
      static const int max_exponent       = std::numeric_limits<double>::max_exponent;
      static const int max_exponent10     = std::numeric_limits<double>::max_exponent10;
      static const bool has_infinity      = std::numeric_limits<double>::has_infinity;
      static const bool has_quiet_NaN     = std::numeric_limits<double>::has_quiet_NaN;
      static const bool has_signaling_NaN = std::numeric_limits<double>::has_signaling_NaN;
      static const bool has_denorm_loss   = std::numeric_limits<double>::has_denorm_loss;
      static const bool is_signed         = std::numeric_limits<double>::is_signed;
      static const bool is_integer        = std::numeric_limits<double>::is_integer;
      static const bool is_exact          = std::numeric_limits<double>::is_exact;
      static const bool is_iec559         = std::numeric_limits<double>::is_iec559;
      static const bool is_bounded        = std::numeric_limits<double>::is_bounded;
      static const bool is_modulo         = std::numeric_limits<double>::is_modulo;
      static const bool traps             = std::numeric_limits<double>::traps;
      static const float_denorm_style has_denorm = std::numeric_limits<double>::has_denorm;
      static const float_round_style round_style = std::numeric_limits<double>::round_style;
   };
}





#endif //!COMPLEX_HPP