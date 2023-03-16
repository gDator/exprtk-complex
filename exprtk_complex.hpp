/*
 **************************************************************
 *         C++ Mathematical Expression Toolkit Library        *
 *                                                            *
 * Custom std type Adaptor                                   *
 * Authors: Arash Partow                                      *
 * URL: http://www.partow.net/programming/exprtk/index.html   *
 *                                                            *
 * Copyright notice:                                          *
 * Free use of the Mathematical Expression Toolkit Library is *
 * permitted under the guidelines and in accordance with the  *
 * most current version of the MIT License.                   *
 * http://www.opensource.org/licenses/MIT                     *
 *                                                            *
 **************************************************************
*/
/*
******************************************************************
*           C++ Mathematical Expression Toolkit Library          *
*                           Complex Extension                    *
*                                                                *
* Author: Daniel Hagen (2022)                                    *
* URL: https://www.partow.net/programming/exprtk/index.html      *
*                                                                *
* Copyright notice:                                              *
* Free use of the C++ Mathematical Expression Toolkit Library is *
* permitted under the guidelines and in accordance with the most *
* current version of the MIT License.                            *
* https://www.opensource.org/licenses/MIT                        *
*                                                                *
* Example expressions:                                           *
* (00) (y + x / y) * (x - y / x)                                 *
* (01) (x^2 / sin(2 * pi / y)) - x / 2                           *
* (02) sqrt(1 - (x^2))                                           *
* (03) 1 - sin(2 * x) + cos(pi / y)                              *
* (04) a * exp(2 * t) + c                                        *
* (05) if(((x + 2) == 3) and ((y + 5) <= 9),1 + w, 2 / z)        *
* (06) (avg(x,y) <= x + y ? x - y : x * y) + 2 * pi / x          *
* (07) z := x + sin(2 * pi / y)                                  *
* (08) u := 2 * (pi * z) / (w := x + cos(y / pi))                *
* (09) clamp(-1,sin(2 * pi * x) + cos(y / 2 * pi),+1)            *
* (10) inrange(-2,m,+2) == if(({-2 <= m} and [m <= +2]),1,0)     *
* (11) (2sin(x)cos(2y)7 + 1) == (2 * sin(x) * cos(2*y) * 7 + 1)  *
* (12) (x ilike 's*ri?g') and [y < (3 z^7 + w)]                  *
*                                                                *
*******************************************************************/
#include <string>
#include "complex.hpp"

#ifndef EXPRTK_COMPLEX_HPP
#define EXPRTK_COMPLEX_HPP
// #ifdef EXPRTK_COMPLEX

typedef ComplexWrapper<double> Complexd;
typedef ComplexWrapper<float> Complexf;


namespace exprtk
{
   namespace details
   {
      namespace numeric { namespace details
      {
         struct complex_type;

         template <typename T> inline T const_pi_impl(complex_type);
         template <typename T> inline T const_e_impl (complex_type);
      }}

      template <typename T> inline bool is_true (const ComplexWrapper<T> v);
      template <typename T> inline bool is_false(const ComplexWrapper<T> v);

      template <typename Iterator, typename T>
      inline bool string_to_real(Iterator& itr_external, const Iterator end, ComplexWrapper<T>& t, details::numeric::details::complex_type);
   }

   namespace rtl { namespace io
   {
      namespace details
      {
         template <typename T> inline void print_type(const std::string& fmt, const T& v, exprtk::details::numeric::details::complex_type);
      }
   }}

   using details::is_true;
}

#include "exprtk.hpp"
namespace exprtk
{
   namespace details
   {
      namespace numeric
      {
         namespace details
         {
            struct complex_type {};

            //needs to be cahnged
            template<typename T> struct number_type< ComplexWrapper<T>> { typedef complex_type type;};

            template <typename T>
            struct epsilon_type<ComplexWrapper<T>>
            {
               static inline ComplexWrapper<T> value()
               {
                  const ComplexWrapper<T> epsilon = ComplexWrapper<T>(0.000000000001, 0.000000000001);
                  return epsilon;
               }
            };
            
            template <typename T>
            inline bool is_nan_impl(const T& v, complex_type)
            {
               return v != v;
            }

            template <typename T> inline ComplexWrapper<T>   abs_impl(const ComplexWrapper<T> v, complex_type) { return abs(v); }
            template <typename T> inline ComplexWrapper<T>  acos_impl(const ComplexWrapper<T> v, complex_type) { return acos (v); }
            template <typename T> inline ComplexWrapper<T> acosh_impl(const ComplexWrapper<T> v, complex_type) { return log(v + sqrt((v * v) - ComplexWrapper<T>(1)));}
            template <typename T> inline ComplexWrapper<T>  asin_impl(const ComplexWrapper<T> v, complex_type) { return asin (v); }
            template <typename T> inline ComplexWrapper<T> asinh_impl(const ComplexWrapper<T> v, complex_type) { return log(v + sqrt((v * v) + ComplexWrapper<T>(1)));}
            template <typename T> inline ComplexWrapper<T>  atan_impl(const ComplexWrapper<T> v, complex_type) { return atan (v); }
            template <typename T> inline ComplexWrapper<T> atanh_impl(const ComplexWrapper<T> v, complex_type) { return atanh(v); }
            template <typename T> inline ComplexWrapper<T>  ceil_impl(const ComplexWrapper<T> v, complex_type) { return ceil (v); }
            template <typename T> inline ComplexWrapper<T>   cos_impl(const ComplexWrapper<T> v, complex_type) { return cos  (v); }
            template <typename T> inline ComplexWrapper<T>  cosh_impl(const ComplexWrapper<T> v, complex_type) { return cosh (v); }
            template <typename T> inline ComplexWrapper<T>   exp_impl(const ComplexWrapper<T> v, complex_type) { return exp  (v); }
            template <typename T> inline ComplexWrapper<T> floor_impl(const ComplexWrapper<T> v, complex_type) { return floor(v); }
            template <typename T> inline ComplexWrapper<T>   log_impl(const ComplexWrapper<T> v, complex_type) { return ln  (v); }
            template <typename T> inline ComplexWrapper<T> log10_impl(const ComplexWrapper<T> v, complex_type) { return log_10(v); }
            template <typename T> inline ComplexWrapper<T>  log2_impl(const ComplexWrapper<T> v, complex_type) { return log_2(v);}
            template <typename T> inline ComplexWrapper<T>   neg_impl(const ComplexWrapper<T> v, complex_type) { return -v;}
            template <typename T> inline ComplexWrapper<T>   pos_impl(const ComplexWrapper<T> v, complex_type) { return  v;}
            template <typename T> inline ComplexWrapper<T>   sin_impl(const ComplexWrapper<T> v, complex_type) { return sin  (v); }
            template <typename T> inline ComplexWrapper<T>  sinh_impl(const ComplexWrapper<T> v, complex_type) { return sinh (v); }
            template <typename T> inline ComplexWrapper<T>  sqrt_impl(const ComplexWrapper<T> v, complex_type) { return sqrt (v); }
            template <typename T> inline ComplexWrapper<T>   tan_impl(const ComplexWrapper<T> v, complex_type) { return tan  (v); }
            template <typename T> inline ComplexWrapper<T>  tanh_impl(const ComplexWrapper<T> v, complex_type) { return tanh (v); }
            template <typename T> inline ComplexWrapper<T>   cot_impl(const ComplexWrapper<T> v, complex_type) { return cot(v); }
            template <typename T> inline ComplexWrapper<T>   sec_impl(const ComplexWrapper<T> v, complex_type) { return sec(v); }
            template <typename T> inline ComplexWrapper<T>   csc_impl(const ComplexWrapper<T> v, complex_type) { return csc(v); }
            template <typename T> inline ComplexWrapper<T>   r2d_impl(const ComplexWrapper<T> v, complex_type) { return r2d(v); }
            template <typename T> inline ComplexWrapper<T>   d2r_impl(const ComplexWrapper<T> v, complex_type) { return d2r(v);  }
            template <typename T> inline ComplexWrapper<T>   d2g_impl(const ComplexWrapper<T> v, complex_type) { return d2g(v); }
            template <typename T> inline ComplexWrapper<T>   g2d_impl(const ComplexWrapper<T> v, complex_type) { return g2d(v); }
            template <typename T> inline ComplexWrapper<T>  notl_impl(const ComplexWrapper<T> v, complex_type) { return notl(v); }
            template <typename T> inline ComplexWrapper<T>  frac_impl(const ComplexWrapper<T> v, complex_type) { return frac(v); }
            template <typename T> inline ComplexWrapper<T> trunc_impl(const ComplexWrapper<T> v, complex_type) { return trunc(v);}
            template <typename T> inline ComplexWrapper<T> min_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type) { return abs(v1)<abs(v0)?v1:v0;}
            template <typename T> inline ComplexWrapper<T> max_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1,complex_type) {return abs(v1)>abs(v0)?v1:v0;}



            template <typename T> inline T const_pi_impl(complex_type) { return complex::details::constant::pi<T>; }
            template <typename T> inline T const_e_impl (complex_type) { return complex::details::constant::e<T>;  }

            template <typename T>
            inline int to_int32_impl(const ComplexWrapper<T> v, complex_type)
            {
               return 0; //// static_cast<int>(abs(v));
            }

            template <typename T>
            inline long long to_int64_impl(const ComplexWrapper<T> v, complex_type)
            {
               return 0; ////static_cast<long long int>(abs(v));
            }

            template <typename T>
            inline bool is_true_impl(const ComplexWrapper<T> v)
            {
               return v != ComplexWrapper<T>(0,0);
            }

            template <typename T>
            inline bool is_false_impl(const ComplexWrapper<T> v)
            {
               return v == ComplexWrapper<T>(0,0);
            }

            template <typename T>
            inline ComplexWrapper<T> expm1_impl(const ComplexWrapper<T> v, complex_type)
            {
               if (abs(v) < ComplexWrapper<T>(0.00001, 0.0))
                  return v + (ComplexWrapper<T>(0.5, 0.0) * v * v);
               else
                  return exp(v) - ComplexWrapper<T>(1, 0.0);
            }

            template <typename T>
            inline ComplexWrapper<T> nequal_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               const ComplexWrapper<T> epsilon  = epsilon_type<ComplexWrapper<T>>::value();
               const ComplexWrapper<T> eps_norm = (std::max(ComplexWrapper<T>(1, 0.0),std::max(abs_impl(v0,complex_type()),abs_impl(v1,complex_type()))) * epsilon);
               return (abs_impl(v0 - v1,complex_type()) > eps_norm) ? ComplexWrapper<T>(1, 0) : ComplexWrapper<T>(0, 0);
            }

            template <typename T>
            inline ComplexWrapper<T> sgn_impl(const ComplexWrapper<T> v, complex_type)
            {
               return modulus(v);
               //      if (v > T(0, 0)) return T(+1);
               // else if (v < T(0, 0)) return T(-1);
               // else               return T( 0);
            }

            template <typename T>
            inline ComplexWrapper<T> log1p_impl(const ComplexWrapper<T> v, complex_type)
            {
               if (v > ComplexWrapper<T>(-1, 0))
               {
                  if (abs_impl(v,complex_type()) > ComplexWrapper<T>(0.0001, 0))
                  {
                     return log_impl(ComplexWrapper<T>(1) + v,complex_type());
                  }
                  else
                     return (ComplexWrapper<T>(-0.5, 0) * v + ComplexWrapper<T>(1, 0)) * v;
               }
               else
                  return ComplexWrapper<T>(std::numeric_limits<T>::quiet_NaN());
            }

            template <typename T>
            inline ComplexWrapper<T> erf_impl(ComplexWrapper<T> v, complex_type)
            {
               const ComplexWrapper<T> t = ComplexWrapper<T>(1, 0) / (ComplexWrapper<T>(1, 0) + ComplexWrapper<T>(0.5, 0) * abs_impl(v,complex_type()));

               static const ComplexWrapper<T> c[] = {
                                      ComplexWrapper<T>( 1.26551223, 0), ComplexWrapper<T>(1.00002368, 0),
                                      ComplexWrapper<T>( 0.37409196, 0), ComplexWrapper<T>(0.09678418, 0),
                                      ComplexWrapper<T>(-0.18628806, 0), ComplexWrapper<T>(0.27886807, 0),
                                      ComplexWrapper<T>(-1.13520398, 0), ComplexWrapper<T>(1.48851587, 0),
                                      ComplexWrapper<T>(-0.82215223, 0), ComplexWrapper<T>(0.17087277, 0)
                                    };

               ComplexWrapper<T> result = ComplexWrapper<T>(1) - t * exp_impl((-v * v) -
                                      c[0] + t * (c[1] + t *
                                     (c[2] + t * (c[3] + t *
                                     (c[4] + t * (c[5] + t *
                                     (c[6] + t * (c[7] + t *
                                     (c[8] + t * (c[9]))))))))),complex_type());

               return (v >= ComplexWrapper<T>(0, 0)) ? result : -result;
            }

            template <typename T>
            inline ComplexWrapper<T> erfc_impl(ComplexWrapper<T> v, complex_type)
            {
               return ComplexWrapper<T>(1, 0) - erf_impl(v,complex_type());
            }

            template <typename T>
            inline ComplexWrapper<T> ncdf_impl(ComplexWrapper<T> v, complex_type)
            {
               ComplexWrapper<T> cnd = ComplexWrapper<T>(0.5, 0) * (ComplexWrapper<T>(1, 0) + erf_impl(
                                           abs_impl(v,complex_type()) /
                                           ComplexWrapper<T>(complex::details::constant::sqrt2<T>),complex_type()));
               return  (v < ComplexWrapper<T>(0, 0)) ? (ComplexWrapper<T>(1, 0) - cnd) : cnd;
            }

            template <typename T>
            inline ComplexWrapper<T> modulus_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               return v0/v1;
            }

            template <typename T>
            inline ComplexWrapper<T> pow_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               return pow(v0, v1); //! Beware
            }

            template <typename T>
            inline ComplexWrapper<T> logn_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               return ln(v0) / ln(v1);
            }

            template <typename T>
            inline ComplexWrapper<T> sinc_impl(ComplexWrapper<T> v, complex_type)
            {
               if (abs_impl(v,complex_type()) >= std::numeric_limits<ComplexWrapper<T>>::epsilon())
                   return(sin_impl(v,complex_type()) / v);
               else
                  return ComplexWrapper<T>(1);
            }

            template <typename T>
            inline ComplexWrapper<T> xor_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               return (is_false_impl(v0) != is_false_impl(v1)) ? ComplexWrapper<T>(1) : ComplexWrapper<T>(0);
            }

            template <typename T>
            inline ComplexWrapper<T> xnor_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               const bool v0_true = is_true_impl(v0);
               const bool v1_true = is_true_impl(v1);
               if ((v0_true &&  v1_true) || (!v0_true && !v1_true))
                  return ComplexWrapper<T>(1, 0);
               else
                  return ComplexWrapper<T>(0, 0);
            }

            template <typename T>
            inline ComplexWrapper<T> equal_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               const ComplexWrapper<T> epsilon = epsilon_type<ComplexWrapper<T>>::value();
               const ComplexWrapper<T> eps_norm = (max(ComplexWrapper<T>(1, 0),max(abs_impl(v0,complex_type()),abs_impl(v1,complex_type()))) * epsilon);
               return (abs_impl(v0 - v1,complex_type()) <= eps_norm) ? ComplexWrapper<T>(1, 0) : ComplexWrapper<T>(0, 0);
            }

            template <typename T>
            inline ComplexWrapper<T> round_impl(const ComplexWrapper<T> v, complex_type)
            {
               return ((v < ComplexWrapper<T>(0, 0)) ? ceil(v - ComplexWrapper<T>(0.5, 0.5)) : floor(v + ComplexWrapper<T>(0.5, 0.5)));
            }

            template <typename T>
            inline ComplexWrapper<T> roundn_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               const int index = std::max<int>(0, std::min<int>(pow10_size - 1, (int)abs(floor(v1))));
               const ComplexWrapper<T> p10 = ComplexWrapper<T>(pow10[index]);
               if (v0 < ComplexWrapper<T>(0, 0))
                  return ComplexWrapper<T>(ceil ((v0 * p10) - ComplexWrapper<T>(0.5, 0.5)) / p10);
               else
                  return ComplexWrapper<T>(floor((v0 * p10) + ComplexWrapper<T>(0.5, 0.5)) / p10);
            }

            template <typename T>
            inline bool is_integer_impl(const ComplexWrapper<T> v, complex_type)
            {
               return false;//(T(0) == modulus_impl(v.real,T(1),complex_type()));
            }

            template <typename T>
            inline ComplexWrapper<T> root_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               return pow(v0, ComplexWrapper<T>(1, 0) / v1);
            }

            template <typename T>
            inline ComplexWrapper<T> hypot_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               return sqrt((v0 * v0) + (v1 * v1));
            }

            template <typename T>
            inline ComplexWrapper<T> atan2_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               return atan(v0+v1); //! Makes no sense always in complex
            }

            template <typename T>
            inline ComplexWrapper<T> shr_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               return v0 * (ComplexWrapper<T>(1, 0) / pow(ComplexWrapper<T>(2, 0),v1));
            }

            template <typename T>
            inline ComplexWrapper<T> shl_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               return v0 * pow(ComplexWrapper<T>(2, 0),v1);
            }

            template <typename T>
            inline ComplexWrapper<T> and_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               return (is_true_impl(v0) && is_true_impl(v1)) ? ComplexWrapper<T>(1, 0) : ComplexWrapper<T>(0, 0);
            }

            template <typename T>
            inline ComplexWrapper<T> nand_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               return (is_false_impl(v0) || is_false_impl(v1)) ? ComplexWrapper<T>(1, 0) : ComplexWrapper<T>(0, 0);
            }

            template <typename T>
            inline ComplexWrapper<T> or_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               return (is_true_impl(v0) || is_true_impl(v1)) ? ComplexWrapper<T>(1, 0) : ComplexWrapper<T>(0, 0);
            }

            template <typename T>
            inline ComplexWrapper<T> nor_impl(const ComplexWrapper<T> v0, const ComplexWrapper<T> v1, complex_type)
            {
               return (is_false_impl(v0) && is_false_impl(v1)) ? ComplexWrapper<T>(1, 0) : ComplexWrapper<T>(0, 0);
            }
         }
      }

      template <typename Iterator, typename T>
      inline bool string_to_real(Iterator& itr_external, const Iterator end, ComplexWrapper<T>& t, details::numeric::details::complex_type)
      {
         std::string f(itr_external, end);
         std::istringstream is('(' + f + ')'); //convert string to complex number
         is >> t;
         return true;
      }

      template <typename T>
      inline bool is_true(const ComplexWrapper<T> v)
      {
         return v != ComplexWrapper<T>(0,0);
      }

      template <typename T>
      inline bool is_false(const ComplexWrapper<T> v)
      {
         return v == ComplexWrapper<T>(0,0);
      }
   }

   namespace rtl { namespace io
   {
      namespace details
      {
         template <typename T>
         inline void print_type(const std::string& fmt, const ComplexWrapper<T> v, exprtk::details::numeric::details::complex_type)
         {
            printf(fmt.c_str(),v.asString());
         }
      }
   }}
}

// #endif /*EXPRTK_COMPLEX*/
#endif