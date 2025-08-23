#ifndef _SCREENED_COULOMB_H_
#define _SCREENED_COULOMB_H_

#include <cmath>
#include <limits>
#include <cassert>
#include <array>

/// Bohr radius [nm]
#define BOHR_RADIUS 0.05291772108

/// Lindhard screening length constant
#define SCREENCONST 0.88534

/// \f$ e^2 / 4 \pi \epsilon_0 = e^2 c^2 \f$ in [eV-nm] */
#define E2C2 1.43996445

/** \file */

/**
 * \brief Enumeration of the supported types of screening
 */
enum class Screening {
    None = 0, /**< Unscreened Coulomb potential */
    Bohr = 1, /**< Bohr */
    KrC = 2, /**< Kr-C */
    Moliere = 3, /**< Moliere */
    ZBL = 4, /**< Ziegler-Biersack-Littmark (ZBL) Universal */
    Invalid = -1 /**< Invalid screening value */
};

/**
 * @brief Definition of the screening function
 *
 * A templated structure that includes the definition of the screening
 * function \f$\Phi(x)\f$ (screening_function::operator()) and the screening length \f$ a \f$
 * (screening_function::screeningLength).
 *
 * Additionally, it returns the name and type of screening as a string or enum, respectively.
 *
 * @tparam ScreeningType (screening_function_t) the type of screening function
 */
template <Screening ScreeningType>
struct screening_function
{
    /**
     * @brief Returns the screening length in [nm]
     * @param Z1 is the atomic number of the projectile
     * @param Z2 is the atomic number of the target
     * @return the screening length [nm]
     */
    static double screeningLength(int Z1, int Z2);

    /**
     * @brief operator () implements the actual screening function
     * @param x is the radial distance in units of the screening length
     * @return the value of \f$ \Phi(x) \f$
     */
    double operator()(const double &x) const;

    /// The name of the screening function
    static const char *screeningName();

    /// The type of screening as a \ref Screening enum type
    static Screening screeningType() { return ScreeningType; }
};

// Explicit specialization for the unscreened Coulomb interaction
template <>
struct screening_function<Screening::None>
{
    static double screeningLength(int Z1, int Z2) { return 1.; } // fake
    double operator()(const double &x) const { return 1.; }
    static const char *screeningName() { return "Unscreened Coulomb"; }
    static Screening screeningType() { return Screening::None; }
};
// Explicit specialization for the Lenz-Jensen potential
template <>
struct screening_function<Screening::Bohr>
{
    static double screeningLength(int Z1, int Z2)
    {
        return BOHR_RADIUS / std::sqrt(std::pow(Z1, 2. / 3) + std::pow(Z2, 2. / 3));
    }
    /* Screening function coefficients */
    constexpr static const int N = 1;
    constexpr static const double C[] = { 1.0 };
    constexpr static const double A[] = { 1.0 };
    double operator()(const double &x) const { return exp(-x); }
    static const char *screeningName() { return "Bohr"; }
    static Screening screeningType() { return Screening::Bohr; }
};
// Explicit specialization for the Kr-C potential
template <>
struct screening_function<Screening::KrC>
{
    static Screening screeningType() { return Screening::KrC; }
    static const char *screeningName() { return "Kr-C"; }
    static double screeningLength(int Z1, int Z2)
    {
        return SCREENCONST * BOHR_RADIUS / (std::pow(Z1, 0.23) + std::pow(Z2, 0.23));
    }

    /* Screening function coefficients */
    constexpr static const int N = 3;
    constexpr static const double C[] = { 0.190945, 0.473674, 0.335381 };
    constexpr static const double A[] = { 0.278544, 0.637174, 1.919249 };

    double operator()(const double &x) const
    {
        return C[0] * exp(-A[0] * x) + C[1] * exp(-A[1] * x) + C[2] * exp(-A[2] * x);
    }
};
// Explicit specialization for the Moliere potential
template <>
struct screening_function<Screening::Moliere>
{
    static Screening screeningType() { return Screening::Moliere; }
    static const char *screeningName() { return "Moliere"; }
    static double screeningLength(int Z1, int Z2)
    {
        return SCREENCONST * BOHR_RADIUS / (std::pow(Z1, 0.23) + std::pow(Z2, 0.23));
    }

    /* Screening function coefficients */
    constexpr static const int N = 3;
    constexpr static const double C[] = { 0.35, 0.55, 0.10 };
    constexpr static const double A[] = { 0.30, 1.20, 6.00 };

    double operator()(const double &x) const
    {
        return C[0] * exp(-A[0] * x) + C[1] * exp(-A[1] * x) + C[2] * exp(-A[2] * x);
    }
};
// Explicit specialization for the Ziegler-Biersack-Littmark (ZBL) potential
template <>
struct screening_function<Screening::ZBL>
{
    static Screening screeningType() { return Screening::ZBL; }
    static const char *screeningName() { return "Ziegler-Biersack-Littmark (ZBL)"; }
    static double screeningLength(int Z1, int Z2)
    {
        return SCREENCONST * BOHR_RADIUS / (std::pow(Z1, 0.23) + std::pow(Z2, 0.23));
    }

    /* Universal screening function coefficients TRIM85 */
    constexpr static const int N = 4;
    constexpr static const double C[] = { 0.18175, 0.50986, 0.28022, 0.028171 };
    constexpr static const double A[] = { 3.19980, 0.94229, 0.40290, 0.201620 };

    double operator()(const double &x) const
    {
        return C[0] * exp(-A[0] * x) + C[1] * exp(-A[1] * x) + C[2] * exp(-A[2] * x)
                + C[3] * exp(-A[3] * x);
    }
};

/**
 * @brief The xs_cms_base class implements basic components for screened Coulomb scattering
 * calculations
 *
 * The calculations refer to the center-of-mass system.
 *
 * The quantities obtained by this class are:
 * - the function F() in the integrant of the scattering integral
 * - the apsis of the trajectory, i.e., the distance of closest approach, apsis()
 * - the scattering angle in the impulse approximation (high enegry and large impact parameter),
 * theta_impulse_approx()
 *
 * @tparam
 *   ScreeningType specifies the type of screening function
 */
template <Screening ScreeningType>
struct xs_cms_base : public screening_function<ScreeningType>
{

public:
    typedef screening_function<ScreeningType> Phi;

    /**
     * @brief The function \f$ F(x) \f$ under the scattering integral
     *
     * \f[
     * F(x) = 1 - \frac{\Phi(x)}{x\,e} - \frac{s^2}{x^2}
     * \f]
     *
     * @param x integral variable
     * @param e reduced energy
     * @param s reduced impact parameter
     * @return the value of the function
     */
    static double F(double x, double e, double s)
    {
        double sx = s / x;
        Phi p;
        return 1. - p(x) / (x * e) - sx * sx;
    }

    /**
     * @brief Finds the apsis of the scattered projectile trajectory
     *
     * The apsis is the distance from the scattering center at the point of closest approach.
     *
     * It is the root of
     * \f$ F(x_0)=0 \f$. In the general case the equation is solved numerically
     * using bisection.
     *
     * For the unscreened Coulomb potential, \f$ x_0 \f$ is given by
     * \f[
     * x_0 = \frac{1}{2\epsilon} \left[ 1 + \sqrt{1 + (2\,\epsilon\,s)^2} \right]
     * \f]
     *
     * @param e reduced energy of the scattered particle
     * @param s reduced impact parameter
     * @return the minimal approach distance (in screening length units)
     */
    static double apsis(double e, double s);

    /**
     * @brief Scattering angle in the Impulse approximation
     *
     * From Lehmann & Leibfried, Z. Physic 172 (1963) 465, the 1st order term
     * in the "momentum approx." for the screened potential \f$ V(r) = (A/r)e^{-r/a} \f$ is
     *
     * \f[
     * \theta_1 = \epsilon^{-1}\,K_1(s)
     * \f]
     *
     * where \f$ K_1(x) \f$ is the modified Hankel function (or modified
     * Bessel function of the 2nd kind).
     *
     * The above formula can be generalized for a screening function
     * expressed as a sum of exponentials, i.e., for the ZBL, Kr-C & Moliere potentials.
     * Namely, for a screening function \f$ \Phi(x) = \sum_i{C_i\, e^{-A_i x}} \f$ the scattering
     * angle in the impulse approx. is
     *
     * \f[
     * \theta_1 = \epsilon^{-1} \sum_i {C_i A_i K_1(A_i s) }
     * \f]
     *
     * For the unscreened Coulomb potential the function returns the
     * exact scattering angle
     *
     * \f[
     * \theta = 2 \arcsin \left[ 1 + (2\epsilon\, s)^2\right]^{-1/2}
     * \f]
     *
     * @param e reduced energy of the scattered particle
     * @param s reduced impact factor
     * @return the scattering angle [rad]
     */
    static double theta_impulse_approx(double e, double s);
};

/**
 * @brief Enumeration of the different quadratures for numerical evaluation of the scattering
 * integral
 *
 */
enum class Quadrature {
    GaussChebyshev, /**< Gauss-Chebushev scheme */
    Lobatto4, /**< 4-th order Lobatto scheme */
    Magic, /**<  MAGIC interpolation formula of Biersack & Haggmark */
    None /**<  No quad needed */
};

namespace detail {

/**
 * @brief Gauss-Chebyshev quadrature in  \f$ [0, 1] \f$
 *
 * A modified Gauss-Chebyshev quadrature scheme to approximate the integral
 * \f[
 * I = \int_0^1 {f(x)\, dx} = \int_0^1 {[f(x)\,w(x)]\, [w(x)]^{-1}\, dx},
 * \f]
 * where \f$ w(x)=\sqrt{1-x^2} \f$.
 *
 * The approximation is given by
 * \f[
 * I \approx \sum_{i=0}^{N-1}{w_i f(x_i)},
 * \f]
 * where
 * \f[
 * x_i = \cos\left[ \frac{\pi}{2N}\,\left( i + \frac{1}{2} \right) \right]
 * \f]
 * \f[
 * w_i = \frac{\pi}{2N} \cdot w(x_i)
 * \f]
 *
 * @tparam Functor the function \f$ f(x) \f$
 * @tparam N the number of terms in the sum
 *
 * \sa https://dlmf.nist.gov/3.5.v
 */
template <class Functor, int N>
class gc_pos_quad
{
public:
    /// @brief Construct a new gc_pos_quad object
    /// @param func A function object
    gc_pos_quad(Functor &func) : func_(func) { }

    /// @brief Return the value of the integral
    double operator()() const
    {
        double s{ 0.0 };
        for (int i = 0; i < N; ++i)
            s += w_[i] * func_(u_[i]);
        return s;
    }

private:
    Functor &func_;
    constexpr static auto u_{ []() constexpr {
        const double d = 0.5 * M_PI / N;
        std::array<double, N> u{};
        for (int i = 0; i < N; ++i) {
            u[i] = std::cos(d * (i + 0.5));
        }
        return u;
    }() };
    constexpr static auto w_{ []() constexpr {
        const double d = 0.5 * M_PI / N;
        std::array<double, N> w{};
        for (int i = 0; i < N; ++i) {
            double u = std::cos(d * (i + 0.5));
            w[i] = d * std::sqrt(1 - u * u);
        }
        return w;
    }() };
};

// The preferred quadrature order
template <Quadrature QuadType>
constexpr int __preferredQuadOrder__()
{
    switch (QuadType) {
    case Quadrature::GaussChebyshev:
        return 16;
    case Quadrature::Lobatto4:
        return 4;
    case Quadrature::Magic:
        return 0;
    case Quadrature::None:
        return 0;
    }
};

/**
 * @brief Numerical evaluation of the scattering angle
 *
 * The center-of-mass scattering angle is \f$ \theta = \pi - 2 s I\f$,
 * where
 * \f[
 * I = \int_{x_0}^{\infty}{x^{-2} \left[ F(x) \right]^{-1/2} dx}
 * \f]
 *
 * Setting \f$ x = x_0/u \f$ the integral becomes
 * \f[
 * I = x_0^{-1} \int_{0}^{1}{\left[ F(x_0/u) \right]^{-1/2} du}
 * \f]
 * which can be evaluated by quadrature. Different schemes are employed depending on the template
 * argument.
 *
 * - Gauss–Chebyshev quadrature using \ref gc_pos_quad . The order is set by the template parameter
 * \a N.
 *
 * - 4th-order Lobatto quadrature from
 * Mendenhall & Weller 2005 (https://dx.doi.org/10.1016/j.nimb.2004.08.014)
 *
 * \f[
 *   I = \frac{\pi}{2\, x_0} \left[
 *   \frac{1 + \lambda_0}{30} +
 *   \sum_{j=1}^{4}{w_j\, \left[ F(x_0/u_j) \right]^{-1/2} }
 *   \right]
 * \f]
 *
 * where
 *
 * \f[
 * \lambda_0 = \left[
 * \frac{1}{2} - \frac{\Phi'(x_0)}{2\,\epsilon} + \frac{s^2}{2x_0^2}
 * \right]^{-1/2}
 * \f]
 *
 * and
 *
 * \f[
 * w = [0.03472124, 0.1476903, 0.23485003, 0.1860249],
 * \f]
 * \f[
 * u = [0.9830235, 0.8465224, 0.5323531, 0.18347974].
 * \f]
 *
 * This method is slightly less accurate but more efficient. N=4.
 *
 * - the MAGIC analytic formula of Biersack & Haggmark (NIM1980). N=0.
 *
 * @tparam ScreeningType The type of screening
 * @tparam QuadType The type of quadrature to use
 * @tparam N the quadrature order
 */
template <Screening ScreeningType, Quadrature QuadType, int N>
struct theta_integrator
{
    /**
     * @brief Returns the scattering angle
     * @param e reduced energy
     * @param s reduced impact parameter
     * @param x0 apsis of the trajectory
     * @return
     */
    static double theta(double e, double s, double x0);

    /// @brief The name of the quadrature scheme
    static const char *quadratureName();

    /// @brief The type of quadrature as a \ref Quadrature enum type
    static Quadrature quadratureType() { return QuadType; }

    /// @brief The quadrature order
    static int quadratureOrder() { return N; };
};

// partial specializations of theta_integrator

template <Screening ScreeningType, int N>
struct theta_integrator<ScreeningType, Quadrature::GaussChebyshev, N>
{
    static double theta(double e, double s, double x0)
    {
        H h(e, s, x0);
        gc_pos_quad<H, N> gc(h);
        return M_PI - 2.0 * s * gc() / x0;
    }
    static const char *quadratureName() { return "Gauss-Chebyshev"; }
    static Quadrature quadratureType() { return Quadrature::GaussChebyshev; }
    static int quadratureOrder() { return N; };

private:
    typedef xs_cms_base<ScreeningType> XS_;

    struct H
    {
        double e, s, x0;
        void set(double ae, double as, double ax0) { e = ae, s = as, x0 = ax0; }
        H(double ae, double as, double ax0) : e(ae), s(as), x0(ax0) { }
        double operator()(double u) { return std::sqrt(1.0 / XS_::F(x0 / u, e, s)); }
    };
};

template <Screening ScreeningType, int N>
struct theta_integrator<ScreeningType, Quadrature::Lobatto4, N>
{
    static_assert(N == 4, "Quadrature::Lobatto4 requires N=4.");

    static double theta(double e, double s, double x0)
    {
        // Lobatto quadrature coefficients
        constexpr static const double w[] = { 0.03472124, 0.1476903, 0.23485003, 0.1860249 };
        constexpr static const double q[] = { 0.9830235, 0.8465224, 0.5323531, 0.18347974 };

        // a from MW2005
        double a = (1.0 + lambda0(e, s, x0)) / 30;
        a += w[0] / std::sqrt(XS_::F(x0 / q[0], e, s));
        a += w[1] / std::sqrt(XS_::F(x0 / q[1], e, s));
        a += w[2] / std::sqrt(XS_::F(x0 / q[2], e, s));
        a += w[3] / std::sqrt(XS_::F(x0 / q[3], e, s));

        return M_PI - s * M_PI / x0 * a;
    }
    static const char *quadratureName() { return "Lobatto"; }
    static Quadrature quadratureType() { return Quadrature::Lobatto4; }
    static int quadratureOrder() { return 4; };

private:
    typedef xs_cms_base<ScreeningType> XS_;

    static double phi_dot(double x)
    {
        const auto &C = XS_::Phi::C;
        const auto &A = XS_::Phi::A;
        double s(0.);
        for (int i = 0; i < XS_::Phi::N; ++i)
            s -= C[i] * A[i] * std::exp(-A[i] * x);
        return s;
    }

    static double lambda0(double e, double s, double x0)
    {
        double s2 = s / x0;
        double l0 = 0.5 * (1. + s2 * s2 - phi_dot(x0) / e);
        return 1.0 / std::sqrt(l0);
    }
};

template <Screening ScreeningType, int N>
struct theta_integrator<ScreeningType, Quadrature::Magic, N>
{
    static_assert(ScreeningType == Screening::ZBL,
                  "Quadrature::Magic analytic formula can be used only for ZBL screening.");

    static_assert(N == 0, "Quadrature::Magic requires N=0.");

    static double theta(double e, double s, double x0) { return 2 * std::acos(cosThetaBy2(e, s)); }
    static const char *quadratureName() { return "Magic"; }
    static Quadrature quadratureType() { return Quadrature::Magic; }
    static int quadratureOrder() { return 0; };

private:
    // Return \f$ \cos(\theta/2) \f$, where \f$ \theta \f$ is the center-of-mass scattering
    // angle
    static double cosThetaBy2(double e, double s);
    // Return the ZBL potential and optionally its derivative
    static double zbl_and_deriv(double R, double *Vprime);
};

template <Screening ScreeningType, int N>
struct theta_integrator<ScreeningType, Quadrature::None, N>
{
    static_assert(ScreeningType == Screening::None,
                  "Quadrature::None can be used only for un-screened potential.");

    static_assert(N == 0, "Quadrature::None requires N=0.");

    static double theta(double e, double s, double x0) { return 0; }
    static const char *quadratureName() { return "None"; }
    static Quadrature quadratureType() { return Quadrature::None; }
    static int quadratureOrder() { return 0; };
};

/**
 * @brief Modified Bessel of the 2nd kind K1(x)
 *
 * For \f$ x>200 \f$ an approximation is used (error below 1%)
 *
 * @param x Real number
 * @return K1(x)
 */
static double bessel_k1(double x)
{
    // sqrt(pi/2)
    constexpr static const double sqrtpihalf = 1.25331413731550012081;
    return (x > 200) ? sqrtpihalf * exp(-x) / std::sqrt(x) : std::cyl_bessel_k(1, x);
}

} // namespace detail

/**
 * @brief The xs_cms class implements screened Coulomb scattering calculations
 * in the center-of-mass system.
 *
 * The class provides functions for calculating all aspects of the scattering
 * process as a function of reduced energy and impact parameter:
 * scattering angle, cross-section, stopping cross-section.
 *
 * The scattering integral is generally evaluated by numerical quadrature according to
 * the method specified by the template parameter Quad.
 *
 * For Screening::ZBL_MAGIC, \f$\theta\f$ is evaluated using the MAGIC approximation
 * of Biersack & Haggmark (NIM1980).
 *
 * For the unscreened Coulomb potential (Screening::None),
 * the exact analytical formulae are used.
 *
 * find_s() inverts numerically the scattering integral to obtain
 * the reduced impact paramter as a function of energy and scattering
 * angle, \f$s = s(\epsilon,\theta)\f$.
 *
 * @tparam
 *   ScreeningType specifies the type of screening function
 * @tparam
 *   QuadType specifies the quadrature method
 * @tparam
 *   N specifies the quadrature order
 *
 * @sa \ref detail::theta_integrator
 */
template <Screening ScreeningType, Quadrature QuadType = Quadrature::GaussChebyshev,
          int N = detail::__preferredQuadOrder__<QuadType>()>
struct xs_cms : public xs_cms_base<ScreeningType>,
                private detail::theta_integrator<ScreeningType, QuadType, N>
{
private:
    typedef xs_cms<ScreeningType, QuadType> My_t_;
    typedef xs_cms_base<ScreeningType> Base_t_;
    typedef detail::theta_integrator<ScreeningType, QuadType, N> QuadInt_t_;

public:
    using QuadInt_t_::quadratureName;
    using QuadInt_t_::quadratureOrder;
    using QuadInt_t_::quadratureType;
    using typename Base_t_::Phi;

    /**
     * @brief Scattering angle in center-of-mass (CM) system
     *
     * Calculated using Gauss-Chebyshev quadrature.
     *
     * For \f$ s\cdot \epsilon^{1/6} > 100 \f$ where the scattering angle is very small
     * (\f$ \theta < 10^{-8} \f$) this function
     * returns the impulse approximation theta_impulse_approx().
     *
     * The region of \f$ (\epsilon,s) \f$ where impulse approximation can
     * be applied has been found empirically.
     *
     * @param e reduced energy of the scattered particle in CM system
     * @param s reduced impact parameter
     * @return
     */
    static double theta(double e, double s)
    {
        // if s > 100/e^(1/6) use impulse approx
        // The limit was found empirically for theta < 1e-9
        double s3 = s * s * s;
        if (e * s3 * s3 > 1.e12)
            return Base_t_::theta_impulse_approx(e, s);

        double x0 = Base_t_::apsis(e, s);
        return QuadInt_t_::theta(e, s, x0);
    }

    /**
     * @brief Return \f$ \sin^2(\theta(\epsilon, s)/2) \f$
     * @param e the reduced energy
     * @param s the reduced impact parameter
     * @return the value of \f$ \sin^2(\theta/2) \f$
     */
    static double sin2Thetaby2(double e, double s)
    {
        double v = std::sin(0.5 * theta(e, s));
        return v * v;
    }

    /**
     * @brief Return the reduced impact parameter \f$ s=s(\epsilon,\theta) \f$
     *
     * Use bisection to invert the function \f$ \theta=\theta(\epsilon,s) \f$
     * and obtain the reduced impact parameter s
     * for given energy and scattering angle
     *
     * @param e reduced energy
     * @param thetaCM scattering angle (rad) in center-of-mass system
     * @param tol abs. tolerance for the value of s. Defaults to the machine epsilon.
     * @return the reduced impact parameter
     */
    static double find_s(double e, double thetaCM,
                         double tol = std::numeric_limits<double>::epsilon());

    /**
     * @brief Differential cross-section in center-of-mass system
     *
     * \f[
     *   \frac{d\sigma}{d\Omega} = \frac{s}{\sin(\theta)}\left| \frac{ds}{d\theta}\right|
     * \f]
     *
     * In units of \f$ a^2 \f$
     *
     * @param e is the reduced energy
     * @param thetaCM is the center-of-mass scattering angle (rad)
     * @param tol abs. tolerance for finding \f$ s=s(\epsilon,\theta) \f$.
     * @return the reduced cross-section
     *
     * \sa find_s()
     */
    static double crossSection(double e, double thetaCM,
                               double tol = std::numeric_limits<double>::epsilon());

    /**
     * @brief Reduced stopping cross-section \f$ s_n(\epsilon) \f$
     *
     * Calculate the reduced stopping cross-section
     *
     * \f[
     * s_n(\epsilon) = \frac{\epsilon}{\pi\, a^2 \gamma E} S_n(E) =
     * \frac{\epsilon}{\pi a^2} \, \frac{1}{4\pi} \cdot
     * \int_0^{4\pi}
     * {\Omega\,\frac{d\sigma}{d\Omega} d\Omega}
     * \f]
     *
     * Optionally, the function calculates the stopping power
     * for scattering angles up to a maximum value.
     * In this case the
     * upper limit in the integral is replaced by \f$\theta_{max}\f$.
     *
     * The integral is evaluated by Gauss-Chebyshev quadrature. The order is specified by the
     * template parameter N.
     *
     * If \f$ x = \Omega / \Omega_{max}\f$ and
     * \f[
     * f(x) = x\,\sqrt{1-x^2}\, \frac{d\sigma}{d\Omega}
     * \f]
     * then
     * \f[
     * s_n(\epsilon) = 4 \epsilon \left( \frac{\Omega_{max}}{4\pi} \right)^2 \int_0^1 {\frac{f(x)
     * dx}{\sqrt{1-x^2}}}
     * \approx 4 \epsilon \left( \frac{\Omega_{max}}{4\pi} \right)^2 w_N \sum_{i=0}^{N-1}{f(x_i)}
     * \f]
     *
     *
     * @param e the reduced energy
     * @param theta_max optional maximum scattering angle, defaults to \f$ \pi \f$
     * @return the energy loss cross-section
     */
    static double sn(double e, double theta_max = M_PI)
    {
        double mu_max = sin(theta_max / 2);
        mu_max = mu_max * mu_max;
        eloss_functor_ F(e, mu_max);
        detail::gc_pos_quad<eloss_functor_, N> gc(F);
        return 8 * e * mu_max * mu_max * gc();
    }

private:
    // integrant for numerical integration to obtain stopping cross-section
    struct eloss_functor_
    {
        double e_, mumax_;
        eloss_functor_(double e, double mumax) : e_(e), mumax_(mumax) { }
        double operator()(double u)
        {
            double u2 = u * u;
            double mu = mumax_ * u2;
            double thetaCM = 2. * std::asin(std::sqrt(mu));
            return crossSection(e_, thetaCM, 1e-10) * u * u2;
        }
    };

    // This function returns π-θ
    static double pi_minus_theta(double e, double s)
    {
        double x0 = Base_t_::apsis(e, s);
        return M_PI - QuadInt_t_::theta(e, s, x0);
    }
};

/**
 * @brief The xs_lab class implements screened Coulomb scattering calculations
 * in the lab system.
 *
 * The class is instantiated for a particular
 * projectile/target combination.
 *
 * It can be used to perform calculations of various
 * scattering parameters.
 *
 * The scatter() function gives the recoil energy and lab scattering angle sine & cosine
 * for the scattering of a projectile with given energy and impact parameter.
 *
 * crossSection() gives the differential cross-section as a function of projectile and recoil
 * energy.
 *
 * Sn() calculates the stopping cross-section as a function of projectile energy.
 *
 * Finally, find_p() returns the impact parameter for given
 * projectile and recoil energy.
 *
 *
 * @tparam ScreeningType The type of screening function
 * @tparam QuadType The quadrature method employed for the scattering integrals
 */
template <Screening ScreeningType, Quadrature QuadType = Quadrature::GaussChebyshev,
          int N = detail::__preferredQuadOrder__<QuadType>()>
class xs_lab : public xs_cms<ScreeningType, QuadType, N>
{
    using xs_cms<ScreeningType, QuadType, N>::screeningLength;
    using xs_cms<ScreeningType, QuadType, N>::sin2Thetaby2;
    using xs_cms<ScreeningType, QuadType, N>::find_s;
    using xs_cms<ScreeningType, QuadType, N>::crossSection;
    using xs_cms<ScreeningType, QuadType, N>::sn;

protected:
    double screening_length_; /* screening length [nm] */
    double mass_ratio_; /* M1/M2 */
    double sqrt_mass_ratio_; /* we will need this occasionally */
    double gamma_; /* 4 M1 M2 / (M1 + M2)^2 */
    double red_E_conv_; /* reduced energy conversion factor */
    double sig0_; /* pi a^2 [nm^2] */

public:
    /**
     * @brief Construct an xs_lab object for a specific projectile-target combination
     * @param Z1 the projectile atomic number
     * @param M1 the projectile mass
     * @param Z2 the target atomic number
     * @param M2 the target mass
     */
    xs_lab(int Z1, double M1, int Z2, double M2)
        : screening_length_(screeningLength(Z1, Z2)),
          mass_ratio_(M1 / M2),
          sqrt_mass_ratio_(std::sqrt(mass_ratio_)),
          gamma_(4 * mass_ratio_ / ((mass_ratio_ + 1) * (mass_ratio_ + 1))),
          red_E_conv_(screening_length_ / ((mass_ratio_ + 1) * Z1 * Z2 * E2C2)),
          sig0_(M_PI * screening_length_ * screening_length_)
    {
        // f = a/(A1/A2+1)/Z1/Z2/e^2
        // e = Ecm*a/Z1/Z2/e^2 = f*E
        // Ecm = E/(A1/A2+1)

        // π*a^2/4 * γΕ/(Ε^2 a^2 )* (A1/A2+1)^2*Z1^2*Z2^2*e^4
        // π (1/E) A1/A2 Z1^2*Z2^2*e^4 /T
    }

    /// Returns the screening length \f$ a(Z_1,Z_2) \f$ in nm
    double screening_length() const { return screening_length_; }
    /// Returns the projectile to target mass ratio \f$ A = M_1 / M_2 \f$
    double mass_ratio() const { return mass_ratio_; }
    /// Returns the precomputed square root of the mass ratio
    double sqrt_mass_ratio() const { return sqrt_mass_ratio_; }
    /// Returns \f$ \gamma = 4 A / (1 + A)^2 \f$
    double gamma() const { return gamma_; }
    /// Return the reduced energy conversion factor
    /// \f[ f = \frac{a} {(1+A)Z_1Z_2e^2}  \f]
    /// such that \f$ \epsilon = f\times E  \f$
    double red_E_conv() const { return red_E_conv_; }

    /**
     * @brief Calculate scattering angle and target recoil energy.
     *
     * Given the initial energy E and impact parameter S of an incoming
     * projectile, this function calculates the target recoil energy and the projectile
     * scattering angle.
     *
     * All quantities refer to the lab system.
     *
     * @param E is the initial projectile energy [eV]
     * @param P is the impact factor [nm]
     * @param recoil_erg is the target recoil energy [eV]
     * @param theta the scattering angle in the lab system [rad]
     */
    void scatter(double E, double P, double &recoil_erg, double &theta) const;

    /**
     * @brief Returns the impact parameter for given initial projectile energy and target recoil
     * energy
     * @param E the projectile initial energy [eV]
     * @param T the target recoil energy [eV]
     * @return the corresponding impact factor [nm]
     */
    double find_p(double E, double T) const;

    /**
     * @brief Differential cross-section \f$ d\sigma(E,T)/dT \f$
     *
     * \f[
     *   \frac{d\sigma}{dT}  = \frac{4\pi a^2}{\gamma E} \frac{d\sigma}{d\Omega_{CM}}
     * \f]
     *
     * @param E projectile energy [eV]
     * @param T recoil energy [eV]
     * @return the cross-section [nm^2/eV]
     */
    double crossSection(double E, double T) const;

    /**
     * @brief Stopping cross-section \f$ S_n(E) \f$
     *
     * This function returns the stopping cross-section
     * \f[
     *   S_n(E) = \int_0^{\gamma E}{T\, d\sigma(E,T)}
     * \f]
     *
     * @param E projectile energy [eV]
     * @return the stopping power [eV-nm^2]
     */
    double Sn(double E) const;

    /**
     * @brief Stopping cross-section \f$ S_n(E,T) \f$
     *
     * This function returns the stopping cross-section for
     * energy tranfers up to \f$ T \f$:
     * \f[
     *   S_n(E,T) = \int_0^{T}{T'\, d\sigma(E,T')}
     * \f]
     *
     * @param E projectile energy [eV]
     * @param T maximum recoil energy [eV]
     * @return the stopping power [eV-nm^2]
     */
    double Sn(double E, double T) const;
};

// xs_cms_base implementation

template <Screening ScreeningType>
double xs_cms_base<ScreeningType>::apsis(double e, double s)
{
    double x2 = 1.0 / (2.0 * e);
    x2 = x2 + sqrt(x2 * x2 + s * s); // inital guesses: Mendenhall & Weller NIMB 58(1991)11, eq. 15

    double x1 = x2 / 10.;
    double f1 = F(x1, e, s); // should be always negative
    double f2 = F(x2, e, s); // should be always positive (starting ~1.0)

    // ensure bracketing
    while (f1 >= 0.) {
        // initial guess for x1 too optimistic
        x1 *= 0.1;
        f1 = F(x1, e, s); // should be always negative
    }
    while (f2 <= 0.) {
        // just in case
        x2 *= 1.001;
        f2 = F(x2, e, s);
    }

    // assert(f1<0.0 && f2>0.0); // values should be on each side of 0

    double xm = 0.5 * (x1 + x2);
    double fm = F(xm, e, s);
    double d = 1.;
    while (std::abs(fm) > std::numeric_limits<double>::epsilon()
           && std::abs(d) > std::numeric_limits<double>::epsilon()) {
        if (fm < 0.)
            x1 = xm;
        else
            x2 = xm;
        double q = 0.5 * (x1 + x2);
        d = (q - xm) / xm;
        xm = q;
        fm = F(xm, e, s);
    }

    return xm;
}

template <Screening ScreeningType>
double xs_cms_base<ScreeningType>::theta_impulse_approx(double e, double s)
{
    const auto &C = Phi::C;
    const auto &A = Phi::A;
    double th(0.0);
    for (int i = 0; i < Phi::N; ++i)
        th += C[i] * A[i] * detail::bessel_k1(A[i] * s);
    return th / e;
}

// Unscreened Coulomb specializations
template <>
inline double xs_cms_base<Screening::None>::apsis(double e, double s)
{
    double x = 1.0 / (2 * e);
    return x + std::sqrt(x * x + s * s);
}

template <>
inline double xs_cms_base<Screening::None>::theta_impulse_approx(double e, double s)
{
    // return the exact analutic theta
    double x = 2 * e * s;
    x = 1. / (1. + x * x);
    return 2 * std::asin(std::sqrt(x));
}

// xs_cms

template <Screening ScreeningType, Quadrature QuadType, int N>
double xs_cms<ScreeningType, QuadType, N>::find_s(double e, double thetaCM, double tol)
{

    if (thetaCM == 0.)
        return std::numeric_limits<double>::infinity();

    // inital guesses: Mendenhall & Weller NIMB 58 (1991) 11, eqs. 23-25
    double d = thetaCM / M_PI;
    double gamma = 1. - d;

    if (gamma == 0.)
        return 0;

    double e0 = (d < 10 * std::numeric_limits<double>::epsilon()) ? 2. * d * e
                                                                  : (1.0 - gamma * gamma) * e;
    double x0 = Base_t_::apsis(e0, 1.e-8);
    double x1, x2;
    if (e >= 1.) {
        x1 = 0.7 * gamma * x0;
        x2 = 1.0 / (2.0 * e * tan(thetaCM / 2.0));
    } else {
        x1 = 0.9 * gamma * x0;
        x2 = 1.4 * gamma * x0;
    }

    if (x2 > 1.e4)
        x2 = 1.e4; // this is as high as it gets. Above causes numerical instability

    // ensure bracketing
    while (thetaCM - theta(e, x1) >= 0.0)
        x1 *= 0.1;
    while (thetaCM - theta(e, x2) <= 0.0 && x2 <= 1.e4)
        x2 *= 1.001;

    if (!(thetaCM - theta(e, x1) < 0.0
          && thetaCM - theta(e, x2) > 0.0)) // values should be on each side of 0
        return std::numeric_limits<double>::quiet_NaN();

    double xm = 0.5 * (x1 + x2);
    double fm = thetaCM - theta(e, xm);
    d = 1.;
    int k = 0;

    while (std::abs(d) > tol && k < 100) {
        if (fm < 0.)
            x1 = xm;
        else
            x2 = xm;
        double q = 0.5 * (x1 + x2);
        d = (q - xm) / xm;
        xm = q;
        fm = thetaCM - theta(e, xm);
        k++;
    }

    return xm;
}

template <Screening ScreeningType, Quadrature QuadType, int N>
double xs_cms<ScreeningType, QuadType, N>::crossSection(double e, double thetaCM, double tol)
{
    // find corresponfding reduced impact parameter
    double s = find_s(e, thetaCM, tol);

    if (s < 1.e-6) { // quasi head on collision

        if (s == 0.) {
            double ds = 1.e-9;
            double dth = pi_minus_theta(e, ds);
            double dsdth = ds / dth;
            return dsdth * dsdth;
        }

        double dth = pi_minus_theta(e, s);
        double ds = s * 0.001;
        double dsdTheta = (12.0 * ds)
                / (-pi_minus_theta(e, s + 2.0 * ds) + 8.0 * pi_minus_theta(e, s + ds)
                   - 8.0 * pi_minus_theta(e, s - ds) + pi_minus_theta(e, s - 2.0 * ds));

        return s / sin(dth) * fabs(dsdTheta);
    }

    // ds/dTheta using five-point stencil
    double ds = s * 0.001;
    double dsdTheta = (12.0 * ds)
            / (-theta(e, s + 2.0 * ds) + 8.0 * theta(e, s + ds) - 8.0 * theta(e, s - ds)
               + theta(e, s - 2.0 * ds));

    return s / sin(thetaCM) * fabs(dsdTheta);
}

// xs_cms template specializations for Unscreened Coulomb (analytic results)
template <>
inline double xs_cms<Screening::None, Quadrature::None, 0>::sin2Thetaby2(double e, double s)
{
    double x = 2 * e * s;
    return 1. / (1. + x * x);
}
template <>
inline double xs_cms<Screening::None, Quadrature::None, 0>::theta(double e, double s)
{
    return 2 * std::asin(std::sqrt(sin2Thetaby2(e, s)));
}
template <>
inline double xs_cms<Screening::None, Quadrature::None, 0>::find_s(double e, double thetaCM, double)
{
    double x = std::sin(thetaCM / 2);
    return 0.5 / e * std::sqrt(1. / (x * x) - 1.);
}
template <>
inline double xs_cms<Screening::None, Quadrature::None, 0>::crossSection(double e, double thetaCM,
                                                                         double)
{
    double x = std::sin(thetaCM / 2);
    x *= x;
    x *= 4 * e;
    x *= x;
    return 1. / x;
}
template <>
inline double xs_cms<Screening::None, Quadrature::None, 0>::sn(double, double)
{
    return std::numeric_limits<double>::infinity();
}

// xs_cms template specializations for ZBL + MAGIC

// this is an approximation for ZBL stopping cross-section
template <>
inline double xs_cms<Screening::ZBL, Quadrature::Magic, 0>::sn(double e, double theta_max)
{
    return 0.5 * std::log(1 + 1.1383 * e)
            / (e + 0.01321 * std::pow(e, 0.21226) + 0.19593 * std::sqrt(e));
}

// detail::theta_integrator implementation for MAGIC

/*
 * Implementation of MAGIC formula
 *
 * Returns \f$ \cos(\theta/2) \f$ where \f$ \theta \f$ is the center-of-mass
 * scattering angle
 *
 * @param e is the reduced energy in center of mass system
 * @param s is the reduced impact parameter
 * @return \f$ \cos(\theta/2) \f$
 */
template <Screening ScreeningType, int N>
inline double detail::theta_integrator<ScreeningType, Quadrature::Magic, N>::cosThetaBy2(double e,
                                                                                         double s)
{
    double cost2; /* cos(theta/2)*/
    double RoC, Delta, R, RR, A, G, alpha, beta, gamma, V, V1, FR, FR1, Q;
    double SQE;

    /* TRIM 85:  */
    static const double C[] = { 0., 0.99229, 0.011615, 0.0071222, 14.813, 9.3066 };

    if (s == 0.0)
        return 0.0;

    /* Initial guess for R: */
    R = s;
    RR = -2.7 * log(e * s);
    if (RR >= s) {
        /*   if(RR<B) calc potential; */
        RR = -2.7 * log(e * RR);
        if (RR >= s) {
            R = RR;
        }
    }
    /* TRIM85: 330 */
    do {
        /* Calculate potential and its derivative */
        V = zbl_and_deriv(R, &V1);
        FR = s * s / R + V * R / e - R;
        FR1 = -s * s / (R * R) + (V + V1 * R) / e - 1.0;
        Q = FR / FR1;
        R = R - Q;
    } while (fabs(Q / R) > 0.001);

    RoC = -2.0 * (e - V) / V1;
    SQE = sqrt(e);

    alpha = 1 + C[1] / SQE;
    beta = (C[2] + SQE) / (C[3] + SQE); /* TRIM85: CC */
    gamma = (C[4] + e) / (C[5] + e);
    A = 2 * alpha * e * pow(s, beta);
    G = gamma / (sqrt((1.0 + A * A)) - A); /* TRIM85: 1/FF */
    Delta = A * (R - s) / (1 + G);

    cost2 = (s + RoC + Delta) / (R + RoC);
    return std::max(0.0, cost2);
}

/*
 * @brief ZBL potential
 *
 * Evaluate the ZBL potential and optionally the 1st derivative dV/dR
 *
 * Implementation taken from ZBL85
 *
 * @param R is the reduced radius
 * @param Vprime if a non-NULL pointer is passed, it receives the value of dV/dR
 * @return the value of the potential
 */
template <Screening ScreeningType, int N>
inline double
detail::theta_integrator<ScreeningType, Quadrature::Magic, N>::zbl_and_deriv(double R,
                                                                             double *Vprime)
{
    const auto &C = screening_function<Screening::ZBL>::C;
    const auto &A = screening_function<Screening::ZBL>::A;

    double EX1 = C[0] * exp(-A[0] * R);
    double EX2 = C[1] * exp(-A[1] * R);
    double EX3 = C[2] * exp(-A[2] * R);
    double EX4 = C[3] * exp(-A[3] * R);

    double V = (EX1 + EX2 + EX3 + EX4) / R;
    if (Vprime)
        *Vprime = -(V + A[0] * EX1 + A[1] * EX2 + A[2] * EX3 + A[3] * EX4) / R;

    return V;
}

// xs_lab implementation

template <Screening ScreeningType, Quadrature QuadType, int N>
void xs_lab<ScreeningType, QuadType, N>::scatter(double E, double S, double &recoil_erg,
                                                 double &theta) const
{
    double e = E * red_E_conv_;
    double sin2thetaby2 = sin2Thetaby2(e, S / screening_length_);
    recoil_erg = E * gamma_ * sin2thetaby2;
    /* convert scattering angle to lab frame of reference: */
    double costheta = 1.f - 2 * sin2thetaby2;
    double sintheta = std::sqrt(1.f - costheta * costheta);
    theta = atan(sintheta / (costheta + mass_ratio_));
}

template <Screening ScreeningType, Quadrature QuadType, int N>
double xs_lab<ScreeningType, QuadType, N>::find_p(double E, double T) const
{
    double thetaCM = 1.0 * T / E / gamma_;
    if (thetaCM > 1.0)
        return std::numeric_limits<double>::quiet_NaN();
    if (thetaCM == 1.0)
        return 0.f;
    thetaCM = 2. * std::asin(std::sqrt(thetaCM));
    return find_s(E * red_E_conv_, thetaCM) * screening_length_;
}

template <Screening ScreeningType, Quadrature QuadType, int N>
double xs_lab<ScreeningType, QuadType, N>::crossSection(double E, double T) const
{
    double thetaCM = T / E / gamma_;
    if (thetaCM > 1.0)
        return std::numeric_limits<double>::quiet_NaN();
    thetaCM = 2. * std::asin(std::sqrt(thetaCM));
    return crossSection(E * red_E_conv_, thetaCM) * 4 * sig0_ / E / gamma_;
}

template <Screening ScreeningType, Quadrature QuadType, int N>
double xs_lab<ScreeningType, QuadType, N>::Sn(double E) const
{
    return sn(E * red_E_conv_) * sig0_ * gamma_ / red_E_conv_;
}

template <Screening ScreeningType, Quadrature QuadType, int N>
double xs_lab<ScreeningType, QuadType, N>::Sn(double E, double T1) const
{
    double theta_max = 1.0 * T1 / E / gamma_;
    if (theta_max >= 1.)
        return Sn(E);
    theta_max = 2. * std::asin(std::sqrt(theta_max));
    return sn(E * red_E_conv_, theta_max) * sig0_ * gamma_ / red_E_conv_;
}

#endif
