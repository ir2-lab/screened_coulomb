#include "screened_coulomb.h"

#include <iomanip>
#include <iostream>
#include <chrono>
#include <random>

using namespace std::chrono;
using namespace std;

#define NGS 8
#define NITER 3000000

template <Quadrature Q_, int N_>
void test1();

template <Screening S_, Quadrature Q_, int N>
void test2();

template <Screening S_, Quadrature Q_, int N>
void test3(size_t M);

int main()
{
    const char *sep = "\n==========================================================\n\n";

    test1<Quadrature::GaussChebyshev, NGS>();

    cout << sep;
    test1<Quadrature::Lobatto4, 4>();

    cout << sep;
    test2<Screening::ZBL, Quadrature::GaussChebyshev, NGS>();

    cout << sep;
    test2<Screening::ZBL, Quadrature::Lobatto4, 4>();

    cout << sep;
    test2<Screening::ZBL, Quadrature::Magic, 0>();

    cout << sep;
    test2<Screening::Moliere, Quadrature::GaussChebyshev, NGS>();

    cout << sep;
    test2<Screening::Bohr, Quadrature::GaussChebyshev, NGS>();

    cout << sep;
    test2<Screening::KrC, Quadrature::GaussChebyshev, NGS>();

    cout << sep;
    test3<Screening::ZBL, Quadrature::GaussChebyshev, NGS>(NITER);

    cout << sep;
    test3<Screening::ZBL, Quadrature::Lobatto4, 4>(NITER);

    cout << sep;
    test3<Screening::ZBL, Quadrature::Magic, 0>(NITER);

    cout << sep;
    test3<Screening::None, Quadrature::None, 0>(NITER);

    return 0;
}

// data from Robinson1970 tables, e=0.01, s=0:1:40

const double apsis_ref[] = { 5.9381104, 5.9530192, 5.9975794, 6.0713055, 6.1734304, 6.3029615,
                             6.4587447, 6.6395283, 6.8440201, 7.0709315, 7.3190089, 7.5870514,
                             7.8739176, 8.178525,  8.4998435, 8.8368866, 9.1887023, 9.5543651,
                             9.93297,   10.323629, 10.725471, 11.137642, 11.559309, 11.989665,
                             12.427932, 12.873369, 13.325273, 13.782983, 14.245885, 14.71341,
                             15.185037, 15.660293, 16.138745, 16.620009, 17.103738, 17.589624,
                             18.077392, 18.566801, 19.057638, 19.549715, 20.04287 };

const double th_ref[] = { 180,       166.50335, 153.18169, 140.19836,  127.69528, 115.78643,
                          104.55506, 94.054367, 84.310654, 75.327989,  67.093278, 59.581085,
                          52.757736, 46.584562, 41.020293, 36.022722,  31.5498,   27.560328,
                          24.014371, 20.87352,  18.101051, 15.66205,   13.523507, 11.654404,
                          10.025784, 8.6108004, 7.3847416, 6.3250273,  5.4111727, 4.6247229,
                          3.9491627, 3.3698044, 2.8736608, 2.4493091,  2.0867502, 1.7772686,
                          1.513297,  1.2882864, 1.0965865, 0.93333417, 0.79435364 };

const double xs_ref[] = { 56.369466, 57.394136, 60.549823, 66.086753, 74.44714,  86.304819,
                          102.62665, 124.762,   154.5703,  194.60113, 248.34784, 320.605,
                          417.97248, 549.5677,  728.03342, 970.96571, 1302.9398, 1758.3872,
                          2385.6841, 3252.9668, 4456.4034, 6131.963,  8472.1601, 11749.873,
                          16352.21,  22828.643, 31959.383, 44852.432, 63081.263, 88879.982,
                          125419.73, 177199.91, 250601.17, 354668.06, 502210.94, 711361.7,
                          1007778,   1427707.5, 2022374.8, 2864069.2, 4054736.6 };

// const double th_ref[] = { 180,        171.14606,  162.29865,  153.46457,   144.65127,  135.86718,
//                           127.12211,  118.42771,  109.79808,  101.25039,   92.805655,  84.4896,
//                           76.333569,  68.375389,  60.660024,  53.239765,   46.17356,   39.525077,
//                           33.359086,  27.736059,  22.70541,   18.298516,   14.523256,  11.361777,
//                           8.7724046,  6.6952387,  5.0598311,  3.7929365,   2.8247631,  2.0930139,
//                           1.5447866,  1.13684,    0.83484556, 0.61214656,  0.44838561,
//                           0.32820827, 0.24013997, 0.17566357, 0.12848826, 0.093984122,
//                           0.068752204 };

// const double xs_ref[] = { 526.11038, 528.4672,  535.62717, 547.86486, 565.66337,    589.75565,
//                           621.19073, 661.43445, 712.5211,  777.28211, 859.69433,    965.4189,
//                           1102.6486, 1283.4673, 1526.0734, 1858.4956, 2324.9253,    2996.7371,
//                           3992.0332, 5510.9246, 7900.2377, 11773.792, 18238.424,    29322.479,
//                           48793.087, 83722.136, 147495.56, 265605.31, 486807.16,    904623.75,
//                           1698761.3, 3214856.6, 6117806.4, 11686785,  22381146,     42925414,
//                           82391622,  158176770, 303624660, 582537600, 1.1168471e+09 };

template <Quadrature Q_, int N_>
void test1()
{
    typedef xs_cms<Screening::Moliere, Q_, N_> XS;

    cout << "TEST-1" << endl;
    cout << "Compare screened_coulomb to refence data from Robinson1970" << endl;
    cout << "Screening:  " << XS::screeningName() << endl;
    cout << "Quadrature: " << XS::quadratureName() << endl;
    cout << "Order: " << XS::quadratureOrder() << endl;
    cout << endl;

    double e = 0.01;
    cout << "e = " << e << endl << endl;

    int w1 = 12;
    cout << setw(6) << "s" << ' ';
    cout << setw(3 * w1 + 2) << "Apsis R (a)" << ' ';
    cout << setw(3 * w1 + 2) << "Scattering Angle Th (Deg)" << ' ';
    cout << setw(3 * w1 + 3) << "Cross-section XS (a^2/4π)" << endl;
    cout << endl;
    cout << setw(6) << " " << ' ';
    cout << setw(w1) << "Rref" << ' ';
    cout << setw(w1) << "R" << ' ';
    cout << setw(w1) << "Rel. Diff" << ' ';
    cout << setw(w1) << "Thref" << ' ';
    cout << setw(w1) << "  Th" << ' ';
    cout << setw(w1) << "Rel. Diff" << ' ';
    cout << setw(w1) << "XSref" << ' ';
    cout << setw(w1) << "XS" << ' ';
    cout << setw(w1) << "Rel. Diff" << endl;

    for (int i = 0; i <= 40; ++i) {
        double s = 0.5 * i;
        double a = XS::apsis(e, s);
        double th = XS::theta(e, s);
        double xs = XS::crossSection(e, th) * 4 * M_PI;
        th *= (180.0 / M_PI);
        cout << setprecision(3) << setw(6) << s << ' ';
        cout << setprecision(9);
        cout << setw(w1) << apsis_ref[i] << ' ';
        cout << setw(w1) << a << ' ';
        cout << setprecision(2) << setw(w1) << (apsis_ref[i] - a) / apsis_ref[i] << ' ';
        cout << setprecision(9) << setw(w1) << th_ref[i] << ' ';
        cout << setw(w1) << th << ' ';
        cout << setprecision(2) << setw(w1) << (th_ref[i] - th) / th_ref[i] << ' ';
        cout << setprecision(9) << setw(w1) << xs_ref[i] << ' ';
        cout << setw(w1) << xs << ' ';
        cout << setprecision(2) << setw(w1) << (xs_ref[i] - xs) / xs_ref[i] << endl;
    }
}

template <Screening S_, Quadrature Q_, int N>
void test2()
{
    double e[] = { 1e-3, 1e-3, 0.1, 0.1, 10., 10. };
    double s[] = { 0.5, 20, 0.2, 8.0, 0.025, 1.0 };

    typedef xs_cms<S_, Quadrature::GaussChebyshev, 128> xs0;
    typedef xs_cms<S_, Q_, N> xs1;

    cout << "TEST-2" << endl;
    cout << "Compare low order quadrature to high-precision Gauss-Cheb. ord. 128" << endl;
    cout << "Screening:  " << xs1::screeningName() << endl;
    cout << "Quadrature: " << xs1::quadratureName() << endl;
    cout << "Order: " << xs1::quadratureOrder() << endl;
    cout << endl;

    int w1 = 14;

    cout << setw(6) << "e" << ' ';
    cout << setw(6) << "s" << ' ';
    cout << setw(3 * w1 + 2) << "Scattering Angle Th (Deg)" << ' ';
    cout << setw(3 * w1 + 2) << "Cross-section XS (a^2/4π)" << ' ';
    cout << setw(3 * w1 + 3) << "stopping XS sn" << endl;
    cout << endl;
    cout << setw(6) << " " << ' ';
    cout << setw(6) << " " << ' ';
    cout << setw(w1) << "TH" << ' ';
    cout << setw(w1) << "TH0" << ' ';
    cout << setw(w1) << "Rel. Diff" << ' ';
    cout << setw(w1) << "XS" << ' ';
    cout << setw(w1) << "XS0" << ' ';
    cout << setw(w1) << "Rel. Diff" << ' ';
    cout << setw(w1) << "SN" << ' ';
    cout << setw(w1) << "SN0" << ' ';
    cout << setw(w1) << "Rel. Diff" << endl;

    for (int i = 0; i < 6; ++i) {
        double th1 = xs1::theta(e[i], s[i]);
        double x1 = xs1::crossSection(e[i], th1);
        double sn1 = xs1::sn(e[i]);
        double th0 = xs0::theta(e[i], s[i]);
        double x0 = xs0::crossSection(e[i], th0);
        double sn0 = xs0::sn(e[i]);
        th1 *= 180 / M_PI;
        th0 *= 180 / M_PI;
        cout << setprecision(2);
        cout << setw(6) << e[i] << ' ';
        cout << setw(6) << s[i] << ' ';
        cout << setprecision(9) << setw(w1) << th1 << ' ';
        cout << setw(w1) << th0 << ' ';
        cout << setprecision(2) << setw(w1) << (th1 - th0) / th0 << ' ';
        cout << setprecision(9) << setw(w1) << x1 << ' ';
        cout << setw(w1) << x0 << ' ';
        cout << setprecision(2) << setw(w1) << (x1 - x0) / x0 << ' ';
        cout << setprecision(9) << setw(w1) << sn1 << ' ';
        cout << setw(w1) << sn0 << ' ';
        cout << setprecision(2) << setw(w1) << (sn1 - sn0) / sn0 << endl;
    }
}

template <Screening S_, Quadrature Q_, int N>
void test3(size_t M)
{
    double e[] = { 1e-3, 1e-3, 0.1, 0.1, 10., 10. };
    double s[] = { 0.5, 20, 0.2, 8.0, 0.025, 1.0 };

    typedef xs_cms<S_, Q_, N> XS;

    cout << "TEST-3" << endl;
    cout << "Timing benchmark" << endl;
    cout << "Screening:   " << XS::screeningName() << endl;
    cout << "Quadrature:  " << XS::quadratureName() << endl;
    cout << "Order:  " << XS::quadratureOrder() << endl;
    cout << "Iterations:  " << M << endl;
    cout << endl;

    random_device rd; // Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<> ue(-3.0, 2.0);
    uniform_real_distribution<> us(0.0, 10.0);

    high_resolution_clock::time_point t1, t2;
    duration<double, std::nano> dt1, dt2;

    cout << "Running ... ";
    cout.flush();

    t1 = high_resolution_clock::now();

    double th = 0.0, se = 0.0, ss = 0.0;

    for (size_t i = 0; i < M; i++) {
        double e = std::pow(10.0, ue(gen));
        double s = us(gen);
        se += e;
        ss += s;
        th += XS::theta(e, s);
    }

    t2 = high_resolution_clock::now();

    dt1 = t2 - t1;

    cout << "done." << endl;

    cout << "Timing rngs ... ";
    cout.flush();

    t1 = high_resolution_clock::now();

    th = 0.0, se = 0.0, ss = 0.0;
    gen.seed();

    for (size_t i = 0; i < M; i++) {
        double e = std::pow(10.0, ue(gen));
        double s = us(gen);
        se += e;
        ss += s;
    }

    t2 = high_resolution_clock::now();

    dt2 = t2 - t1;

    cout << "done." << endl << endl;

    cout << setprecision(3);
    cout << "dt/theta: " << (dt1.count() - dt2.count()) / M << " ns" << std::endl;
}
