#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>

#define _USE_MATH_DEFINES

using namespace std;
using ldouble = long double;
using lcomplex = complex<long double>;

ldouble theta = 0.05L;
ldouble q = 1.L;

ldouble om = 0.057L;
ldouble c = 137.03599L;

long int n = 21000000001;
ldouble a_0 = 2.L;
ldouble a_1 = 2600.L;
ldouble a_2 = 0.001L;

ldouble q2 = q * q;
ldouble step_1 = a_1 / n;
ldouble kp = 0.5L * step_1;

lcomplex integrate_0(lcomplex f(ldouble x)) {
    ldouble step = a_0 / n;
    lcomplex result = 0.L;

    for (long int l = 0; l <= n; l++) {
        result += f(l * step) * step;
    }
    return 2.L * result;
}
lcomplex integrate_1(lcomplex f(ldouble x)) {
    lcomplex result = 0.L;

    for (long int l = 0; l <= n; l++) {
        /*if (abs(kp - l * step) >= step) {
            result += f(l * step) * step;
        } else if (kp < l * step) {
            result += f(kp + step) * step;
        } else {
            result += f(kp - step) * step;
        }*/
        result += f(l * step_1) * step_1;
    }
    return 2.L * result;
}
lcomplex integrate_2(lcomplex f(ldouble x)) {
    ldouble step_2 = 2.L * a_2 / n;
    lcomplex result = 0.L;

    for (long int l = 0; l <= n; l++) {
        result += f(kp - a_2 + l * step_2) * step_2;
    }
    return 2.L * result;
}
lcomplex z_0(ldouble k) {
    return 1.L / (q2 * q2 * exp(-2.il * theta) * c * c * k * k / 2.L + 2.L * M_PI * M_PI * c * c * (om * om - exp(-2.il * theta) * c * c * k * k) * (om * om - exp(-2.il * theta) * c * c * k * k));
}
lcomplex z_1(ldouble k) {
    return k * k * c * (1.L / (kp - k) + 1.L / (kp + k)) / (q2 * q2 * exp(-2.il * theta) * k * k + 4.L * M_PI * M_PI * (om * om - exp(-2.il * theta) * c * c * k * k) * (om * om - exp(-2.il * theta) * c * c * k * k));
}
int main() {
    lcomplex Omega = sqrt(om * om - q2 * q2 / (16.L * M_PI * M_PI * c * c)) - 1il * q2 / (4.L * M_PI * c);
    /*
        lcomplex K_0 = 1.L / (1.L - 1.il * q2 * Omega / (4.L * M_PI * c * om * om)) + exp(-1.il * theta) * c * c * om * om * q2 * integrate_0(z_0);
        cout << "(1, 0) = " << setprecision(10) << K_0 << endl;
    */
    lcomplex K_1 = Omega * Omega * c * kp / ((Omega * Omega - exp(-2.il * theta) * c * c * kp * kp) * (om * om - 1.il * q2 * Omega / (4.L * M_PI * c))) + c * kp * (exp(-2.il * theta) * c * c * kp * kp - om * om) / (q2 * q2 * exp(-2.il * theta) * kp * kp / (4.L * M_PI * M_PI) + (om * om - exp(-2.il * theta) * c * c * kp * kp) * (om * om - exp(-2.il * theta) * c * c * kp * kp)) + q2 * exp(-1.il * theta) * integrate_1(z_1);
    cout << "k = " << kp << " : (0, 0) = " << setprecision(10) << K_1 << endl;

    kp = step_1 * 105000.5L;
    K_1 = Omega * Omega * c * kp / ((Omega * Omega - exp(-2.il * theta) * c * c * kp * kp) * (om * om - 1.il * q2 * Omega / (4.L * M_PI * c))) + c * kp * (exp(-2.il * theta) * c * c * kp * kp - om * om) / (q2 * q2 * exp(-2.il * theta) * kp * kp / (4.L * M_PI * M_PI) + (om * om - exp(-2.il * theta) * c * c * kp * kp) * (om * om - exp(-2.il * theta) * c * c * kp * kp)) + q2 * exp(-1.il * theta) * integrate_1(z_1);
    cout << "k = " << kp << " : (0, 0) = " << setprecision(10) << K_1 << endl;

    kp = step_1 * 10500000.5L;
    K_1 = Omega * Omega * c * kp / ((Omega * Omega - exp(-2.il * theta) * c * c * kp * kp) * (om * om - 1.il * q2 * Omega / (4.L * M_PI * c))) + c * kp * (exp(-2.il * theta) * c * c * kp * kp - om * om) / (q2 * q2 * exp(-2.il * theta) * kp * kp / (4.L * M_PI * M_PI) + (om * om - exp(-2.il * theta) * c * c * kp * kp) * (om * om - exp(-2.il * theta) * c * c * kp * kp)) + q2 * exp(-1.il * theta) * integrate_1(z_1);
    cout << "k = " << kp << " : (0, 0) = " << setprecision(10) << K_1 << endl;

    kp = step_1 * 1050000000.5L;
    K_1 = Omega * Omega * c * kp / ((Omega * Omega - exp(-2.il * theta) * c * c * kp * kp) * (om * om - 1.il * q2 * Omega / (4.L * M_PI * c))) + c * kp * (exp(-2.il * theta) * c * c * kp * kp - om * om) / (q2 * q2 * exp(-2.il * theta) * kp * kp / (4.L * M_PI * M_PI) + (om * om - exp(-2.il * theta) * c * c * kp * kp) * (om * om - exp(-2.il * theta) * c * c * kp * kp)) + q2 * exp(-1.il * theta) * integrate_1(z_1);
    cout << "k = " << kp << " : (0, 0) = " << setprecision(10) << K_1 << endl;

    kp = step_1 * 10500000000.5L;
    K_1 = Omega * Omega * c * kp / ((Omega * Omega - exp(-2.il * theta) * c * c * kp * kp) * (om * om - 1.il * q2 * Omega / (4.L * M_PI * c))) + c * kp * (exp(-2.il * theta) * c * c * kp * kp - om * om) / (q2 * q2 * exp(-2.il * theta) * kp * kp / (4.L * M_PI * M_PI) + (om * om - exp(-2.il * theta) * c * c * kp * kp) * (om * om - exp(-2.il * theta) * c * c * kp * kp)) + q2 * exp(-1.il * theta) * integrate_1(z_1);
    cout << "k = " << kp << " : (0, 0) = " << setprecision(10) << K_1 << endl;

    kp = step_1 * 105000.5L;
    K_1 = Omega * Omega * c * kp / ((Omega * Omega - exp(-2.il * theta) * c * c * kp * kp) * (om * om - 1.il * q2 * Omega / (4.L * M_PI * c))) + c * kp * (exp(-2.il * theta) * c * c * kp * kp - om * om) / (q2 * q2 * exp(-2.il * theta) * kp * kp / (4.L * M_PI * M_PI) + (om * om - exp(-2.il * theta) * c * c * kp * kp) * (om * om - exp(-2.il * theta) * c * c * kp * kp)) + q2 * exp(-1.il * theta) * integrate_2(z_1);
    cout << endl
         << "Around k = " << kp << " : " << setprecision(10) << K_1 << endl;

    n = 2100000001;
    kp = step_1 * 105000.5L;
    K_1 = Omega * Omega * c * kp / ((Omega * Omega - exp(-2.il * theta) * c * c * kp * kp) * (om * om - 1.il * q2 * Omega / (4.L * M_PI * c))) + c * kp * (exp(-2.il * theta) * c * c * kp * kp - om * om) / (q2 * q2 * exp(-2.il * theta) * kp * kp / (4.L * M_PI * M_PI) + (om * om - exp(-2.il * theta) * c * c * kp * kp) * (om * om - exp(-2.il * theta) * c * c * kp * kp)) + q2 * exp(-1.il * theta) * integrate_2(z_1);
    cout << "10x less points, Around k = " << kp << " : " << setprecision(10) << K_1 << endl;

    n = 210000001;
    kp = step_1 * 105000.5L;
    K_1 = Omega * Omega * c * kp / ((Omega * Omega - exp(-2.il * theta) * c * c * kp * kp) * (om * om - 1.il * q2 * Omega / (4.L * M_PI * c))) + c * kp * (exp(-2.il * theta) * c * c * kp * kp - om * om) / (q2 * q2 * exp(-2.il * theta) * kp * kp / (4.L * M_PI * M_PI) + (om * om - exp(-2.il * theta) * c * c * kp * kp) * (om * om - exp(-2.il * theta) * c * c * kp * kp)) + q2 * exp(-1.il * theta) * integrate_2(z_1);
    cout << "100x less points, Around k = " << kp << " : " << setprecision(10) << K_1 << endl;

    n = 21000001;
    kp = step_1 * 105000.5L;
    K_1 = Omega * Omega * c * kp / ((Omega * Omega - exp(-2.il * theta) * c * c * kp * kp) * (om * om - 1.il * q2 * Omega / (4.L * M_PI * c))) + c * kp * (exp(-2.il * theta) * c * c * kp * kp - om * om) / (q2 * q2 * exp(-2.il * theta) * kp * kp / (4.L * M_PI * M_PI) + (om * om - exp(-2.il * theta) * c * c * kp * kp) * (om * om - exp(-2.il * theta) * c * c * kp * kp)) + q2 * exp(-1.il * theta) * integrate_2(z_1);
    cout << "1000x less points, Around k = " << kp << " : " << setprecision(10) << K_1 << endl;

    return 0;
}