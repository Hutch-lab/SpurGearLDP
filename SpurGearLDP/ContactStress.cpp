#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "FunctDefs.h"
#include <tuple>
#include <fstream>
using namespace std;
// This is a static model that doesn't make physical sense. Going to create a dynamic model.
double calcContactStress(const SpurGear* pinion, const SpurGear* gear, double Wt, double V, double Ko, int Qv, int gearQuality, double Ks, double ZR) {
    /* Shigley's Mechanical Design Equation 14-16.
     sigma_c = Z_e * sqrt(Wt * Ko * Kv * Ks * KH / (dw1 * b) * ZR/ZI)

     Definitions
        sigma_c: contact stress
        Ze: elastic coefficient [sqrt(N/mm^2) or sqrt(MPa)] - eq 14-13
        Wt: tangential transmitted load [N] - fig 14-1
        Ko: overload factor - Figure 14-17/14-18
        Kv: dynamic factor - eq 14-27
        Ks: size factor - 14-10 a
        KH: load-distribution factor - eq 14-30
        dw1: pitch diameter of pinion [mm]
        b: facewidth of the narrower member [mm]
        ZR: surface condition factor - greater than 1 - eq 14-16
        Zl: geometry factor for pitting resistance - eq 14-23
     */

     // Input Validation
    if (!pinion || !gear || Wt < 0 || V < 0 || Ko < 1.0 || Ks <= 0 || ZR <= 0 || pinion->pitchDiameter <= 0 || pinion->faceWidth <= 0 || gear->faceWidth <= 0 || pinion->youngsModulus <= 0 || gear->youngsModulus <= 0 || pinion->gearRatio <= 0) { return -1.0; }

    double dp = pinion->pitchDiameter; // Pitch diameter of pinion
    double F = min(pinion->faceWidth, gear->faceWidth); // Face width of the narrower member

    // --- Convert units for AGMA formulas ---
// V needs to be in ft/min for the Kv formula used
    double V_fpm = V * 196.850; // Convert m/s to ft/min (1 m/s = 196.85 ft/min)
    // F needs to be in inches for the KH formula used
    double F_in = F / 25.4; // Convert mm to inches

    // Calculate Factors
    double Ze = calculateZe(pinion, gear); // Elastic Coefficient
    double ZI = calculateZI(pinion, gear); // Geometry Factor (dimensionless)
    double Kv = calculateKv(V_fpm, Qv);    // Dynamic Factor (dimensionless)
    double KH = calculateKH(pinion, gear, F_in, gearQuality); // Load Distribution Factor (dimensionless)

    if (Ze <= 0 || ZI <= 0 || Kv <= 0 || KH <= 0) {
        cerr << "Error: Calculated factors (Ze, ZI, Kv, KH) must be positive." << endl;
        return -1.0;
    }


    return Ze * sqrt((Wt * Ko * Kv * Ks * KH * ZR) / (dp * F * ZI));
}

double calculateZe(const SpurGear* pinion, const SpurGear* gear) {
    double nu_p = pinion->poissonsRatio;
    double E_p = pinion->youngsModulus;
    double nu_g = gear->poissonsRatio;
    double E_g = gear->youngsModulus;

    if (E_p <= 0 || E_g <= 0) {
        cerr << "Error: Young's Modulus must be positive." << endl;
        return 0.0;
    }

    double term_p = (1.0 - nu_p * nu_p) / E_p;
    double term_g = (1.0 - nu_g * nu_g) / E_g;

    return sqrt(1.0 / (M_PI * (term_p + term_g)));
}

double calculateZI(const SpurGear* pinion, const SpurGear* gear) {
    double phi_rad = pinion->pressureAngle * M_PI / 180.0;
    double mg = gear->gearRatio;

    if (mg <= 0) {
        cerr << "Error: Invalid gear ratio." << endl;
        return 0.0;
    }

    return (cos(phi_rad) * sin(phi_rad) / 2.0) * (mg / (mg + 1.0));
}

double calculateKv(double V, int Qv) {
    // V: Pitch line velocity (ft/min)
    // Qv: Transmission accuracy level (e.g., 6-11, lower is better)

    if (Qv < 3 || Qv > 11) {
        cerr << "Warning: Qv = " << Qv << " is outside the typical range (3-11). Using Qv=11." << endl;
        Qv = 11; // Clamp value
    }
    if (V < 0) {
        cerr << "Error: Pitch line velocity cannot be negative." << endl;
        return 1.0;
    }
    if (V == 0) { return 1.0; } // Static case

    // AGMA formula (approximate, check standard for specifics)
    double B = 0.25 * pow((12.0 - Qv), (2.0 / 3.0));
    double A = 50.0 + 56.0 * (1.0 - B);

    // Ensure A is not zero and pow is non-negative
    if (A == 0) {
        cerr << "Error: Calculated A is zero in Kv calculation." << endl;
    }
    if (V < 0) V = 0; // Velocity shouldn't be negative

    double Kv = pow(((A + sqrt(V)) / A), B);

    return Kv; // generally be >= 1.0
}

double calculateKH(const SpurGear* pinion, const SpurGear* /*gear*/, double F, int gearQuality) {
    // F: Face width of narrower member (in)
    // gearQuality: AGMA quality number or similar index affecting Cma

    double dp = pinion->pitchDiameter; // Pinion pitch diameter

    // Cmc: Lead correction factor (1.0 for no crowning, < 1.0 for crowning)
    double Cmc = 1.0; // Assume no lead correction unless specified

    // Cpf: Pinion proportion factor
    double Cpf = 0.0;
    double F_inches = F / 25.4; // Convert F (assumed mm) to inches for AGMA formulas
    double dp_inches = dp / 25.4; // Convert dp (assumed mm) to inches

    double f_div_10d = F_inches / dp_inches; // F/d ratio

    if (F_inches <= 1.0) {
        Cpf = (f_div_10d / 10.0) - 0.025;
    }
    else if (F_inches <= 17.0) {
        Cpf = (f_div_10d / 10.0) - 0.0375 + 0.00469 * F_inches;
    }
    else {
        Cpf = (f_div_10d / 10.0) - 0.0375 + 0.00469 * 17.0;
        cerr << "Warning: Face width F > 17 inches, Cpf calculation inaccurate?. Capping at F=17 value." << endl;
    }
    // Ensure Cpf is not negative 
    if (Cpf < 0) Cpf = 0;

    // Cpm: Pinion proportion modifier (1.0 if pinion centered, >1.0 otherwise)
    double Cpm = 1.0; // Assume pinion is mounted symmetrically unless specified

    // Cma: Mesh alignment factor
    double A_ma, B_ma, C_ma_coeff;
    if (gearQuality >= 11) { A_ma = 0.0380; B_ma = 0.0102; C_ma_coeff = -0.822e-4; }
    else if (gearQuality >= 8) { A_ma = 0.0675; B_ma = 0.0128; C_ma_coeff = -0.926e-4; }
    else if (gearQuality >= 6) { A_ma = 0.127; B_ma = 0.0158; C_ma_coeff = -0.930e-4; }
    else { A_ma = 0.247; B_ma = 0.0167; C_ma_coeff = -0.765e-4; }

    // Using F in inches for these standard AGMA formulas
    double Cma = A_ma + B_ma * F_inches + C_ma_coeff * pow(F_inches, 2);

    // Ce: Mesh alignment correction factor (1.0 if unmodified, 0.8 if adjusted for alignment)
    double Ce = 1.0; // Assume no special alignment adjustment

    // Calculate KH
    double KH = 1.0 + Cmc * (Cpf * Cpm + Cma * Ce);

    return max(1.0, KH); // KH should be >= 1.0
}
