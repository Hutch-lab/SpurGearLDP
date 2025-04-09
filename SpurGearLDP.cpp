/* Physics 5810 - Computational Physics Final Project - Nick Hutchison */
/*  Spur Gear Load Distribution Program */
#define _USE_MATH_DEFINES
#include <cmath>
#include "FunctDefs2.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <tuple>
using namespace std;

/*
f1 = m_pinion * x_pinion_doubledot
f2 = m_pinion * y_pinion_doubledot
f3 = momentInertia_pinion * phi_pinion_doubledot
f4 = m_gear * x_gear_doubledot
f5 = m_gear * y_gear_doubledot
f6 = momentInertia_gear * phi_gear_doubledot
*/
#define f1(F_pinion_sx, F_pinion_1x, F_pinion_2x, F_pinion_3x, F_pinion_4x) (F_pinion_sx - F_pinion_1x - F_pinion_2x - F_pinion_3x - F_pinion_4x)
#define f2(F_pinion_sy, F_pinion_1y, F_pinion_2y, F_pinion_3y, F_pinion_4y) (F_pinion_sy - F_pinion_1y - F_pinion_2y - F_pinion_3y - F_pinion_4y) - m_pinion * g
#define f3(T_pinion_1, T_pinion_2, T_pinion_3, T_pinion_4, T_d) (T_pinion_1 + T_pinion_2 + T_pinion_3 + T_pinion_4 + T_d)
#define f4(F_gear_sx, F_gear_1x, F_gear_2x, F_gear_3x, F_gear_4x) (F_gear_sx - F_gear_1x - F_gear_2x - F_gear_3x - F_gear_4x)
#define f5(F_gear_sy, F_gear_1y, F_gear_2y, F_gear_3y, F_gear_4y) (F_gear_sy - F_gear_1y - F_gear_2y - F_gear_3y - F_gear_4y - m_gear * g)
#define f6(T_gear_1, T_gear_2, T_gear_3, T_gear_4, T_l) (T_gear_1 + T_gear_2 + T_gear_3 + T_gear_4 + T_l)

int main() {
    SpurGear pinion;
    SpurGear gear;
    DynamicsVars dynvars;

    startupmsg();
    cout << endl;

    // Get input method
    int method;
    cout << " Enter 0 to manually enter parameters or 1 to load parameters from file";
    cin >> method;

    if (method == 0) {
        int cont = getUserInput(&pinion, &gear);
    }
    else if (method == 1) {
        // add read from file
    }
    else {
        cout << "Input Error";
    }
    // Setup involute shape for the pinion and gear
    setupProfiles(&pinion, &gear, &dynvars);
    initialToothPinion(&pinion, &dynvars);
    initialToothGear(&gear, &dynvars);

    // Copy the tooth profile for every tooth
    completePinionProfile(&pinion, &dynvars);
    completeGearProfile(&gear, &dynvars);

    // Compute meshing line and contact points
    meshingLine(&pinion, &gear, &dynvars);

    State state;
    // Initialize state to zeros
    for (int i = 0; i < 6; i++) {
        state.pos[i] = 0.0;
        state.vel[i] = 0.0;
    }
    double t = 0.0;
    double dt = 0.001;  // Time step (adjust as needed)
    int steps = 100;   // Number of simulation steps

    cout << "\nRunning dynamic simulation:" << endl;
    for (int i = 0; i < steps; i++) {
        // Advance state using Runge-Kutta integration.
        state = rungeKuttaStep(state, t, dt, pinion, gear);
        t += dt;
        
        moveGearPosition(state, &pinion, &gear, &dynVars); // Update the current tooth profiles based on the new state.

        // Output some state information.
        cout << "Time: " << t
            << " | Pinion pos: (" << state.pos[0] << ", " << state.pos[1] << ")"
            << " | Gear pos: (" << state.pos[3] << ", " << state.pos[4] << ")" << endl;
    }

    cout << "\nSimulation complete." << endl;
    return 0;
}


void startupmsg() {
    cout << " Physics 5810 - Computational Physics Final Project ";
    cout << " ----------------- Nick Hutchison ----------------- ";
    cout << " ------ Spur Gear Load Distribution Program ------- ";
}

int getUserInput(SpurGear *pinion, SpurGear *gear){
    /*

    Returns 0 when the user wants to exit and 1 if the user enters valid input data.

    INPUT:  
    OUTPUT: 
    CALLS: 	
     */
    cout << "Enter number of teeth for the pinion and gear: ";
    cin >> pinion->numTeeth >> gear->numTeeth;
    cout << "Enter module: ";
    cin >> pinion->module;
    gear->module = pinion->module;
    cout << "Enter pressure angle: ";
    cin >> pinion->pressureAngle;
    gear->pressureAngle = pinion->pressureAngle;
    cout << "Enter effective outside diameter for the pinion and gear: ";
    cin >> pinion->outsideDiameter >> gear->outsideDiameter;
    cout << "Enter root diameter for the pinion and gear: ";
    cin >> pinion->rootDiameter >> gear->rootDiameter;
    cout << "Enter center distance between the gear and pinion (enter 0 for auto-calc): ";
    cin >> pinion->centerDistance;
    gear->centerDistance = pinion->centerDistance;
    cout << "Enter faceWidth for the pinion and gear: ";
    cin >> pinion->faceWidth >> gear->faceWidth;
    cout << "Enter face width offset: ";
    cin >> gear->offsetFaceWidth;
    cout << "Enter transverse tooth thickness and measured diameter for the pinion: ";
    cin >> pinion->transverseToothThickness >> pinion->diameterToothThickness;
    cout << "Enter transverse tooth thickness and measured diameter for the gear: ";
    cin >> gear->transverseToothThickness >> gear->diameterToothThickness;
    cout << "Enter Youngs modulus for the pinion and gear: ";
    cin >> pinion->youngsModulus >> gear->youngsModulus;
    cout << "Enter the Poissons ratio of the pinion and gear: ";
    cin >> pinion->poissonsRatio >> gear->poissonsRatio;

    if (pinion->centerDistance == 0) {
        double cd = (gear->pitchDiameter + pinion->pitchDiameter) / 2;
        gear->centerDistance = cd;
        pinion->centerDistance = cd;
    }
    pinion->gearRatio = gear->numTeeth / pinion->numTeeth;
    gear->gearRatio = pinion->gearRatio;
    pinion->pitchDiameter = (2 * pinion->centerDistance) / (pinion->gearRatio + 1);
    gear->pitchDiameter = (2 * gear->centerDistance * gear->gearRatio) / (gear->gearRatio + 1);
    pinion->baseDiameter = pinion->pitchDiameter * pinion->pressureAngle;
    gear->baseDiameter = gear->pitchDiameter * gear->pressureAngle;

    return TRUE;
}

int setupProfiles(SpurGear *pinion, SpurGear *gear, DynamicsVars *dynamicvars){
    dynamicvars->n = 100;
    pinion->theta_e = M_PI / pinion->numTeeth;
    gear->theta_e = M_PI / gear->numTeeth;
    pinion->theta_a = computeThetaA(pinion->pressureAngle);
    gear->theta_a = computeThetaA(gear->pressureAngle);

    dynamicvars->trig_pinion = -1.0 * (pinion->theta_e + 2.0 * pinion->theta_a);
    pinion->theta_r = M_PI_2 + pinion->theta_a - pinion->theta_e / 2.;
    dynamicvars->trig_gear = -1.0 * (gear->theta_e + 2.0 * gear->theta_a);
    gear->theta_r = M_PI_2 + gear->theta_a - gear->theta_e / 2.;
    double stepsize_pinion = ((pinion->outsideDiameter / 2.0) - (pinion->baseDiameter / 2.0)) / n;
    double stepsize_gear = ((gear->outsideDiameter / 2.0) - (gear->baseDiameter / 2.0)) / n;

    pinion->position_x_left.resize(n);
    pinion->position_y_left.resize(n);
    pinion->position_x_right.resize(n);
    pinion->position_y_right.resize(n);

    gear->position_x_left.resize(n);
    gear->position_y_left.resize(n);
    gear->position_x_right.resize(n);
    gear->position_y_right.resize(n);

    int i = 0;
    //double mu_k = gear->pressureAngle + theta_k; // need to define theta_k
    for (double r = pinion->baseDiameter / 2.0; r <= (pinion->outsideDiameter / 2); r+= stepsize_pinion) {
        double mu_k = compute_mu_k_from_r(r, pinion->baseDiameter / 2);
        pinion->position_x_left[i] = (pinion->baseDiameter) / 2.0 * (sin(mu_k) - mu_k * cos(mu_k));
        pinion->position_y_left[i] = (pinion->baseDiameter) / 2.0 * (cos(mu_k) + mu_k * sin(mu_k));
        pinion->position_x_right[i] = (-1. * pinion->position_x_left[i]) * cos(dynamicvars->trig_pinion) + pinion->position_y_left[i] * -1.0 * sin(dynamicvars->trig_pinion);
        pinion->position_y_right[i] = (-1. * pinion->position_x_left[i]) * sin(dynamicvars->trig_pinion) + pinion->position_y_left[i] * -1.0 * cos(dynamicvars->trig_pinion);
        i = i + 1;
    }
    i = 0;
    for (double r = gear->baseDiameter / 2.0; r <= (gear->outsideDiameter / 2); r+= stepsize_gear) {
        double mu_k = compute_mu_k_from_r(r, gear->baseDiameter / 2);
        gear->position_x_left[i] = (gear->baseDiameter) / 2.0 * (sin(mu_k) - mu_k * cos(mu_k));
        gear->position_y_left[i] = (gear->baseDiameter) / 2.0 * (cos(mu_k) + mu_k * sin(mu_k));
        gear->position_x_right[i] = (-1. * gear->position_x_left[i]) * cos(dynamicvars->trig_gear) + gear->position_y_left[i] * -1.0 * sin(dynamicvars->trig_gear);
        gear->position_y_right[i] = (-1. * gear->position_x_left[i]) * sin(dynamicvars->trig_gear) + gear->position_y_left[i] * -1.0 * cos(dynamicvars->trig_gear);
        i = i + 1;
    }
    return 1;
}

int dynamicsDiffEqSolver(SpurGear *pinion, SpurGear *gear, DynamicsVars *dynamicvars){
   // Section 3
    

}

int initialToothPinion(SpurGear* pinion, DynamicsVars* dynamicvars){
    // Eq (1)
    pinion.initialToothProfileXLeft.resize(pinion.numTeeth, vector<double>(dynamicVars.n, 0.0));
    pinion.initialToothProfileYLeft.resize(pinion.numTeeth, vector<double>(dynamicVars.n, 0.0));
    pinion.initialToothProfileXRight.resize(pinion.numTeeth, vector<double>(dynamicVars.n, 0.0));
    pinion.initialToothProfileYRight.resize(pinion.numTeeth, vector<double>(dynamicVars.n, 0.0));


    for (int i = 0; i < dynamicvars->n; i++) {
        pinion->initial_position_x_left[i] = cos(-1. * pinion->theta_r) * pinion->position_x_left[i] - sin(-1. * pinion->theta_r) * pinion->position_y_left[i];
        pinion->initial_position_y_left[i] = sin(-1. * pinion->theta_r) * pinion->position_x_left[i] + cos(-1. * pinion->theta_r) * pinion->position_y_left[i];
        pinion->initial_position_x_right[i] = cos(-1. * pinion->theta_r) * pinion->position_x_right[i] - sin(-1. * pinion->theta_r) * pinion->position_y_right[i];
        pinion->initial_position_y_right[i] = sin(-1. * pinion->theta_r) * pinion->position_x_right[i] + cos(-1. * pinion->theta_r) * pinion->position_y_right[i];
    }
}

int initialToothGear(SpurGear* gear, DynamicsVars* dynamicvars) {
    // Eq (2)
    gear.initialToothProfileXLeft.resize(gear.numTeeth, vector<double>(dynamicVars.n, 0.0));
    gear.initialToothProfileYLeft.resize(gear.numTeeth, vector<double>(dynamicVars.n, 0.0));
    gear.initialToothProfileXRight.resize(gear.numTeeth, vector<double>(dynamicVars.n, 0.0));
    gear.initialToothProfileYRight.resize(gear.numTeeth, vector<double>(dynamicVars.n, 0.0));


    for (int i = 0; i < dynamicvars->n; i++) {
        gear->initial_position_x_left[i] = cos(gear->theta_r) * gear->position_x_left[i] - sin(gear->theta_r) * gear->position_y_left[i];
        gear->initial_position_y_left[i] = sin(gear->theta_r) * gear->position_x_left[i] + cos(gear->theta_r) * gear->position_y_left[i];
        gear->initial_position_x_right[i] = cos(gear->theta_r) * gear->position_x_right[i] - sin(gear->theta_r) * gear->position_y_right[i];
        gear->initial_position_y_right[i] = sin(gear->theta_r) * gear->position_x_right[i] + cos(gear->theta_r) * gear->position_y_right[i];
    }
}

int completePinionProfile(SpurGear* pinion, DynamicsVars* dynamicvars) {
    // Eq (5)
   // pinion->toothProfileXLeft.resize(dynamicvars->n);
    //pinion->toothProfileYLeft.resize(dynamicvars->n);
    //pinion->toothProfileXRight.resize(dynamicvars->n);
    //pinion->toothProfileYRight.resize(dynamicvars->n);
    for (int i = 0; i < pinion->numTeeth; i++) {
        double z = 2 * i * pinion->theta_e;
        for (int j = 0; j < dynamicvars->n; j++) {
            pinion->initialToothProfileXLeft[i][j] = cos(z) * pinion->initial_position_x_left[j] - sin(z) * pinion->initial_position_y_left[j];
            pinion->initialToothProfileYLeft[i][j] = sin(z) * pinion->initial_position_x_left[j] + cos(z) * pinion->initial_position_y_left[j];
            pinion->initialToothProfileXRight[i][j] = cos(z) * pinion->initial_position_x_right[j] - sin(z) * pinion->initial_position_y_right[j];
            pinion->initialToothProfileYRight[i][j] = sin(z) * pinion->initial_position_x_right[j] + cos(z) * pinion->initial_position_y_right[j];
        }
    }
    return 1;
}

int completeGearProfile(SpurGear* gear, DynamicsVars* dynamicvars) {
    // Eq (6)
    //gear->toothProfileXLeft.resize(dynamicvars->n);
    //gear->toothProfileYLeft.resize(dynamicvars->n);
    //gear->toothProfileXRight.resize(dynamicvars->n);
    //gear->toothProfileYRight.resize(dynamicvars->n);
    for (int i = 0; i < gear->numTeeth; i++) {
        double z = 2.0 * i * gear->theta_e;
        for (int j = 0; j < dynamicvars->n; j++) {
            gear->initialToothProfileXLeft[i][j] = cos(z) * gear->initial_position_x_left[j] - sin(z) * gear->initial_position_y_left[j] + gear->centerDistance;
            gear->initialToothProfileYLeft[i][j] = sin(z) * gear->initial_position_x_left[j] + cos(z) * gear->initial_position_y_left[j];
            gear->initialToothProfileXRight[i][j] = cos(z) * gear->initial_position_x_right[j] - sin(z) * gear->initial_position_y_right[j] + gear->centerDistance;
            gear->initialToothProfileYRight[i][j] = sin(z) * gear->initial_position_x_right[j] + cos(z) * gear->initial_position_y_right[j];
        }
    }
    return 1;
}

double computeThetaA(double pressureAngle) {
    double pressureAngle_rad = pressureAngle * M_PI / 180.0;
    // Compute constant C = (1 - cos(pressureAngle)) / cos(pressureAngle)
    double C = (1.0 - cos(pressureAngle_rad)) / cos(pressureAngle_rad);

    // Initial guess for theta_a (in radians). For small angles, theta^3/3 ~ C.
    double theta = cbrt(3.0 * C);

    // Newton's method iteration
    const int maxIter = 100;
    const double tol = 1e-6;
    for (int i = 0; i < maxIter; i++) {
        double f = tan(theta) - theta - C;
        // f'(theta) = sec^2(theta) - 1 = 1/cos^2(theta) - 1
        double fp = (1.0 / (cos(theta) * cos(theta))) - 1.0;
        if (fabs(fp) < tol) break; // avoid division by zero
        double theta_new = theta - f / fp;
        if (fabs(theta_new - theta) < tol) {
            theta = theta_new;
            break;
        }
        theta = theta_new;
    }
    return theta;
}

void computeInvolutePoint(double r_base, double mu_k, double &x, double &y){
    // Compute the coordinates on the involute given parameter u and base circle radius r_base.
    // mu_k is the sum of theta_k and alpha_k.
    x = r_base * (sin(mu_k) - u * cos(mu_k));
    y = r_base * (cos(mu_k) + u * sin(mu_k));
    return x, y;
}

double compute_mu_k_from_r(double r, double r_base) {
    // Given a desired point on the involute (for example, determined by a specific radius r)
    // compute mu_k via r = r_base / cos(mu_k)
    return acos(r_base / r); // since cos(mu_k) = r_base / r, so u = arccos(r_base/r)
}

State computeDerivatives(const State& state, double t, const SpurGear* pinion, const SpurGear& gear) {

    State dState;

    // Derivative of positions are the velocities.
    for (int i = 0; i < 6; i++) {
        dState.pos[i] = state.vel[i];
    }

    // Compute the forces and moments for the pinion and gear.
    // These would be computed based on your model (contact forces, etc.).
    // For example, let’s assume you have computed:
    double F_pinion_sx, F_pinion_1x, F_pinion_2x, F_pinion_3x, F_pinion_4x;
    double F_pinion_sy, F_pinion_1y, F_pinion_2y, F_pinion_3y, F_pinion_4y;
    double T_pinion_1, T_pinion_2, T_pinion_3, T_pinion_4, T_d;
    double F_gear_sx, F_gear_1x, F_gear_2x, F_gear_3x, F_gear_4x;
    double F_gear_sy, F_gear_1y, F_gear_2y, F_gear_3y, F_gear_4y;
    double T_gear_1, T_gear_2, T_gear_3, T_gear_4, T_l;

    // ... (Compute these forces and moments from your contact model)
    double F_pinion_sx  = -1.0 * (pinion->elasticCoeff * state.pos[0] + pinion->dampingCoeff * state.vel[0]);
    double F_pinion_sy  = -1.0 * (pinion->elasticCoeff * state.pos[1] + pinion->dampingCoeff * state.vel[1]);
    double F_gear_sx    = -1.0 * (gear->elasticCoeff * (state.pos[4] - gear->centerDistance) + gear->dampingCoeff * state.vel[4]);
    double F_gear_sy    = -1.0 * (gear->elasticCoeff * state.pos[5] + gear->dampingCoeff * state.vel[5]);

    // Use the macros to compute accelerations:
    double x_double_dot_p = f1(F_pinion_sx, F_pinion_1x, F_pinion_2x, F_pinion_3x, F_pinion_4x) / pinion.mass;
    double y_double_dot_p = f2(F_pinion_sy, F_pinion_1y, F_pinion_2y, F_pinion_3y, F_pinion_4y) / pinion.mass;
    double phi_double_dot_p = f3(T_pinion_1, T_pinion_2, T_pinion_3, T_pinion_4, T_d) / pinion.momentInertia;

    double x_double_dot_g = f4(F_gear_sx, F_gear_1x, F_gear_2x, F_gear_3x, F_gear_4x) / gear.mass;
    double y_double_dot_g = f5(F_gear_sy, F_gear_1y, F_gear_2y, F_gear_3y, F_gear_4y) / gear.mass;
    double phi_double_dot_g = f6(T_gear_1, T_gear_2, T_gear_3, T_gear_4, T_l) / gear.momentInertia;

    // Fill in the derivative of the velocities:
    dState.vel[0] = x_double_dot_p;
    dState.vel[1] = y_double_dot_p;
    dState.vel[2] = phi_double_dot_p;
    dState.vel[3] = x_double_dot_g;
    dState.vel[4] = y_double_dot_g;
    dState.vel[5] = phi_double_dot_g;

    return dState;

}

State rungeKuttaStep(const State& state, double t, double dt, const SpurGear& pinion, const SpurGear& gear) {
    State k1 = computeDerivatives(state, t, pinion, gear);
    State state2;
    for (int i = 0; i < 6; i++) {
        state2.pos[i] = state.pos[i] + 0.5 * dt * k1.pos[i];
        state2.vel[i] = state.vel[i] + 0.5 * dt * k1.vel[i];
    }

    State k2 = computeDerivatives(state2, t + 0.5 * dt, pinion, gear);
    State state3;
    for (int i = 0; i < 6; i++) {
        state3.pos[i] = state.pos[i] + 0.5 * dt * k2.pos[i];
        state3.vel[i] = state.vel[i] + 0.5 * dt * k2.vel[i];
    }

    State k3 = computeDerivatives(state3, t + 0.5 * dt, pinion, gear);
    State state4;
    for (int i = 0; i < 6; i++) {
        state4.pos[i] = state.pos[i] + dt * k3.pos[i];
        state4.vel[i] = state.vel[i] + dt * k3.vel[i];
    }

    State k4 = computeDerivatives(state4, t + dt, pinion, gear);

    State newState;
    for (int i = 0; i < 6; i++) {
        newState.pos[i] = state.pos[i] + (dt / 6.0) * (k1.pos[i] + 2 * k2.pos[i] + 2 * k3.pos[i] + k4.pos[i]);
        newState.vel[i] = state.vel[i] + (dt / 6.0) * (k1.vel[i] + 2 * k2.vel[i] + 2 * k3.vel[i] + k4.vel[i]);
    }
    return newState;
}

void moveGearPosition(const State& state, SpurGear* pinion, SpurGear* gear, DynamicsVars* dynamicvars) {
    pinion.currentToothProfileXLeft.resize(pinion.numTeeth, vector<double>(dynamicVars.n, 0.0));
    pinion.currentToothProfileYLeft.resize(pinion.numTeeth, vector<double>(dynamicVars.n, 0.0));
    pinion.currentToothProfileXRight.resize(pinion.numTeeth, vector<double>(dynamicVars.n, 0.0));
    pinion.currentToothProfileYRight.resize(pinion.numTeeth, vector<double>(dynamicVars.n, 0.0));

    gear.currentToothProfileXLeft.resize(gear.numTeeth, vector<double>(dynamicVars.n, 0.0));
    gear.currentToothProfileYLeft.resize(gear.numTeeth, vector<double>(dynamicVars.n, 0.0));
    gear.currentToothProfileXRight.resize(gear.numTeeth, vector<double>(dynamicVars.n, 0.0));
    gear.currentToothProfileYRight.resize(gear.numTeeth, vector<double>(dynamicVars.n, 0.0));

    for (int i = 0; i < gear->numTeeth; i++) {
        for (int j = 0; j < dynamicvars->n; j++) {
            pinion->currentToothProfileXLeft[i][j] = cos(dynamicvars->phi_k) * pinion->initialToothProfileXLeft[i][j] - sin(dynamicvars->phi_k) * pinion->initialToothProfileYLeft[i][j] + state.pos[0];
            pinion->currentToothProfileYLeft[i][j] = sin(dynamicvars->phi_k) * pinion->initialToothProfileXLeft[i][j] + cos(dynamicvars->phi_k) * pinion->initialToothProfileYLeft[i][j] + state.pos[1];
            pinion->currentToothProfileXRight[i][j] = cos(dynamicvars->phi_k) * pinion->initialToothProfileXRight[i][j] - sin(dynamicvars->phi_k) * pinion->initialToothProfileYRight[i][j] + state.pos[0];
            pinion->currentToothProfileYRight[i][j] = sin(dynamicvars->phi_k) * pinion->initialToothProfileXRight[i][j] + cos(dynamicvars->phi_k) * pinion->initialToothProfileYRight[i][j] + state.pos[1];

            gear->currentToothProfileXLeft[i][j] = cos(dynamicvars->phi_k) * (gear->initialToothProfileXLeft[i][j] - gear->centerDistance) - sin(dynamicvars->phi_k) * gear->initialToothProfileYLeft[i][j] + state.pos[3];
            gear->currentToothProfileYLeft[i][j] = sin(dynamicvars->phi_k) * (gear->initialToothProfileYLeft[i][j] - gear->centerDistance) + cos(dynamicvars->phi_k) * gear->initialToothProfileYLeft[i][j] + state.pos[4];
            gear->currentToothProfileXRight[i][j] = cos(dynamicvars->phi_k) * (gear->initialToothProfileXRight[i][j] - gear->centerDistance) - sin(dynamicvars->phi_k) * gear->initialToothProfileYRight[i][j] + state.pos[3];
            gear->currentToothProfileYRight[i][j] = sin(dynamicvars->phi_k) * (gear->initialToothProfileYRight[i][j] - gear->centerDistance) + cos(dynamicvars->phi_k) * gear->initialToothProfileYRight[i][j] + state.pos[4];
        }
    }
}

void meshingLine(SpurGear* pinion, SpurGear* gear, DynamicsVars* dynamicvars) {
    double xcenterpinion = 0;
    double ycenterpinion = 0;

    double xcentergear = gear->centerDistance;
    double ycentergear = 0;

    // Theta_pk is the angle between the line connecting the gear and pinion's rotating axis and the x-axis.
    double theta_pk = atan2(ycentergear - ycenterpinion, xcentergear - xcenterpinion);
    
    double transmissionratio = pinion->gearRatio;
    double a_k_Pinion = sqrt(pow(xcentergear - xcenterpinion, 2) + pow(ycentergear - ycenterpinion, 2));
    double pitchCircleRadiusPinion = a_k_Pinion / (1 + transmissionratio);

    // Eq (16)
    double pitchPointX, pitchPointY;
    if (ycentergear >= 0) {
        pitchPointX = pitchCircleRadiusPinion * cos(theta_pk) + xcenterpinion;
        pitchPointY = pitchCircleRadiusPinion * sin(theta_pk) + ycenterpinion;
    }
    else {
        pitchPointX = pitchCircleRadiusPinion * cos(theta_pk) + xcenterpinion;
        pitchPointY = -1.0 * pitchCircleRadiusPinion * sin(theta_pk) + ycenterpinion;
    }
    double k1 = computeInvoluteSlope((pinion->baseDiameter / 2.0), (pinion->pitchDiameter / 2.0));
    double k2 = computeInvoluteSlope((gear->baseDiameter / 2.0), (gear->pitchDiameter / 2.0));

    // Line of mesh point - Eq (17) and Eq (18)
    // y = k1 * (x - pitchPointX) + pitchPointY;
    // y = k2 * (x - pitchPointX) + pitchPointY;

    double u_mesh1 = solveForU((pinion->baseDiameter / 2.0), k1, pitchPointX, pitchPointY, 0.1);
    double x_mesh1, y_mesh1;
    computeInvolutePoint((pinion->baseDiameter / 2.0), u_mesh1, x_mesh1, y_mesh1);

    double u_mesh2 = solveForU((pinion->baseDiameter / 2.0), k1, pitchPointX, pitchPointY, 0.1);
    double x_mesh2, y_mesh2;
    computeInvolutePoint((pinion->baseDiameter / 2.0), u_mesh2, x_mesh2, y_mesh2);

    cout << "Mesh point coordinates: (" << x_mesh << ", " << y_mesh << ")" << endl;

    // Calculate distance between two candidate mesh points.
    double delta = sqrt(pow(x_mesh1 - x_mesh2, 2) + pow(y_mesh1 - y_mesh2, 2));
    if (delta <= dynamicvars->deltaJudgmentThreshold) {
        ContactPoint cp;
        cp.x = x_mesh1;
        cp.y = y_mesh1;
        pinion->contactPoints.push_back(cp);
        gear->contactPoints.push_back(cp);
        cout << "Single meshing point detected: (" << cp.x << ", " << cp.y << ")" << endl;
    }
    else {
        // The candidate points are distinct -> there is two separate contact pairs
        ContactPoint cp1, cp2;
        cp1.x = x_mesh1;
        cp1.y = y_mesh1;
        cp2.x = x_mesh2;
        cp2.y = y_mesh2;
        pinion->contactPoints.push_back(cp1);
        pinion->contactPoints.push_back(cp2);
        gear->contactPoints.push_back(cp1);
        gear->contactPoints.push_back(cp2);
        cout << "Two distinct meshing points detected:" << endl;
        cout << "Point 1: (" << cp1.x << ", " << cp1.y << ")" << endl;
        cout << "Point 2: (" << cp2.x << ", " << cp2.y << ")" << endl;
    }
}

// Helper functions for meshingLine();

double computeInvoluteSlope(double baseRadius, double pitchRadius) {
    // Ensure pitchRadius is greater than baseRadius to avoid invalid input.
    if (pitchRadius <= baseRadius) {
        cout << "Error: pitchRadius must be greater than baseRadius." << endl;
        return 0.0;
    }
    double u = acos(baseRadius / pitchRadius);
    // Compute slope: k = cot(u) = cos(u) / sin(u)
    double sin_u = sin(u);
    if (fabs(sin_u) < 1e-8) { // Avoid division by zero
        cout << "Error: sin(acos(baseRadius / pitchRadius)) is < 1e-8" << endl;
        return 0.0;
    }
    return (cos(u) / sin_u);
}

void computeInvolutePoint(double r_b, double u, double& x, double& y) {
    x = r_b * (sin(u) - u * cos(u));
    y = r_b * (cos(u) + u * sin(u));
}

double F(double u, double r_b, double k, double x_J, double y_J) {
    double x_u, y_u;
    computeInvolutePoint(r_b, u, x_u, y_u);
    return y_u - (k * (x_u - x_J) + y_J);
}

double dFdu(double u, double r_b, double k, double x_J, double y_J, double h = 1e-6) {
    // numerical derivative for F(u)
    return (F(u + h, r_b, k, x_J, y_J) - F(u - h, r_b, k, x_J, y_J)) / (2 * h);
}

double solveForU(double r_b, double k, double x_J, double y_J, double initialGuess = 0.1, int maxIter = 100, double tol = 1e-6) {
    // Newton's method to solve F(u) = 0
    double u = initialGuess;
    for (int i = 0; i < maxIter; i++) {
        double fVal = F(u, r_b, k, x_J, y_J);
        double fDeriv = dFdu(u, r_b, k, x_J, y_J);
        if (fabs(fDeriv) < tol) {  // Avoid division by near-zero.
            cout << "Derivative too small in Newton's method." << endl;
            break;
        }
        double u_next = u - fVal / fDeriv;
        if (fabs(u_next - u) < tol)
            return u_next;
        u = u_next;
    }
    return u;
}

// End of helper functions

void computeContactVectors(const Vec2& p1, const Vec2& p2, const SpurGear& gear, Vec2& n_g, Vec2& n1, Vec2& n1_t, Vec2& n2, Vec2& n2_t) {
    // p1: contact point on pinion, p2: contact point on gear
    // Eq (20): unit normal vector in the deformation direction.
    n_g = computeUnitVector(p1, p2);

    // Eq (21): unit vector from pinion center (0,0) to contact point p1.
    Vec2 O1 = { 0.0, 0.0 };
    n1 = computeUnitVector(p1, O1);

    // Eq (22): unit vector from gear center (centerDistance, 0) to contact point p2.
    Vec2 O2 = { gear.centerDistance, 0.0 };
    n2 = computeUnitVector(p2, O2);

    // Eq (23) and Eq (24): rotate the above unit vectors by 90 degrees.
    n1_t = rotate90(n1);
    n2_t = rotate90(n2);

    // debugging output
    cout << "Deformation unit vector n_g = (" << n_g.x << ", " << n_g.y << ")\n";
    cout << "Pinion unit vector n1 = (" << n1.x << ", " << n1.y << ")\n";
    cout << "Pinion tangential unit vector n1_t = (" << n1_t.x << ", " << n1_t.y << ")\n";
    cout << "Gear unit vector n2 = (" << n2.x << ", " << n2.y << ")\n";
    cout << "Gear tangential unit vector n2_t = (" << n2_t.x << ", " << n2_t.y << ")\n";
}

// Helper functions for computerContactVectors()

Vec2D computeUnitVector(const Vec2& p1, const Vec2& p2) {
    Vec2 diff;
    diff.x = p1.x - p2.x;
    diff.y = p1.y - p2.y;
    double norm = sqrt(diff.x * diff.x + diff.y * diff.y);
    return { diff.x / norm, diff.y / norm };
}

Vec2D rotate90(const Vec2& v) {
    // Rotation 90° CCW
    return { -v.y, v.x };
}

// End of helper functions

Vec2D computeContactForce(double delta, double delta_dot, double kv, double cd, double mu, double fF, const Vec2D& v_t, const Vec2D& n, const Vec2D& n_t) {
    // Computes the net contact force at a mesh point
    // Inputs:
    //   delta    : deformation
    //   delta_dot: time derivative of deformation
    //   kv       : meshing stiffness
    //   cd       : damping coefficient
    //   mu       : friction coefficient
    //   fF       : friction parameter
    //   v_t      : relative tangential velocity vector
    //   n        : unit normal (deformation) vector
    //   n_t      : unit tangential vector
    // Returns:
    //   A 2D vector representing the net contact force
    
    double v0 = 0.1; // m/s, below which friction is not fully activated
    double vL = 0.5; // m/s, above which full sliding friction is assumed

    // Compute the normal contact force
    double F_cg = kv * delta + cd * delta_dot;

    // Compute the friction force
    Vec2D F_friction = { 0.0, 0.0 };
    double vt_mag = sqrt(v_t.x * v_t.x + v_t.y * v_t.y);
    double fF = computeFrictionParameter(vt_mag, v0, vL);
    if (vt_mag > 1e-6) { // Avoid divide by zero.
        // Friction force 
        F_friction.x = -fF * mu * F_cg * (v_t.x / vt_mag);
        F_friction.y = -fF * mu * F_cg * (v_t.y / vt_mag);
    }
    // Compute the net contact force
    Vec2D F;
    F.x = F_cg * n.x + F_friction.x * n_t.x;
    F.y = F_cg * n.y + F_friction.y * n_t.y;
    return F;
}

Vec2D normalize(const Vector2D& v) {
    double norm = sqrt(v.x * v.x + v.y * v.y);
    if (norm < 1e-8) {
        cout << "Error: norm too small in normalization." << endl;
        return { 0.0, 0.0 };
    }
    return { v.x / norm, v.y / norm };
}

double computeFrictionParameter(double vtMagnitude, double v0, double vL) {
    // Eq (34)
    if (vtMagnitude <= v0)
        return 0.0;
    else if (vtMagnitude <= vL)
        return (vtMagnitude - v0) / (vL - v0);
    else
        return 1.0;
}

double computeMoment(const double contactX, const double contactY, const double centerX, const double centerY, const double forceX, const double forceY) {
    // Eq (36) and Eq (37)
    return (contactY - centerY) * forceX - (contactX - centerX) * forceY;
}