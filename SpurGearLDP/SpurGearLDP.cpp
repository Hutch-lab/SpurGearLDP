/* Physics 5810 - Computational Physics Final Project - Nick Hutchison */
/*  Spur Gear Load Distribution Program */
#define _USE_MATH_DEFINES
#include <cmath>
#include "FunctDefs2.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <fstream>
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
#define f2(F_pinion_sy, F_pinion_1y, F_pinion_2y, F_pinion_3y, F_pinion_4y) (F_pinion_sy - F_pinion_1y - F_pinion_2y - F_pinion_3y - F_pinion_4y) - pinion->mass * g
#define f3(T_pinion_1, T_pinion_2, T_pinion_3, T_pinion_4, T_d) (T_pinion_1 + T_pinion_2 + T_pinion_3 + T_pinion_4 + T_d)
#define f4(F_gear_sx, F_gear_1x, F_gear_2x, F_gear_3x, F_gear_4x) (F_gear_sx - F_gear_1x - F_gear_2x - F_gear_3x - F_gear_4x)
#define f5(F_gear_sy, F_gear_1y, F_gear_2y, F_gear_3y, F_gear_4y) (F_gear_sy - F_gear_1y - F_gear_2y - F_gear_3y - F_gear_4y - gear->mass * g)
#define f6(T_gear_1, T_gear_2, T_gear_3, T_gear_4, T_l) (T_gear_1 + T_gear_2 + T_gear_3 + T_gear_4 + T_l)

int main() {
    SpurGear pinion;
    SpurGear gear;
    DynamicsVars dynvars;
    dynvars.deltaJudgmentThreshold = 1e-5;
    startupmsg();
    cout << endl;

    // Get input method
    int method;
    cout << " Enter 0 to manually enter parameters or 1 to load parameters from file" << endl;
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
    outputGearProfile(&pinion, &dynvars, "pinionProfile.csv");
    outputGearProfile(&gear, &dynvars, "gearProfile.csv");

    // Compute meshing line and contact points
    meshingLine(&pinion, &gear, &dynvars);
    /*
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

        moveGearPosition(state, &pinion, &gear, &dynvars); // Update the current tooth profiles based on the new state.

        // Output some state information.
        cout << "Time: " << t
            << " | Pinion pos: (" << state.pos[0] << ", " << state.pos[1] << ")"
            << " | Gear pos: (" << state.pos[3] << ", " << state.pos[4] << ")" << endl;
    }

    cout << "\nSimulation complete." << endl;
    */
    return 0;
}


void startupmsg() {
    cout << " Physics 5810 - Computational Physics Final Project " << endl;
    cout << " ----------------- Nick Hutchison ----------------- " << endl;
    cout << " ------ Spur Gear Load Distribution Program ------- " << endl;
}

int getUserInput(SpurGear* pinion, SpurGear* gear) {
    /*
    Returns 0 when the user wants to exit and 1 if the user enters valid input data.

    INPUT: pinion gear struct, spur gear struct
    OUTPUT: 1 if no error, 0 otherwise
    CALLS: none
     */

    /*cout << "Enter number of teeth for the pinion and gear: ";
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
    cin >> pinion->poissonsRatio >> gear->poissonsRatio;*/

    pinion->numTeeth = 30.0;
    gear->numTeeth = 50.0;
    pinion->pressureAngle = 20.0;
    gear->pressureAngle = pinion->pressureAngle;
    pinion->rootDiameter = 40.0;
    gear->rootDiameter = 70.0;
    pinion->centerDistance = 0.0;
    pinion->faceWidth = 20.0;
    gear->faceWidth = 20.0;
    gear->offsetFaceWidth = 0.0;
    pinion->module = 2.5;
    gear->module = pinion->module;
   
    pinion->gearRatio = gear->numTeeth / pinion->numTeeth;
    gear->gearRatio = pinion->gearRatio;
    pinion->pitchDiameter = pinion->numTeeth * pinion->module;
    gear->pitchDiameter = gear->numTeeth * gear->module;
    pinion->baseDiameter = pinion->pitchDiameter * cos(M_PI / 180.0 * pinion->pressureAngle);
    gear->baseDiameter = gear->pitchDiameter * cos(M_PI / 180.0 * gear->pressureAngle);
    pinion->outsideDiameter = pinion->module * 2.0 + pinion->pitchDiameter;
    gear->outsideDiameter = gear->module * 2.0 + gear->pitchDiameter;
    pinion->pitchRadius = pinion->pitchDiameter / 2.0;
    pinion->baseRadius = pinion->baseDiameter / 2.0;
    pinion->outsideRadius = pinion->outsideDiameter / 2.0;
    gear->pitchRadius = gear->pitchDiameter / 2.0;
    gear->baseRadius = gear->baseDiameter / 2.0;
    gear->outsideRadius = gear->outsideDiameter / 2.0;

    if (pinion->centerDistance == 0.0) {
        // Calculate the nominal center distance between gears
        double cd = (gear->pitchDiameter + pinion->pitchDiameter) / 2.0;
        gear->centerDistance = cd;
        pinion->centerDistance = cd;
    }

    cout << "pinion->centerDistance: " << pinion->centerDistance << endl;
    cout << "gear->centerDistance: " << gear->centerDistance << endl;
    cout << "pinion->gearRatio: " << pinion->gearRatio << endl;
    cout << "pinion->pitchDiameter: " << pinion->pitchDiameter << endl;
    cout << "gear->pitchDiameter: " << gear->pitchDiameter << endl;
    cout << "pinion->baseDiameter: " << pinion->baseDiameter << endl;
    cout << "gear->baseDiameter: " << gear->baseDiameter << endl;
    cout << "pinion->outsideDiameter: " << pinion->outsideDiameter << endl;
    cout << "gear->outsideDiameter: " << gear->outsideDiameter << endl;

    cout << "parameters loaded!" << endl;
    return 1;
}

int setupProfiles(SpurGear* pinion, SpurGear* gear, DynamicsVars* dynamicvars) {
    dynamicvars->n = 1000;
    pinion->theta_e = M_PI / pinion->numTeeth;
    gear->theta_e   = M_PI / gear->numTeeth;
    pinion->theta_a = computeThetaA(pinion->pressureAngle);
    gear->theta_a   = computeThetaA(gear->pressureAngle);

    dynamicvars->trig_pinion = -pinion->theta_e - 2.0 * pinion->theta_a;
    dynamicvars->trig_gear   = -gear->theta_e   - 2.0 * gear->theta_a;
    pinion->theta_r = M_PI_2 - pinion->theta_a - pinion->theta_e / 2.0;
    gear->theta_r   = M_PI_2 + gear->theta_a   - gear->theta_e   / 2.0;
    double stepsize_pinion = (pinion->outsideRadius - pinion->baseRadius) / dynamicvars->n;
    double stepsize_gear = (gear->outsideRadius - gear->baseRadius) / dynamicvars->n;
    cout << "pinion->theta_e: " << pinion->theta_e << endl;
    cout << "gear->theta_e: " << gear->theta_e << endl;
    cout << "pinion->theta_a: " << pinion->theta_a << endl;
    cout << "gear->theta_a: " << gear->theta_a << endl;
    cout << "dynamicvars->trig_pinion: " << dynamicvars->trig_pinion << endl;
    cout << "dynamicvars->trig_gear: " << dynamicvars->trig_gear << endl;
    cout << "pinion->theta_r: " << pinion->theta_r << endl;
    cout << "gear->theta_r: " << gear->theta_r << endl;
    pinion->position_x_left.resize(dynamicvars->n);
    pinion->position_y_left.resize(dynamicvars->n);
    pinion->position_x_right.resize(dynamicvars->n);
    pinion->position_y_right.resize(dynamicvars->n);

    gear->position_x_left.resize(dynamicvars->n);
    gear->position_y_left.resize(dynamicvars->n);
    gear->position_x_right.resize(dynamicvars->n);
    gear->position_y_right.resize(dynamicvars->n);

    
    double mu_max = sqrt((pinion->outsideRadius * pinion->outsideRadius) / (pinion->baseRadius * pinion->baseRadius) - 1.0);
    double mu_step = mu_max / dynamicvars->n;
    int i = 0;
    for (double mu = 0; mu <= mu_max && i < (dynamicvars->n); mu += mu_step, i++) {
        double x = pinion->baseRadius * (sin(mu) - mu * cos(mu));
        double y = pinion->baseRadius * (cos(mu) + mu * sin(mu));
        // For the left side we assume these are the coordinates.
        pinion->position_x_left[i] = x;
        pinion->position_y_left[i] = y;
        // mirror reflection about the y-axis for the right profile
        //pinion->position_x_right[i] = -x;
        //pinion->position_y_right[i] = y;
        pinion->position_x_right[i] = (-1.0 * pinion->position_x_left[i]) * cos(dynamicvars->trig_pinion) - pinion->position_y_left[i] * sin(dynamicvars->trig_pinion);
        pinion->position_y_right[i] = (-1.0 * pinion->position_x_left[i]) * sin(dynamicvars->trig_pinion) + pinion->position_y_left[i] * cos(dynamicvars->trig_pinion);
    }

    mu_max = sqrt((gear->outsideRadius * gear->outsideRadius) / (gear->baseRadius * gear->baseRadius) - 1.0);
    mu_step = mu_max / dynamicvars->n;
    i = 0;
    for (double mu = 0; mu <= mu_max && i < (dynamicvars->n); mu += mu_step, i++) {
        double x = gear->baseRadius * (sin(mu) - mu * cos(mu));
        double y = gear->baseRadius * (cos(mu) + mu * sin(mu));
        // For the left side we assume these are the coordinates.
        gear->position_x_left[i] = x;
        gear->position_y_left[i] = y;
        // mirror reflection about the y-axis for the right profile
        //gear->position_x_right[i] = -x;
        //gear->position_y_right[i] = y;
        gear->position_x_right[i] = (-1.0 * gear->position_x_left[i]) * cos(dynamicvars->trig_gear) - gear->position_y_left[i] * sin(dynamicvars->trig_gear);
        gear->position_y_right[i] = (-1.0 * gear->position_x_left[i]) * sin(dynamicvars->trig_gear) + gear->position_y_left[i] * cos(dynamicvars->trig_gear);
    }


    //for (double r = pinion->baseRadius; r <= pinion->outsideRadius; r += stepsize_pinion) {
        /*double mu_k = acos(pinion->baseRadius/r);
        if (mu_k < 0) {
            cout << "mu_k less than 0: " << mu_k << endl;
        }
        pinion->position_x_left[i]  = pinion->baseRadius * (sin(mu_k) - mu_k * cos(mu_k));
        pinion->position_y_left[i]  = pinion->baseRadius * (cos(mu_k) + mu_k * sin(mu_k));
        pinion->position_x_right[i] = (-1.0 * pinion->position_x_left[i]) * cos(dynamicvars->trig_pinion) - pinion->position_y_left[i] * sin(dynamicvars->trig_pinion);
        pinion->position_y_right[i] = (-1.0 * pinion->position_x_left[i]) * sin(dynamicvars->trig_pinion) + pinion->position_y_left[i] * cos(dynamicvars->trig_pinion);
        i = i + 1;
        */
    //}
    //for (double r = gear->baseRadius; r <= gear->outsideRadius; r += stepsize_gear) {
    //    //if (i == dynamicvars->n) {
    //        //break;
    //    //}
    //    double mu_k = -acos(gear->baseRadius/r);
    //    gear->position_x_left[i] = gear->baseRadius * (sin(mu_k) - mu_k * cos(mu_k));
    //    gear->position_y_left[i] = gear->baseRadius * (cos(mu_k) + mu_k * sin(mu_k));
    //    gear->position_x_right[i] = (-1.0 * gear->position_x_left[i]) * cos(dynamicvars->trig_gear) - gear->position_y_left[i] * sin(dynamicvars->trig_gear);
    //    gear->position_y_right[i] = (-1.0 * gear->position_x_left[i]) * sin(dynamicvars->trig_gear) + gear->position_y_left[i] * cos(dynamicvars->trig_gear);
    //}
    return 1;
}

int initialToothPinion(SpurGear* pinion, DynamicsVars* dynamicvars) {
    // Resize the dynamic arrays
    pinion->initial_position_x_left.resize(dynamicvars->n);
    pinion->initial_position_y_left.resize(dynamicvars->n);
    pinion->initial_position_x_right.resize(dynamicvars->n);
    pinion->initial_position_y_right.resize(dynamicvars->n);

    pinion->initialToothProfileXLeft.resize(pinion->numTeeth, vector<double>(dynamicvars->n, 0.0));
    pinion->initialToothProfileYLeft.resize(pinion->numTeeth, vector<double>(dynamicvars->n, 0.0));
    pinion->initialToothProfileXRight.resize(pinion->numTeeth, vector<double>(dynamicvars->n, 0.0));
    pinion->initialToothProfileYRight.resize(pinion->numTeeth, vector<double>(dynamicvars->n, 0.0));
    
    // Eq (1)
    for (int i = 0; i < dynamicvars->n; i++) {
        pinion->initial_position_x_left[i]  = cos(-1.0 * pinion->theta_r) * pinion->position_x_left[i]  - sin(-1.0 * pinion->theta_r) * pinion->position_y_left[i];
        pinion->initial_position_y_left[i]  = sin(-1.0 * pinion->theta_r) * pinion->position_x_left[i]  + cos(-1.0 * pinion->theta_r) * pinion->position_y_left[i];
        pinion->initial_position_x_right[i] = cos(-1.0 * pinion->theta_r) * pinion->position_x_right[i] - sin(-1.0 * pinion->theta_r) * pinion->position_y_right[i];
        pinion->initial_position_y_right[i] = sin(-1.0 * pinion->theta_r) * pinion->position_x_right[i] + cos(-1.0 * pinion->theta_r) * pinion->position_y_right[i];
    }
    return 1;
}

int initialToothGear(SpurGear* gear, DynamicsVars* dynamicvars) {
    // Resize the dynamic arrays
    gear->initial_position_x_left.resize(dynamicvars->n);
    gear->initial_position_y_left.resize(dynamicvars->n);
    gear->initial_position_x_right.resize(dynamicvars->n);
    gear->initial_position_y_right.resize(dynamicvars->n);

    gear->initialToothProfileXLeft.resize(gear->numTeeth, vector<double>(dynamicvars->n, 0.0));
    gear->initialToothProfileYLeft.resize(gear->numTeeth, vector<double>(dynamicvars->n, 0.0));
    gear->initialToothProfileXRight.resize(gear->numTeeth, vector<double>(dynamicvars->n, 0.0));
    gear->initialToothProfileYRight.resize(gear->numTeeth, vector<double>(dynamicvars->n, 0.0));

    // Eq (2)
    for (int i = 0; i < dynamicvars->n; i++) {
        gear->initial_position_x_right[i]  = cos(gear->theta_r) * gear->position_x_left[i]  - sin(gear->theta_r) * gear->position_y_left[i];
        gear->initial_position_y_right[i]  = sin(gear->theta_r) * gear->position_x_left[i]  + cos(gear->theta_r) * gear->position_y_left[i];
        gear->initial_position_x_left[i] = cos(gear->theta_r) * gear->position_x_right[i] - sin(gear->theta_r) * gear->position_y_right[i];
        gear->initial_position_y_left[i] = sin(gear->theta_r) * gear->position_x_right[i] + cos(gear->theta_r) * gear->position_y_right[i];
    }
    return 1;
}

int completePinionProfile(SpurGear* pinion, DynamicsVars* dynamicvars) {
    // Eq (5)
    // Apply rotation matrix to copy tooth profile for every tooth location
    for (int i = 1; i <= pinion->numTeeth; i++) {
        double z = 2.0 * (i - 1) * pinion->theta_e;
        for (int j = 0; j < dynamicvars->n; j++) {
            pinion->initialToothProfileXLeft[i - 1][j]  = cos(z) * pinion->initial_position_x_left[j]  - sin(z) * pinion->initial_position_y_left[j];
            pinion->initialToothProfileYLeft[i - 1][j]  = sin(z) * pinion->initial_position_x_left[j]  + cos(z) * pinion->initial_position_y_left[j];
            pinion->initialToothProfileXRight[i - 1][j] = cos(z) * pinion->initial_position_x_right[j] - sin(z) * pinion->initial_position_y_right[j];
            pinion->initialToothProfileYRight[i - 1][j] = sin(z) * pinion->initial_position_x_right[j] + cos(z) * pinion->initial_position_y_right[j];
        }
    }
    return 1;
}

int completeGearProfile(SpurGear* gear, DynamicsVars* dynamicvars) {
    // Eq (6)
    // Apply rotation matrix to copy tooth profile for every tooth location
    for (int i = 1; i <= gear->numTeeth; i++) {
        double z = 2.0 * (i-1) * gear->theta_e;
        for (int j = 0; j < dynamicvars->n; j++) {
            // Offset the x values by the center distance b/c the pinion is centered at (0,0)
            gear->initialToothProfileXLeft[i-1][j] = cos(z) * gear->initial_position_x_left[j] - sin(z) * gear->initial_position_y_left[j] + gear->centerDistance;
            gear->initialToothProfileYLeft[i-1][j] = sin(z) * gear->initial_position_x_left[j] + cos(z) * gear->initial_position_y_left[j];
            gear->initialToothProfileXRight[i-1][j] = cos(z) * gear->initial_position_x_right[j] - sin(z) * gear->initial_position_y_right[j] + gear->centerDistance;
            gear->initialToothProfileYRight[i-1][j] = sin(z) * gear->initial_position_x_right[j] + cos(z) * gear->initial_position_y_right[j];
        }
    }
    return 1;
}

void outputGearProfile(const SpurGear* gear, const DynamicsVars* dynVars, const string& filename) {
    ofstream outFile(filename);
    outFile << "Tooth,ProfilePoint,x_left,y_left,x_right,y_right\n"; // Write file header

    for (int tooth = 0; tooth < gear->numTeeth; tooth++) { // Loop over each tooth (column) and each profile point (row)
        for (int i = 0; i < dynVars->n; i++) {
            // Assuming your complete gear profile is stored in initialToothProfileXLeft, etc.
            outFile << tooth << "," << i << ","
                << gear->initialToothProfileXLeft[tooth][i] << ","
                << gear->initialToothProfileYLeft[tooth][i] << ","
                << gear->initialToothProfileXRight[tooth][i] << ","
                << gear->initialToothProfileYRight[tooth][i] << "\n";
        }
    }
    outFile.close();
    cout << "Gear profile data written to " << filename << endl;
}

double computeThetaA(double pressureAngle) {
    double alpha_rad = pressureAngle * M_PI / 180.0; // Convert the pressure angle from degrees to radians.
    return tan(alpha_rad) - alpha_rad; // Return the involute function for alpha.
}

State rungeKuttaStep(State& state, double t, double dt, SpurGear& pinion, SpurGear& gear, DynamicsVars* dynVars) {
    State k1 = computeDerivatives(state, t, &pinion, &gear, dynVars);
    State state2;
    for (int i = 0; i < 6; i++) {
        state2.pos[i] = state.pos[i] + 0.5 * dt * k1.pos[i];
        state2.vel[i] = state.vel[i] + 0.5 * dt * k1.vel[i];
    }

    State k2 = computeDerivatives(state2, t + 0.5 * dt, &pinion, &gear, dynVars);
    State state3;
    for (int i = 0; i < 6; i++) {
        state3.pos[i] = state.pos[i] + 0.5 * dt * k2.pos[i];
        state3.vel[i] = state.vel[i] + 0.5 * dt * k2.vel[i];
    }

    State k3 = computeDerivatives(state3, t + 0.5 * dt, &pinion, &gear, dynVars);
    State state4;
    for (int i = 0; i < 6; i++) {
        state4.pos[i] = state.pos[i] + dt * k3.pos[i];
        state4.vel[i] = state.vel[i] + dt * k3.vel[i];
    }

    State k4 = computeDerivatives(state4, t + dt, &pinion, &gear, dynVars);

    State newState;
    for (int i = 0; i < 6; i++) {
        newState.pos[i] = state.pos[i] + (dt / 6.0) * (k1.pos[i] + 2 * k2.pos[i] + 2 * k3.pos[i] + k4.pos[i]);
        newState.vel[i] = state.vel[i] + (dt / 6.0) * (k1.vel[i] + 2 * k2.vel[i] + 2 * k3.vel[i] + k4.vel[i]);
    }
    return newState;
}

State computeDerivatives(State& state, double t, SpurGear* pinion, SpurGear* gear, DynamicsVars* dynVars) {

    State dState;
    double g = 9.80665; // gravitational acceleration

    // Derivative of positions are the velocities.
    for (int i = 0; i < 6; i++) {
        dState.pos[i] = state.vel[i];
    }
    double F_pinion_1x = 0, F_pinion_2x = 0, F_pinion_3x = 0, F_pinion_4x = 0;
    double F_pinion_1y = 0, F_pinion_2y = 0, F_pinion_3y = 0, F_pinion_4y = 0;
    double T_pinion_1 = 0, T_pinion_2 = 0, T_pinion_3 = 0, T_pinion_4 = 0, T_d = 0;
    double F_gear_1x = 0, F_gear_2x = 0, F_gear_3x = 0, F_gear_4x = 0;
    double F_gear_1y = 0, F_gear_2y = 0, F_gear_3y = 0, F_gear_4y = 0;
    double T_gear_1 = 0, T_gear_2 = 0, T_gear_3 = 0, T_gear_4 = 0, T_l = 0;
    double F_pinion_sx = solveF_Pinion_sx(pinion, state);
    double F_pinion_sy = solveF_Pinion_sy(pinion, state);
    double F_gear_sx = solveF_Gear_sx(gear, state);
    double F_gear_sy = solveF_Gear_sy(gear, state);

    // Use the macros to compute accelerations:
    double x_double_dot_p = f1(F_pinion_sx, F_pinion_1x, F_pinion_2x, F_pinion_3x, F_pinion_4x) / pinion->mass;
    double y_double_dot_p = f2(F_pinion_sy, F_pinion_1y, F_pinion_2y, F_pinion_3y, F_pinion_4y) / pinion->mass;
    double phi_double_dot_p = f3(T_pinion_1, T_pinion_2, T_pinion_3, T_pinion_4, T_d) / pinion->momentInertia;

    double x_double_dot_g = f4(F_gear_sx, F_gear_1x, F_gear_2x, F_gear_3x, F_gear_4x) / gear->mass;
    double y_double_dot_g = f5(F_gear_sy, F_gear_1y, F_gear_2y, F_gear_3y, F_gear_4y) / gear->mass;
    double phi_double_dot_g = f6(T_gear_1, T_gear_2, T_gear_3, T_gear_4, T_l) / gear->momentInertia;

    // Fill in the derivative of the velocities:
    dState.vel[0] = x_double_dot_p;
    dState.vel[1] = y_double_dot_p;
    dState.vel[2] = phi_double_dot_p;
    dState.vel[3] = x_double_dot_g;
    dState.vel[4] = y_double_dot_g;
    dState.vel[5] = phi_double_dot_g;

    return dState;
}

// ---------------------------------------------- //
// -- Solving for Resultant Forces and Moments -- //
// ---------------------------------------------- //

double solveF_Pinion_sx(SpurGear* pinion, State& state) {
    double F_pinion_sx = -1.0 * (pinion->elasticCoeff * state.pos[0] + pinion->dampingCoeff * state.vel[0]);
    return F_pinion_sx;
}
double solveF_Pinion_sy(SpurGear* pinion, State& state) {
    double F_pinion_sy = -1.0 * (pinion->elasticCoeff * state.pos[1] + pinion->dampingCoeff * state.vel[1]);
    return F_pinion_sy;
}
double solveF_Gear_sx(SpurGear* gear, State& state) {
    double F_gear_sx = -1.0 * (gear->elasticCoeff * (state.pos[4] - gear->centerDistance) + gear->dampingCoeff * state.vel[4]);
    return F_gear_sx;
}
double solveF_Gear_sy(SpurGear* gear, State& state) {
    double F_gear_sy = -1.0 * (gear->elasticCoeff * state.pos[5] + gear->dampingCoeff * state.vel[5]);
    return F_gear_sy;
}



void moveGearPosition(const State& state, SpurGear* pinion, SpurGear* gear, DynamicsVars* dynamicvars) {
    pinion->currentToothProfileXLeft.resize(pinion->numTeeth, vector<double>(dynamicvars->n, 0.0));
    pinion->currentToothProfileYLeft.resize(pinion->numTeeth, vector<double>(dynamicvars->n, 0.0));
    pinion->currentToothProfileXRight.resize(pinion->numTeeth, vector<double>(dynamicvars->n, 0.0));
    pinion->currentToothProfileYRight.resize(pinion->numTeeth, vector<double>(dynamicvars->n, 0.0));

    gear->currentToothProfileXLeft.resize(gear->numTeeth, vector<double>(dynamicvars->n, 0.0));
    gear->currentToothProfileYLeft.resize(gear->numTeeth, vector<double>(dynamicvars->n, 0.0));
    gear->currentToothProfileXRight.resize(gear->numTeeth, vector<double>(dynamicvars->n, 0.0));
    gear->currentToothProfileYRight.resize(gear->numTeeth, vector<double>(dynamicvars->n, 0.0));

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
    cout << "gear->centerDistance: " << gear->centerDistance << endl;

    // Theta_pk is the angle between the line connecting the gear and pinion's rotating axis and the x-axis
    double theta_pk = atan2(ycentergear - ycenterpinion, xcentergear - xcenterpinion); // Will always be 0 since both gears center on x-axis
    cout << "theta_pk: " << theta_pk << endl;

    double transmissionratio = pinion->gearRatio;
    double a_k_Pinion = sqrt(pow(xcentergear - xcenterpinion, 2) + pow(ycentergear - ycenterpinion, 2));
    double pitchCircleRadiusPinion = a_k_Pinion / (1 + transmissionratio);
    cout << "transmissionratio: " << transmissionratio << endl;
    cout << "a_k_Pinion: " << a_k_Pinion << endl;
    cout << "pitchCircleRadiusPinion: " << pitchCircleRadiusPinion << endl;

    // Eq (16)
    double pitchPointX, pitchPointY;
    if (ycentergear >= 0) {
        pitchPointX = pitchCircleRadiusPinion * cos(theta_pk) + xcenterpinion;
        pitchPointY = pitchCircleRadiusPinion * sin(theta_pk) + ycenterpinion;
        cout << "ycentergear >= 0" << endl;
    }
    else {
        pitchPointX = pitchCircleRadiusPinion * cos(theta_pk) + xcenterpinion;
        pitchPointY = -1.0 * pitchCircleRadiusPinion * sin(theta_pk) + ycenterpinion;
        cout << "ycentergear < 0" << endl;
    }

    cout << "pitchPointX: " << pitchPointX << endl;
    cout << "pitchPointY: " << pitchPointY << endl;
    double k1 = computeInvoluteSlope(pinion->baseRadius, pinion->pitchRadius);
    double k2 = computeInvoluteSlope(gear->baseRadius, gear->pitchRadius);

    cout << "k1: " << k1 << endl;
    cout << "k2: " << k2 << endl;

    // Line of mesh point - Eq (17) and Eq (18)
    // y = k1 * (x - pitchPointX) + pitchPointY;
    // y = k2 * (x - pitchPointX) + pitchPointY;

    double initialguess = 1.0;

    double u_mesh1 = solveForU(pinion->baseRadius, k1, pitchPointX, pitchPointY, initialguess);
    double x_mesh1, y_mesh1;
    computeInvolutePoint(pinion->baseRadius, u_mesh1, x_mesh1, y_mesh1);

    double u_mesh2 = solveForU(gear->baseRadius, k2, pitchPointX, pitchPointY, initialguess);
    double x_mesh2, y_mesh2;
    computeInvolutePoint(pinion->baseRadius, u_mesh2, x_mesh2, y_mesh2);
    cout << "u_mesh1: " << u_mesh1 << endl;
    cout << "u_mesh2: " << u_mesh2 << endl;
    cout << "Mesh point 1 coordinates: (" << x_mesh1 << ", " << y_mesh1 << ")" << endl;
    cout << "Mesh point 2 coordinates: (" << x_mesh2 << ", " << y_mesh2 << ")" << endl;

    // Calculate distance between two candidate mesh points.
    double delta = sqrt(pow(x_mesh1 - x_mesh2, 2) + pow(y_mesh1 - y_mesh2, 2));
    cout << "delta: " << delta << endl;
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
        cout << endl;
        cout << "Error: pitchRadius must be greater than baseRadius." << endl;
        cout << "pitchRadius: " << pitchRadius << " baseRadius: " << baseRadius << endl;
        cout << endl;
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

void computeInvolutePoint(double r_base, double mu_k, double& x, double& y) {
    // Compute the coordinates on the involute given parameter u and base circle radius r_base.
    // mu_k is the sum of theta_k and alpha_k.
    x = r_base * (sin(mu_k) - mu_k * cos(mu_k));
    y = r_base * (cos(mu_k) + mu_k * sin(mu_k));
}

//void computeInvolutePoint(double r_b, double u, double& x, double& y) {
//    x = r_b * (sin(u) - u * cos(u));
//    y = r_b * (cos(u) + u * sin(u));
//}

double F(double mu_k, double r_b, double k, double x_J, double y_J) {
    //double x_u, y_u;
    //computeInvolutePoint(r_b, u, x_u, y_u);
    double x_u = r_b * (sin(mu_k) - mu_k * cos(mu_k));
    double y_u = r_b * (cos(mu_k) + mu_k * sin(mu_k));
    return y_u - (k * (x_u - x_J) + y_J);
}

double dFdu(double u, double r_b, double k, double x_J, double y_J, double h) {
    // numerical derivative for F(u)
    return (F(u + h, r_b, k, x_J, y_J) - F(u - h, r_b, k, x_J, y_J)) / (2 * h);
}

double solveForU(double r_b, double k, double x_J, double y_J, double initialGuess, int maxIter, double tol) {
    // Newton's method to solve F(u) = 0
    double u = initialGuess;
    for (int i = 0; i < maxIter; i++) {
        double fVal = F(u, r_b, k, x_J, y_J);
        double fDeriv = dFdu(u, r_b, k, x_J, y_J);
        if (fabs(fDeriv) < tol) {  // Avoid division by near zero.
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

void computeContactVectors(const Vec2D& p1, const Vec2D& p2, const SpurGear& gear, Vec2D& n_g, Vec2D& n1, Vec2D& n1_t, Vec2D& n2, Vec2D& n2_t) {
    // p1: contact point on pinion, p2: contact point on gear
    // Eq (20): unit normal vector in the deformation direction.
    n_g = computeUnitVector(p1, p2);

    // Eq (21): unit vector from pinion center (0,0) to contact point p1.
    Vec2D O1 = { 0.0, 0.0 };
    n1 = computeUnitVector(p1, O1);

    // Eq (22): unit vector from gear center (centerDistance, 0) to contact point p2.
    Vec2D O2 = { gear.centerDistance, 0.0 };
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

Vec2D computeUnitVector(const Vec2D& p1, const Vec2D& p2) {
    Vec2D diff;
    diff.x = p1.x - p2.x;
    diff.y = p1.y - p2.y;
    double norm = sqrt(diff.x * diff.x + diff.y * diff.y);
    return { diff.x / norm, diff.y / norm };
}

Vec2D rotate90(const Vec2D& v) {
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
    fF = computeFrictionParameter(vt_mag, v0, vL);
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

Vec2D normalize(const Vec2D& v) {
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