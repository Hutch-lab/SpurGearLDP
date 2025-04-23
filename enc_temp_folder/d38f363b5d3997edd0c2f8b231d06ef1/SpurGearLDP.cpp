/* Physics 5810 - Computational Physics Final Project - Nick Hutchison */
/*  Spur Gear Load Distribution Program */
#include "FunctDefs.h"
#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <tuple>
#include <fstream>
//#include <gsl/gsl_interp.h>
//#include <gsl/gsl_spline.h>
using namespace std;

/*
f1 = m_pinion * x_pinion_doubledot
f2 = m_pinion * y_pinion_doubledot
f3 = momentInertia_pinion * phi_pinion_doubledot
f4 = m_gear * x_gear_doubledot
f5 = m_gear * y_gear_doubledot
f6 = momentInertia_gear * phi_gear_doubledot
*/

int main() {
    SpurGear pinion;
    SpurGear gear;
    DynamicsVars dynvars;
    startupmsg();
    cout << endl;
    dynvars.g = 9.80665;
    dynvars.frictionCoeff = FRICTION_COEFF;
    dynvars.v0_friction = V0_FRICTION_SI;
    dynvars.vL_friction = VL_FRICTION_SI;
    dynvars.deltaJudgmentThreshold = CONTACT_THRESHOLD_SI * 1000.0;
    dynvars.stepCounter = 0; // Initialize debug counter

    int method;
    cout << " Enter 0 to manually enter parameters or 1 to load parameters from file: ";
    cin >> method;
    if (method == 0) { if (!getUserInput(&pinion, &gear)) { cerr << "Failed input.\n"; return 1; } }
    else if (method == 1) { cout << "Read from file not implemented.\n"; return 1; }
    else { cout << "Invalid input method.\n"; return 1; }
    if (!setupProfiles(&pinion, &gear, &dynvars)) { cerr << "Failed setup.\n"; return 1; }
    initialToothPinion(&pinion, &dynvars); initialToothGear(&gear, &dynvars);
    completePinionProfile(&pinion, &dynvars); completeGearProfile(&gear, &dynvars);

    double density = 7.85 / 1.0e6; // g/cm^3 -> kg/mm^3

    double pinion_volume = M_PI * pow(pinion.pitchRadius, 2) * pinion.faceWidth; // mm^3
    pinion.mass = pinion_volume * density; // kg
    pinion.momentInertia = 0.5 * pinion.mass * pow(pinion.pitchRadius / 1000.0, 2); // kg*m^2

    double gear_volume = M_PI * pow(gear.pitchRadius, 2) * gear.faceWidth;
    gear.mass = gear_volume * density;
    gear.momentInertia = 0.5 * gear.mass * pow(gear.pitchRadius / 1000.0, 2);

    cout << "  Estimated Pinion Mass: " << pinion.mass << " kg, Inertia: " << pinion.momentInertia << " kg*m^2" << endl;
    cout << "  Estimated Gear Mass:   " << gear.mass << " kg, Inertia: " << gear.momentInertia << " kg*m^2" << endl;

    pinion.elasticCoeff = 1.0e7;
    pinion.dampingCoeff = 5.0e3;
    gear.elasticCoeff = 1.0e7;
    gear.dampingCoeff = 5.0e3;
    cout << "  Using Support K=" << pinion.elasticCoeff << " N/m, C=" << pinion.dampingCoeff << " Ns/m" << endl;
    dynvars.meshStiffness = MESH_STIFFNESS_SI;
    dynvars.meshDampingCoeff = getMeshDampingCoeff(dynvars.meshStiffness, MESH_DAMPING_RATIO, pinion.mass, gear.mass);
    cout << "  Using Mesh Kv=" << dynvars.meshStiffness << " N/m, Cd=" << dynvars.meshDampingCoeff << " Ns/m (zeta=" << MESH_DAMPING_RATIO << ")" << endl;

    // Initialize dynamic stress profile vector
    if (dynvars.n > 0) {
        pinion.contactStressProfile.assign(dynvars.n, 0.0); // Initialize with zeros
    }
    else {
        cerr << "Error: dynvars.n is zero or negative, cannot initialize stress profile." << endl;
        return 1;
    }
    cout << "\nInitialized Peak Dynamic Stress Profile vector (size " << dynvars.n << ") with zeros." << endl;

    // Static Stress Calc
    cout << "\nCalculating representative contact stress profile..." << endl;
    double pitchRadius_m = pinion.pitchRadius / 1000.0;
    double Wt_static = 0.0;
    if (abs(pitchRadius_m) > 1e-9) { Wt_static = pinion.torque / pitchRadius_m; }
    double omega_pinion = pinion.rpmSpeed * (2.0 * M_PI) / 60.0;
    double V_static = omega_pinion * pitchRadius_m;
    cout << "  Representative Static Tangential Load (Wt): " << Wt_static << " N" << endl;
    cout << "  Representative Pitch Line Velocity (V): " << V_static << " m/s" << endl;
    double Ko = 1.25;
    int Qv = 8;
    int gearQuality = 8;
    double Ks = 1.0;
    double ZR = 1.0;
    cout << "  Using Factors: Ko=" << Ko << ", Qv=" << Qv << ", Quality=" << gearQuality << ", Ks=" << Ks << ", ZR=" << ZR << endl;
    if (dynvars.n <= 0) { return 1; }
    pinion.contactStressProfile.resize(dynvars.n);
    cout << "  Calculating stress for " << dynvars.n << " points on the profile..." << endl;

    //for (int j = 0; j < dynvars.n; ++j) {
    //    pinion.contactStressProfile[j] = calcContactStress(&pinion, &gear, Wt_static, V_static, Ko, Qv, gearQuality, Ks, ZR);
    //    if (pinion.contactStressProfile[j] < 0) {
    //        pinion.contactStressProfile[j] = 0;
    //    }
    //}
    //cout << "  Calculated Contact Stress Profile (MPa): [";
    //for (int j = 0; j < dynvars.n; ++j) {
    //    cout << fixed << setprecision(2) << pinion.contactStressProfile[j] << (j == dynvars.n - 1 ? "" : ", ");
    //}
    //cout << "]" << endl;
	// End Static Stress Calc

    // Dynamic Simulation
    State state;

    // Initialize position and velocity states
    
    state.pos[0] = 0;
    state.pos[1] = 0;
    state.pos[2] = 0;
    state.pos[3] = pinion.centerDistance;
    state.pos[4] = 0;
    state.pos[5] = 0;
    double omega_pinion_init = pinion.rpmSpeed * (2.0 * M_PI) / 60.0;
    state.vel[0] = 0;
    state.vel[1] = 0;
    state.vel[2] = omega_pinion_init;
    state.vel[3] = 0;
    state.vel[4] = 0;
    state.vel[5] = -omega_pinion_init / pinion.gearRatio;
    
    double t = 0.0;
    double dt = 5e-2;
    double pinion_period = 1.0;

    if (abs(pinion.rpmSpeed) > 1e-6) {
        pinion_period = 1.0 / (pinion.rpmSpeed / 60.0);
    }
    else { pinion_period = 1.0; }

    double t_final = 2 * pinion_period;

    int steps = static_cast<int>(ceil(t_final / dt));
    cout << "\nRunning dynamic simulation for " << t_final << " s (" << steps << " steps with dt=" << dt << ")..." << endl;
    ofstream outfile("dynamics_output.txt"); outfile << "Time\t\tPinX\t\tPinY\t\tPinPhi\t\tGearX\t\tGearY\t\tGearPhi\t\tPinVelX\t\tPinVelY\t\tPinVelPhi\tGearVelX\tGearVelY\tGearVelPhi" << endl; outfile << fixed << setprecision(8);
    moveGearPosition(state, &pinion, &gear, &dynvars);

    for (int i = 0; i < steps; i++) {
        state = rungeKuttaStep(state, t, dt, pinion, gear, dynvars);
        
        t += dt;
        if (i % (steps > 1000 ? steps / 100 : 10) == 0 || i == steps - 1) { // Adjust output frequency
            cout << "Time: " << fixed << setprecision(7) << t << " | PinPhi(deg): " << fixed << setprecision(2) << state.pos[2] * 180.0 / M_PI << " | GearPhi(deg): " << fixed << setprecision(2) << state.pos[5] * 180.0 / M_PI << " | PinVelPhi: " << fixed << setprecision(2) << state.vel[2] << " | GearVelPhi: " << fixed << setprecision(2) << state.vel[5] << endl;
            outfile << t; for (int k = 0; k < 6; ++k) outfile << "\t" << state.pos[k]; for (int k = 0; k < 6; ++k) outfile << "\t" << state.vel[k]; outfile << endl;

            // Check for NaN/Inf and stop if unstable
            bool unstable = false;
            for (int k = 0; k < 6; ++k) { if (isnan(state.pos[k]) || isinf(state.pos[k]) || isnan(state.vel[k]) || isinf(state.vel[k])) { unstable = true; break; } }
            if (unstable) {
                cerr << "Error: Simulation became unstable at t = " << t << endl;
                break;
            }
        }
    }
    outfile.close();

    // Print Final Peak Dynamic Stress Profile
    cout << "\nFinal Peak Dynamic Contact Stress Profile (MPa):" << endl;
    cout << "Point Index | Stress (MPa)" << endl;
    cout << "---------------------------" << endl;
    for (size_t j = 0; j < pinion.contactStressProfile.size(); ++j) {
        cout << setw(11) << j << " | " << fixed << setprecision(2) << setw(12) << pinion.contactStressProfile[j] << endl;
    }
    cout << "---------------------------" << endl;
    return 0;
}

/*
int main() {
    SpurGear pinion;
    SpurGear gear;
    DynamicsVars dynvars;

    startupmsg();
    cout << endl;

    // Initialize dynamic parameters in dynvars
    dynvars.g = 9.80665; // m/s^2
    dynvars.frictionCoeff = FRICTION_COEFF;
    dynvars.v0_friction = V0_FRICTION_SI;
    dynvars.vL_friction = VL_FRICTION_SI;
    dynvars.deltaJudgmentThreshold = CONTACT_THRESHOLD_SI * 1000.0;
    dynvars.meshDampingCoeff = getMeshDampingCoeff(MESH_STIFFNESS_SI, MESH_DAMPING_RATIO, 1.0, 2.0) / 1000.0;
    dynvars.stepCounter = 0;

    // Get input method
    int method;
    cout << " Enter 0 to manually enter parameters or 1 to load parameters from file";
    cin >> method;

    if (method == 0) {
        int cont = getUserInput(&pinion, &gear);
    }
    else if (method == 1) {
        // add read from file
        cout << "Read from file not implemented yet." << endl;
        return 1;
    }
    else {
        cout << "Input Error";
        return 1;
    }
    // Setup involute shape for the pinion and gear
    setupProfiles(&pinion, &gear, &dynvars);
    cout << "setupProfiles worked!" << endl;
    initialToothPinion(&pinion, &dynvars);
    cout << "initialToothPinion worked!" << endl;
    initialToothGear(&gear, &dynvars);
    cout << "initialToothGear worked!" << endl;

    // Copy the tooth profile for every tooth
    completePinionProfile(&pinion, &dynvars);
    completeGearProfile(&gear, &dynvars);

    // Best guess at gear weight
    double density = 7.85e-9; // kg/mm^3 (8620 steel density)
    double pinion_volume = M_PI * pow(pinion.pitchRadius, 2) * pinion.faceWidth; // mm^3
    pinion.mass = pinion_volume * density; // kg
    pinion.momentInertia = 0.5 * pinion.mass * pow(pinion.pitchRadius / 1000.0, 2); // kg*m^2

    double gear_volume = M_PI * pow(gear.pitchRadius, 2) * gear.faceWidth; // mm^3
    gear.mass = gear_volume * density; // kg
    gear.momentInertia = 0.5 * gear.mass * pow(gear.pitchRadius / 1000.0, 2); // kg*m^2

    // Calculate Contact Stress @ static
    cout << "\nCalculating contact stress at static cond." << endl;
    double pitchRadius_m = pinion.pitchRadius / 1000.0; // Convert mm to m
    double Wt_static = pinion.torque / pitchRadius_m; // N
    cout << "  Static Tangential Load (Wt): " << Wt_static << " N" << endl;
    double omega_pinion = pinion.rpmSpeed * (2.0 * M_PI) / 60.0; // rad/s
    double V_static = omega_pinion * pitchRadius_m; // m/s
    cout << "  Pitch Line Velocity (V): " << V_static << " m/s" << endl;

    // Define stiffness / damping parameters
    pinion.elasticCoeff = 1.0e7; // N/m
    pinion.dampingCoeff = 5.0e3; // Ns/m
    gear.elasticCoeff = 1.0e7; // N/m
    gear.dampingCoeff = 5.0e3; // Ns/m
    dynvars.meshDampingCoeff = getMeshDampingCoeff(MESH_STIFFNESS_SI, MESH_DAMPING_RATIO, pinion.mass, gear.mass); // Ns/m
    dynvars.meshDampingCoeff /= 1000.0; // Convert Ns/m to Ns/mm for computeContactForce if it expects mm/s velocity

    // Define parameters for contact stress calculation
    double Ko = 1.25; // Overload factor (1.0=uniform, 1.25=light shock, 1.5=moderate, 1.75+=heavy)
    int Qv = 8;       // Transmission Accuracy Level (e.g., 6-11, lower is better)
    int gearQuality = 8; // AGMA Quality Number (e.g., 6-12)
    double Ks = 1.0;  // Size Factor (often 1.0 for contact stress)
    double ZR = 1.0;  // Surface Condition Factor (1.0 for typical surface finish)

    cout << "  Using Factors: Ko=" << Ko << ", Qv=" << Qv << ", Quality=" << gearQuality << ", Ks=" << Ks << ", ZR=" << ZR << endl;

    // This static calculation isn't representative of what actually occurs, kinda a gross oversimplification
    pinion.contactStressProfile.resize(dynvars.n);
    cout << "  Calculating stress for " << dynvars.n << " points on the profile..." << endl;
    for (int j = 0; j < dynvars.n; ++j) {
        // Call the stress function using the representative static values
        pinion.contactStressProfile[j] = calcContactStress(&pinion, &gear, Wt_static, V_static, Ko, Qv, gearQuality, Ks, ZR);
    }

    // Dynamic Simulation
    State state;
    // Initialize position and velocity states
    state.pos[0] = 0.0;
    state.pos[1] = 0.0;
    state.pos[2] = 0.0;
    state.pos[3] = pinion.centerDistance;
    state.pos[4] = 0.0;
    state.pos[5] = 0.0;
    state.vel[0] = 0.0;
    state.vel[1] = 0.0;
    state.vel[2] = omega_pinion;
    state.vel[3] = 0.0;
    state.vel[4] = 0.0;
    state.vel[5] = -omega_pinion / pinion.gearRatio;

    double t = 0.0;
    double dt = 5.0e-7;  // Time step (adjust as needed)
    double pinion_period = 1.0 / (pinion.rpmSpeed / 60.0);
    double t_final = 5.0 * pinion_period;
    int steps = static_cast<int>(ceil(t_final / dt));

    cout << "\nRunning dynamic simulation:" << endl;
    ofstream outfile("dynamics_output.txt");
    outfile << "Time\tPinX\tPinY\tPinPhi\tGearX\tGearY\tGearPhi\tPinVelX\tPinVelY\tPinVelPhi\tGearVelX\tGearVelY\tGearVelPhi" << endl;
    outfile << fixed << setprecision(8);
    for (int i = 0; i < steps; i++) {
            t += dt;
        //moveGearPosition(state, &pinion, &gear, &dynvars); // update tooth coords
        // Advance state using Runge-Kutta integration.
        state = rungeKuttaStep(state, t, dt, pinion, gear, dynvars);
        dynvars.phi_k = state.pos[2];

        if (i % (steps / 100) == 0 || i == steps - 1) {
            // Output full state to file
            outfile << t;
            for (int k = 0; k < 6; ++k) outfile << "\t" << state.pos[k];
            for (int k = 0; k < 6; ++k) outfile << "\t" << state.vel[k];
            outfile << endl;
        }
    }

    cout << "Simulation complete, Output saved to dynamics_output.txt" << endl;
    return 0;
}
*/

void outputGearProfile(const SpurGear& gear, const DynamicsVars& dynVars, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open file " << filename << " for writing." << std::endl;
        return;
    }

    // Write header: tooth index, profile point index, x_left, y_left, x_right, y_right
    outFile << "Tooth,ProfilePoint,x_left,y_left,x_right,y_right\n";

    // Loop over each tooth (column) and each profile point (row)
    for (int tooth = 0; tooth < gear.numTeeth; tooth++) {
        for (int i = 0; i < dynVars.n; i++) {
            // Assuming your complete gear profile is stored in initialToothProfileXLeft, etc.
            outFile << tooth << "," << i << ","
                << gear.initialToothProfileXLeft[tooth][i] << ","
                << gear.initialToothProfileYLeft[tooth][i] << ","
                << gear.initialToothProfileXRight[tooth][i] << ","
                << gear.initialToothProfileYRight[tooth][i] << "\n";
        }
    }
    outFile.close();
    std::cout << "Gear profile data written to " << filename << std::endl;
}

void startupmsg() {
    cout << " Physics 5810 - Computational Physics Final Project " << endl;
    cout << " ----------------- Nick Hutchison ----------------- " << endl;
    cout << " ------ Spur Gear Load Distribution Program ------- " << endl;
}

int getUserInput(SpurGear* pinion, SpurGear* gear) {
    /*
    Returns 0 when the user wants to exit and 1 if the user enters valid input data.

    INPUT:
    OUTPUT:
    CALLS:
     */
     /*
    cout << "Enter speed of pinion (rpm): ";
    cin >> pinion->rpmSpeed;
    cout << "Enter torque applied to pinion (Nm): ";
    cin >> pinion->torque;
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
        */
    pinion->rpmSpeed = 2000;
    pinion->torque = 50;
    pinion->numTeeth = 30.0;
    gear->numTeeth = 50.0;
    pinion->pressureAngle = 20.0;
    gear->pressureAngle = pinion->pressureAngle;
    //pinion->rootDiameter = 40.0;
    //gear->rootDiameter = 70.0;
    pinion->centerDistance = 0.0;
    pinion->faceWidth = 20.0;
    gear->faceWidth = 20.0;
    gear->offsetFaceWidth = 0.0;
    pinion->module = 2.5;
    gear->module = pinion->module;
    pinion->gearQuality = 8;
    pinion->Ko = 1.25;
    pinion->Ks = 1.0;
    pinion->ZR = 1.0;
    pinion->Qv = 8;

    pinion->addendum = 1.0 * pinion->module;
    gear->addendum = 1.0 * gear->module;

    pinion->dedendum = 1.25 * pinion->module;
    gear->dedendum = 1.25 * gear->module;

    pinion->gearRatio = gear->numTeeth / pinion->numTeeth;
    gear->gearRatio = pinion->gearRatio;
    pinion->pitchDiameter = pinion->numTeeth * pinion->module;
    gear->pitchDiameter = gear->numTeeth * gear->module;
    pinion->baseDiameter = pinion->pitchDiameter * cos(M_PI / 180.0 * pinion->pressureAngle);
    gear->baseDiameter = gear->pitchDiameter * cos(M_PI / 180.0 * gear->pressureAngle);
    pinion->outsideDiameter = pinion->pitchDiameter + 2 * pinion->addendum;
    gear->outsideDiameter = gear->pitchDiameter + 2 * gear->addendum;
    pinion->pitchRadius = pinion->pitchDiameter / 2.0;
    gear->pitchRadius = gear->pitchDiameter / 2.0;
    pinion->baseRadius = pinion->baseDiameter / 2.0;
    gear->baseRadius = gear->baseDiameter / 2.0;
    pinion->outsideRadius = pinion->outsideDiameter / 2.0;
    gear->outsideRadius = gear->outsideDiameter / 2.0;

    pinion->poissonsRatio = 0.28;
    gear->poissonsRatio = pinion->poissonsRatio;
    pinion->youngsModulus = 530.0e3; // MPa
    gear->youngsModulus = pinion->youngsModulus;

    if (pinion->centerDistance == 0.0) {
        // Calculate the nominal center distance between gears
        double cd = (gear->pitchDiameter + pinion->pitchDiameter) / 2.0;
        gear->centerDistance = cd;
        pinion->centerDistance = cd;
    }

    // Print calculated parameters
    cout << "\nCalculated Gear Parameters" << endl;
    cout << fixed << setprecision(3);
    cout << "Pinion:" << endl;
    cout << "  Pitch Diameter: " << pinion->pitchDiameter << " mm" << endl;
    cout << "  Base Diameter:  " << pinion->baseDiameter << " mm" << endl;
    cout << "  Outside Dia:    " << pinion->outsideDiameter << " mm" << endl;
    cout << "  Addendum:       " << pinion->addendum << " mm" << endl;
    cout << "  Dedendum:       " << pinion->dedendum << " mm" << endl;
    cout << "Gear:" << endl;
    cout << "  Pitch Diameter: " << gear->pitchDiameter << " mm" << endl;
    cout << "  Base Diameter:  " << gear->baseDiameter << " mm" << endl;
    cout << "  Outside Dia:    " << gear->outsideDiameter << " mm" << endl;
    cout << "  Addendum:       " << gear->addendum << " mm" << endl;
    cout << "  Dedendum:       " << gear->dedendum << " mm" << endl;
    cout << "Pair:" << endl;
    cout << "  Gear Ratio:     " << pinion->gearRatio << endl;
    cout << "  Center Dist:    " << pinion->centerDistance << " mm" << endl;
    cout << "----------------------------------" << endl;

    cout << "Parameters loaded!" << endl;
    return 1; // Indicate success
}

int setupProfiles(SpurGear* pinion, SpurGear* gear, DynamicsVars* dynamicvars) {
    //#define f1(F_pinion_sx, F_pinion_1x, F_pinion_2x, F_pinion_3x, F_pinion_4x) (F_pinion_sx - F_pinion_1x - F_pinion_2x - F_pinion_3x - F_pinion_4x)
    //#define f2(F_pinion_sy, F_pinion_1y, F_pinion_2y, F_pinion_3y, F_pinion_4y) (F_pinion_sy - F_pinion_1y - F_pinion_2y - F_pinion_3y - F_pinion_4y) - m_pinion * dynamicvars->g
    //#define f3(T_pinion_1, T_pinion_2, T_pinion_3, T_pinion_4, T_d) (T_pinion_1 + T_pinion_2 + T_pinion_3 + T_pinion_4 + T_d)
    //#define f4(F_gear_sx, F_gear_1x, F_gear_2x, F_gear_3x, F_gear_4x) (F_gear_sx - F_gear_1x - F_gear_2x - F_gear_3x - F_gear_4x)
    //#define f5(F_gear_sy, F_gear_1y, F_gear_2y, F_gear_3y, F_gear_4y) (F_gear_sy - F_gear_1y - F_gear_2y - F_gear_3y - F_gear_4y - m_gear * dynamicvars->g)
    //#define f6(T_gear_1, T_gear_2, T_gear_3, T_gear_4, T_l) (T_gear_1 + T_gear_2 + T_gear_3 + T_gear_4 + T_l)
    dynamicvars->n = 100;			// number points per profile
    pinion->theta_e = M_PI / pinion->numTeeth; // Angular pitch
    gear->theta_e = M_PI / gear->numTeeth;
    pinion->theta_a = computeThetaA(pinion->pressureAngle); // Roll angle to the pitch point
    gear->theta_a = computeThetaA(gear->pressureAngle);


    // Angle between left and right flanks on a tooth + transformations
    dynamicvars->trig_pinion = -1.0 * (pinion->theta_e + 2.0 * pinion->theta_a);
    pinion->theta_r = M_PI_2 + pinion->theta_a - pinion->theta_e / 2.;
    dynamicvars->trig_gear = -1.0 * (gear->theta_e + 2.0 * gear->theta_a);
    gear->theta_r = M_PI_2 + gear->theta_a - gear->theta_e / 2.;

    // Stepsize for generating points from base circle -> outside diameter
    double stepsize_pinion = (pinion->outsideRadius - pinion->baseRadius) / (dynamicvars->n);
    double stepsize_gear = (gear->outsideRadius - gear->baseRadius) / (dynamicvars->n);

    cout << "Resizing pinion vectors to " << dynamicvars->n << endl;

    pinion->position_x_left.resize(dynamicvars->n);
    pinion->position_y_left.resize(dynamicvars->n);
    pinion->position_x_right.resize(dynamicvars->n);
    pinion->position_y_right.resize(dynamicvars->n);

    pinion->initial_position_x_left.resize(dynamicvars->n);
    pinion->initial_position_y_left.resize(dynamicvars->n);
    pinion->initial_position_x_right.resize(dynamicvars->n);
    pinion->initial_position_y_right.resize(dynamicvars->n);

    cout << "  Pinion vector sizes after resize:" << endl;
    cout << "    pos_x_l: " << pinion->position_x_left.size() << endl;
    cout << "    init_pos_x_l: " << pinion->initial_position_x_left.size() << endl;

    cout << "Resizing gear vectors to " << dynamicvars->n << endl;

    gear->position_x_left.resize(dynamicvars->n);
    gear->position_y_left.resize(dynamicvars->n);
    gear->position_x_right.resize(dynamicvars->n);
    gear->position_y_right.resize(dynamicvars->n);

    gear->initial_position_x_left.resize(dynamicvars->n);
    gear->initial_position_y_left.resize(dynamicvars->n);
    gear->initial_position_x_right.resize(dynamicvars->n);
    gear->initial_position_y_right.resize(dynamicvars->n);

    cout << "  Gear vector sizes after resize:" << endl;
    cout << "    pos_x_l: " << gear->position_x_left.size() << endl;
    cout << "    init_pos_x_l: " << gear->initial_position_x_left.size() << endl;

    int i = 0;
    for (double r = pinion->baseRadius; r <= pinion->outsideRadius; r += stepsize_pinion) {
        double mu_k = compute_mu_k_from_r(r, pinion->baseRadius);
        pinion->position_x_left[i] = pinion->baseRadius * (sin(mu_k) - mu_k * cos(mu_k));
        pinion->position_y_left[i] = pinion->baseRadius * (cos(mu_k) + mu_k * sin(mu_k));
        pinion->position_x_right[i] = (-1.0 * pinion->position_x_left[i]) * cos(dynamicvars->trig_pinion) + pinion->position_y_left[i] * -1.0 * sin(dynamicvars->trig_pinion);
        pinion->position_y_right[i] = (-1.0 * pinion->position_x_left[i]) * sin(dynamicvars->trig_pinion) + pinion->position_y_left[i] * -1.0 * cos(dynamicvars->trig_pinion);
        cout << "Pinion Initial point # " << i << " [XY Left] [XY Right]: [" << pinion->position_x_left[i] << pinion->position_y_left[i] << "] [" << pinion->position_x_right[i] << pinion->position_y_right[i] << "]" << endl;
        if( i >= dynamicvars->n - 1){
            cout << "Pinion profile setup broke out!" << endl;
            break;
        }
        i = i + 1;
    }
    i = 0;
    for (double r = gear->baseRadius; r <= gear->outsideRadius; r += stepsize_gear) {
        double mu_k = compute_mu_k_from_r(r, gear->baseRadius);
        gear->position_x_left[i] = gear->baseRadius * (sin(mu_k) - mu_k * cos(mu_k));
        gear->position_y_left[i] = gear->baseRadius * (cos(mu_k) + mu_k * sin(mu_k));
        gear->position_x_right[i] = (-1.0 * gear->position_x_left[i]) * cos(dynamicvars->trig_gear) + gear->position_y_left[i] * -1.0 * sin(dynamicvars->trig_gear);
        gear->position_y_right[i] = (-1.0 * gear->position_x_left[i]) * sin(dynamicvars->trig_gear) + gear->position_y_left[i] * -1.0 * cos(dynamicvars->trig_gear);
        if (i >= dynamicvars->n - 1) {
            cout << "Gear profile setup broke out!" << endl;
            break;
        }
        i = i + 1;
    }
    return 1;
}

void initialToothPinion(SpurGear* pinion, DynamicsVars* dynamicvars) {
    // Eq (1) from paper
    // Create the initial left + right involute profiles for the pinion with dynamicvars->n number of points
    
    // Check for size or pointer errors
    if (!pinion || !dynamicvars) { cerr << "Error: Null pointer to initialToothPinion.\n"; return; }
    int n_points_int = dynamicvars->n;
    if (n_points_int <= 0) { cerr << "Error: Invalid n_points in initialToothPinion.\n"; return; }
    size_t n_points = static_cast<size_t>(n_points_int);
    if (pinion->position_x_left.size() != n_points || pinion->position_y_left.size() != n_points ||
        pinion->position_x_right.size() != n_points || pinion->position_y_right.size() != n_points ||
        pinion->initial_position_x_left.size() != n_points || pinion->initial_position_y_left.size() != n_points ||
        pinion->initial_position_x_right.size() != n_points || pinion->initial_position_y_right.size() != n_points) {
        cerr << "Error: Vector size mismatch in initialToothPinion. Expected " << n_points << endl;
        // cerr << "Sizes: pos_l=" << pinion->position_x_left.size() << " init_l=" << pinion->initial_position_x_left.size() << endl;
        return;
    }


    for (int i = 0; i < dynamicvars->n; i++) {
        pinion->initial_position_x_left[i]  = cos(-1.0 * pinion->theta_r) * pinion->position_x_left[i]  - sin(-1.0 * pinion->theta_r) * pinion->position_y_left[i];
        pinion->initial_position_y_left[i]  = sin(-1.0 * pinion->theta_r) * pinion->position_x_left[i]  + cos(-1.0 * pinion->theta_r) * pinion->position_y_left[i];
        pinion->initial_position_x_right[i] = cos(-1.0 * pinion->theta_r) * pinion->position_x_right[i] - sin(-1.0 * pinion->theta_r) * pinion->position_y_right[i];
        pinion->initial_position_y_right[i] = sin(-1.0 * pinion->theta_r) * pinion->position_x_right[i] + cos(-1.0 * pinion->theta_r) * pinion->position_y_right[i];
    }
    return;
}

void initialToothGear(SpurGear* gear, DynamicsVars* dynamicvars) {
    // Eq (2) from paper
    // Create the initial left + right involute profiles for the gear with dynamicvars->n number of points

    // Check for size or pointer errors
    if (!gear || !dynamicvars) { cerr << "Error: Null pointer to initialToothGear.\n"; return; }
    int n_points_int = dynamicvars->n;
    if (n_points_int <= 0) { cerr << "Error: Invalid n_points in initialToothGear.\n"; return; }
    size_t n_points = static_cast<size_t>(n_points_int);
    if (gear->position_x_left.size() != n_points || gear->position_y_left.size() != n_points ||
        gear->position_x_right.size() != n_points || gear->position_y_right.size() != n_points ||
        gear->initial_position_x_left.size() != n_points || gear->initial_position_y_left.size() != n_points ||
        gear->initial_position_x_right.size() != n_points || gear->initial_position_y_right.size() != n_points) {
        cerr << "Error: Vector size mismatch in initialToothGear. Expected " << n_points << endl;
        return;
    }

    for (int i = 0; i < dynamicvars->n; i++) {
        gear->initial_position_x_left[i]  = cos(gear->theta_r) * gear->position_x_left[i]  - sin(gear->theta_r) * gear->position_y_left[i];
        gear->initial_position_y_left[i]  = sin(gear->theta_r) * gear->position_x_left[i]  + cos(gear->theta_r) * gear->position_y_left[i];
        gear->initial_position_x_right[i] = cos(gear->theta_r) * gear->position_x_right[i] - sin(gear->theta_r) * gear->position_y_right[i];
        gear->initial_position_y_right[i] = sin(gear->theta_r) * gear->position_x_right[i] + cos(gear->theta_r) * gear->position_y_right[i];
    }
    return;
}

void completePinionProfile(SpurGear* pinion, DynamicsVars* dynamicvars) {
    // Eq (5)
    pinion->initialToothProfileXLeft.resize(pinion->numTeeth, vector<double>(dynamicvars->n));
    pinion->initialToothProfileYLeft.resize(pinion->numTeeth, vector<double>(dynamicvars->n));
    pinion->initialToothProfileXRight.resize(pinion->numTeeth, vector<double>(dynamicvars->n));
    pinion->initialToothProfileYRight.resize(pinion->numTeeth, vector<double>(dynamicvars->n));
    for (int i = 0; i < pinion->numTeeth; i++) {
        double angle = 2 * i * pinion->theta_e;
        for (int j = 0; j < dynamicvars->n; j++) {
            pinion->initialToothProfileXLeft[i][j] = cos(angle) * pinion->initial_position_x_left[j] - sin(angle) * pinion->initial_position_y_left[j];
            pinion->initialToothProfileYLeft[i][j] = sin(angle) * pinion->initial_position_x_left[j] + cos(angle) * pinion->initial_position_y_left[j];
            pinion->initialToothProfileXRight[i][j] = cos(angle) * pinion->initial_position_x_right[j] - sin(angle) * pinion->initial_position_y_right[j];
            pinion->initialToothProfileYRight[i][j] = sin(angle) * pinion->initial_position_x_right[j] + cos(angle) * pinion->initial_position_y_right[j];
        }
    }
    outputGearProfile(*pinion, *dynamicvars, "pinionProfile.csv");
    return;
}

void completeGearProfile(SpurGear* gear, DynamicsVars* dynamicvars) {
    // Eq (6)
    gear->initialToothProfileXLeft.resize(gear->numTeeth, vector<double>(dynamicvars->n));
    gear->initialToothProfileYLeft.resize(gear->numTeeth, vector<double>(dynamicvars->n));
    gear->initialToothProfileXRight.resize(gear->numTeeth, vector<double>(dynamicvars->n));
    gear->initialToothProfileYRight.resize(gear->numTeeth, vector<double>(dynamicvars->n));
    for (int i = 0; i < gear->numTeeth; i++) {
        double angle = 2.0 * i * gear->theta_e;
        for (int j = 0; j < dynamicvars->n; j++) {
            gear->initialToothProfileXLeft[i][j] = cos(angle) * gear->initial_position_x_left[j] - sin(angle) * gear->initial_position_y_left[j] + gear->centerDistance;
            gear->initialToothProfileYLeft[i][j] = sin(angle) * gear->initial_position_x_left[j] + cos(angle) * gear->initial_position_y_left[j];
            gear->initialToothProfileXRight[i][j] = cos(angle) * gear->initial_position_x_right[j] - sin(angle) * gear->initial_position_y_right[j] + gear->centerDistance;
            gear->initialToothProfileYRight[i][j] = sin(angle) * gear->initial_position_x_right[j] + cos(angle) * gear->initial_position_y_right[j];
        }
    }
    outputGearProfile(*gear, *dynamicvars, "gearProfile.csv");
    return;
}

double computeThetaA(double pressureAngle) {
    double pressureAngle_rad = pressureAngle * M_PI / 180.0;
    if (pressureAngle <= 0 || pressureAngle >= 90) return 0;
    double cos_pa = cos(pressureAngle_rad);
    if (abs(cos_pa) < 1e-9) return 0; // Avoid division by zero
    return tan(pressureAngle_rad);
}

void computeInvolutePoint(double r_b, double mu_k, double& x, double& y) {
    // r_b: base radius (mm)
    // mu_k: roll angle (radians)
    // x, y: output coordinates (mm)
    x = r_b * (sin(mu_k) - mu_k * cos(mu_k));
    y = r_b * (cos(mu_k) + mu_k * sin(mu_k));
}

double compute_mu_k_from_r(double r, double r_base) {
    // Given a desired point on the involute (for example, determined by a specific radius r)
    // compute mu_k via r = r_base / cos(mu_k)
    return acos(r_base / r); // since cos(mu_k) = r_base / r, so u = arccos(r_base/r)
}

State computeDerivatives(const State& state, double /*t*/, SpurGear* pinion, SpurGear* gear, DynamicsVars* dynvars) {
    State dState;

    // Calculate Dynamic Contact Forces/Torques
    Vec2D Fp_contact_total, Fg_contact_total; // N
    double Tp_contact_total, Tg_contact_total; // Nm
    dynamicMeshingUpdate(state, pinion, gear, dynvars, Fp_contact_total, Fg_contact_total, Tp_contact_total, Tg_contact_total);

    double pos_m[6], vel_m[6];
    for (int i = 0; i < 6; ++i) { pos_m[i] = state.pos[i] / 1000.0; vel_m[i] = state.vel[i] / 1000.0; }
    pos_m[2] = state.pos[2]; vel_m[2] = state.vel[2];
    pos_m[5] = state.pos[5]; vel_m[5] = state.vel[5];

    // Compute forces
    double F_pinion_sx = -1.0 * (pinion->elasticCoeff * pos_m[0] + pinion->dampingCoeff * vel_m[0]); // N
    double F_pinion_sy = -1.0 * (pinion->elasticCoeff * pos_m[1] + pinion->dampingCoeff * vel_m[1]); // N
    double center_dist_m = pinion->centerDistance / 1000.0;
    double gear_dx_m = pos_m[3] - center_dist_m;
    double gear_dy_m = pos_m[4];
    double F_gear_sx = -1.0 * (gear->elasticCoeff * gear_dx_m + gear->dampingCoeff * vel_m[3]); // N
    double F_gear_sy = -1.0 * (gear->elasticCoeff * gear_dy_m + gear->dampingCoeff * vel_m[4]); // N

    // External Torques
    double T_d = pinion->torque; // Driving torque on pinion (Nm)
    // Load torque on gear - needs sign convention
    double T_l = -pinion->torque / pinion->gearRatio; // Load torque on gear (Nm) - simplified static reaction

    // Gravitational force
    double grav_force_p = pinion->mass * dynvars->g;
    double grav_force_g = gear->mass * dynvars->g;

    // Derivative of positions are the velocities.
    for (int i = 0; i < 6; i++) {
        dState.pos[i] = state.vel[i];
    }

    // Calculate accelerations using Newton-Euler equations (f1-f6 macros)
    if (pinion->mass <= 0 || pinion->momentInertia <= 0 || gear->mass <= 0 || gear->momentInertia <= 0) {
        cerr << "Error: Mass or Moment of Inertia is zero or negative in computeDerivatives." << endl;
        for (int i = 0; i < 6; ++i) dState.vel[i] = 0;
        return dState;
    }

    // Pinion accelerations
    dState.vel[0] = (F_pinion_sx - Fp_contact_total.x) / pinion->mass * 1000.0;
    dState.vel[1] = (F_pinion_sy - Fp_contact_total.y - grav_force_p) / pinion->mass * 1000.0;
    dState.vel[2] = (Tp_contact_total + T_d) / pinion->momentInertia;

    // Gear accelerations
    dState.vel[3] = (F_gear_sx - Fg_contact_total.x) / gear->mass * 1000.0;
    dState.vel[4] = (F_gear_sy - Fg_contact_total.y - grav_force_g) / gear->mass * 1000.0;
    dState.vel[5] = (Tg_contact_total + T_l) / gear->momentInertia;
    return dState;
}

State rungeKuttaStep(const State& state, double t, double dt, SpurGear& pinion, SpurGear& gear, DynamicsVars& dynvars) {
    State k1 = computeDerivatives(state, t, &pinion, &gear, &dynvars);
    State state2 = state;
    for (int i = 0; i < 6; i++) {
        state2.pos[i] += 0.5 * dt * k1.pos[i];
        state2.vel[i] += 0.5 * dt * k1.vel[i];
    }

    State k2 = computeDerivatives(state2, t + 0.5 * dt, &pinion, &gear, &dynvars);
    State state3 = state;
    for (int i = 0; i < 6; i++) {
        state3.pos[i] += 0.5 * dt * k2.pos[i];
        state3.vel[i] += 0.5 * dt * k2.vel[i];
    }

    State k3 = computeDerivatives(state3, t + 0.5 * dt, &pinion, &gear, &dynvars);
    State state4 = state;
    for (int i = 0; i < 6; i++) {
        state4.pos[i] += dt * k3.pos[i];
        state4.vel[i] += dt * k3.vel[i];
    }
    moveGearPosition(state4, &pinion, &gear, &dynvars);
    State k4 = computeDerivatives(state4, t + dt, &pinion, &gear, &dynvars);

    State newState = state;
    for (int i = 0; i < 6; i++) {
        newState.pos[i] += (dt / 6.0) * (k1.pos[i] + 2.0 * k2.pos[i] + 2.0 * k3.pos[i] + k4.pos[i]);
        newState.vel[i] += (dt / 6.0) * (k1.vel[i] + 2.0 * k2.vel[i] + 2.0 * k3.vel[i] + k4.vel[i]);
    }
    return newState;
}

void moveGearPosition(const State& state, SpurGear* pinion, SpurGear* gear, DynamicsVars* /*dynvars*/) {
    int num_pinion_teeth = static_cast<int>(pinion->numTeeth);
    int num_gear_teeth = static_cast<int>(gear->numTeeth);
    size_t n_points = pinion->initialToothProfileXLeft[0].size();
    size_t st_num_pinion = static_cast<size_t>(num_pinion_teeth);
    n_points = static_cast<size_t>(n_points);

    // Ensure current profile vectors are allocated
    if (pinion->currentToothProfileXLeft.size() != num_pinion_teeth || pinion->currentToothProfileXLeft[0].size() != n_points) {
        pinion->currentToothProfileXLeft.resize(num_pinion_teeth, vector<double>(n_points));
        pinion->currentToothProfileYLeft.resize(num_pinion_teeth, vector<double>(n_points));
        pinion->currentToothProfileXRight.resize(num_pinion_teeth, vector<double>(n_points));
        pinion->currentToothProfileYRight.resize(num_pinion_teeth, vector<double>(n_points));
    }
    if (gear->currentToothProfileXLeft.size() != num_gear_teeth || gear->currentToothProfileXLeft[0].size() != n_points) {
        gear->currentToothProfileXLeft.resize(num_gear_teeth, vector<double>(n_points));
        gear->currentToothProfileYLeft.resize(num_gear_teeth, vector<double>(n_points));
        gear->currentToothProfileXRight.resize(num_gear_teeth, vector<double>(n_points));
        gear->currentToothProfileYRight.resize(num_gear_teeth, vector<double>(n_points));
    }

    double phi_p = state.pos[2]; // Pinion angle (rad)
    double cos_p = cos(phi_p);
    double sin_p = sin(phi_p);
    double x_p = state.pos[0]; // Pinion center x (mm)
    double y_p = state.pos[1]; // Pinion center y (mm)

    double phi_g = state.pos[5]; // Gear angle (rad)
    double cos_g = cos(phi_g);
    double sin_g = sin(phi_g);
    double x_g = state.pos[3]; // Gear center x (mm)
    double y_g = state.pos[4]; // Gear center y (mm)
    double nominal_center_dist = pinion->centerDistance; // Assuming this is nominal

    // Update Pinion Profiles
    for (int i = 0; i < num_pinion_teeth; ++i) {
        for (int j = 0; j < n_points; ++j) {
            double x0_l = pinion->initialToothProfileXLeft[i][j];
            double y0_l = pinion->initialToothProfileYLeft[i][j];
            pinion->currentToothProfileXLeft[i][j] = x0_l * cos_p - y0_l * sin_p + x_p;
            pinion->currentToothProfileYLeft[i][j] = x0_l * sin_p + y0_l * cos_p + y_p;

            double x0_r = pinion->initialToothProfileXRight[i][j];
            double y0_r = pinion->initialToothProfileYRight[i][j];
            pinion->currentToothProfileXRight[i][j] = x0_r * cos_p - y0_r * sin_p + x_p;
            pinion->currentToothProfileYRight[i][j] = x0_r * sin_p + y0_r * cos_p + y_p;
            cout << "Pinion Tooth "<< i << " point # " << j << " [XY Left] [XY Right]: [" << pinion->currentToothProfileXLeft[i][j] << pinion->currentToothProfileYLeft[i][j] << "] [" << pinion->currentToothProfileXRight[i][j] << pinion->currentToothProfileYRight[i][j] << "]" << endl;
        }
    }
    // Update Gear Profiles (Rotate around nominal center, then translate along x-axis)
    for (int i = 0; i < num_gear_teeth; ++i) {
        for (int j = 0; j < n_points; ++j) {
            // Point relative to gear's nominal center (nominal_center_dist, 0)
            double x0_l_rel = gear->initialToothProfileXLeft[i][j] - nominal_center_dist;
            double y0_l_rel = gear->initialToothProfileYLeft[i][j];
            // Rotate around nominal center
            double x_rot_l = x0_l_rel * cos_g - y0_l_rel * sin_g;
            double y_rot_l = x0_l_rel * sin_g + y0_l_rel * cos_g;
            // Translate to current gear center (x_g, y_g)
            gear->currentToothProfileXLeft[i][j] = x_rot_l + x_g;
            gear->currentToothProfileYLeft[i][j] = y_rot_l + y_g;

            double x0_r_rel = gear->initialToothProfileXRight[i][j] - nominal_center_dist;
            double y0_r_rel = gear->initialToothProfileYRight[i][j];
            double x_rot_r = x0_r_rel * cos_g - y0_r_rel * sin_g;
            double y_rot_r = x0_r_rel * sin_g + y0_r_rel * cos_g;
            gear->currentToothProfileXRight[i][j] = x_rot_r + x_g;
            gear->currentToothProfileYRight[i][j] = y_rot_r + y_g;
        }
    }
}

void dynamicMeshingUpdate(const State& state, SpurGear* pinion, SpurGear* gear, DynamicsVars* dynvars, Vec2D& total_pinion_force, Vec2D& total_gear_force, double& total_pinion_torque, double& total_gear_torque) {

    total_pinion_force = Vec2D(0, 0);
    total_gear_force = Vec2D(0, 0);
    total_pinion_torque = 0;
    total_gear_torque = 0;
    pinion->activeContactPairs.clear();
    Vec2D O1(state.pos[0] / 1000.0, state.pos[1] / 1000.0);
    Vec2D O2(state.pos[3] / 1000.0, state.pos[4] / 1000.0);
    double phi1 = state.pos[2];
    double phi2 = state.pos[5];
    double omega1 = state.vel[2];
    double omega2 = state.vel[5];
    Vec2D V1(state.vel[0] / 1000.0, state.vel[1] / 1000.0); Vec2D V2(state.vel[3] / 1000.0, state.vel[4] / 1000.0);
    double pressureAngle_rad = pinion->pressureAngle * M_PI / 180.0;
    int num_pinion_teeth = static_cast<int>(pinion->numTeeth);
    int num_gear_teeth = static_cast<int>(gear->numTeeth);
    size_t n_points = static_cast<size_t>(dynvars->n);
    if (n_points == 0 || num_pinion_teeth <= 0 || num_gear_teeth <= 0) { return; }
    if (pinion->currentToothProfileXLeft.size() != static_cast<size_t>(num_pinion_teeth) || pinion->currentToothProfileXLeft[0].size() != n_points ||
        gear->currentToothProfileXLeft.size() != static_cast<size_t>(num_gear_teeth) || gear->currentToothProfileXLeft[0].size() != n_points) {
        return;
    }

    double Ze = calculateZe(pinion, gear);
    double ZI = calculateZI(pinion, gear); // Using static version
    double F_mm = min(pinion->faceWidth, gear->faceWidth);
    double F_in = F_mm / 25.4;
    int Qv = pinion->Qv;
    int gearQuality = pinion->gearQuality;
    double Ko = pinion->Ko;
    double Ks = pinion->Ks;
    double ZR = pinion->ZR;
    double KH = calculateKH(pinion, gear, F_in, gearQuality);

    double pinion_angle_tooth0 = pinion->theta_r;
    double gear_angle_tooth0 = gear->theta_r;
    double center_line_angle = atan2(O2.y - O1.y, O2.x - O1.x);
    double target_angle_p = center_line_angle - pressureAngle_rad;
    double target_angle_g = center_line_angle + pressureAngle_rad - M_PI;

    target_angle_p = fmod(target_angle_p, 2.0 * M_PI);
    if (target_angle_p < 0) target_angle_p += 2.0 * M_PI;
    target_angle_g = fmod(target_angle_g, 2.0 * M_PI);
    if (target_angle_g < 0) target_angle_g += 2.0 * M_PI;

    double current_tooth0_angle_p = fmod(phi1 + pinion_angle_tooth0, 2.0 * M_PI);
    if (current_tooth0_angle_p < 0) current_tooth0_angle_p += 2.0 * M_PI;
    double current_tooth0_angle_g = fmod(phi2 + gear_angle_tooth0, 2.0 * M_PI);
    if (current_tooth0_angle_g < 0) current_tooth0_angle_g += 2.0 * M_PI;

    double angle_diff_p = target_angle_p - current_tooth0_angle_p;
    double angle_diff_g = target_angle_g - current_tooth0_angle_g;

    angle_diff_p = fmod(angle_diff_p, 2.0 * M_PI);
    if (angle_diff_p < 0) angle_diff_p += 2.0 * M_PI;
    angle_diff_g = fmod(angle_diff_g, 2.0 * M_PI);
    if (angle_diff_g < 0) angle_diff_g += 2.0 * M_PI;

    int p_idx_center = static_cast<int>(round(angle_diff_p / (2.0 * pinion->theta_e))) % num_pinion_teeth;
    int g_idx_center = static_cast<int>(round(angle_diff_g / (2.0 * gear->theta_e))) % num_gear_teeth;
    if (p_idx_center < 0) p_idx_center += num_pinion_teeth;
    if (g_idx_center < 0) g_idx_center += num_gear_teeth;

    bool contact_found_this_step = false; // Flag for debug print

    for (int p_offset = -1; p_offset <= 1; ++p_offset) {
        int p_idx = (p_idx_center + p_offset + num_pinion_teeth) % num_pinion_teeth;
        for (int g_offset = -1; g_offset <= 1; ++g_offset) {
            int g_idx = (g_idx_center + g_offset + num_gear_teeth) % num_gear_teeth;
            for (size_t j = 0; j < n_points; ++j) {
                Vec2D p1(pinion->currentToothProfileXLeft[p_idx][j] / 1000.0, pinion->currentToothProfileYLeft[p_idx][j] / 1000.0);
                Vec2D p2(gear->currentToothProfileXLeft[g_idx][j] / 1000.0, gear->currentToothProfileYLeft[g_idx][j] / 1000.0);
                Vec2D diff = p1 - p2;
                double dist = mag(diff);
                double dist_sq = dot(diff, diff);

                // Debug - Distance Check
                // Print only if points are kinda close
                bool print_debug = (dynvars->stepCounter < 50); // Limit debug prints
                if (print_debug && dist < CONTACT_THRESHOLD_SI * 20.0) { // Check if within 20x threshold
                    cout << fixed << setprecision(9);
                    cout << "  DEBUG Step " << dynvars->stepCounter << ": P" << p_idx << "[" << j << "] G" << g_idx << "[" << j << "] Dist=" << dist * 1e6 << "um";
                }

                if (dist_sq < pow(CONTACT_THRESHOLD_SI * 2.0, 2)) {
                    double normal_angle = center_line_angle - pressureAngle_rad;
                    Vec2D n_g = Vec2D(cos(normal_angle), sin(normal_angle));
                    double penetration = dot(diff, n_g);

                    if (print_debug) { // Print if close enough
                        cout << " | Pen=" << penetration * 1e6 << "um";
                    }
                    if (penetration > PENETRATION_THRESHOLD_SI) {
                        contact_found_this_step = true; // Mark contact found

                        // Debug - Contact detected
                        if (print_debug) {
                            cout << " | CONTACT DETECTED " << endl;
                        }

                        double delta = penetration;
                        Vec2D r1 = p1 - O1;
                        Vec2D r2 = p2 - O2;
                        Vec2D v1 = V1 + Vec2D(-omega1 * r1.y, omega1 * r1.x);
                        Vec2D v2 = V2 + Vec2D(-omega2 * r2.y, omega2 * r2.x);
                        Vec2D v_rel = v1 - v2;
                        double delta_dot = dot(v_rel, n_g);
                        Vec2D v_t = v_rel - n_g * delta_dot;
                        double vt_mag = mag(v_t);
                        double kv = getMeshStiffness(pinion, gear, p1);
                        double cd = dynvars->meshDampingCoeff;
                        double mu = dynvars->frictionCoeff;
                        double fF = computeFrictionParameter(vt_mag, dynvars->v0_friction, dynvars->vL_friction);
                        Vec2D F_contact = computeContactForce(delta, delta_dot, kv, cd, mu, fF, v_t, n_g);

                        // having issue with stress values being all zeros. printing statements to the terminal to debug issue.
                        if (dynvars->stepCounter < 100) {
                            cout << "Debug - contact detected at step " << dynvars->stepCounter << " between P" << p_idx << "/G" << g_idx << " point " << j << " (Penetration: " << delta * 1e6 << " um" << endl;
                        }
                        // Accumulate forces/torques
                        total_pinion_force = total_pinion_force + F_contact;
                        total_gear_force = total_gear_force - F_contact;
                        total_pinion_torque += computeMoment(p1, O1, F_contact);
                        total_gear_torque += computeMoment(p2, O2, -F_contact);
                        pinion->activeContactPairs.push_back({ p1, p2 });

                        // Calculate Dynamic Stress & Update Peak Profile
                        double Fn_mag = dot(F_contact, n_g); // Normal component of total contact force
                        if (Fn_mag > 0) { // Only calculate stress if there's positive normal force
                            // Approximate Wt from Normal Force
                            double Wt_dyn = Fn_mag / cos(pressureAngle_rad); // N

                            // Dynamic Velocity for Kv (use instantaneous pitch line velocity)
                            double V_dyn_mps = abs(omega1 * (pinion->pitchRadius / 1000.0)); // m/s
                            double V_dyn_fpm = V_dyn_mps * 196.850; // ft/min
                            double Kv_dyn = calculateKv(V_dyn_fpm, Qv);

                            // Calculate instantaneous stress using dynamic Wt, Kv and static factors
                            double sigma_c_inst = calcContactStress(pinion, gear, Wt_dyn, V_dyn_mps, Ko, Qv, gearQuality, Ks, ZR);

                            if (dynvars->stepCounter < 100) {
                                cout << "  Debug - Wt_dyn=" << Wt_dyn << " V_dyn=" << V_dyn_mps << " Kv_dyn=" << Kv_dyn << " -> sigma_c_inst=" << sigma_c_inst << endl;
                            }

                            if (sigma_c_inst >= 0) { // Check for valid stress calculation
                                // Transform p1 (global, m) to local initial frame (mm)
                                Vec2D p1_rel_center_m = p1 - O1;
                                // Rotate back by current angle phi1
                                Vec2D p1_local_rotated_m = rotate(p1_rel_center_m, -phi1);
                                Vec2D p1_local_initial_mm = p1_local_rotated_m * 1000.0;
                                size_t closest_j = findClosestProfileIndex(p1_local_initial_mm, pinion);
                                if (closest_j < pinion->contactStressProfile.size()) {
                                    pinion->contactStressProfile[closest_j] = max(pinion->contactStressProfile[closest_j], sigma_c_inst);
                                }
                            }
                        }
                    }
                    else {
                        // Debug - Close but No Penetration
                        if (print_debug) {
                            cout << " | No Penetration" << endl;
                        }
                    }
                }
                else {
                    // Debug - Too Far
                    if (print_debug && dist < CONTACT_THRESHOLD_SI * 20.0) { // Only print if kinda close
                        cout << " | Too Far" << endl;
                    }
                }
            } // end loop j
        } // end loop g_offset
    }// end loop p_offset
   // Increment step counter (moved outside loops)
    dynvars->stepCounter++;
}

/*
void dynamicMeshingUpdate(const State& state, SpurGear* pinion, SpurGear* gear, DynamicsVars* dynvars, Vec2D& total_pinion_force, Vec2D& total_gear_force, double& total_pinion_torque, double& total_gear_torque) {
    // Reset totals
    total_pinion_force = Vec2D(0, 0);
    total_gear_force = Vec2D(0, 0);
    total_pinion_torque = 0;
    total_gear_torque = 0;
    pinion->activeContactPairs.clear(); // Clear previous points

    // Current centers (m)
    Vec2D O1(state.pos[0] / 1000.0, state.pos[1] / 1000.0);
    Vec2D O2(state.pos[3] / 1000.0, state.pos[4] / 1000.0);

    // Current angular positions (rad) and velocities (rad/s)
    double phi1 = state.pos[2];
    double phi2 = state.pos[5];
    double omega1 = state.vel[2];
    double omega2 = state.vel[5];

    // Center velocities (m/s)
    Vec2D V1(state.vel[0] / 1000.0, state.vel[1] / 1000.0);
    Vec2D V2(state.vel[3] / 1000.0, state.vel[4] / 1000.0);

    // Gear parameters in SI (m)
    double pressureAngle_rad = pinion->pressureAngle * M_PI / 180.0; // Assuming same

    // 1. Identify potentially meshing teeth pairs
    int num_pinion_teeth = static_cast<int>(pinion->numTeeth);
    int num_gear_teeth = static_cast<int>(gear->numTeeth);
    int n_points = dynvars->n; // Points per profile flank

    if (n_points <= 0 || num_pinion_teeth <= 0 || num_gear_teeth <= 0) { return; }
    if (pinion->currentToothProfileXLeft.size() != static_cast<size_t>(num_pinion_teeth) ||
        pinion->currentToothProfileXLeft[0].size() != static_cast<size_t>(n_points) ||
        gear->currentToothProfileXLeft.size() != static_cast<size_t>(num_gear_teeth) ||
        gear->currentToothProfileXLeft[0].size() != static_cast<size_t>(n_points)) {
        return; // Vectors not sized correctly
    }

    double pinion_angle_tooth0 = pinion->theta_r;
    double gear_angle_tooth0 = gear->theta_r;
    double center_line_angle = atan2(O2.y - O1.y, O2.x - O1.x);
    double target_angle_p = center_line_angle - pressureAngle_rad;
    double target_angle_g = center_line_angle + pressureAngle_rad - M_PI;
    target_angle_p = fmod(target_angle_p, 2.0 * M_PI); if (target_angle_p < 0) target_angle_p += 2.0 * M_PI;
    target_angle_g = fmod(target_angle_g, 2.0 * M_PI); if (target_angle_g < 0) target_angle_g += 2.0 * M_PI;
    double current_tooth0_angle_p = fmod(phi1 + pinion_angle_tooth0, 2.0 * M_PI); if (current_tooth0_angle_p < 0) current_tooth0_angle_p += 2.0 * M_PI;
    double current_tooth0_angle_g = fmod(phi2 + gear_angle_tooth0, 2.0 * M_PI); if (current_tooth0_angle_g < 0) current_tooth0_angle_g += 2.0 * M_PI;
    double angle_diff_p = target_angle_p - current_tooth0_angle_p;
    double angle_diff_g = target_angle_g - current_tooth0_angle_g;
    angle_diff_p = fmod(angle_diff_p, 2.0 * M_PI); if (angle_diff_p < 0) angle_diff_p += 2.0 * M_PI;
    angle_diff_g = fmod(angle_diff_g, 2.0 * M_PI); if (angle_diff_g < 0) angle_diff_g += 2.0 * M_PI;

    int p_idx_center = static_cast<int>(round(angle_diff_p / (2.0 * pinion->theta_e))) % num_pinion_teeth;
    int g_idx_center = static_cast<int>(round(angle_diff_g / (2.0 * gear->theta_e))) % num_gear_teeth;
    if (p_idx_center < 0) p_idx_center += num_pinion_teeth;
    if (g_idx_center < 0) g_idx_center += num_gear_teeth;

    // Check tooth pairs around these central indices (+/- 1 or 2)
    for (int p_offset = -1; p_offset <= 1; ++p_offset) {
        int p_idx = (p_idx_center + p_offset + num_pinion_teeth) % num_pinion_teeth;
        for (int g_offset = -1; g_offset <= 1; ++g_offset) {
            int g_idx = (g_idx_center + g_offset + num_gear_teeth) % num_gear_teeth;

            // Check for contact between points on these teeth (Driving: Pinion Left vs Gear Left)
            for (int j = 0; j < n_points; ++j) {
                Vec2D p1(pinion->currentToothProfileXLeft[p_idx][j] / 1000.0, pinion->currentToothProfileYLeft[p_idx][j] / 1000.0); // m
                Vec2D p2(gear->currentToothProfileXLeft[g_idx][j] / 1000.0, gear->currentToothProfileYLeft[g_idx][j] / 1000.0); // m
                Vec2D diff = p1 - p2;
                double dist_sq = dot(diff, diff);

                if (dist_sq < pow(CONTACT_THRESHOLD_SI * 2.0, 2)) {
                    double normal_angle = center_line_angle - pressureAngle_rad;
                    Vec2D n_g = Vec2D(cos(normal_angle), sin(normal_angle));
                    double penetration = dot(diff, n_g);

                    if (penetration > PENETRATION_THRESHOLD_SI) {
                        double delta = penetration;
                        Vec2D r1 = p1 - O1;
                        Vec2D r2 = p2 - O2;
                        Vec2D v1 = V1 + Vec2D(-omega1 * r1.y, omega1 * r1.x);
                        Vec2D v2 = V2 + Vec2D(-omega2 * r2.y, omega2 * r2.x);
                        Vec2D v_rel = v1 - v2;
                        double delta_dot = dot(v_rel, n_g);
                        Vec2D v_t = v_rel - n_g * delta_dot;
                        double vt_mag = mag(v_t);
                        double kv = getMeshStiffness(pinion, gear, p1);
                        double cd = dynvars->meshDampingCoeff; // Use SI value from dynvars
                        double mu = dynvars->frictionCoeff;
                        double fF = computeFrictionParameter(vt_mag, dynvars->v0_friction, dynvars->vL_friction);
                        Vec2D F_contact = computeContactForce(delta, delta_dot, kv, cd, mu, fF, v_t, n_g); // N

                        total_pinion_force = total_pinion_force + F_contact;
                        total_gear_force = total_gear_force - F_contact;
                        total_pinion_torque += computeMoment(p1, O1, F_contact);
                        total_gear_torque += computeMoment(p2, O2, -F_contact);
                        pinion->activeContactPairs.push_back({ p1, p2 });
                    }
                }
            }
            // TODO: Add check for backlash contact (Pinion Right vs Gear Right)
        }
    }
}
*/

void computeContactVectors(const Vec2D& p1, const Vec2D& p2, const Vec2D& O1, const Vec2D& O2, Vec2D& n_g, Vec2D& n1, Vec2D& n1_t, Vec2D& n2, Vec2D& n2_t) {
    // p1: contact point on pinion, p2: contact point on gear
    // Eq (20): unit normal vector in the deformation direction.
    n_g = normalize(p1 - p2); // Vector from p2 to p1, normalized

    // Eq (21): unit vector from pinion center O1 to contact point p1.
    n1 = normalize(p1 - O1);

    // Eq (22): unit vector from gear center O2 to contact point p2.
    n2 = normalize(p2 - O2);

    // Eq (23) & (24): Tangential vectors (rotate 90 deg CCW assuming standard axes)
    n1_t = rotate90(n1);
    n2_t = rotate90(n2);
}

// Helper functions for computeContactVectors()

Vec2D computeUnitVector(const Vec2D& p1, const Vec2D& p2) {
    Vec2D diff;
    diff.x = p1.x - p2.x;
    diff.y = p1.y - p2.y;
    double norm = sqrt(diff.x * diff.x + diff.y * diff.y);
    return { diff.x / norm, diff.y / norm };
}

// End of helper functions

size_t findClosestProfileIndex(const Vec2D& contact_point_local_mm, const SpurGear* gear) {
    if (!gear || gear->initial_position_x_left.empty()) {
        return 0; // Return 0 if error or empty profile
    }

    size_t n_points = gear->initial_position_x_left.size();
    size_t closest_j = 0;
    double min_dist_sq = numeric_limits<double>::max();

    for (size_t j = 0; j < n_points; ++j) {
        // Get point j from the initial profile (left flank)
        Vec2D profile_point_j(gear->initial_position_x_left[j], gear->initial_position_y_left[j]);
        // Calculate squared distance
        Vec2D diff = contact_point_local_mm - profile_point_j;
        double dist_sq = dot(diff, diff);

        if (dist_sq < min_dist_sq) {
            min_dist_sq = dist_sq;
            closest_j = j;
        }
    }
    return closest_j;
}

Vec2D computeContactForce(double delta, double delta_dot, double kv, double cd, double mu, double fF, const Vec2D& v_t, const Vec2D& n_g) {
    // Inputs:
    // delta:     Deformation (penetration) (m)
    // delta_dot: Normal relative velocity (m/s)
    // kv:        Mesh stiffness (N/m)
    // cd:        Mesh damping coefficient (Ns/m)
    // mu:        Coefficient of friction
    // fF:        Friction parameter (from computeFrictionParameter)
    // v_t:       Tangential relative velocity vector (m/s)
    // n_g:       Unit normal vector
    // t_hat:     Unit tangential vector
    // Returns:   Force vector ON PINION (N)

    // Normal Force Calculation
    double Fn_mag = kv * delta + cd * delta_dot;
    // Ensure contact force is not adhesive (cannot pull)
    if (Fn_mag < 0) {
        Fn_mag = 0;
    }
    Vec2D Fn_vec = n_g * Fn_mag; // Normal force ON PINION (along n_g)

    // Friction Force Calculation
    // Friction opposes relative tangential motion (v_t)
    Vec2D t_hat = normalize(v_t);
    Vec2D Ff_vec = t_hat * (-fF * mu * Fn_mag); // Friction force ON PINION

    // Total Force ON PINION
    return Fn_vec + Ff_vec; // N
}

double computeFrictionParameter(double vtMagnitude, double v0, double vL) {
    // Eq (34)
    vtMagnitude = abs(vtMagnitude);
    if (vtMagnitude <= v0)
        return 0.0;
    else if (vtMagnitude <= vL)
        return (vtMagnitude - v0) / (vL - v0);
    else
        return 1.0;
}

double getMeshStiffness(const SpurGear* pinion, const SpurGear* gear, const Vec2D& contact_point_pinion) {
    // TODO: Implement Time-Varying Mesh Stiffness calculation based on contact point position
    // For now, return a constant average value in SI units
    return MESH_STIFFNESS_SI;
}

double getMeshDampingCoeff(double kv, double zeta, double m_pinion, double m_gear) {
    double m_eff = 1.0 / (1.0 / m_pinion + 1.0 / m_gear);
    return 2.0 * zeta * sqrt(m_eff * kv); // Units: Ns/m if kv is N/m, m_eff is kg
}

double computeMoment(const Vec2D& contactPoint, const Vec2D& centerPoint, const Vec2D& forceVector) {
    Vec2D r = contactPoint - centerPoint; // Lever arm vector
    return r.x * forceVector.y - r.y * forceVector.x; // Units: m * N = Nm
}
