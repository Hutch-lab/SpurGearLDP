#ifndef FUNCTDEFS_H
#define FUNCTDEFS_H
#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <tuple>
#include <fstream>
using namespace std;

// Simulation Constants
const double MESH_STIFFNESS_SI = 2.0e6; // N/m (replace with TVMS later)
const double MESH_DAMPING_RATIO = 0.15; // Dimensionless damping ratio (zeta) - typical range 0.03-0.1
const double FRICTION_COEFF = 0.06;     // Coefficient of friction (mu)
const double V0_FRICTION_SI = 0.001;    // m/s (Threshold velocity for friction activation)
const double VL_FRICTION_SI = 0.01;     // m/s (Velocity for full friction)
const double CONTACT_THRESHOLD_SI = 0.05 / 1000.0; // m (Max distance for contact check, 5 microns)
const double PENETRATION_THRESHOLD_SI = 0.0; // m (Max allowed penetration before force acts)

struct ContactPoint {
    double x;
    double y;
};

struct Vec2D { // I must admit, ChatGPT came up with this struct and it's been quite helpful
    double x, y;
    Vec2D() : x(0), y(0) {}
    Vec2D(double _x, double _y) :x(_x), y(_y) {}
    Vec2D operator+(const Vec2D& o)const { return{ x + o.x, y + o.y }; }
    Vec2D operator-(const Vec2D& o)const { return{ x - o.x, y - o.y }; }
    Vec2D operator*(double s)      const { return{ x * s,   y * s }; }
    Vec2D operator-() const { return { -x, -y }; }
};
static double dot(const Vec2D& a, const Vec2D& b) { return a.x * b.x + a.y * b.y; }
static double mag(const Vec2D& v) { return std::sqrt(dot(v, v)); }
static Vec2D normalize(const Vec2D& v) {
    double m = mag(v);
    return m > 1e-12 ? Vec2D{ v.x / m,v.y / m } : Vec2D{ 1,0 };
}
static Vec2D rotate(const Vec2D& v, double angle_rad) {
    double cos_a = cos(angle_rad);
    double sin_a = sin(angle_rad);
    return { v.x * cos_a - v.y * sin_a, v.x * sin_a + v.y * cos_a };
}
static Vec2D rotate90(const Vec2D& v) {
    return { -v.y, v.x };
}

typedef struct {
    double torque; // Input torque (Nm)
    double rpmSpeed; // Input speed (rpm)
    double numTeeth;
    double gearRatio;
    double module; // (mm)
    double pressureAngle; // (degrees)
    double outsideDiameter; // (mm)
    double rootDiameter; // (mm)
    double pitchDiameter; // (mm)
    double baseDiameter; // (mm)
    double centerDistance; // (mm)
    double pitchRadius; // (mm)
    double baseRadius; // (mm)
    double outsideRadius; // (mm)
    double faceWidth; // (mm)
    double offsetFaceWidth; // (mm)
    double transverseToothThickness; // (mm)
    double diameterToothThickness; // (mm)
    double addendum; // (mm)
    double dedendum; // (mm)

    double youngsModulus; // (MPa)
    double poissonsRatio;
    double momentInertia; // (kg*m^2)
    double mass; // (kg)
    double elasticCoeff; // Support stiffness
    double dampingCoeff; // Support damping

    double theta_e; // Angular pitch (radians)
    double theta_a; // Involute roll angle to pitch point (radians)
    double theta_r; // Rotation angle for initial tooth setup (radians)

    // Vectors defining the geometry of ONE tooth in its base orientation
    vector<double>  position_x_left;
    vector<double>  position_x_right;
    vector<double>  position_y_left;
    vector<double>  position_y_right;
    // Vectors defining the geometry of ONE tooth in its initial rotated orientation
    vector<double>  initial_position_x_left;
    vector<double>  initial_position_x_right;
    vector<double>  initial_position_y_left;
    vector<double>  initial_position_y_right;
    // Vectors defining the geometry of ALL teeth in their initial positions
    vector<vector<double>> initialToothProfileXLeft;
    vector<vector<double>> initialToothProfileYLeft;
    vector<vector<double>> initialToothProfileXRight;
    vector<vector<double>> initialToothProfileYRight;
    // Vectors defining the geometry of ALL teeth in their CURRENT dynamic positions
    vector<vector<double>> currentToothProfileXLeft;
    vector<vector<double>> currentToothProfileYLeft;
    vector<vector<double>> currentToothProfileXRight;
    vector<vector<double>> currentToothProfileYRight;

    // Parameters needed for dynamic stress calculation
    int Qv;             // Transmission Accuracy level (for Kv)
    int gearQuality;    // Gear Quality (for KH)
    double Ko;          // Overload Factor
    double Ks;   				// Size Factor
    double ZR;   				// Surface Condition Factor


    vector<double> contactStressProfile;

    vector<pair<Vec2D, Vec2D>> activeContactPairs; // Store pairs (pinion_pt, gear_pt)
} SpurGear;

typedef struct DynamicsVars {
    int n; // Number of points per involute profile
    double phi_k; // Current rotation angle
    double deltaJudgmentThreshold; // Threshold for detecting contact points
    double trig_pinion; // Rotation angle for generating right flank from left
    double trig_gear; // Rotation angle for generating right flank from left
    double g;
    double meshDampingCoeff; // Ns/m
    double meshStiffness;
    double frictionCoeff; // mu
    double v0_friction; // m/s
    double vL_friction; // m/s
    long long stepCounter;
} DynamicsVars;

typedef struct {
    // Positions: [x_p, y_p, phi_p, x_g, y_g, phi_g]
    // Velocities: [x_dot_p, y_dot_p, phi_dot_p, x_dot_g, y_dot_g, phi_dot_g]
    double pos[6];
    double vel[6];
} State;

void startupmsg();
int getUserInput(SpurGear* pinion, SpurGear* gear);
double computeThetaA(double pressureAngle);
void computeInvolutePoint(double r_b, double mu_k, double& x, double& y);
double compute_mu_k_from_r(double r, double r_base);
double computeInvoluteSlope(double baseRadius, double pitchRadius);
int setupProfiles(SpurGear* pinion, SpurGear* gear, DynamicsVars* dynamicvars);
void initialToothPinion(SpurGear* pinion, DynamicsVars* dynamicvars);
void initialToothGear(SpurGear* gear, DynamicsVars* dynamicvars);
void completePinionProfile(SpurGear* pinion, DynamicsVars* dynamicvars);
void completeGearProfile(SpurGear* gear, DynamicsVars* dynamicvars);
void moveGearPosition(const State& state, SpurGear* pinion, SpurGear* gear, DynamicsVars* dynamicvars);
//void meshingLine(SpurGear* pinion, SpurGear* gear, DynamicsVars* dynamicvars);
State rungeKuttaStep(const State& state, double t, double dt, SpurGear& pinion, SpurGear& gear, DynamicsVars& dynvars);
State computeDerivatives(const State& state, double t, SpurGear* pinion, SpurGear* gear, DynamicsVars* dynvars);
//double solveF_Pinion_sx(SpurGear* pinion, State& state);
//double solveF_Pinion_sy(SpurGear* pinion, State& state);
//double solveF_Gear_sx(SpurGear* gear, State& state);
//double solveF_Gear_sy(SpurGear* gear, State& state);
//double F(double u, double r_b, double k, double x_J, double y_J);
//double dFdu(double u, double r_b, double k, double x_J, double y_J, double h = 1e-6);
//double solveForU(double r_b, double k, double x_J, double y_J, double initialGuess = 0.1, int maxIter = 100, double tol = 1e-6);
double computeFrictionParameter(double vtMagnitude, double v0, double vL);
Vec2D computeContactForce(double delta, double delta_dot, double kv, double cd, double mu, double fF, const Vec2D& v_t, const Vec2D& n_g);
double computeMoment(const Vec2D& contactPoint, const Vec2D& centerPoint, const Vec2D& forceVector);
void computeContactVectors(const Vec2D& p1, const Vec2D& p2, const Vec2D& O1, const Vec2D& O2, Vec2D& n_g, Vec2D& n1, Vec2D& n1_t, Vec2D& n2, Vec2D& n2_t);
Vec2D computeUnitVector(const Vec2D& p1, const Vec2D& p2);
// void outputGearProfile(const SpurGear* gear, const DynamicsVars* dynVars, const std::string& filename);

// Dynamics simulation functions
void dynamicMeshingUpdate(const State& state, SpurGear* pinion, SpurGear* gear, DynamicsVars* dynvars, Vec2D& total_pinion_force, Vec2D& total_gear_force, double& total_pinion_torque, double& total_gear_torque);
double getMeshStiffness(const SpurGear* pinion, const SpurGear* gear, const Vec2D& contact_point_pinion);
double getMeshDampingCoeff(double kv, double zeta, double m_pinion, double m_gear);
size_t findClosestProfileIndex(const Vec2D& contact_point_local_mm, const SpurGear* gear);

double calculateZe(const SpurGear* pinion, const SpurGear* gear);
double calculateZI(const SpurGear* pinion, const SpurGear* gear);
double calculateKv(double V, int Qv);
double calculateKH(const SpurGear* pinion, const SpurGear* gear, double F, int gearQuality);
double calcContactStress(const SpurGear* pinion, const SpurGear* gear, double Wt, double V, double Ko, int Qv, int gearQuality, double Ks, double ZR);

#endif
