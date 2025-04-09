#ifndef FUNCTDEFS_H
#define FUNCTDEFS_H

#define _USE_MATH_DEFINES
#include <vector>
#include <string>
using namespace std;

struct ContactPoint {
    double x;
    double y;
};

typedef struct {
    double numTeeth;
    double gearRatio;
    double module;
    double pressureAngle;
    double outsideDiameter;
    double rootDiameter;
    double pitchDiameter;
    double baseDiameter;
    double centerDistance;
    double pitchRadius;
    double baseRadius;
    double outsideRadius;
    double faceWidth;
    double offsetFaceWidth;
    double transverseToothThickness;
    double diameterToothThickness;

    double youngsModulus;
    double poissonsRatio;
    double momentInertia;
    double mass;
    double elasticCoeff;
    double dampingCoeff;

    double theta_e;
    double theta_a;
    double theta_r;

    vector<double>  position_x_left;
    vector<double>  position_x_right;
    vector<double>  position_y_left;
    vector<double>  position_y_right;

    vector<double>  initial_position_x_left;
    vector<double>  initial_position_x_right;
    vector<double>  initial_position_y_left;
    vector<double>  initial_position_y_right;

    vector<vector<double>> initialToothProfileXLeft;
    vector<vector<double>> initialToothProfileYLeft;
    vector<vector<double>> initialToothProfileXRight;
    vector<vector<double>> initialToothProfileYRight;

    vector<vector<double>> currentToothProfileXLeft;
    vector<vector<double>> currentToothProfileYLeft;
    vector<vector<double>> currentToothProfileXRight;
    vector<vector<double>> currentToothProfileYRight;

    vector<ContactPoint> contactPoints;
} SpurGear;

typedef struct {
    int n;
    double phi_k;
    double deltaJudgmentThreshold;
    double trig_pinion;
    double trig_gear;
    double g;
} DynamicsVars;

typedef struct {
    // Positions: [x_p, y_p, phi_p, x_g, y_g, phi_g]
    // Velocities: [x_dot_p, y_dot_p, phi_dot_p, x_dot_g, y_dot_g, phi_dot_g]
    double pos[6];
    double vel[6];
} State;

struct Vec2D {
    double x;
    double y;
};

void startupmsg();
int getUserInput(SpurGear* pinion, SpurGear* gear);
double computeThetaA(double pressureAngle);
void computeInvolutePoint(double r_b, double u, double& x, double& y);
double compute_mu_k_from_r(double r, double r_base);
double computeInvoluteSlope(double baseRadius, double pitchRadius);
int setupProfiles(SpurGear* pinion, SpurGear* gear, DynamicsVars* dynamicvars);
int initialToothPinion(SpurGear* pinion, DynamicsVars* dynamicvars);
int initialToothGear(SpurGear* gear, DynamicsVars* dynamicvars);
int completePinionProfile(SpurGear* pinion, DynamicsVars* dynamicvars);
int completeGearProfile(SpurGear* gear, DynamicsVars* dynamicvars);
void moveGearPosition(const State& state, SpurGear* pinion, SpurGear* gear, DynamicsVars* dynamicvars);
void meshingLine(SpurGear* pinion, SpurGear* gear, DynamicsVars* dynamicvars);
State rungeKuttaStep( State& state, double t, double dt,  SpurGear& pinion,  SpurGear& gear);
State computeDerivatives(State& state, double t, SpurGear* pinion, SpurGear* gear, DynamicsVars* dynVars);
double solveF_Pinion_sx(SpurGear* pinion, State& state);
double solveF_Pinion_sy(SpurGear* pinion, State& state);
double solveF_Gear_sx(SpurGear* gear, State& state);
double solveF_Gear_sy(SpurGear* gear, State& state);
double F(double u, double r_b, double k, double x_J, double y_J);
double dFdu(double u, double r_b, double k, double x_J, double y_J, double h = 1e-6);
double solveForU(double r_b, double k, double x_J, double y_J, double initialGuess = 0.1, int maxIter = 100, double tol = 1e-6);
double computeFrictionParameter(double vtMagnitude, double v0, double vL);
Vec2D computeContactForce(double delta, double delta_dot, double kv, double cd, double mu, double fF, const Vec2D& v_t, const Vec2D& n, const Vec2D& n_t);
double computeMoment(const double contactX, const double contactY, const double centerX, const double centerY, const double forceX, const double forceY);
void computeContactVectors(const Vec2D& p1, const Vec2D& p2, const SpurGear& gear, Vec2D& n_g, Vec2D& n1, Vec2D& n1_t, Vec2D& n2, Vec2D& n2_t);
Vec2D computeUnitVector(const Vec2D& p1, const Vec2D& p2);
Vec2D rotate90(const Vec2D& v);
void outputGearProfile(const SpurGear* gear, const DynamicsVars* dynVars, const std::string& filename);



#endif
