#ifndef functdefs.h
#define functdefs.h
#include <vector>
using namespace std;

typedef struct {
    int numTeeth;
    double gearRatio;
    double module;
    double pressureAngle;
    double outsideDiameter;
    double rootDiameter;
    double pitchDiameter;
    double baseDiameter
    double centerDistance;
    double faceWidth;
    double offsetFaceWidth;
    double transverseToothThickness;
    double diameterToothThickness;

    double youngsModulus;
    double poissonsRatio;
    double momentInertia;
    double mass;
    double elasticCoeff;
    double dampingCoeff

    double theta_e;
    double theta_a;
    double theta_r;

    vector<ContactPoint> contactPoints;

    vector<double>  position_x_left[];
    vector<double>  position_x_right[];
    vector<double>  position_y_left[];
    vector<double>  position_y_right[];

    vector<double>  initial_position_x_left[];
    vector<double>  initial_position_x_right[];
    vector<double>  initial_position_y_left[];
    vector<double>  initial_position_y_right[];

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
} DynamicsVars;

typedef struct {
    // Positions: [x_p, y_p, phi_p, x_g, y_g, phi_g]
    // Velocities: [x_dot_p, y_dot_p, phi_dot_p, x_dot_g, y_dot_g, phi_dot_g]
    double pos[6];
    double vel[6];
} State;
int getUserInput(SpurGear* pinion, SpurGear* gear);

struct Vec2D {
    double x;
    double y;
};

#endif
