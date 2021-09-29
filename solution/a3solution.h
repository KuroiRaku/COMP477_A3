#ifndef A2SOLUTION_H
#define A2SOLUTION_H

#include <vector>
#include <algorithm>

#include "OpenGL/elements/joint2D.h"
#include "OpenGL/elements/spring2d.h"
#include "dependencies/Eigen/Dense"

using Eigen::Vector2f;
using Eigen::VectorXf;
using Eigen::MatrixXf;

class A3Solution
{
public:
    A3Solution(std::vector<Joint2D*>& joints, std::vector<Spring2D*>& springs, float& gravity, float& positional_damping, float& mass, float& timestep, bool& implicit, float& stiffness);

    // OpenGL members (these are updated for you in the default solution)
    std::vector<Joint2D*>& m_joints;
    std::vector<Spring2D*>& m_links;
    float& m_gravity;
    float& m_positional_damping;
    float& m_mass;
    float& m_timestep;
    bool& m_implicit;
    float& m_stiffness;

    // Separate tracking of positions and velocities
    void update();
    void update(Joint2D* selected, QVector2D mouse_pos);
    void explicitEuler(VectorXf& vectorYK, VectorXf& vectorYKPrime);

    static void test_eigen_library();

    // Data
    std::vector<Joint2D*> m_moving_joints;
    VectorXf m_yk;
    VectorXf m_yk_prime;
    bool isUptoDate = false;

    // Have to use 2f since there are two directions x,y
    Vector2f totalForce;
    Vector2f dampingForce;


    std::vector<Joint2D*> allCurrentJoints;
    int jointCount = -1;

    Joint2D* selected;
    QVector2D mousePos;
    int selectedIndex = -1;

private:
    bool checkIfSystemUpdated();

    Vector2f calculateSpringForce(Joint2D* joint, std::vector<Spring2D*> springs);

    Vector2f calculateDampingForce(Vector2f velocity);

    void initialiseYK();

    // Therer are four numbers in the constant are used to indicate the position in Eigen:VectorXF
    // In one vectorXF they would have xposition, yPosition, xVector, yVector
    // Constants
    int xPosition = 0;
    int yPosition = 1;
    int xVelocity = 2;
    int yVelocity = 3;

    // Although we can reuse the constant above but to prevent confusion, I decided to create another one.
    int xVelocityPrime = 0;
    int yVelocityPrime = 1;
    int xAccelerationPrime = 2;
    int yAccelerationPrime = 3;
};

#endif // A2SOLUTION_H
