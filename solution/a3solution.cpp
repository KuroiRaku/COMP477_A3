#include "a3solution.h"

#include "dependencies/Eigen/Dense"
#include "QDebug"
#include "QElapsedTimer"

using Eigen::Vector2f;
using Eigen::Vector4f;
using Eigen::VectorXf;
using Eigen::MatrixXd;
using Eigen::MatrixXf;

A3Solution::A3Solution(std::vector<Joint2D*>& joints, std::vector<Spring2D*>& springs, float& gravity, float& positional_damping, float& mass, float& timestep, bool& implicit, float& stiffness)
    :m_joints(joints),
    m_links(springs),
    m_gravity(gravity),
    m_positional_damping(positional_damping),
    m_mass(mass),
    m_timestep(timestep),
    m_implicit(implicit),
    m_stiffness(stiffness)
    {

    this->isInitialized = false;
}


void A3Solution::update(Joint2D* selected, QVector2D mouse_pos){
    if (selected->is_locked()) {
        return;
    }

    this->selected = selected;
    this->mousePos = mouse_pos;
    selected->set_position(mouse_pos);
}

void A3Solution::update(){

    if (this->checkIfSystemUpdated()) {
        this->isInitialized = false;
    }

    if (!isInitialized) {
        this->initializeYk();
    }

    for (int i=0; i<m_moving_joints.size(); i++) {
        Joint2D* joint = m_moving_joints[i];

        // if the joint is the selected joint, skip and set selectedIndex = i
        if (joint == this->selected) {
            this->selectedIndex = i;
            continue;
        }

        // F = Kx, look at the function to understand,
        Vector2f springForces = calculateSpringForces(joint, joint->get_springs());

        // Damping forces changes every second.
        // However, positional damping forces were given so merely need to multiple with the velocity
        Vector2f dampingForces = -(m_positional_damping * Vector2f(m_yk[i*4+xVelocity], m_yk[i*4+yVelocity]));

        // compute forces for all joints and set
        // don't forget gravity forces acting downwards
        // gravity force = m*g
        Vector2f gravityForces = - (m_mass * m_gravity * Vector2f::UnitY());

        Vector2f totalForces = gravityForces + springForces + dampingForces;
        // f = ma
        // a = f/m
        Vector2f totalAcceleration = totalForces / this->m_mass ;

        m_ykPrime[i*4 +xAccelerationPrime] = totalAcceleration.x();
        m_ykPrime[i*4 +yAccelerationPrime] = totalAcceleration.y();
    }

    // do Explict Euler
    explicitEuler(m_yk, m_ykPrime);

    // update every positions in Ui
    for (int i=0; i<m_moving_joints.size(); ++i) {

        // on the joints that selected by mouse, continue
        if (i >= 0 && i == selectedIndex) {
            selected = nullptr;
            selectedIndex = -1;
            continue;
        }

        m_moving_joints[i]->set_position(QVector2D(m_yk[i*4+xPosition],-m_yk[i*4+yPosition]));
    }
}


bool A3Solution::checkIfSystemUpdated() {

    // The jointCount will be the newest join count while m_joints is the one holding the previous joints
    if (this->jointCount != this->m_joints.size() || this->allCurrentJoints != this->m_joints) {
        this->allCurrentJoints.clear();
        for (Joint2D* joint : this->m_joints) {
            this->allCurrentJoints.push_back(joint);
        }

        this->jointCount = this->m_joints.size();
        return true;
    }

    return false;
}

void A3Solution::explicitEuler(VectorXf& yk, VectorXf& ykPrime){

    if (!m_implicit){
        // Pretty much explicit Euler lol.
        yk = yk + (ykPrime * this->m_timestep);

        // need to update ykPrime too
        for (int i=0; i<this->m_moving_joints.size(); ++i) {
            if (i == selectedIndex) {
                continue;
            }

            ykPrime[i*4+xVelocityPrime] = yk[i*4+xVelocity];
            ykPrime[i*4+yVelocityPrime] = yk[i*4+yVelocity];
        }
    }
    // For Implicit
    else
    {

    }

    // selection of mouse case
    if (selectedIndex >= 0) {
        int index = selectedIndex * 4;

        // We only setting the position because when we are moving the joints
        // velocity and acceleration of that spring supposed to be zero
        yk[index + xPosition] = this->mousePos.x();
        yk[index + yPosition] = -(this->mousePos.y());

        yk[index +xVelocity] = 0.0f;
        yk[index +yVelocity] = 0.0f;
        ykPrime[index +xVelocityPrime] = 0.0f;
        ykPrime[index +yVelocityPrime] = 0.0f;

        ykPrime[index +xAccelerationPrime] = 0.0f;
        ykPrime[index +yAccelerationPrime] = 0.0f;
    }
}

// Initialize everything needed for VK and joints
void A3Solution::initializeYk(){
    // set all moveable joints
    this->m_moving_joints.clear();

    for (Joint2D* joint : this->m_joints) {
        if (!joint->is_locked()) {
            this->m_moving_joints.push_back(joint);
        }
    }

    // now set yk
    int ykSize = 4 * this->m_moving_joints.size();
    VectorXf yk(ykSize);
    VectorXf ykPrime(ykSize);

    for (int i=0; i<ykSize/4; i++) {
        int index = i*4;
        Joint2D* joint = m_moving_joints[i];

        yk[index+xPosition] = joint->get_position().x();
        yk[index+yPosition] = -joint->get_position().y();

        if (!isInitialized) {
            yk[index+xVelocity] = 0.0f;
            yk[index+yVelocity] = 0.0f;
        }

        // 0, 1
        ykPrime[index+xVelocityPrime] = yk[index+xVelocity];
        ykPrime[index+yVelocityPrime] = yk[index+yVelocity];

        // 2, 3
        ykPrime[index+xAccelerationPrime] = 0.0f;
        ykPrime[index+yAccelerationPrime] = 0.0f;
    }

    this->m_yk = yk;
    this->m_ykPrime = ykPrime;
    isInitialized = true;
}

Vector2f A3Solution::calculateSpringForces(Joint2D* joint, std::vector<Spring2D*> springs) {
    // initializing the variable needed.
    Vector2f totalSpringForce = Vector2f(0.0f,0.0f);

    for(Spring2D* spring : springs) {
        Joint2D* other = spring->get_other_joint(joint);
        // vector to the joint being pulled
        Vector2f vecToJoint = Vector2f(joint->get_position().x() - other->get_position().x(),-(joint->get_position().y() - other->get_position().y()) );
        // Hooke's Law
        // F = kx
        // X = current length - resting length
        Vector2f force = -1 * this->m_stiffness * (spring->get_length()-spring->get_rest_length()) * vecToJoint.normalized();

        totalSpringForce.x() += force.x();
        totalSpringForce.y() += force.y();
    }

    return totalSpringForce;
}

void A3Solution::test_eigen_library(){

    // create a simple matrix 5 by 6
    MatrixXd mat(5,6);

    // Fills in matrix
    // Important Note: Eigen matrices are row major
    // so mat(0,1) references the 0-th column and 1-th row
    for(unsigned int row=0;row<mat.rows();row++){
        for(unsigned int col=0;col<mat.cols();col++){
            mat(row,col) = row+col;
        }
    }

    // create the pseudoinverse
    MatrixXd pseudo_inv = mat.completeOrthogonalDecomposition().pseudoInverse();

    // print the pseudoinverse
    std::cout << "--------------------------" << std::endl;
    std::cout << pseudo_inv << std::endl;
    std::cout << "--------------------------" << std::endl;

}
