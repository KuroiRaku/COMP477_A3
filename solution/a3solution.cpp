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
    isUptoDate = false;
    }


void A3Solution::update(Joint2D* selected, QVector2D mouse_pos){

    if (selected->is_locked())
    {
        return;
    }

    this->mousePos = mouse_pos;
    this->selected = selected;
    selected->set_position(mouse_pos);
}

// This function will only get called when in FileMode = 4 which is the Demo.
void A3Solution::update(){

    if (checkIfSystemUpdated())
    {
        isUptoDate = false;
    }

    if (!isUptoDate)
    {
        initialiseYK();
    }

    qInfo()<< "SelectedIndex: " <<selected;


    for(int i=0 ; i > this->m_moving_joints.size(); ++i)
    {
        Joint2D* joint = m_moving_joints[i];

        if (joint == this->selected) {
            this->selectedIndex = i;
            continue;
        }

        // compute forces for all joints and set
        Vector2f GravityForce = - (Vector2f::UnitY() * this->m_mass * m_gravity);

        Vector2f springForce = this->calculateSpringForce(joint, joint->get_springs());

        // positional damping are calculated in other places.
        Vector2f dampForce = this->calculateDampingForce(Vector2f(m_yk[i*4 + xVelocity], m_yk[i*4 + yVelocity]));

        Vector2f totalForce = GravityForce + springForce + dampForce;
        // F = ma
        // a = F/m
        Vector2f totalAccel = totalForce / m_mass ;

        m_yk_prime[i*4 + xAccelerationPrime ] = totalAccel.x();
        m_yk_prime[i*4 + yAccelerationPrime ] = totalAccel.y();
    }

    this->explicitEuler(m_yk,m_yk_prime);

    // Updating everything in UI
    for (int i=0; i<m_moving_joints.size(); ++i) {

        if (i >= 0 && i == this->selectedIndex) {
            this->selected = nullptr;
            this->selectedIndex = -1;
            continue;
        }

        // qInfo()<< "Update: " <<m_yk[i*4 + xPosition] << " " <<m_yk[i*4 + yPosition];
        m_moving_joints[i]->set_position(QVector2D(m_yk[i*4 + xPosition], - (m_yk[i*4 + yPosition])));
    }
}

// Initializing the VK
void A3Solution::initialiseYK(){
    // set all moveable joints
    this->m_moving_joints.clear();

    qInfo()<< "m_joints: " << m_joints.size();
    for (Joint2D* joint : this->m_joints) {
        if (!joint->is_locked()) {
            this->m_moving_joints.push_back(joint);
        }
    }

    // now set yk
    int yk_length = 4*this->m_moving_joints.size();

    VectorXf yk(yk_length);
    VectorXf yk_prime(yk_length);

    for (int i=0; i<yk_length; i=i+4) {
        int jointIndex = i/4;
        Joint2D* joint = m_moving_joints[jointIndex];
        Vector2f jointMathPos = Vector2f(joint->get_position().x(),-(joint->get_position().y()));

        yk[i+xPosition] = jointMathPos.x();
        yk[i+yPosition] = jointMathPos.y();

        if (!isUptoDate) {
            yk[i+xVelocity] = 0.0f;
            yk[i+yVelocity] = 0.0f;
        }

        yk_prime[i+xVelocityPrime] = yk[i+xVelocity];
        yk_prime[i+yVelocityPrime] = yk[i+yVelocity];

        yk_prime[i+xAccelerationPrime] = 0.0f;
        yk_prime[i+yAccelerationPrime] = 0.0f;
    }

    this->m_yk = yk;
    this->m_yk_prime = yk_prime;
    isUptoDate = true;
}

Vector2f A3Solution::calculateDampingForce(Vector2f velocity) {
    return -1 * m_positional_damping * velocity;
}

Vector2f A3Solution::calculateSpringForce(Joint2D* joint, std::vector<Spring2D*> springs)
{
    std::vector<Vector2f> forces;
    Vector2f force_spring_total = Vector2f(0.0f,0.0f);

    for(Spring2D* spring : springs) {
        float currL = spring->get_length();
        float restL = spring->get_rest_length();

        Joint2D* other = spring->get_other_joint(joint);
        QVector2D temp = joint->get_position() - other->get_position();
        Vector2f vec_to_joint = Vector2f(temp.x(), - temp.y());

        // Hooke's law
        // F = kx
        Vector2f force = -1 * this->m_stiffness * (currL-restL) * vec_to_joint.normalized();
            forces.push_back(force);
        }

        for (Vector2f sForce : forces) {
            force_spring_total.x() += sForce.x();
            force_spring_total.y() += sForce.y();
        }

        return force_spring_total;
}

bool A3Solution::checkIfSystemUpdated(){

    if (!isUptoDate) {
        this->allCurrentJoints.clear();
        for (Joint2D* joint : this->m_joints) {
            this->allCurrentJoints.push_back(joint);
        }
        this->jointCount = this->m_joints.size();
        return true;
    }

    // if the size is different or the joints are different, remove everything and add everything again.
    // Not optimized but yeah
    if (jointCount != m_joints.size() || this->allCurrentJoints != this->m_joints) {
        this->allCurrentJoints.clear();
        for (Joint2D* joint : this->m_joints) {
            this->allCurrentJoints.push_back(joint);
        }
        this->jointCount = this->m_joints.size();

        return true;
    }

    return false;
}

void A3Solution::explicitEuler(VectorXf& vectorYK, VectorXf& vectorYKPrime)
{
    qInfo()<< "Explcit Euler " << m_moving_joints.size();

    vectorYK = vectorYK + (vectorYKPrime * this->m_timestep);

    // update VK prime
    for (int i=0; i<this->m_moving_joints.size(); ++i) {
        if (i == selectedIndex) {
            continue;
        }

        int index = i*4 ;
        vectorYKPrime[index + xVelocityPrime]=vectorYK[index+xVelocity];
        vectorYKPrime[index + yVelocityPrime]=vectorYK[index+yVelocity];
    }

    // if there is selection of the mouse.
    if (selectedIndex >= 0) {

        qInfo()<< "Mouse selected joints";

        // *4 because each vector has 4 elements
        int index = selectedIndex * 4;
        // Set Position

        vectorYK[index+xPosition] = mousePos.x();

        // qt and eigen math are different
        vectorYK[index+yPosition] = - mousePos.y();

        // The velocity and the acceleration have to be recalculated
        // Setting velocity to 0 while during selection
        vectorYK[index+xVelocity] = 0.0f;
        vectorYK[index+yVelocity] = 0.0f;
        vectorYKPrime[index + xVelocityPrime]= 0.0f;
        vectorYKPrime[index + yVelocityPrime]= 0.0f;

        // setting acceleration to 0 while during selection
        vectorYKPrime[index + xAccelerationPrime] = 0.0f;
        vectorYKPrime[index + yAccelerationPrime] = 0.0f;
    }
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
