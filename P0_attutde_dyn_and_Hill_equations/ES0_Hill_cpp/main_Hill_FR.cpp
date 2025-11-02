#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>

class LTIPlant {
private:
    Eigen::MatrixXd A, B;   
    Eigen::VectorXd x;     

public:
    LTIPlant(const Eigen::MatrixXd& A_, const Eigen::MatrixXd& B_, const Eigen::VectorXd& x0)
        : A(A_), B(B_), x(x0) {}

    Eigen::VectorXd getState() const { return x; }

    void setState(const Eigen::VectorXd& x0) { x = x0; }

    Eigen::VectorXd xdot(const Eigen::VectorXd& x_local, const Eigen::VectorXd& u) const {
        return A * x_local + B * u;
    }

    void step_rk4(double dt, const Eigen::VectorXd& u) {
        Eigen::VectorXd k1 = xdot(x, u);
        Eigen::VectorXd k2 = xdot(x + 0.5 * dt * k1, u);
        Eigen::VectorXd k3 = xdot(x + 0.5 * dt * k2, u);
        Eigen::VectorXd k4 = xdot(x + dt * k3, u);
        x += (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
    }
};




class Simulation {
private:
    LTIPlant& plant;
    double T;
    double dt;
    std::string filename;

public:
    Simulation(LTIPlant& plant_, double T_, double dt_, const std::string& file_)
        : plant(plant_), T(T_), dt(dt_), filename(file_) {}

    void run() {
        int steps = static_cast<int>(T / dt);
        int nx = plant.getState().size();

        std::ofstream fout(filename);
        if (!fout.is_open()) {
            std::cerr << "Error: cannot open file " << filename << "\n";
            return;
        }

        Eigen::VectorXd u = Eigen::VectorXd::Zero(3);
        int nu = u.size();

        fout << "t";
        for (int i = 0; i < nx; ++i) fout << ",x" << i;
        for (int i = 0; i < nu; ++i) fout << ",u" << i;

        fout << "\n";

        // Simulation loop
        for (int k = 0; k <= steps; ++k) {
            double t = k * dt;
            Eigen::VectorXd x = plant.getState();



            fout << t;
            for (int i = 0; i < nx; ++i) fout << "," << x(i);
            for (int i = 0; i < nu; ++i) fout << "," << u(i);
            fout << "\n";

            plant.step_rk4(dt, u);
        }

        fout.close();
        std::cout << "Simulation done. Data written to " << filename << "\n";
    }
};


int main() {
    //Hill equation simulation (free response)
    const double omega_orb_radps = 0.001259; 
    const double m_chaser_kg = 100; 
    // for now let's assume the mass is constant (no propulsion in free response)
    
    Eigen::MatrixXd A(6,6);
    A.setZero();
    A(0,3) = 1.0;  A(1,4) = 1.0;  A(2,5) = 1.0;
    A(3,5) =  2.0 * omega_orb_radps;
    A(4,1) = -omega_orb_radps * omega_orb_radps;
    A(5,3) =  -2.0 * omega_orb_radps;
    A(5,2) = 3.0 * omega_orb_radps * omega_orb_radps;

    // ------------------ B matrix (6x3) ------------------
    Eigen::MatrixXd B(6,3);
    B.setZero();
    B(3,0) = 1.0 / m_chaser_kg;
    B(4,1) = 1.0 / m_chaser_kg;
    B(5,2) = 1.0 / m_chaser_kg;

    // Initial condition
    Eigen::Vector3d r0_rel_LVLH_m(50.0, -20.0, 10.0);      
    Eigen::Vector3d v0_rel_LVLH_mps(0.02, 0.00, -0.01);    

    Eigen::VectorXd x0(6);
    x0 << r0_rel_LVLH_m, v0_rel_LVLH_mps;

    LTIPlant plant(A, B, x0);

    Simulation sim(plant, 15000.0, 0.1, "Hill_eq_Free_Response.csv");
    sim.run();

    return 0;
}
