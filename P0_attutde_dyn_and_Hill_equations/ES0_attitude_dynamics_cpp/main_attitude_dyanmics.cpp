#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>

class ScAttitudePlant {
private:
    Eigen::Matrix3d I_B_kgm2, I_B_inv_kgm2;
    Eigen::VectorXd x; // [qscalar, qv_x, qv_y, qv_z, omega_x, omega_y, omega_z]

public:
    ScAttitudePlant(const Eigen::Matrix3d& I_in, const Eigen::VectorXd& x0)
        : I_B_kgm2(I_in), I_B_inv_kgm2(I_in.inverse()), x(x0) {}

    Eigen::VectorXd getState() const { return x; }

    Eigen::Matrix3d quatToDCM_NB(const Eigen::Quaterniond& q) const {
        // Rotation matrix from the body frame (B) to the inertial frame (N)
        return q.toRotationMatrix();
    }

    Eigen::VectorXd dynamics(const Eigen::VectorXd& x_local, const Eigen::Vector3d& M_ext_B_Nm) const {
        Eigen::Quaterniond q_NB(x_local(0), x_local(1), x_local(2), x_local(3)); // (w, x, y, z)
        Eigen::Vector3d omega_BN_B_radps = x_local.segment<3>(4);

        // Euler rotational dynamics
        Eigen::Vector3d omega_dot_BN_B_radps2 =
            I_B_inv_kgm2 * (M_ext_B_Nm - omega_BN_B_radps.cross(I_B_kgm2 * omega_BN_B_radps));

        // Quaternion Kinematics
        Eigen::Vector4d qvec(x_local(0), x_local(1), x_local(2), x_local(3));
        const double wx = omega_BN_B_radps.x();
        const double wy = omega_BN_B_radps.y();
        const double wz = omega_BN_B_radps.z();

        Eigen::Matrix4d W;
        W <<  0.0,  -wx,  -wy,  -wz,
              wx,   0.0,   wz,  -wy,
              wy,  -wz,   0.0,   wx,
              wz,   wy,  -wx,   0.0;

        Eigen::Vector4d qdot_vec = 0.5 * W * qvec;

        Eigen::VectorXd xdot(7);
        xdot(0) = qdot_vec(0);
        xdot(1) = qdot_vec(1);
        xdot(2) = qdot_vec(2);
        xdot(3) = qdot_vec(3);
        xdot.segment<3>(4) = omega_dot_BN_B_radps2;

        return xdot;
    }

    void step_rk4(double dt, const Eigen::Vector3d& M_ext_B_Nm) {
        Eigen::VectorXd k1 = dynamics(x, M_ext_B_Nm);
        Eigen::VectorXd k2 = dynamics(x + 0.5 * dt * k1, M_ext_B_Nm);
        Eigen::VectorXd k3 = dynamics(x + 0.5 * dt * k2, M_ext_B_Nm);
        Eigen::VectorXd k4 = dynamics(x + dt * k3, M_ext_B_Nm);
        x += (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);

        Eigen::Quaterniond q_NB(x(0), x(1), x(2), x(3));
        q_NB.normalize();
        x(0) = q_NB.w();
        x(1) = q_NB.x();
        x(2) = q_NB.y();
        x(3) = q_NB.z();
    }

    Eigen::Vector3d getAngularMomentum_B() const {
        Eigen::Vector3d omega_BN_B_radps = x.segment<3>(4);
        return I_B_kgm2 * omega_BN_B_radps;
    }

    Eigen::Vector3d getAngularMomentum_I() const {
        Eigen::Quaterniond q_NB(x(0), x(1), x(2), x(3));
        Eigen::Matrix3d R_NB = quatToDCM_NB(q_NB);
        return R_NB * getAngularMomentum_B();
    }

    double getKineticEnergy_J() const {
        Eigen::Vector3d omega_BN_B_radps = x.segment<3>(4);
        return 0.5 * omega_BN_B_radps.transpose() * I_B_kgm2 * omega_BN_B_radps;
    }
};


class Simulation {
private:
    ScAttitudePlant& plant;
    double T_sim_s;
    double dt_s;
    std::string filename;

public:
    Simulation(ScAttitudePlant& plant_, double T_sim_s_, double dt_s_, const std::string& file_)
        : plant(plant_), T_sim_s(T_sim_s_), dt_s(dt_s_), filename(file_) {}

    void run() {
        int steps = static_cast<int>(T_sim_s / dt_s);
        int nx = plant.getState().size();

        std::ofstream fout(filename);
        if (!fout.is_open()) {
            std::cerr << "Error: cannot open file " << filename << "\n";
            return;
        }

        fout << "t_s,"
             << "q_NB_w,q_NB_x,q_NB_y,q_NB_z,"
             << "omega_BN_B_x_radps,omega_BN_B_y_radps,omega_BN_B_z_radps,"
             << "H_B_x_Nms,H_B_y_Nms,H_B_z_Nms,"
             << "H_I_x_Nms,H_I_y_Nms,H_I_z_Nms,"
             << "H_I_mag_Nms,KE_rot_J\n";

        for (int k = 0; k <= steps; ++k) {
            double t_s = k * dt_s;

            Eigen::Vector3d M_ext_B_Nm = Eigen::Vector3d::Zero();

            Eigen::VectorXd x = plant.getState();
            Eigen::Vector3d H_B_Nms = plant.getAngularMomentum_B();
            Eigen::Vector3d H_I_Nms = plant.getAngularMomentum_I();
            double H_I_mag_Nms = H_I_Nms.norm();
            double KE_rot_J = plant.getKineticEnergy_J();

            fout << t_s;
            for (int i = 0; i < nx; ++i) fout << "," << x(i);
            fout << "," << H_B_Nms(0) << "," << H_B_Nms(1) << "," << H_B_Nms(2);
            fout << "," << H_I_Nms(0) << "," << H_I_Nms(1) << "," << H_I_Nms(2);
            fout << "," << H_I_mag_Nms << "," << KE_rot_J << "\n";

            plant.step_rk4(dt_s, M_ext_B_Nm);
        }

        fout.close();
        std::cout << "Simulation done. Data written to " << filename << "\n";
    }
};


int main() {

    Eigen::Matrix3d I_B_kgm2 = Eigen::Matrix3d::Zero();
    I_B_kgm2(0,0) = 10.0;
    I_B_kgm2(1,1) = 15.0;
    I_B_kgm2(2,2) = 20.0;

    Eigen::Vector4d q_NB_Init = 0.5 * Eigen::Vector4d::Ones();  // [qs, qv_x, qv_y, qv_z]
    Eigen::Vector3d omega_BN_B_Init_radps(0.0, 0.1, 0.2);

    Eigen::VectorXd x0(7);
    x0 << q_NB_Init, omega_BN_B_Init_radps;

    ScAttitudePlant plant(I_B_kgm2, x0);

    double T_sim_s = 60.0;
    double dt_s = 0.01;

    Simulation sim(plant, T_sim_s, dt_s, "Attitude_Dynamics_Full_Sim.csv");
    sim.run();

    return 0;
}
