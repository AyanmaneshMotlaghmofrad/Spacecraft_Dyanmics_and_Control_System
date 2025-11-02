# Spacecraft_Dynamics_and_Control_System

Projects and assignments from the course *Dinamica e Controllo di Veicoli Spaziali* at Politecnico di Torino.  
This is an ongoing project that will be progressively developed and completed throughout the semester.

---

## P0: Simulation of Spacecraft Attitude Dynamics and Hill’s Equations

This project focuses on simulating:
1. The **rotational dynamics** of a spacecraft (attitude)  
2. The **relative translational motion** using **Hill’s (Clohessy–Wiltshire) equations** for rendezvous

### Tools & Languages
- **MATLAB / Simulink** — modeling and simulation of **attitude dynamics** and **Hill (LVLH) equations**
- **C++ (Eigen3)** — implementation of the **Hill equations** using an LTI model with RK4 integration.  
  *Note:* the **C++ simulation of Hill** is available; attitude is currently done in **MATLAB/Simulink**.
- **Python (pandas, matplotlib)** — post-processing and visualization of the C++ simulation results (plots of \(x,y,z\) and \(v_x,v_y,v_z\)).
