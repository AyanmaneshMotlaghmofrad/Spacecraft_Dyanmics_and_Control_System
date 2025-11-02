import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("Hill_eq_Free_Response.csv")

# Relative Position of the chaser w.r.t the target in LVLH
plt.figure(figsize=(10, 5))
plt.plot(data["t"], data["x0"], label="x [m]")
plt.plot(data["t"], data["x1"], label="y [m]")
plt.plot(data["t"], data["x2"], label="z [m]")
plt.xlabel("Time [s]")
plt.ylabel("Relative Position [m]")
plt.title("Relative Position (Hill Equations - Free Response)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("Hill_Free_Response_Position.png", dpi=600)

# Relative velocity of the chaser w.r.t the target in LVLH frame
plt.figure(figsize=(10, 5))
plt.plot(data["t"], data["x3"], label="vx [m/s]")
plt.plot(data["t"], data["x4"], label="vy [m/s]")
plt.plot(data["t"], data["x5"], label="vz [m/s]")
plt.xlabel("Time [s]")
plt.ylabel("Relative Velocity [m/s]")
plt.title("Relative Velocity (Hill Equations - Free Response)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("Hill_Free_Response_Velocity.png", dpi=600)

plt.show()
