import numpy as np
import matplotlib.pyplot as plt

type = ["explicit", "implicit", "improved", "crank"]
step = "100"
data = []

colors = ['#1f77b4', '#ff7f0e', '#d62728', '#2ca02c']


data.append(np.loadtxt(f"data/circuit_explicit.txt", usecols=(0, 1, 2)))

for i in range(len(data)):
    plt.plot(data[i][:,0], data[i][:,1], label=f'U_0')
    plt.plot(data[i][:,0], data[i][:,2], label=f'U_C')
plt.xlabel('time')
plt.ylabel('Velocity')
plt.title('Voltage change in time')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.grid()
plt.savefig(f"./plots/ElectricCircuit_{step}.png")
plt.show()