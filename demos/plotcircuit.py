import numpy as np
import matplotlib.pyplot as plt

types = ["explicit", "implicit", "improved", "crank"]
type = types[3]
# steps = ["100", "1000", "10000"]
steps = ["1000", "5000", "10000"]
# steps = ["5000", "10000", "15000"]
# steps = ["100", "1000", "10000", "20000"]
step = steps[0]
set = "2"
data1, data2, data3, data4 = [], [], [], []
data = [data1, data2, data3, data4]

colors = ['#1f77b4', '#ff7f0e', '#d62728', '#2ca02c', '#9467bd']

for i in range(len(steps)):
    data[i].append(np.loadtxt(f"demos/data/electric_circuit/set{set}/circuit_{type}{steps[i]}.txt", usecols=(0, 1, 2)))

for i in range(len(data1)):
    plt.plot(data3[i][:,0], np.cos(100*np.pi*data3[i][:,0]), label=f'U_0', color=colors[0])
    plt.plot(data1[i][:,0], data1[i][:,1], label=f'U_C steps = {steps[0]}', color=colors[1])
    plt.plot(data2[i][:,0], data2[i][:,1], label=f'U_C steps = {steps[1]}', color=colors[2])
    plt.plot(data3[i][:,0], data3[i][:,1], '--',label=f'U_C steps = {steps[2]}', color=colors[3])
    # plt.plot(data4[i][:,0], data4[i][:,1], label=f'U_C steps = {steps[2]}', color=colors[4])


plt.xlabel('Time [s]')
plt.ylabel('Voltage [V]')
plt.xlim(0, 0.04)
# plt.ylim(-1.5, 1.75)
plt.title(f'Electric Circuit - Crank-Nicolson Method')
plt.legend(loc='upper right')
plt.grid()
plt.savefig(f"demos/plots/ElectricCircuit_set2_crank2.png")
plt.show()
