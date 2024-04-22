import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("energy.dat", delim_whitespace=True, header=None)

walkers = data[0]
local_energy = data[1]

plt.figure()
plt.plot(walkers, local_energy)
plt.xlabel("Walkers")
plt.ylabel("Local Energy (a. u.)")
plt.title("Local Energy for H3+")
plt.xlim(0, 500)

plt.savefig("local_energy.png")
plt.show()

