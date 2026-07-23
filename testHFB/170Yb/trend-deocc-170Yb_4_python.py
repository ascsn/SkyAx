import numpy as np
import matplotlib.pyplot as plt
import scienceplots
import matplotlib as mpl

# Apply the 'science' style for the plot design
plt.style.use('science')

# --- Set global font sizes for all text elements ---
mpl.rcParams['axes.labelsize'] = 25    # font size of the x and y labels
mpl.rcParams['xtick.labelsize'] = 28   # font size of the x tick labels
mpl.rcParams['ytick.labelsize'] = 28   # font size of the y tick labels
mpl.rcParams['legend.fontsize'] = 20   # legend font size
mpl.rcParams['axes.titlesize'] = 22    # font size of the axes title
mpl.rcParams['lines.linewidth'] = 4.0

# Hardcoded data
data_20 = [
    [20, 2, 408, -1374.529, 5.240], [20, 4, 776, -1374.537, 5.240],
    [20, 6, 1322, -1374.535, 5.240], [20, 8, 2068, -1374.535, 5.240],
    [20, 10, 3050, -1374.535, 5.240], [20, 12, 4306, -1374.523, 5.240],
    [20, 14, 5844, -1374.530, 5.240]
]

data_40 = [
    [40, 2, 408, -1374.2818, 5.2902], [40, 4, 776, -1374.5733, 5.2916],
    [40, 6, 1322, -1374.5795, 5.2916], [40, 8, 2068, -1374.5838, 5.2915],
    [40, 10, 3050, -1374.5778, 5.2915], [40, 12, 4306, -1374.5721, 5.2914],
    [40, 14, 5844, -1374.5608, 5.2914]
]

data_60 = [
    [60, 2, 408, -1374.2187, 5.2890], [60, 4, 776, -1374.5150, 5.2917],
    [60, 6, 1322, -1374.6366, 5.2921], [60, 8, 2068, -1374.6449, 5.2920],
    [60, 10, 3050, -1374.6171, 5.2919], [60, 12, 4306, -1374.5585, 5.2918],
    [60, 14, 5844, -1374.5528, 5.2918]
]

# Convert to numpy arrays
datasets = {
    20: np.array(data_20),
    40: np.array(data_40),
    60: np.array(data_60)
}

# Colors for different cuts
colors = {20: "tab:blue", 40: "tab:orange", 60: "tab:green"}

# Create two-panel plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 9), sharex=True)

for cut, data in datasets.items():
    N = data[:, 2]
    EB = data[:, 3]
    rC = data[:, 4]
    x = N ** (1/3)

    ax1.plot(x, EB, "-", label=f"cut={cut}", color=colors[cut])
    ax2.plot(x, rC, "-", label=f"cut={cut}", color=colors[cut])

# Formatting
ax1.set_ylabel(r"total energy (MeV)")
ax1.legend()
ax1.grid(True)

ax2.set_xlabel(r"$N_{st}^{1/3}$")
ax2.set_ylabel(r"charge radius ($10^{-3}$ fm)")
# ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.savefig("varydeocc_170Yb_4_python.pdf", dpi=300)
# plt.show()
