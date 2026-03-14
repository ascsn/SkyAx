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
mpl.rcParams["text.usetex"] = True
mpl.rcParams["axes.linewidth"] = 1.5
mpl.rcParams["axes.formatter.useoffset"] = False

# Increased font sizes and using numeric fontweight for maximum boldness
# 900 typically corresponds to 'black' or 'heavy'
fontDict = {"family": "helvetica", "fontsize": 32, "fontweight": 900}
fontDictTitle = {"family": "helvetica", "fontsize": 32, "fontweight": 900}

# Mapping from shell numbers to approximate N states based on the gnuplot file
# shell  N_p/2  N_n/2  N
# 2     87   117 408
# 4    174   214 776
# 6    303   358 1322
# 8    479   555 2068
# 10   710  815  3050
# 12  1009  1144 4306
# 14  1378  1544 5844
shell_to_N = {2: 408, 4: 776, 6: 1322, 8: 2068, 10: 3050, 12: 4306, 13: 5000, 14: 5844}

# Data from summary_170Yb_Fayans.txt
# Format: [Cutoff, Deocc, N, Energy, RMSC, beta2]
# All deocc values from summary file

# data_5 = [
#     [5, 2, 408, -1374.300777, 5.26540, 0.333],
#     [5, 4, 776, -1374.300777, 5.26540, 0.333],
#     [5, 6, 1322, -1374.300777, 5.26540, 0.333],
#     [5, 8, 2068, -1374.300777, 5.26540, 0.333],
#     [5, 10, 3050, -1374.300777, 5.26540, 0.333],
#     [5, 12, 4306, -1374.300777, 5.26540, 0.333],
#     [5, 13, 5000, -1374.300777, 5.26540, 0.333],
#     [5, 14, 5844, -1374.300777, 5.26540, 0.333]
# ]

# data_10 = [
#     [10, 2, 408, -1374.621157, 5.27645, 0.329],
#     [10, 4, 776, -1374.621157, 5.27645, 0.329],
#     [10, 6, 1322, -1374.621157, 5.27645, 0.329],
#     [10, 8, 2068, -1374.621157, 5.27645, 0.329],
#     [10, 10, 3050, -1374.621157, 5.27645, 0.329],
#     [10, 12, 4306, -1374.621157, 5.27645, 0.329],
#     [10, 13, 5000, -1374.621157, 5.27645, 0.329],
#     [10, 14, 5844, -1374.621157, 5.27645, 0.329]
# ]

data_20 = [
    [20, 2, 408, -1376.049953, 5.34111, 0.134],
    [20, 4, 776, -1378.075074, 5.36095, 0.079],
    [20, 6, 1322, -1378.075074, 5.36095, 0.079],
    [20, 8, 2068, -1378.075074, 5.36095, 0.079],
    [20, 10, 3050, -1378.075074, 5.36095, 0.079],
    [20, 12, 4306, -1378.075074, 5.36095, 0.079],
    [20, 13, 5000, -1378.075074, 5.36095, 0.079],
    [20, 14, 5844, -1378.075074, 5.36095, 0.079]
]

data_40 = [
    [40, 2, 408, -1377.726714, 5.35458, 0.091],
    [40, 4, 776, -1399.643809, 5.46593, 0.022],
    [40, 6, 1322, -1404.670300, 5.48143, 0.024],
    [40, 8, 2068, -1405.217095, 5.48314, 0.025],
    [40, 10, 3050, -1406.710413, 5.48677, 0.028],
    [40, 12, 4306, -1408.002041, 5.48935, 0.027],
    [40, 13, 5000, -1408.258817, 5.49012, 0.026],
    [40, 14, 5844, -1408.480379, 5.49064, 0.026]
]

data_60 = [
    [60, 2, 408, -1377.726714, 5.35458, 0.091],
    [60, 4, 776, -1401.570997, 5.47173, 0.025],
    [60, 6, 1322, -1416.548919, 5.51809, 0.032],
    [60, 8, 2068, -1421.197496, 5.53227, 0.033],
    [60, 10, 3050, -1422.450166, 5.53595, 0.034],
    [60, 12, 4306, -1423.575354, 5.54106, 0.036],
    [60, 13, 5000, -1424.099735, 5.54273, 0.036],
    [60, 14, 5844, -1424.397472, 5.54369, 0.038]
]

# Convert to numpy arrays
datasets = {
    # 5: np.array(data_5),
    # 10: np.array(data_10),
    20: np.array(data_20),
    40: np.array(data_40),
    60: np.array(data_60)
}

# Colors for different cuts
colors = {5: "tab:purple", 10: "tab:red", 20: "tab:blue", 40: "tab:orange", 60: "tab:green"}

# Create two-panel plot (matching gnuplot style)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 9), sharex=True)

for cut, data in datasets.items():
    N = data[:, 2]
    EB = data[:, 3]
    rC = data[:, 4]
    beta2 = data[:, 5]
    x = N ** (1/3)
    
    # Subtract the last value to get relative change (matching gnuplot)
    EB_relative = EB - EB[-1]  # E_B - E_B(last)
    rC_relative = (rC - rC[-1]) * 1000  # (r_C - r_C(last)) * 1000 to get 10^-3 fm

    ax2.plot(x, EB_relative, "o-", label=f"$\\varepsilon_{{cut}}$ = {cut} MeV", color=colors[cut], markersize=8)
    ax1.plot(x, rC_relative, "o-", label=f"$\\varepsilon_{{cut}}$ = {cut} MeV", color=colors[cut], markersize=8)

# Formatting (top panel - charge radius)
ax1.set_ylabel(r"charge radius  ($10^{-3}$ fm)")
ax1.legend(loc='best', frameon=True)
# ax1.grid(True, alpha=0.3)
ax1.set_title(r"$^{170}$Yb, FyDr (Fayans), HFB", pad=10)

# Formatting (bottom panel - energy)
ax2.set_xlabel(r"$N_{\mathrm{states}}^{1/3}$")
ax2.set_ylabel(r"total energy  (MeV)")
# ax2.grid(True, alpha=0.3)

# Apply tick and label settings to both axes
for ax in [ax1, ax2]:
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=28,
        labelright=False,
        direction="in",
        right=True,
        top=True,
        size=8,
    )
    ax.tick_params(
        axis="both",
        which="minor",
        labelsize=28,
        labelright=False,
        direction="in",
        right=True,
        top=True,
        size=4,
    )

# Adjust plot layout
fig.align_ylabels([ax1, ax2])


plt.tight_layout()
plt.savefig("trend-deocc-170Yb_Fayans.pdf", dpi=300)
print("Plot saved as trend-deocc-170Yb_Fayans.pdf")
# plt.show()
