"""
This code is used to generate panel plots by reading data from the collect_cat_ files.
The code reads data using read_data function, and then plots the data using make_plots
function. The preprocess_line function is used to preprocess data that is not separated 
by space.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, re
import scienceplots


dataPath = "./collect_cat_120Sn"
cutoff = np.arange(15, 101)


plt.style.use("science")
mpl.rcParams["text.usetex"] = True
mpl.rcParams["axes.linewidth"] = 1.5

colors = [
    "r",
    "b",
    "g",
    "k",
    "purple",
    "m",
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "orange",
    "y",
]


def preprocess_line(line):
    # Add a space before negative signs that are not part of scientific notation
    line = re.sub(r"(\d)(-)", r"\1 -", line)
    # Ensure that six decimal places are not split by adding spaces correctly
    line = re.sub(r"(-?\d\.\d{6})(?=-?\d)", r"\1 ", line)
    return line


def read_data(pairingType=None, functional=None):
    # Read data from collect_cat files
    fileList = []
    hfbwidths = []
    totalEnergyArray = []
    chargeRadiusArray = []
    pairingEnergyArray = []
    neutronGapArray = []
    for file in os.listdir(dataPath):
        fileList.append(file)
    # Sort the files in increasing order of either hfbwidth or box size.
    fileList = sorted(fileList)
    for file in fileList:
        totalEnergy = []
        chargeRadius = []
        pairingEnergy = []
        neutronGap = []
        # Checks the name of the file to contain pairing type, functional name,
        # element name and whether a cutoff is set
        print(pairingType, functional)
        if (
            pairingType in file
            and functional in file
            and "120Sn" in file
            and "cut100" in file
        ):
            print("Pair type", pairingType)
            # Create an array of elements separated by "_" in a filename.
            fileName = file.rstrip().split("_")
            print(fileName)
            # Split the element in the fileName array separated by "d"
            hfbwidths.append(fileName[4].split("d")[0])
            # Split the element in the fileName array separated by "." and round
            # to the nearest integer
            # hfbwidths.append(round(float(fileName[-1].split(".")[0]) * 0.3))
            with open(dataPath + "/" + file) as f:
                data = f.read()
            for line in data.splitlines():
                if line.startswith("Nmass"):
                    dataLine = preprocess_line(line).split(" ")
                    # print(dataLine)
                else:
                    dataLine = line.rstrip().split()
                # Check if size of array dataLine > 1
                if len(dataLine) > 1:
                    print(dataLine[5])
                    totalEnergy.append(float(dataLine[5]))
                    chargeRadius.append(float(dataLine[12]))
                    pairingEnergy.append(float(dataLine[10]))
                    neutronGap.append(float(dataLine[18]))
            totalEnergy = np.array(
                totalEnergy
            )  # - np.array(totalEnergy)[-1] * np.ones(
            #     len(totalEnergy)
            # )
            chargeRadius = np.array(chargeRadius)  # - np.array(chargeRadius)[
            #     -1
            # ] * np.ones(len(chargeRadius))
            pairingEnergy = np.array(pairingEnergy)  # - np.array(pairingEnergy)[
            #     -1
            # ] * np.ones(len(pairingEnergy))
            neutronGap = np.array(neutronGap)  # - np.array(neutronGap)[-1] * np.ones(
            #     len(neutronGap)
            # )
        # Collect all arrays into one big array
        if len(totalEnergy) > 0:
            totalEnergyArray.append(totalEnergy)
            chargeRadiusArray.append(chargeRadius)
            pairingEnergyArray.append(pairingEnergy)
            neutronGapArray.append(neutronGap)
    totalEnergyArray = np.array(totalEnergyArray)
    chargeRadiusArray = np.array(chargeRadiusArray)
    pairingEnergyArray = np.array(pairingEnergyArray)
    neutronGapArray = np.array(neutronGapArray)
    print(totalEnergyArray.shape)
    # totalEnergyArray = totalEnergyArray[totalEnergyArray.astype(bool)]
    return (
        hfbwidths,
        totalEnergyArray,
        chargeRadiusArray,
        pairingEnergyArray,
        neutronGapArray,
    )


def make_plots(data, labels, functional):
    (
        hfbwidth,
        totalEnergyArray,
        chargeRadiusArray,
        pairingEnergyArray,
        neutronGapArray,
    ) = data
    labelpadding = [5, 21, 20, 30]
    # yrange = [[-1025.2, -1015.5], [4.617, 4.654], [-92, -2], [0.51, 2.95]]
    fig, ax = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(6, 12), dpi=600)
    for j in range(4):
        for i, val in enumerate(data[j + 1]):
            ax[j].plot(
                cutoff,
                val,
                label=f"width = {hfbwidth[i]} fm",
                linewidth=2,
                color=colors[i],
            )
            # ax[j].axvline(x=60, color="teal", ls="--")
        dict1 = labels[j]
        xl = dict1["xlabel"]
        yl = dict1["ylabel"]
        fl = dict1["filename"]
        ttl = dict1["title"]
        # ax[j].set_ylim(yrange[j])
        ax[j].set_ylabel(yl, labelpad=labelpadding[j], fontsize=16)
        ax[j].tick_params(
            axis="both",
            which="major",
            labelsize=16,
            labelright=False,
            direction="in",
            right=True,
            top=True,
            size=8,
        )
        ax[j].tick_params(
            axis="both",
            which="minor",
            labelsize=16,
            labelright=False,
            direction="in",
            right=True,
            top=True,
            size=4,
        )
    fig.text(0.5, 0.07, r"$\epsilon_{{cut}}$ [MeV]", ha="center", fontsize=16)
    plt.subplots_adjust(hspace=0)
    ax[0].set_title(rf"$^{{{120}}}$Sn - {functional} - bcs", fontsize=18)
    plt.legend(loc="lower right", prop={"size": 7})
    fileName = f"{functional}_120Sn_bcs_props_Cut100.pdf"
    plt.savefig(fileName)
    plt.close()


def main():
    pairingType = ["bcs"]
    eltName = r"$^{120}$Sn"
    functional = "SVbas"
    for pairtype in pairingType:
        data = read_data(pairtype, functional)
        # Total Energy plot
        labels1 = {
            "xlabel": r"$\epsilon_{{cut}}$",
            "ylabel": r"Total Energy [MeV]",
            "filename": f"TotalEnergy_{pairtype}_{eltName}_{functional}_Cut100.pdf",
            "title": f"{eltName}",
        }

        # Radii plot
        labels2 = {
            "xlabel": r"$\epsilon_{{cut}}$",
            "ylabel": r"Charge Radius [fm]",
            "filename": rf"Radius_{pairtype}_{eltName}_{functional}_Cut100.pdf",
            "title": f"{eltName}",
        }

        # Neutron pairing energy plot
        labels3 = {
            "xlabel": r"$\epsilon_{{cut}}$",
            "ylabel": r"Pairing Energy [MeV]",
            "filename": rf"PairingEnergy_{pairtype}_{eltName}_{functional}_Cut100.pdf",
            "title": eltName,
        }

        # Average neutron pairing gap plot
        labels4 = {
            "xlabel": r"$\epsilon_{{cut}}$",
            "ylabel": r"Pairing Gap [MeV]",
            "filename": rf"PairingGap_{pairtype}_{eltName}_{functional}_Cut100.pdf",
            "title": eltName,
        }

        labels = [labels1, labels2, labels3, labels4]
        make_plots(data, labels, functional)

        os.system(rf"mv *_{pairtype}_*.pdf plots_paper/")


if __name__ == "__main__":
    main()
