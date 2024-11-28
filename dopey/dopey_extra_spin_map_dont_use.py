"""
This is an extra module for manipulating spinmap data.

Methods:

concatSpinmaps(D1, D2)
    D1 and D2 are two spinmaps, measured with the same settings except for one deflector axis. Returns a new dopey spinmap dict.

Double check the methods below. Not ready.

    merge(D1, D2, coil, rotator)
        D1 and D2 are the data loaded by dopey.load() for each of the data sets to be merged.
        coil is either 1 or 2, rotator is either -1 or 1.

    plotMerged(D, intensity)
        D is the data set from merge().
        intensity is a string telling the method what to plot. It can be:
            intensity_minus, intensity_plus, asymmetry, component_intensity_minus, component_intensity_plus, or component_intensity

    polarization(D, polar, sherman)
        D is a list of data sets from merge(), i.e [c2rm, c2rp, c1].
        polar is the manipulator polar angle (default 0).
        sherman is... well, the sherman function (default 0.3).

    plotPolarization(D, plot_type, vmin, vmax)
        D is the output from polarization().
        plot_type is a string telling the method what to plot. It can be p, c, px, py, pz, cx, cy, and cz.


"""



__version__ = "24.11.26"
__author__  = "Mats Leandersson"


import numpy as np
from colorama import Fore
from copy import deepcopy
import matplotlib.pyplot as plt
from copy import deepcopy

print(f"{__name__}, {__version__}, beta")
print(Fore.LIGHTBLACK_EX + "Note: This module is intended for merging spin maps with fixed energies, nothing else." + Fore.RESET + "\n")



def concatSpinmaps(D1 = {}, D2 = {}, shup = True):
    """
    This method contatenates two spin maps that for instance has been reccorded for deflector range -8 to 1 and for 2 to 8.
    It returns a new dopey dict of type spin_map.
    """
    DD = {}
    #
    try: typ1, typ2 = D1.get("type", ""), D2.get("type", "")
    except: typ1, typ2 = "", ""
    if not typ1 == typ2:
        print(Fore.RED + "concatSpinmaps(): Arguments D1 and D2 must both be dopey spin map dicts." + Fore.RESET); return DD
    if not typ1 == "spin_map":
        print(Fore.RED + "concatSpinmaps(): Arguments D1 and D2 must be dopey spin map dicts." + Fore.RESET); return DD
    #
    energy1, energy2 = D1["x"], D2["x"]
    if not np.array_equal(energy1, energy2):
        print(Fore.RED + "concatSpinmaps(): Arguments D1 and D2 does not have the same energy (axis)." + Fore.RESET); return DD
    
    # ---------- check for what axis we are going to concatenate
    def11, def12 = D1["y"], D2["y"]
    def21, def22 = D1["z"], D2["z"]
    if np.array_equal(def11, def12) and not np.array_equal(def21, def22): conc_ax = "def2"
    elif not np.array_equal(def11, def12) and np.array(def21, def22): conc_ax = "def1"
    else: conc_ax = ""
    if not conc_ax in ["def1", "def2"]:
        print(Fore.RED + "concatSpinmaps(): Arguments D1 and D2 must have one common deflector axis and one that will be concatenated." + Fore.RESET); return DD
    
    # intensity: [polarity][scan][delf][defl][data]

    if conc_ax == "def2": axis1, axis2, axis = np.copy(D1["y"]), np.append(D1["z"], D2["z"]), 0
    else: axis1, axis2, axis = np.append(D1["y"], D2["y"]), np.copy(D1["z"]), 1
    DM, DP = [], []
    for dm1, dm2, dp1, dp2 in zip(D1["intensity"][0], D2["intensity"][0], D1["intensity"][1], D2["intensity"][1]):
        dm, dp = np.append(dm1, dm2, axis = axis), np.append(dp1, dp2, axis = axis)
        DM.append(dm)
        DP.append(dp)
    DM, DP = np.array(DM), np.array(DP)

    DD = deepcopy(D1)
    experiment = DD.get("experiment"); experiment.update({"Spectrum_ID": 999999}); DD.update({"experiment": experiment})
    DD.update({"y": axis1, "z": axis2, "intensity": np.array([DM, DP])})
    
    mean_m, mean_p = np.zeros([len(axis1), len(axis2), 1]), np.zeros([len(axis1), len(axis2), 1])
    for curve in DM: mean_m += curve
    for curve in DP: mean_p += curve
    mean_m, mean_p = mean_m / len(DM), mean_p / len(DP)
    DD.update({"intensity_mean": np.array([mean_m, mean_p])})

    return DD






# ====================================== EVERYTHING BELOW IS UNTESTED OR SHOULD BE REWRITTEN ====================================
# ===============================================================================================================================
# ===============================================================================================================================
# ===============================================================================================================================
# ===============================================================================================================================



def mergeSpinmaps(D1 = {}, D2 = {}, coil = 1, rotator = 1, shup = False):
    ret_dict = {"axis1": np.array([]), 
                "axis2": np.array([]), 
                "intensity_minus": np.array([]), 
                "intensity_plus": np.array([]), 
                "asymmetry": np.array([]),
                "component_intensity_minus": np.array([]),
                "component_intensity_plus": np.array([]),
                "coil": 0,
                "rotator": 0,
                "axis1_label": "", 
                "axis2_label": ""}

    ok = True
    typ1, typ2 = "", ""
    try: typ1 = D1.get("type", "")
    except: ok = False
    try: typ2 = D2.get("type", "")
    except: ok = False
    if not typ1 == "spin_map":
        print(Fore.RED + "mergeSpinmaps(): The first argument must be a dopey spin-map dict." + Fore.RESET); ok = False
    if not typ2 == "spin_map":
        print(Fore.RED + "mergeSpinmaps(): The second argument must be a dopey spin-map dict." + Fore.RESET); ok = False
    if not ok: return ret_dict
    #
    try: coil = int(coil)
    except: coil = 1
    try: rotator = int(rotator)
    except: rotator = 1
    if not coil in [1,2]:
        coil = 1
        print(Fore.MAGENTA + "mergeSpinmaps(): The coil has been set to 1. The argument coil should be 1 or 2." + Fore.RESET)
    if not rotator in [-1,1]:
        print(Fore.MAGENTA + "mergeSpinmaps(): The rotator has been set to 1. The argument rotator should be -1 or 1." + Fore.RESET)
    ret_dict.update({"coil": coil, "rotator": rotator})
    #
    x1, x2 = D1["x"], D2["x"]
    y1, y2 = D1["y"], D2["y"]
    z1, z2 = D1["z"], D2["z"]
    #
    if not len(x1) == 1 and not len(x2) == 1:
        print(Fore.RED + "mergeSpinmaps(): This method is right now only dealing with data collected in FE scan mode." + Fore.RESET); return ret_dict
    if not x1[0] == x2[0]: 
        print(Fore.RED + "mergeSpinmaps(): The kinetic energy in the two data sets does not match." + Fore.RESET); return ret_dict
    #
    merge_type = 0
    yaxes_same, zaxes_same = False, False
    if len(y1) == len(y2):
        if (y1 == y2).all(): yaxes_same = True
    else:
        yaxes_same = False
    if len(z1) == len(z2):
        if (z1 == z2).all(): zaxes_same = True
    else:
        zaxes_same = False
    if yaxes_same and zaxes_same:
        print(Fore.RED + "mergeSpinmaps(): The two data sets share the same deflector ranges so they can not be merged by this method." + Fore.RESET); return ret_dict
    elif yaxes_same and not zaxes_same:
        merge_type = "z"
        if not shup: print(Fore.BLUE + f"mergeSpinmaps(): Will try to merge the data sets along {D1['labels']['z']}." + Fore.RESET)
    elif not yaxes_same and zaxes_same:
        merge_type = "y"
        if not shup: print(Fore.BLUE + f"mergeSpinmaps(): Will try to merge the data sets along {D1['labels']['y']}." + Fore.RESET)
    else:
        print(Fore.RED + "mergeSpinmaps(): The two data sets do not share any deflector range so they can not be merged by this method." + Fore.RESET); return ret_dict
    #
    if merge_type == "z":
        yaxis = np.copy(y1)
        zaxis = np.append(z1, z2)
        int1m, int1p = D1["intensity_mean"][0], D1["intensity_mean"][1]
        int2m, int2p = D2["intensity_mean"][0], D2["intensity_mean"][1]

        intm, intp = [], []
        for curve in int1m: intm.append(curve)
        for curve in int2m: intm.append(curve)
        for curve in int1p: intp.append(curve)
        for curve in int2p: intp.append(curve)
        intm, intp = np.array(intm), np.array(intp)
        asymmetry = (intp - intm) / (intp + intm)
        cintm = (intm + intp) / 2 * (1 - asymmetry)
        cintp = (intm + intp) / 2 * (1 + asymmetry)

        ret_dict.update({"axis1": yaxis, "axis2": zaxis, "intensity_minus": intm, "intensity_plus": intp, 
                        "asymmetry": asymmetry, "component_intensity_minus": cintm, "component_intensity_plus": cintp, 
                        "axis1_label": D1["labels"]["y"], "axis2_label": D1["labels"]["z"]})
        return ret_dict
    
    elif merge_type == "y":
        print(Fore.MAGENTA + "mergeSpinmaps(): Sorry, have not finished the method for merging data set along this deflector axis. Easy fix, let me knonw." + Fore.RESET)



def plotMerged(D = None, intensity = "", ax = None):
    if not type(D) is dict:
        print(Fore.RED + "plotMerged(): Argument D must be a dict from merge()." + Fore.RESET); return ax
    req_keys = ["axis1", "axis2", "intensity_minus", "intensity_plus", "asymmetry", "component_intensity_minus", "component_intensity_plus", "coil", "rotator", "axis1_label", "axis2_label"]
    if D == {}:
        print(Fore.RED + "plotMerged(): Argument D must be a dict from merge()." + Fore.RESET); return ax
    for key in D.keys():
        if not key in req_keys:
            print(Fore.RED + "plotMerged(): The dict in argument D is missing expected keys." + Fore.RESET); return ax
    #
    if not type(intensity) is str: intensity = ""
    allowed_intensities = ["intensity_minus", "intensity_plus", "asymmetry", "component_intensity_minus", "component_intensity_plus", "component_intensity"]
    if intensity == "":
        fig, ax = plt.subplots(nrows = 2, ncols = 3, figsize = (10, 6))
        ax = ax.flatten()
        extent = [D["axis1"][0], D["axis1"][-1], D["axis2"][-1], D["axis2"][0]]
        ax[0].imshow(D["intensity_minus"], extent = extent, aspect = "equal", cmap = "bone_r")
        ax[1].imshow(D["intensity_plus"], extent = extent, aspect = "equal", cmap = "bone_r")
        ax[2].imshow(D["asymmetry"], extent = extent, aspect = "equal", cmap = None)
        ax[3].imshow(D["component_intensity_minus"], extent = extent, aspect = "equal", cmap = "bone_r")
        ax[4].imshow(D["component_intensity_plus"], extent = extent, aspect = "equal", cmap = "bone_r")
        ax[5].imshow(D["component_intensity_plus"] - D["component_intensity_minus"], extent = extent, aspect = "equal", cmap = "bwr")
        for i, txt in enumerate(["Scattered M-", "Scattered M+", "Asymmetry", "Component intensity +", "Component intensity -", ""]):
            ax[i].invert_yaxis()
            ax[i].set_title(txt, fontsize = 12)
            if i in [0, 3]: ax[i].set_ylabel(D["axis2_label"])
            if i in [3, 4, 5]: ax[i].set_xlabel(D["axis1_label"])
        fig.tight_layout()
    else:
        fig = None
        if type(ax) is type(None):
            fig, ax = plt.subplots(figsize = (6,3))
            if intensity == "component_intensity":
                map = D["component_intensity_plus"] - D["component_intensity_minus"]
                cmap = "bwr"
            else:
                map = D[intensity]
                cmap = "bone_r"
            extent = [D["axis1"][0], D["axis1"][-1], D["axis2"][-1], D["axis2"][0]]
            ax.imshow(map, extent = extent, aspect = "equal", cmap = cmap)
            ax.invert_yaxis()
            ax.set_ylabel(D["axis2_label"])
            ax.set_xlabel(D["axis1_label"])
            if not type(fig) is type(None): fig.tight_layout()
    #
    return ax


def polarization(D = [{}, {}, {}], polar = 0., sherman = 0.3, shup = False):
    """
    """
    ret_dict = {"axis1": np.array([]), "axis2": np.array([]),
                "p": "", "px": np.array([]), "py": np.array([]), "pz": np.array([]),
                "cxm": np.array([]), "cxp": np.array([]), "cym": np.array([]), "cyp": np.array([]), "czm": np.array([]), "czp": np.array([]),
                "polar": np.NaN, "axis1_label": "", "axis2_label": ""}
    #
    if not type(D) is list:
        print(Fore.RED + "polarization(): Argument D must be a list of dicts from merge()." + Fore.RESET); return ret_dict
    if not len(D) in [1,2,3]:
        print(Fore.RED + "polarization(): Argument D must be a list of length 1, 2, or 3." + Fore.RESET); return ret_dict
    for d in D:
        if not "coil" in d:
            print(Fore.RED + "polarization(): Argument D (a list) must cointain dicts from merge()." + Fore.RESET); return ret_dict
    #
    DD, ptyp = [None, None, None], ""
    if len(D) == 1:
        if D[0]["coil"] == 1:
            ptyp = "z"
            DD[3] = deepcopy(D[0])
        else:
            print(Fore.RED + "polarization(): Argument D contains one data set. It must be coil = 1. It is not." + Fore.RESET)
            return ret_dict
    elif len(D) == 2:
        if D[0]["coil"] == 2 and D[1]["coil"] == 2:
            if D[0]["rotator"] == -1 and D[1]["rotator"] == 1:
                ptyp = "xy"
                DD[0], DD[1] = deepcopy(D[0]), deepcopy(D[1])
            elif D[0]["rotator"] == 1 and D[1]["rotator"] == -1:
                ptyp = "xy"
                DD[0], DD[1] = deepcopy(D[1]), deepcopy(D[0])
            else:
                print(Fore.RED + "polarization(): Argument D contains two data set. They both must be coil = 2, and rotators -1 and +1. They are not." + Fore.RESET)
                return ret_dict
    elif len(D) == 3:
        c1 = 0
        if D[0]["coil"] == 1:
            DD[2] = deepcopy(D[0]); c1 += 1
        elif D[1]["coil"] == 1: 
            DD[2] = deepcopy(D[1]); c1 += 1
        elif D[2]["coil"] == 1: 
            DD[2] = deepcopy(D[2]); c1 += 1
        if c1 == 1: # -----
            if D[0]["coil"] == 2:
                if D[0]["rotator"] == -1: DD[0] = deepcopy(D[0])
                elif D[0]["rotator"] == 1: DD[1] = deepcopy(D[0])
            if D[1]["coil"] == 2:
                if D[1]["rotator"] == -1: DD[0] = deepcopy(D[1])
                elif D[1]["rotator"] == 1: DD[1] = deepcopy(D[1])
            if D[2]["coil"] == 2:
                if D[2]["rotator"] == -1: DD[0] = deepcopy(D[2])
                elif D[2]["rotator"] == 1: DD[1] = deepcopy(D[2])
            #
            if not (type(DD[0]) is type(None) or type(DD[1]) is type(None) or type(DD[2]) is type(None)):
                ptyp = "xyz"
    else:
        print(Fore.RED + "polarization(): Argument D (list) must contain 1, 2, or 3 data sets. It is not." + Fore.RESET)
        return ret_dict
    #
    #print(Fore.GREEN + "--- debug ---")
    #print(f"{ptyp}")
    #for i, dd in enumerate(DD):
    #    try: print(f"{i}, {dd['coil'] = }, {dd['rotator'] = }")
    #    except: pass
    #print("--- ----- ---" + Fore.RESET)
    #
    if not ptyp in ["xy", "z", "xyz"]:
        print(Fore.RED + "polarization(): " + Fore.RESET)
        print(Fore.RED + "                " + Fore.RESET)
        return ret_dict
    ret_dict.update({"p": ptyp})
    #
    if ptyp in ["xy", "xyz"]:
        ret_dict.update({"axis1": D[0]["axis1"], "axis2": D[0]["axis2"]})
        ret_dict.update({"axis1_label": D[0]["axis1_label"], "axis2_label": D[0]["axis2_label"]})
    else:
        ret_dict.update({"axis1": D[2]["axis1"], "axis2": D[2]["axis2"]})
        ret_dict.update({"axis1_label": D[2]["axis1_label"], "axis2_label": D[2]["axis2_label"]})
    #
    if ptyp == "xyz":
        a1 = D[0]["asymmetry"]
        a2 = D[1]["asymmetry"]
        a3 = D[2]["asymmetry"]
        px = (a2 - a1) / np.sqrt(2) / sherman
        py = -(a2 + a1) / np.sqrt(2) / sherman
        pz = a3 / sherman
        tot_int = (D[0]["intensity_minus"] + D[0]["intensity_plus"] + D[1]["intensity_minus"] + D[1]["intensity_plus"] + D[2]["intensity_minus"] + D[2]["intensity_plus"]) / 6
        cxm = tot_int * (1 - px)
        cxp = tot_int * (1 + px)
        cym = tot_int * (1 - py)
        cyp = tot_int * (1 + py)
        czm = tot_int * (1 - pz)
        czp = tot_int * (1 + pz)
        ret_dict.update({"px": px, "py": py, "pz": pz, "cxm": cxm, "cxp": cxp, "cym": cym, "cyp": cyp, "czm": czm, "czp": czp})
    #
    elif ptyp == "xy":
        a1 = D[0]["asymmetry"]
        a2 = D[1]["asymmetry"]
        px = (a2 - a1) / np.sqrt(2) / sherman
        py = -(a2 + a1) / np.sqrt(2) / sherman
        tot_int = (D[0]["intensity_minus"] + D[0]["intensity_plus"] + D[1]["intensity_minus"] + D[1]["intensity_plus"]) / 4
        cxm = tot_int * (1 - px)
        cxp = tot_int * (1 + px)
        cym = tot_int * (1 - py)
        cyp = tot_int * (1 + py)
        ret_dict.update({"px": px, "py": py, "cxm": cxm, "cxp": cxp, "cym": cym, "cyp": cyp})
    #
    elif ptyp == "z":
        a3 = D[2]["asymmetry"]
        pz = a3 / sherman
        tot_int = (D[2]["intensity_minus"] + D[2]["intensity_plus"]) / 2
        czm = tot_int * (1 - pz)
        czp = tot_int * (1 + pz)
        ret_dict.update({"pz": pz, "czm": czm, "czp": czp})
    #
    else:
        print(Fore.RED + "polarization(): Something went wrong, most likely due to sloppy coding." + Fore.RESET)
        return ret_dict
    #

    try: polar = float(polar)
    except:
        print(Fore.RED + "polarization(): The polar argument must be a number. Setting it to polar = 0." + Fore.RESET)
        polar = 0
    if abs(polar) >= 90:
        print(Fore.RED + "polarization(): The polar argument is too large so ignoring it! Setting polar = 0." + Fore.RESET)
        polar = 0
    ret_dict.update({"polar": polar})
    #
    if not shup:
        if ptyp == "xyz": txt = "polarizations Px, Py, and Pz"
        elif ptyp == "xy": txt = "polarizations Px and Py"
        elif ptyp == "z": txt = "polarization Pz"
        print(Fore.BLUE + f"polarization(): Calculated {txt} in the analyzer coordinate system." + Fore.RESET)
    #
    if not polar == 0:
        if ptyp in ["xy", "xyz"]: ishp = np.shape(DD[0]["intensity_minus"])
        else: ishp = np.shape(DD[2]["intensity_minus"])
        px = ret_dict.get("px", np.zeros(ishp))
        py = ret_dict.get("py", np.zeros(ishp))
        pz = ret_dict.get("pz", np.zeros(ishp))
        angle = np.deg2rad(polar)
        #
        p1, p2, p3 = np.copy(pz), np.copy(px), np.copy(py)
        P1 = p1 * np.cos(angle) - p2 * np.sin(angle)
        P2 = p1 * np.sin(angle) + p2 * np.cos(angle)
        P3 = p3 # not rotated by the polar angle
        Px = np.copy(P2)
        Py = np.copy(P3)
        Pz = np.copy(P1)

        if ptyp == "xyz":
            tot_int = (D[0]["intensity_minus"] + D[0]["intensity_plus"] + D[1]["intensity_minus"] + D[1]["intensity_plus"] + D[2]["intensity_minus"] + D[2]["intensity_plus"]) / 6
        elif ptyp == "xy":
            tot_int = (D[0]["intensity_minus"] + D[0]["intensity_plus"] + D[1]["intensity_minus"] + D[1]["intensity_plus"]) / 4
        elif ptyp == "z":
            tot_int = (D[2]["intensity_minus"] + D[2]["intensity_plus"]) / 2
        
        cxm = tot_int * (1 - Px)
        cxp = tot_int * (1 + Px)
        cym = tot_int * (1 - Py)
        cyp = tot_int * (1 + Py)
        czm = tot_int * (1 - Pz)
        czp = tot_int * (1 + Pz)
        ret_dict.update({"px": Px, "py": Py, "pz": Pz, "cxm": cxm, "cxp": cxp, "cym": cym, "cyp": cyp, "czm": czm, "czp": czp})

        if not shup:
            print(Fore.BLUE + f"polarization(): Rotated the polarization with an angle of {polar} degrees." + Fore.RESET)

    return ret_dict


def plotPolarization(D = {}, plot_type = "p", vmin = None, vmax = None, shup = False):
    """
    """
    plot_types = ["p", "c", "px", "py", "pz", "cx", "cy", "cz"]
    if not shup:
        print(Fore.BLUE + "plotPolarization(): Pass argument plot_type as 'p' for polarization, 'c' for component intensity, ... ")
        print("    plot_type = p              all polarizations")
        print("                c              all component intensities")
        print("                px, py, pz     polarization")
        print("                cx, cy, cz     component intensity")

    if not plot_type in plot_types:
        print(Fore.RED + "plotPolarization(): Argument plot_type must be 'p' or 'c'." + Fore.RESET); return
    #
    if not type(D) is dict:
        print(Fore.RED + f"plotPolarization(): Argument D must be a dict from polarization()." + Fore.RESET); return ax
    if not "px" in D and not "py" in D and not "pz" in D:
        print(Fore.RED + f"plotPolarization(): Argument D must be a dict from polarization()." + Fore.RESET); return ax
    #
    extent = [D["axis1"][0], D["axis1"][-1], D["axis2"][-1], D["axis1"][0]]
    #
    if plot_type in ["p", "c"]:
        fig, ax = plt.subplots(ncols = 3, figsize = (10,4))
        if plot_type == "p":
            ax[0].imshow(D["px"], extent = extent, aspect = "equal", cmap = "bone_r", vmin = vmin, vmax = vmax)
            ax[1].imshow(D["py"], extent = extent, aspect = "equal", cmap = "bone_r", vmin = vmin, vmax = vmax)
            ax[2].imshow(D["pz"], extent = extent, aspect = "equal", cmap = "bone_r", vmin = vmin, vmax = vmax)
        elif plot_type == "c":
            ax[0].imshow(D["cxp"] - D["cxm"], extent = extent, aspect = "equal", cmap = "bwr", vmin = vmin, vmax = vmax)
            ax[1].imshow(D["cyp"] - D["cym"], extent = extent, aspect = "equal", cmap = "bwr", vmin = vmin, vmax = vmax)
            ax[2].imshow(D["czp"] - D["czm"], extent = extent, aspect = "equal", cmap = "bwr", vmin = vmin, vmax = vmax)
        for a in ax: a.set_xlabel(D["axis1_label"], fontsize = 10)
        for a in ax: a.invert_yaxis()
        for i, ttl in enumerate(["Px", "Py", "Pz"]): ax[i].set_title(ttl, fontsize = 12)
        ax[0].set_ylabel(D["axis2_label"], fontsize = 10)
        fig.tight_layout()
    #
    elif plot_type in ["px", "py", "pz"]:
        fig, ax = plt.subplots(figsize = (4,4))
        ax.imshow(D[plot_type], extent = extent, cmap = "bone_r", aspect = "equal", vmin = vmin, vmax = vmax)
        ax.invert_yaxis()
        ax.set_xlabel(D["axis1_label"], fontsize = 10)
        ax.set_ylabel(D["axis2_label"], fontsize = 10)
        ttl = f"P{plot_type[1]}"
        ax.set_title(ttl, fontsize = 12)
        fig.tight_layout()
        return ax
    #
    elif plot_type in ["cx", "cy", "cz"]:
        fig, ax = plt.subplots(figsize = (4,4))
        km, kp = f"{plot_type}m", f"{plot_type}p"
        map = D[kp] - D[km]
        ax.imshow(map, extent = extent, cmap = "bwr", aspect = "equal", vmin = vmin, vmax = vmax)
        ax.invert_yaxis()
        ax.set_xlabel(D["axis1_label"], fontsize = 10)
        ax.set_ylabel(D["axis2_label"], fontsize = 10)
        ttl = f"P{plot_type[1]}"
        ax.set_title(ttl, fontsize = 12)
        fig.tight_layout()
        return ax











# ====================

def merge2(D1 = {}, D2 = {}):
    """
    """

    # check the dicts

    
    





            




        







