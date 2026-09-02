"""Order-of-magnitude X-ray attenuation of lead-bearing Herculaneum ink at the
Vesuvius fragment beam energies, from the only published lead loadings
(Brun et al. 2016 PNAS; Tack et al. 2016 Sci Rep: 16 +/- 5 and 84 +/- 5 ug/cm2
on two Paris fragments; their estimated writing thickness 50 um).

Mass attenuation coefficients (cm2/g) from NIST Hubbell & Seltzer tables,
log-log interpolated. Pb K-edge is 88.0045 keV: the '88 keV' fragment scans
sit ON the edge, so both sides are reported. Everything here is an estimate
with two unknowns stated in the accompanying note: the lead loading of OUR
fragments (never measured) and the depth distribution of the ink.
"""
import numpy as np, json
E  = np.array([50, 60, 80, 100.0])
Pb = np.array([8.041, 5.021, 2.419, 5.549])   # 100 keV is above the K-edge
C  = np.array([0.1871, 0.1753, 0.1610, 0.1514])
def interp(Ei, Eo, mu):  # log-log between bracketing table points
    i = np.searchsorted(E, Ei) - 1
    return np.exp(np.interp(np.log(Ei), np.log(E[i:i+2]), np.log(mu[i:i+2])))
mu_pb_54 = interp(54, E, Pb); mu_c_54 = interp(54, E, C)
mu_pb_88_below, mu_pb_88_above = 1.910, 7.683  # NIST edge values
mu_c_88 = np.exp(np.interp(np.log(88), np.log(E[2:4]), np.log(C[2:4])))
out = {"mu_over_rho_cm2_per_g": {"Pb_54keV": round(float(mu_pb_54),3), "C_54keV": round(float(mu_c_54),4),
        "Pb_88keV_below_edge": mu_pb_88_below, "Pb_88keV_above_edge": mu_pb_88_above, "C_88keV": round(float(mu_c_88),4)},
       "rows": []}
for load in (16, 84):                       # ug/cm2
    for t_um in (5, 50, 200):               # ink confined to a skin / Tack's 50 um / soaked through the sheet
        rho_pb = load*1e-6/(t_um*1e-4)      # g/cm3 of Pb averaged over that layer
        out["rows"].append({"Pb_ug_cm2": load, "layer_um": t_um,
            "mu_ink_Pb_54keV_per_cm": round(float(mu_pb_54*rho_pb),4),
            "mu_ink_Pb_88_below": round(mu_pb_88_below*rho_pb,4), "mu_ink_Pb_88_above": round(mu_pb_88_above*rho_pb,4)})
# papyrus: carbon; fibre-scale density ~1.2 g/cm3 (voxels resolve fibres at 3.24 um), sheet-average ~0.4 g/cm3
out["mu_papyrus_54keV_per_cm"] = {"fibre_1.2g": round(float(mu_c_54*1.2),4), "sheet_avg_0.4g": round(float(mu_c_54*0.4),4)}
out["mu_papyrus_88keV_per_cm"] = {"fibre_1.2g": round(float(mu_c_88*1.2),4), "sheet_avg_0.4g": round(float(mu_c_88*0.4),4)}
print(json.dumps(out, indent=1))
json.dump(out, open(__file__.replace('.py','.json'),'w'), indent=1)
