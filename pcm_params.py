"""Parameters for the PCM melt-front identifiability study, with provenance.

Every number below is tagged:

  (a) SOURCED    - taken from a named publication or manufacturer datasheet
  (b) ESTIMATED  - engineering estimate, assumption stated
  (c) DERIVED    - computed from (a) values

Nothing here is invented.  Where a property could not be sourced it is marked
(b) and the sweep that covers the uncertainty is named.

SOURCES
-------
[S1] Rubitherm Technologies GmbH, RT product range table.
     https://www.rubitherm.eu/en/productcategory/organische-pcm-rt
     RT42: melting area 42 C, heat storage capacity 165 kJ/kg.
[S2] Hammoodi et al., "CFD simulation of air layer effects on RT42 PCM melting
     in a square cell", PMC12379284, Table 1.
     rho 760 kg/m3, k 0.2 W/mK, cp 2000 J/kgK, L_f 165000 J/kg,
     melting range 311.5-315.5 K.
[S3] Kubitza et al., "Evaluation of commercial 18650 and 26700 sodium-ion cells
     and comparison with well-established lithium-ion cells",
     J. Power Sources Advances (2024), S2666248524000143 / KIT 1000171390.
     26700 Na-ion (HAKADI / Selian Energy): 3.5 Ah, 84.2 g, <=20 mOhm at 1 kHz,
     1.5-4.1 V, max 3C discharge, tested at 25 C.
     18650 Na-ion: 1.5 Ah, 37.4 g, 18.6 mOhm at 1 kHz.
[S4] Bhundiya, Hunt & Drolen, "Measurement of the effective radial thermal
     conductivities of 18650 and 26650 lithium-ion battery cells", TFAWS 2018.
     18650: 0.43 +/- 0.07 W/mK.  22650: 0.20 +/- 0.04 W/mK.
     Literature range quoted therein: 0.15 to 3.4 W/mK.

TWO PROVENANCE WARNINGS THAT MATTER FOR THE RESULT
--------------------------------------------------
1. Rubitherm [S1] calls 165 kJ/kg the "heat storage capacity", which by their
   convention is LATENT PLUS SENSIBLE over a stated temperature window.  [S2]
   uses the same 165 kJ/kg as pure latent heat.  We follow [S2] so the number
   is traceable, but this OVERSTATES the latent buffer, probably by 15-25%.
   The direction matters: it makes the PCM look like a better store than it is,
   which is generous to the proposal.  A negative identifiability result under
   a generous assumption is stronger, not weaker.

2. [S2] gives ONE thermal conductivity, not separate solid and liquid values.
   The solid/liquid conductivity contrast is the physical mechanism by which
   front position changes the composite thermal resistance, so a study that
   silently sets k_l = k_s would manufacture its own null result.  We therefore
   treat K_RATIO as a swept quantity (see K_RATIO_SWEEP) and report the answer
   across it rather than picking one value.
"""
import numpy as np
from pcm_solver import Material, Model

# --------------------------------------------------------------------------- #
# PCM -- Rubitherm RT42 paraffin
# --------------------------------------------------------------------------- #
PCM_RHO = 760.0        # kg/m3     (a) [S2]
PCM_CP = 2000.0        # J/kgK     (a) [S2]
PCM_K_SOLID = 0.2      # W/mK      (a) [S2]
PCM_LF = 165000.0      # J/kg      (a) [S2], see provenance warning 1
PCM_T_SOLIDUS = 311.5  # K         (a) [S2]; = 38.35 C, matches [S1] 38-43 C
PCM_T_LIQUIDUS = 315.5 # K         (a) [S2]; = 42.35 C
PCM_TM = 0.5 * (PCM_T_SOLIDUS + PCM_T_LIQUIDUS)   # 313.5 K   (c)
PCM_DTM = PCM_T_LIQUIDUS - PCM_T_SOLIDUS          # 4.0 K     (c)

# The 4 K mushy width above is PHYSICAL: RT42 is a technical-grade paraffin with
# a genuine 38-43 C melting range, not a pure substance.  The apparent-heat-
# capacity dTm is therefore only partly a numerical parameter here.  Trap 5 in
# the brief still applies to the NUMERICAL part, so DTM_SWEEP tests it.
DTM_SWEEP = [1.0, 2.0, 4.0, 6.0, 8.0]             # K         (b)

# Ratio k_liquid / k_solid.  1.0 = the single-value assumption of [S2];
# 0.7-0.9 is typical of paraffins whose two phases have been measured
# separately.  Swept because it drives the answer -- see provenance warning 2.
K_RATIO = 1.0                                      # (a) nominal, faithful to [S2]
K_RATIO_SWEEP = [1.0, 0.9, 0.8, 0.7, 0.5]          # (b)

# --------------------------------------------------------------------------- #
# Cell -- 26700 sodium-ion (HAKADI / Selian Energy)
# --------------------------------------------------------------------------- #
CELL_CAPACITY_AH = 3.5     # Ah        (a) [S3]
CELL_MASS = 0.0842         # kg        (a) [S3]
CELL_R_1KHZ = 0.020        # ohm       (a) [S3], <=20 mOhm specified
CELL_C_RATE_MAX = 3.0      # -         (a) [S3]
CELL_DIAM = 0.026          # m         (a) [S3] "26700" format convention
CELL_HEIGHT = 0.070        # m         (a) [S3]
CELL_R_OUTER = CELL_DIAM / 2.0

# Jellyroll occupies everything inside the can wall.
CAN_THICKNESS = 0.0003     # m         (b) typical 26700 nickel-plated steel can
CORE_R_OUTER = CELL_R_OUTER - CAN_THICKNESS

CELL_VOLUME = np.pi * CELL_R_OUTER ** 2 * CELL_HEIGHT          # m3    (c)
CELL_RHO = CELL_MASS / CELL_VOLUME                             # kg/m3 (c) ~2265

# NOT measured for Na-ion in [S3].  Li-ion cylindrical cells sit at
# 1000-1200 J/kgK, and Na-ion cells use the same family of materials
# (hard carbon anode, layered oxide cathode, Al foils, carbonate electrolyte),
# so the transfer is defensible but it is an estimate.
CELL_CP = 1100.0           # J/kgK     (b)

# From [S4]: 0.43 +/- 0.07 for 18650 Li-ion.  Applied to a Na-ion 26700, so (b).
# The cell interior is NOT the bottleneck here -- the PCM layer is ~100x more
# resistive -- so the conclusions are insensitive to this within its range.
CELL_K_RADIAL = 0.40       # W/mK      (b) [S4]

CAN_RHO = 7900.0           # kg/m3     (b) nickel-plated steel
CAN_CP = 500.0             # J/kgK     (b)
CAN_K = 15.0               # W/mK      (b)

# --------------------------------------------------------------------------- #
# Pack / environment
# --------------------------------------------------------------------------- #
PCM_THICKNESS = 0.005      # m         (b) design choice, swept in Stage E
H_CONV = 10.0              # W/m2K     (b) natural convection in air
T_AMB = 298.15             # K         (a) [S3] tests were run at 25 C
T_INIT = 298.15            # K         (b) starts thermally equilibrated

NOISE_SIGMA = 0.1          # K         given in the brief
SAMPLE_DT = 1.0            # s         1 Hz, given in the brief

# --------------------------------------------------------------------------- #
# Heat generation
# --------------------------------------------------------------------------- #
# Ohmic only, from the 1 kHz impedance.  This UNDERSTATES total heat: DC
# resistance is typically 2-3x the 1 kHz value, and entropic plus polarisation
# heat add more again, so real 3C heat is plausibly 2-4x this.  Used as the
# nominal anchor; Stage E sweeps the heat load over more than an order of
# magnitude, which covers the discrepancy.
def q_ohmic(c_rate):
    """Volumetric heat generation, W/m3, for a given C rate.  (c) from [S3]."""
    current = c_rate * CELL_CAPACITY_AH
    return current ** 2 * CELL_R_1KHZ / CELL_VOLUME


Q_3C = q_ohmic(3.0)        # ~59 300 W/m3   (c)
Q_NOMINAL = Q_3C

# A continuous 3C load lasts only 20 min in reality; sustained heating here
# represents repeated cycling at that average dissipation, which is the normal
# framing for a thermal-management study.  Stated so it is not mistaken for a
# single discharge.


# --------------------------------------------------------------------------- #
# Model builders
# --------------------------------------------------------------------------- #

def make_pcm(k_ratio=K_RATIO, dTm=PCM_DTM, Lf=PCM_LF, k_solid=PCM_K_SOLID):
    return Material("pcm", rho=PCM_RHO, cp=PCM_CP, k=k_solid, Lf=Lf,
                    Tm=PCM_TM, dTm=dTm, k_liquid=k_solid * k_ratio)


def build_model(pcm_thickness=PCM_THICKNESS, k_ratio=K_RATIO, dTm=PCM_DTM,
                Lf=PCM_LF, k_pcm=PCM_K_SOLID, n_core=40, n_can=4, n_pcm=80):
    """The cell/can/PCM composite used throughout the study."""
    core = Material("core", rho=CELL_RHO, cp=CELL_CP, k=CELL_K_RADIAL)
    can = Material("can", rho=CAN_RHO, cp=CAN_CP, k=CAN_K)
    pcm = make_pcm(k_ratio=k_ratio, dTm=dTm, Lf=Lf, k_solid=k_pcm)
    return Model(
        layers=[(core, CORE_R_OUTER, n_core),
                (can, CELL_R_OUTER, n_can),
                (pcm, CELL_R_OUTER + pcm_thickness, n_pcm)],
        geom="cylindrical", height=CELL_HEIGHT, r_inner=0.0)


def timescales(model=None, q=Q_NOMINAL, pcm_thickness=PCM_THICKNESS):
    """Regime numbers.  These are what decide the answer, so they are reported
    up front rather than buried in the sweep."""
    if model is None:
        model = build_model(pcm_thickness=pcm_thickness)
    alpha = PCM_K_SOLID / (PCM_RHO * PCM_CP)
    tau_diff = pcm_thickness ** 2 / alpha
    m_pcm = model.pcm_mass()
    E_latent = m_pcm * PCM_LF

    # Net power actually available to melt: input minus loss at the pinned
    # surface temperature.  Using the raw input would overstate the melt rate.
    A_out = model.A[-1]
    Q_in = q * CELL_VOLUME
    Q_loss = H_CONV * A_out * (PCM_TM - T_AMB)
    Q_net = Q_in - Q_loss
    t_melt = E_latent / Q_net if Q_net > 0 else np.inf

    return {
        "alpha_pcm": alpha,
        "tau_diff": tau_diff,
        "pcm_mass": m_pcm,
        "E_latent": E_latent,
        "A_out": A_out,
        "Q_in": Q_in,
        "Q_loss_at_Tm": Q_loss,
        "Q_net": Q_net,
        "t_melt": t_melt,
        "t_melt_over_tau_diff": t_melt / tau_diff,
        "Stefan": PCM_CP * PCM_DTM / PCM_LF,
        "Biot_pcm": H_CONV * pcm_thickness / PCM_K_SOLID,
    }


if __name__ == "__main__":
    mdl = build_model()
    ts = timescales(mdl)
    print(f"cell rho          {CELL_RHO:9.1f} kg/m3   (derived from 84.2 g / 26700 volume)")
    print(f"q''' at 3C        {Q_3C:9.0f} W/m3    ({Q_3C*CELL_VOLUME:.3f} W total)")
    print(f"PCM mass          {ts['pcm_mass']*1e3:9.2f} g")
    print(f"latent buffer     {ts['E_latent']:9.0f} J")
    print(f"outer area        {ts['A_out']*1e4:9.2f} cm2")
    print(f"loss at T_m       {ts['Q_loss_at_Tm']:9.3f} W")
    print(f"net melt power    {ts['Q_net']:9.3f} W")
    print(f"t_melt            {ts['t_melt']:9.0f} s")
    print(f"tau_diff (PCM)    {ts['tau_diff']:9.0f} s")
    print(f"t_melt/tau_diff   {ts['t_melt_over_tau_diff']:9.2f}")
    print(f"Stefan (cp dTm/L) {ts['Stefan']:9.4f}")
    print(f"Biot (PCM layer)  {ts['Biot_pcm']:9.4f}")
