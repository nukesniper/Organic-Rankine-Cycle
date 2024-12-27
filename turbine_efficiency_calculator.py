import numpy as np

def calculate_efficiency(VH, VR, stage):
    """
    Calculate the efficiency (eta) for different stages based on input VH (SP), VR, and stage number.

    Parameters:
    VH (float): Input value for SP (logarithmic terms applied).
    VR (float): Input value for Vr (logarithmic terms applied).
    stage (int): Stage number (1, 2, or 3) to select the corresponding coefficients.

    Returns:
    float: Calculated efficiency (eta).
    """
    # Coefficients for each stage (columns from the table)
    coefficients = {
        1: [
            0.90831500, -0.05248690, -0.04799080, -0.01710380, -0.00244002, 0.0,
            0.04961780, -0.04894860, 0.01171650, -0.00100473, 0.05645970,
            -0.01859440, 0.01288860, 0.00178187, -0.00021196, 0.00078667
        ],
        2: [
            0.923406, -0.0221021, -0.0233814, -0.00844961, -0.0012978, -0.00069293,
            0.0146911, -0.0102795, 0.0, 0.000317241, 0.0163959, -0.00515265,
            0.00358361, 0.000554726, 0.0, 0.000293607
        ],
        3: [
            0.932274, -0.01243, -0.018, -0.00716, -0.00118, -0.00044,
            0.0, 0.0, -0.0016, 0.000298, 0.005959, -0.00163,
            0.001946, 0.000163, 0.0, 0.000211
        ]
    }

    # Calculate terms Fi
    Fi = [
        1,                                # F0: Constant term
        np.log(VH),                       # F1: ln(SP)
        np.log(VH)**2,                    # F2: ln(SP)^2
        np.log(VH)**3,                    # F3: ln(SP)^3
        np.log(VH)**4,                    # F4: ln(SP)^4
        VR,                               # F5: Vr
        np.log(VR),                       # F6: ln(Vr)
        np.log(VR)**2,                    # F7: ln(Vr)^2
        np.log(VR)**3,                    # F8: ln(Vr)^3
        np.log(VR)**4,                    # F9: ln(Vr)^4
        np.log(VR) * np.log(VH),          # F10: ln(Vr) * ln(SP)
        (np.log(VR)**2) * np.log(VH),     # F11: ln(Vr)^2 * ln(SP)
        np.log(VR) * (np.log(VH)**2),     # F12: ln(Vr) * ln(SP)^2
        (np.log(VR)**3) * np.log(VH),     # F13: ln(Vr)^3 * ln(SP)
        (np.log(VR)**3) * (np.log(VH)**2),# F14: ln(Vr)^3 * ln(SP)^2
        (np.log(VR)**2) * (np.log(VH)**3) # F15: ln(Vr)^2 * ln(SP)^3
    ]

    # Calculate efficiency for the given stage
    eta = sum(Ai * Fi_i for Ai, Fi_i in zip(coefficients[stage], Fi))
    return eta

def get_best_efficiency(VH, VR, turbine_isentropic_efficiency):
    """
    Compare turbine isentropic efficiency with the calculated maximum efficiency.

    Parameters:
    VH (float): Input value for SP (VH).
    VR (float): Input value for Vr (VR).
    turbine_isentropic_efficiency (float): Current turbine isentropic efficiency.

    Returns:
    None: Prints results and suggestions.
    """
    # Calculate efficiency for all stages
    efficiencies = {stage: calculate_efficiency(VH, VR, stage) for stage in [1, 2, 3]}

    # Find the stage with the highest efficiency
    best_stage = max(efficiencies, key=efficiencies.get)
    max_eta = efficiencies[best_stage]

    # Print results
    for stage, eta in efficiencies.items():
        print(f"Stage {stage}: Efficiency (eta) = {eta:.6f}")

    # Compare turbine efficiency
    if turbine_isentropic_efficiency < max_eta:
        print(f"\nYou are currently using a turbine isentropic efficiency of {turbine_isentropic_efficiency:.6f}. "
              f"The suggested isentropic efficiency for the turbine is higher: {max_eta:.6f}, valid for a turbine of {best_stage} stages."
              f" Please consider changing the isentropic efficiency to avoid underestimating the thermal efficiency.")
    elif turbine_isentropic_efficiency > max_eta:
        print(f"\nYou are currently using a turbine isentropic efficiency of {turbine_isentropic_efficiency:.6f}. "
              f"The suggested isentropic efficiency for the turbine is lower: {max_eta:.6f}, valid for a turbine of {best_stage} stages."
              f" Please consider changing the isentropic efficiency to avoid overestimating the thermal efficiency.")
    else:
        print(f"\nYou are currently using the suggested turbine isentropic efficiency, valid for a turbine of {best_stage} expansion stages.")

    return max_eta, best_stage

