import streamlit as st
from A_toluene import toluene_orc
import matplotlib.pyplot as plt
import os

# Streamlit app
st.markdown("""
# Welcome to the **Organic Rankine Cycle Simulator by Lucas**!
This app is designed to simulate and analyze the performance of an ORC.
""")
st.write("This application simulates the performance of an Organic Rankine Cycle on the secondary side of a microreactor based on your input data. The reference cycle is the following: ")
# Ensure the path is properly formatted
image_path = r"C:\Users\lucas\Downloads\Organic RAnkine Cycle Streamlit.png"
if not os.path.exists(image_path):
    st.error(f"Image file not found at: {image_path}")
else:
    st.image(image_path, caption="T-s diagram and component representation of the ORC", use_container_width=True)
st.write("Please consider that T1' and T2' are fixed at 310 °C and 190 °C respectively and that the input heat is 15 MW.")

# Input fields
st.header("Input Data")
st.subheader("ORC information")
# Initial variables definition
pressure_loss = st.number_input("Enter the pressure losses within 2-3, 3-4, 5-6 and 6-1 in %:", value=1.0, max_value=5.0)
pressure_loss = pressure_loss / 100
p2 = st.number_input("Enter the pump's output pressure (p2) in bar:", value=10.0, max_value=41.264)
T1 = st.number_input("Enter T1, the condenser's outlet temperature (°C):", value=40.0)
T4 = st.number_input("Enter T4, the maximum temperature of the cycle (°C):", value=240.0, max_value=310.0)

st.subheader("Cycle efficiencies")
st.write("The code will suggest a valid efficiency for the turbine once the code will be run. Default data is based on organic rankine cycles literature.")
# Efficiencies
pump_isentropic_efficiency = st.number_input("Pump Isentropic Efficiency:", value=0.8, max_value=1.0)
pump_electric_motor_efficiency = st.number_input("Pump Electric Motor Efficiency:", value=0.9, max_value=1.0)
pump_mechanical_efficiency = st.number_input("Pump Mechanical Efficiency:", value=0.95, max_value=1.0)
turbine_isentropic_efficiency = st.number_input("Turbine Isentropic Efficiency:", value=0.9, max_value=1.0)
turbine_mechanical_efficiency = st.number_input("Turbine Mechanical Efficiency:", value=0.92, max_value=1.0)
electric_generator_efficiency = st.number_input("Electric Generator Efficiency:", value=0.96, max_value=1.0)
recuperator_efficiency = st.number_input("Recuperator Efficiency:", value=0.73, max_value=1.0)

st.subheader("Technical data valid for a shell and tube HX")
st.write("Please use default data (based on real ORC HX) unless other designs are available.")
# Input data for shell and tube HX
D_t = st.number_input("Internal Tube Diameter (D_t in m):", value=0.019, format="%.3f")
D_s = st.number_input("Shell Diameter (D_s in m):", value=0.203, format="%.3f")
N_tubes = st.number_input("Number of Tubes (N_t):", value=24, min_value=1)
N_t_passes = st.number_input("Number of Tube Passes (N_p):", value=2, min_value=1)
P_T = st.number_input("Square Pitch Length (P_T in m):", value=0.238, format="%.3f")
C = st.number_input("Tube Clearance (C in m):", value=0.005, format="%.3f")
B = st.number_input("Baffle Spacing (B in m):", value=0.305, format="%.3f")
k_SS_shellandtube = st.number_input("Thermal Conductivity Tube (in W/m K):", value=16)

st.subheader("Technical data valid for a plate HX")
st.write("Please use default data (based on real ORC HX) unless other designs are available.")
# Input data for plate HX
L = st.number_input("Effective Length (L in m):", value=0.4)
W = st.number_input("Effective Width (W in m):", value=0.2)
t_plate = st.number_input("Plate Thickness (t in m):", value=0.001, format="%.3f")
Dp = st.number_input("Port Diameter (Dp in m):", value=0.018)
Lv = st.number_input("Vertical Port Distance (Lv in m):", value=0.01)
Lh = st.number_input("Horizontal Port Distance (Lh in m):", value=0.02)
lambda_ = st.number_input("Corrugation Wavelength:", value=0.009, format="%.3f")
beta_angle = st.number_input("Corrugation Angle:", value=60)
Lc = st.number_input("Compressed Plate Pack Length (Lc in m):", value=2)
Phi = st.number_input("Surface Enlargement Factor (Phi):", value=1.25)
k_SS_plate = st.number_input("Stainless Steel Thermal Conductivity (k_SS):", value=16)
N_plates = st.number_input("Number of Plates (N_t):", value=30, min_value=1)
N_p_passes = st.number_input("Number of Plate Passes (N_p_passes):", value=1, min_value=1)

# Run simulation
if st.button("Run Simulation"):

    # Run the simulation
   (efficiency, max_eta, best_stage, pinchpoint, cold_T, cold_Q, hot_T, hot_Q, Area_plate, Area_shelltube, Cost_plate,
    Cost_shellandtube, Cost_turbine, Cost_pump, Cost_generator, Cost_cooling_tower) = toluene_orc(pressure_loss, p2, T1, T4,
    pump_isentropic_efficiency, pump_electric_motor_efficiency,
    pump_mechanical_efficiency, turbine_isentropic_efficiency,
    turbine_mechanical_efficiency, electric_generator_efficiency,
    recuperator_efficiency,
    D_t, D_s, N_tubes, N_t_passes, P_T, C, B, k_SS_shellandtube,
    L, W, t_plate, Dp, Lv, Lh, lambda_, beta_angle,
    Lc, Phi, k_SS_plate, N_plates, N_p_passes)

if pinchpoint < 0:
    st.error(f"Error: the pinch point is negative: please change the system conditions.")
# Display results
if turbine_isentropic_efficiency < max_eta:
    st.warning(f"\nYou are currently using a turbine isentropic efficiency of {turbine_isentropic_efficiency:.6f}. "
          f"The suggested isentropic efficiency for the turbine is higher: {max_eta:.6f}, valid for a turbine of {best_stage} stages."
          f" Please consider changing the isentropic efficiency to avoid underestimating the thermal efficiency.")
elif turbine_isentropic_efficiency > max_eta:
    st.warning(f"\nYou are currently using a turbine isentropic efficiency of {turbine_isentropic_efficiency:.6f}. "
          f"The suggested isentropic efficiency for the turbine is lower: {max_eta:.6f}, valid for a turbine of {best_stage} stages."
          f" Please consider changing the isentropic efficiency to avoid overestimating the thermal efficiency.")
else:
    st.success(f"\nYou are currently using the suggested turbine isentropic efficiency, valid for a turbine of {best_stage} expansion stages.")

st.write(f"The cycle thermal efficiency is {efficiency}, and your pinchpoint is {pinchpoint}. Here is a graph of the temperature range in the evaporator:")

# Plot the data
fig, ax = plt.subplots()  # Create a figure and axis for Streamlit
ax.plot(hot_Q, hot_T, linestyle='-', label='Hot Fluid')
ax.plot(cold_Q, cold_T, linestyle='-', label='Organic Cycle')

# Add labels and title
ax.set_xlabel('Energy Exchanged')
ax.set_ylabel('Temperature')
ax.set_title('Temperature vs. Energy Exchanged')

# Add legend and grid
ax.legend()
ax.grid(True)

# Display the plot in the Streamlit app
st.pyplot(fig)

st.write(f"The total heat exchange area required for a plate HX is: {Area_plate}.")
st.write(f"The total heat exchange area required for a shell and tube HX is: {Area_shelltube}.")
if Cost_shellandtube > Cost_plate:
    st.write(f"The most economically convenient heat exchanger is the plate HX with a cost of {Cost_plate}.")
else:
    st.write(f"The most economically convenient heat exchanger is the shell and tube HX with a cost of {Cost_shellandtube}.")
st.write(f"The cost of the turbine is {Cost_turbine}, the cost of the centrifugal pump (incl. motor) is {Cost_pump}, the cost of the electric generator is {Cost_generator}."
          f" and the cost of the cooling tower is {Cost_cooling_tower}.")