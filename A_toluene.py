import sys
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.service import Service
import time
from bs4 import BeautifulSoup
from turbine_efficiency_calculator import get_best_efficiency
from shellandtube_calculator import calculate_area_s
from plate_calculator import calculate_area_p
import shutil

# Detect paths dynamically
chromedriver_path = shutil.which("chromedriver")
chromium_path = shutil.which("chromium-browser")

print(f"ChromeDriver path: {chromedriver_path}")
print(f"Chromium path: {chromium_path}")

if not chromedriver_path:
    raise FileNotFoundError("ChromeDriver not found. Ensure it is installed and in the PATH.")

if not chromium_path:
    raise FileNotFoundError("Chromium not found. Ensure it is installed and in the PATH.")

# Configure Chrome options
options = Options()
options.add_argument("--headless")
options.add_argument("--no-sandbox")
options.add_argument("--disable-dev-shm-usage")
options.binary_location = chromium_path

# Initialize WebDriver
service = Service(chromedriver_path)

"""
# Initial variables definition
# Constants (remaining the same)
p2 = 10
T1 = 40
T4 = 240
# Efficiencies
# Pump
pump_isentropic_efficiency = 0.8
pump_electric_motor_efficiency = 0.9
pump_mechanical_efficiency = 0.95
# Turbine
turbine_isentropic_efficiency = 0.9
turbine_mechanical_efficiency = 0.92
# Others
electric_generator_efficiency = 0.96
recuperator_efficiency = 0.73
Q_th = 15e6

# Input data for the shell and tube HX
# Geometrical Parameters
D_t = 0.019  # Internal Tube Diameter in meters
D_s = 0.203  # Shell Diameter in meters
N_tubes = 24     # Number of Tubes
N_p = 2      # Number of Tube Passes
P_T = 0.238  # Square Pitch Length in meters
C = 0.005    # Tube Clearance in meters
B = 0.305    # Baffle Spacing in meters
k_SS = 16   # Thermal conductivity tube
Model = 10 # Number of segments

# Input data for the plate HX
L = 0.4  # Effective Length
W = 0.2  # Effective Width
t = 0.001  # Plate Thickness
Dp = 0.018  # Port Diameter
Lv = 0.01 # Vertical Port Distance
Lh = 0.02  # Horizontal Port Distance
lambda_ = 0.009  # Corrugation Wavelength
beta = 60  # Corrugation Angle
Lc = 2  # Compressed Plate Pack Length
Phi = 1.25  # Surface Enlargement Factor
k_SS = 16  # Stainless Steel Thermal Conductivity
N_plates = 30 # Number of plates
N_passes = 1 # Number of passes
beta = 60 # Beta angle
Model = 10 # Number of segments
# Economic parameters
K1 = 4.0336
K2 = 0.2341
K3 = 0.0497
B1 = 0.96
B2 = 1.21
FM = 2.9
FS = 1.7
C1 = -0.1250
C2 = 0.15361
C3 = -0.02861
"""
def toluene_orc(pressure_loss, p2, T1, T4,
    pump_isentropic_efficiency, pump_electric_motor_efficiency,
    pump_mechanical_efficiency, turbine_isentropic_efficiency,
    turbine_mechanical_efficiency, electric_generator_efficiency,
    recuperator_efficiency,
    D_t, D_s, N_tubes, N_t_passes, P_T, C, B, k_SS_shellandtube,
    L, W, t, Dp, Lv, Lh, lambda_, beta,
    Lc, Phi, k_SS_plate, N_plates, N_p_passes):


    #Definition of the variables outside the choice of the user
    Q_th = 15e6
    T1nd = 15
    T2nd = 25
    T1st = 190
    T2st = 310
    Model = 10 # Number of segments
    # Economic parameters
    K1 = 4.0336
    K2 = 0.2341
    K3 = 0.0497
    B1 = 0.96
    B2 = 1.21
    FM = 2.9
    FS = 1.7
    C1 = -0.1250
    C2 = 0.15361
    C3 = -0.02861

    #Links
    Link1 = "https://webbook.nist.gov/cgi/fluid.cgi?TUnit=C&PUnit=bar&DUnit=mol%2Fl&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm&Type=IsoTherm&RefState=DEF&Action=Page&ID=C108883"
    Link2 = "https://webbook.nist.gov/cgi/fluid.cgi?ID=C108883&TUnit=C&PUnit=bar&DUnit=mol%2Fl&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm&Type=IsoBar&RefState=DEF&Action=Page"

    # Set up the Selenium WebDriver in headless mode
    driver = webdriver.Chrome(service=service, options=options)

    # Open the NIST Fluid Properties web page
    driver.get(Link1)

    # Wait for the page to load fully
    WebDriverWait(driver, 60).until(
        EC.presence_of_element_located((By.NAME, "T"))
    )

    # Locate the pressure input field, clear it, and input the value of p2
    pressure_input = driver.find_element(By.NAME, "T")
    pressure_input.clear()
    pressure_input.send_keys(str(T1))

    # Submit the form
    pressure_input.send_keys(Keys.RETURN)

    # Wait for the page to load after form submission
    time.sleep(5)  # Adjust sleep time as needed or use WebDriverWait

    # Get the page source after form submission
    page_source = driver.page_source

    # Use BeautifulSoup to parse the HTML page source
    soup = BeautifulSoup(page_source, 'html.parser')

    # Find the table using BeautifulSoup
    table = soup.find('table', {'class': 'small', 'border': '1'})

    for row in table.find_all('tr')[1:]:  # Skip the header row
        cols = row.find_all('td')
        entropy = float(cols[6].text.strip())
        enthalpy = float(cols[5].text.strip())
        pressure = float(cols[1].text.strip())
        status = cols[13].text.strip()

        if status == "liquid":
            s1 = entropy
            h1 = enthalpy
            p1 = pressure
            break  # No need to continue after finding the required saturation data

    driver.close()  # Closes the current active tab

    # Save the variables considering a 1% pressure loss apart from turbine and pump
    p3 = p2 * (1 - pressure_loss)
    p4 = p3 * (1 - pressure_loss)
    p6 = (1 + pressure_loss) * p1
    p5 = (1 + pressure_loss) * p6
    s2is = s1

    # Set up the driver again to calculate properties at point 2
    ddriver = webdriver.Chrome(service=service, options=options)

    # Open the NIST Fluid Properties web page
    driver.get(Link2)

    # Wait for the page to load fully
    WebDriverWait(driver, 60).until(
        EC.presence_of_element_located((By.NAME, "P"))
    )

    # Locate the pressure input field, clear it, and input the value of p2
    pressure_input = driver.find_element(By.NAME, "P")
    pressure_input.clear()
    pressure_input.send_keys(str(p2))

    # Submit the form
    pressure_input.send_keys(Keys.RETURN)

    # Wait for the page to load after form submission
    time.sleep(10)  # Adjust sleep time as needed or use WebDriverWait

    # Get the page source after form submission
    page_source = driver.page_source

    # Use BeautifulSoup to parse the HTML page source
    soup = BeautifulSoup(page_source, 'html.parser')

    # Find the table using BeautifulSoup
    table = soup.find('table', {'class': 'small', 'border': '1'})

    # Extract the columns of the highest isobar
    data_s_high = []
    data_h_high = []
    data_T_high = []

    for row in table.find_all('tr')[1:]:  # Skip the header row
        cols = row.find_all('td')
        data_s_high.append(float(cols[6].text.strip()))
        data_h_high.append(float(cols[5].text.strip()))
        data_T_high.append(float(cols[0].text.strip()))

    driver.close()  # Closes the current active tab

    # Extract the data for the interpolation of 4
    i = None
    for idx1 in range(1, len(data_s_high)):
        if (s2is - data_s_high[idx1]) * (s2is - data_s_high[idx1 - 1]) < 0:
            i = idx1
            break

    #Calculate h2is
    h2is = (s2is - data_s_high[i - 1]) / (data_s_high[i] - data_s_high[i - 1]) * (data_h_high[i] - data_h_high[i - 1]) + data_h_high[i - 1]
    T2is = (s2is - data_s_high[i - 1]) / (data_s_high[i] - data_s_high[i - 1]) * (data_T_high[i] - data_T_high[i - 1]) + data_T_high[i - 1]

    # Calculate h2 using the definition
    h2 = h1 + (h2is - h1) / pump_isentropic_efficiency

    # Set up the driver again to calculate properties at point 4
    driver = webdriver.Chrome(service=service, options=options)

    # Open the NIST Fluid Properties web page
    driver.get(Link2)

    # Wait for the page to load fully
    WebDriverWait(driver, 60).until(
        EC.presence_of_element_located((By.NAME, "P"))
    )

    # Locate the pressure input field, clear it, and input the value of p2
    pressure_input = driver.find_element(By.NAME, "P")
    pressure_input.clear()
    pressure_input.send_keys(str(p4))

    # Submit the form
    pressure_input.send_keys(Keys.RETURN)

    # Wait for the page to load after form submission
    time.sleep(5)  # Adjust sleep time as needed or use WebDriverWait

    # Get the page source after form submission
    page_source = driver.page_source

    # Use BeautifulSoup to parse the HTML page source
    soup = BeautifulSoup(page_source, 'html.parser')

    # Find the table using BeautifulSoup
    table = soup.find('table', {'class': 'small', 'border': '1'})

    # Extract the columns of the evaporation isobar
    data_s_evap = []
    data_h_evap = []
    data_T_evap = []
    data_saturation_evap = []
    data_d_evap = []

    for row in table.find_all('tr')[1:]:  # Skip the header row
        cols = row.find_all('td')
        data_T_evap.append(float(cols[0].text.strip()))
        data_d_evap.append(float(cols[2].text.strip()))
        data_h_evap.append(float(cols[5].text.strip()))
        data_s_evap.append(float(cols[6].text.strip()))
        data_saturation_evap.append(cols[13].text.strip())

    driver.close()  # Closes the current active tab

    # Extract the data for the interpolation of 4
    i = None
    j = None
    for idx2 in range(1, len(data_T_evap)):
        if data_saturation_evap[idx2] == "vapor" and data_saturation_evap[idx2 - 1] == "liquid":
            i = idx2
        if (T4 - data_T_evap[idx2]) * (T4 - data_T_evap[idx2- 1]) < 0:
            j = idx2

    h4 = (T4 - data_T_evap[j - 1]) / (data_T_evap[j] - data_T_evap[j - 1]) * (data_h_evap[j] - data_h_evap[j - 1]) + data_h_evap[j - 1]
    s4 = (T4 - data_T_evap[j - 1]) / (data_T_evap[j] - data_T_evap[j - 1]) * (data_s_evap[j] - data_s_evap[j - 1]) + data_s_evap[j - 1]
    s5is = s4
    d4 = (T4 - data_T_evap[j - 1]) / (data_T_evap[j] - data_T_evap[j - 1]) * (data_d_evap[j] - data_d_evap[j - 1]) + data_d_evap[j - 1]
    saturation_enthalpy_liquid = data_h_evap[i - 1]
    saturation_enthalpy_gas = data_h_evap[i]
    saturation_temperature = data_T_evap[i]

    # Set up the driver again to calculate properties at point 5
    driver = webdriver.Chrome(service=service, options=options)

    # Open the NIST Fluid Properties web page
    driver.get(Link2)

    # Wait for the page to load fully
    WebDriverWait(driver, 60).until(
        EC.presence_of_element_located((By.NAME, "P"))
    )

    # Locate the pressure input field, clear it, and input the value of p5
    pressure_input = driver.find_element(By.NAME, "P")
    pressure_input.clear()
    pressure_input.send_keys(str(p5))

    # Submit the form
    pressure_input.send_keys(Keys.RETURN)

    # Wait for the page to load after form submission
    time.sleep(5)  # Adjust sleep time as needed or use WebDriverWait

    # Get the page source after form submission
    page_source = driver.page_source

    # Use BeautifulSoup to parse the HTML page source
    soup = BeautifulSoup(page_source, 'html.parser')

    # Find the table using BeautifulSoup
    table = soup.find('table', {'class': 'small', 'border': '1'})
    prev_entropy = 0
    prev_enthalpy = 0
    prev_temperature = 0
    prev_density = 0

    for row in table.find_all('tr')[1:]:  # Skip the header row
        cols = row.find_all('td')
        entropy = float(cols[6].text.strip())
        enthalpy = float(cols[5].text.strip())
        temperature = float(cols[0].text.strip())
        density = float(cols[2].text.strip())

        if s5is < entropy:
            h5is = (s5is - prev_entropy) / (entropy - prev_entropy) * (enthalpy - prev_enthalpy) + prev_enthalpy
            T5is = (s5is - prev_entropy) / (entropy - prev_entropy) * (temperature - prev_temperature) + prev_temperature
            d5is = (s5is - prev_entropy) / (entropy - prev_entropy) * (density - prev_density) + prev_density
            break

        prev_entropy = entropy
        prev_enthalpy = enthalpy
        prev_temperature = temperature
        prev_density = density

    h5 = h4 - turbine_isentropic_efficiency * (h4 - h5is)

    prev_entropy = 0
    prev_enthalpy = 0
    prev_temperature = 0

    for row in table.find_all('tr')[1:]:  # Skip the header row
        cols = row.find_all('td')
        entropy = float(cols[6].text.strip())
        enthalpy = float(cols[5].text.strip())
        temperature = float(cols[0].text.strip())

        if h5 < enthalpy:
            T5 = (h5 - prev_enthalpy) / (enthalpy - prev_enthalpy) * (temperature - prev_temperature) + prev_temperature
            s5 = (h5 - prev_enthalpy) / (enthalpy - prev_enthalpy) * (entropy - prev_entropy) + prev_entropy
            break

        prev_entropy = entropy
        prev_enthalpy = enthalpy
        prev_temperature = temperature

    driver.close()  # Closes the current active tab

    # Compare T5 with enthalpy values of the higher isobar to find hx
    i = None
    for idx3 in range(1, len(data_T_high)):
        if i is None and T5 < data_T_high[idx3]:
            i = idx3
            break  # No need to continue after finding i

    hx = (T5 - data_T_high[i - 1]) / (data_T_high[i] - data_T_high[i - 1]) * (data_h_high[i] - data_h_high[i - 1]) + data_h_high[i - 1]
    h6 = h5 - recuperator_efficiency * (hx - h2)
    h3 = h5 - h6 + h2

    # Set up the driver again to calculate properties at point 3
    driver = webdriver.Chrome(service=service, options=options)

    # Open the NIST Fluid Properties web page
    driver.get(Link2)

    # Wait for the page to load fully
    WebDriverWait(driver, 60).until(
        EC.presence_of_element_located((By.NAME, "P"))
    )

    # Locate the pressure input field, clear it, and input the value of p2
    pressure_input = driver.find_element(By.NAME, "P")
    pressure_input.clear()
    pressure_input.send_keys(str(p3))

    # Submit the form
    pressure_input.send_keys(Keys.RETURN)

    # Wait for the page to load after form submission
    time.sleep(5)  # Adjust sleep time as needed or use WebDriverWait

    # Get the page source after form submission
    page_source = driver.page_source

    # Use BeautifulSoup to parse the HTML page source
    soup = BeautifulSoup(page_source, 'html.parser')

    # Find the table using BeautifulSoup
    table = soup.find('table', {'class': 'small', 'border': '1'})
    prev_entropy = 0
    prev_enthalpy = 0
    prev_temperature = 0

    for row in table.find_all('tr')[1:]:  # Skip the header row
        cols = row.find_all('td')
        entropy = float(cols[6].text.strip())
        enthalpy = float(cols[5].text.strip())
        temperature = float(cols[0].text.strip())

        if h3 < enthalpy:
            T3 = (h3 - prev_enthalpy) / (enthalpy - prev_enthalpy) * (temperature - prev_temperature) + prev_temperature
            s3 = (h3 - prev_enthalpy) / (enthalpy - prev_enthalpy) * (entropy - prev_entropy) + prev_entropy
            break

        prev_entropy = entropy
        prev_enthalpy = enthalpy
        prev_temperature = temperature

    driver.close()  # Closes the current active tab

    mass_flow_rate = Q_th / (h4 - h3) / 1000
    VR = d4 / d5is
    d4_1 = 92.14 * d4
    VH = math.sqrt(mass_flow_rate/d4_1) / ((h4 - h5is) ** 0.25)

    # Call the turbine efficiency function
    max_eta, best_stage = get_best_efficiency(VH, VR, turbine_isentropic_efficiency)

    # Efficiency definition with related constraints
    if h3 > saturation_enthalpy_liquid and h3 < saturation_enthalpy_gas or T4 < saturation_temperature:
        efficiency = None
    else:
        turbine_power = (h4 - h5) * turbine_mechanical_efficiency * electric_generator_efficiency
        pump_power = (h2 - h1) / pump_mechanical_efficiency / pump_electric_motor_efficiency
        efficiency = ((mass_flow_rate * (turbine_power - pump_power) * 1000) / Q_th) * 100

    i = None
    j = None
    # Creation of the plots
    for i in range(1, len(data_T_evap)):
        if (T3 - data_T_evap[i-1]) * (T3 - data_T_evap[i]) == 0:
            j = i
        else:
            if (T3 - data_T_evap[i-1]) * (T3 - data_T_evap[i]) < 0:
                j = i
        if T4 - data_T_evap[i] == 0:
            w = i
        else:
            if (T4 - data_T_evap[i-1]) * (T4 - data_T_evap[i]) < 0:
                w = i - 1

    print(f"The cycle thermal efficiency is {efficiency}")

    Area_plate, pinchpoint, cold_T, cold_Q, hot_T, hot_Q = calculate_area_p(L, W, t, Dp, Lv, Lh, lambda_, beta, Lc, Phi, k_SS_plate, N_plates, N_p_passes, Model,
                                              saturation_temperature, p4, T3, T4, T1st, T2st,
                                              mass_flow_rate, Q_th, saturation_enthalpy_liquid, saturation_enthalpy_gas, Link2)
    if pinchpoint < 0:
        print(f"Error: the pinch point is negative: please change the system conditions.")
        efficiency, max_eta, best_stage, cold_T, cold_Q, hot_T, hot_Q, Area_plate, Area_shelltube, Cost_plate, Cost_shellandtube, Cost_turbine, Cost_pump, Cost_generator, Cost_cooling_tower = 0
        return efficiency, max_eta, best_stage, pinchpoint, cold_T, cold_Q, hot_T, hot_Q, Area_plate, Area_shelltube, Cost_plate, Cost_shellandtube, Cost_turbine, Cost_pump, Cost_generator, Cost_cooling_tower

    print(f"The total heat exchange area required for a plate HX is: {Area_plate}")
    print(f"Your pinch point is: {pinchpoint}")

    Area_shelltube = calculate_area_s(D_t, D_s, N_tubes, N_t_passes, P_T, C, B, k_SS_shellandtube, Model,
                                saturation_temperature, p4, T3, T4, T1st, T2st,
                                mass_flow_rate, Q_th, saturation_enthalpy_liquid, saturation_enthalpy_gas, Link2)
    print(f"The total heat exchange area required for a shell and tube HX is: {Area_shelltube}")

    # Cost of the HXs
    if Area_plate < 10000:
        C_HX_0 = 10 ** (K1 + K2 * math.log10(Area_plate) + K3 * (math.log10(Area_plate) ** 2))
    else:
        C_HX_0 = (Area_plate / 10000) * 10 ** (K1 + K2 * math.log10(10000) + K3 * (math.log10(10000) ** 2))

    log_P = math.log10(p4)
    FP = 10 ** (C1 + C2 * log_P + C3 * (log_P ** 2))
    Cost_plate = (525.7 / 397) * (B1 + B2 * FM * FP) * FS * C_HX_0
    Cost_shellandtube = 3.28e4 * (Area_shelltube / 4)**0.68

    if Cost_shellandtube > Cost_plate:
        print(f"The most economically convenient heat exchanger is the plate HX with a cost of {Cost_plate}")
    else:
        print(f"The most economically convenient heat exchanger is the shell and tube HX with a cost of {Cost_shellandtube}")

    # Cost of the rest
    Cost_turbine = -1.66e4 + 716 * (turbine_power)**0.8
    Cost_pump = 9.84e3 * (pump_power / 4)**0.55
    Cost_generator = 2447 * (turbine_power)**0.49
    condenser_power = mass_flow_rate * (h6 - h1)
    Cost_cooling_tower = 22.582 * condenser_power + 1924.6

    print(f"The cost of the turbine is {Cost_turbine}, the cost of the centrifugal pump (incl. motor) is {Cost_pump}, the cost of the electric generator is {Cost_generator}"
          f"and the cost of the cooling tower is {Cost_cooling_tower}")

    return efficiency, max_eta, best_stage, pinchpoint, cold_T, cold_Q, hot_T, hot_Q, Area_plate, Area_shelltube, Cost_plate, Cost_shellandtube, Cost_turbine, Cost_pump, Cost_generator, Cost_cooling_tower






