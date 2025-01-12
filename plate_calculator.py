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
from plate_hx_correlation_calculator import get_heat_transfer_coefficients
import shutil

# Detect paths dynamically
chromedriver_path = shutil.which("chromedriver")
chromium_path = shutil.which("chromium-browser")

if not chromedriver_path:
    raise FileNotFoundError("ChromeDriver not found. Ensure it is installed and in the PATH.")

if not chromium_path:
    raise FileNotFoundError("Chromium not found. Ensure it is installed and in the PATH.")

# Configure Chrome options
options = Options()
options.add_argument("--headless")
options.add_argument("--no-sandbox")
options.add_argument("--disable-dev-shm-usage")
options.add_argument("--disable-gpu")
options.binary_location = chromium_path

# Initialize WebDriver
service = Service(chromedriver_path)

def calculate_area_p(
    L, W, t, Dp, Lv, Lh, lambda_, beta, Lc, Phi, k_SS, N_t, N_p, Model,
    saturation_temperature, p4, T3, T4, T1st, T2st, mass_flow_rate_cold, Q_th, saturation_enthalpy_l, saturation_enthalpy_g, Link):

    saturation_enthalpy_l = saturation_enthalpy_l * 1000
    saturation_enthalpy_g = saturation_enthalpy_g * 1000

    # Initialize the driver
    driver = webdriver.Chrome(service=service, options=options)

    # Open the NIST Fluid Properties web page
    driver.get(Link)

    # Wait for the page to load fully
    WebDriverWait(driver, 60).until(
        EC.presence_of_element_located((By.NAME, "P"))
    )

    # Locate the temperature and pressure input fields
    pressure_input = driver.find_element(By.NAME, "P")
    # Locate the temperature input fields by their 'name' attributes
    TLow_input = driver.find_element(By.NAME, "TLow")
    THigh_input = driver.find_element(By.NAME, "THigh")

    # Clear and input values for temperatures and pressure
    TLow_input.clear()
    TLow_input.send_keys(str(T3))

    THigh_input.clear()
    THigh_input.send_keys(str(T4))

    pressure_input.clear()
    pressure_input.send_keys(str(p4))

    # Submit the form by pressing RETURN in one of the input fields
    pressure_input.send_keys(Keys.RETURN)

    # Wait for the page to load after form submission
    time.sleep(5)  # Adjust sleep time as needed or use WebDriverWait

    # Get the page source after form submission
    page_source = driver.page_source

    # Use BeautifulSoup to parse the HTML page source
    soup = BeautifulSoup(page_source, 'html.parser')

    # Find the table using BeautifulSoup
    table = soup.find('table', {'class': 'small', 'border': '1'})

    # Extract the columns of the highest isobar
    data_s_cold = []
    data_h_cold = []
    data_T_cold = []
    data_mu_cold = []
    data_k_cold = []
    data_cp_cold = []
    data_V_cold = []

    for row in table.find_all('tr')[1:]:  # Skip the header row
        cols = row.find_all('td')
        data_T_cold.append(float(cols[0].text.strip()))
        data_V_cold.append(float(cols[3].text.strip()))
        data_h_cold.append(float(cols[5].text.strip()))
        data_s_cold.append(float(cols[6].text.strip()))
        data_cp_cold.append(float(cols[8].text.strip()))
        data_mu_cold.append(float(cols[11].text.strip()))
        data_k_cold.append(float(cols[12].text.strip()))

    driver.close()  # Closes the current active tab

    # identification of the saturation variables
    i = None
    for i in range(0, len(data_cp_cold)):
        if saturation_temperature == data_T_cold[i]:
            liquid_k = data_k_cold[i]
            gas_k = data_k_cold[i + 1]
            liquid_mu = data_mu_cold[i]
            gas_mu = data_mu_cold[i + 1]
            liquid_cp = data_cp_cold[i]
            gas_cp = data_cp_cold[i + 1]
            break

    # Save the Heat values
    data_h_cold = np.array(data_h_cold)
    data_h_cold = data_h_cold * 1000
    data_Q_cold = (data_h_cold - data_h_cold[0]) * mass_flow_rate_cold
    saturation_enthalpy_l = (saturation_enthalpy_l - data_h_cold[0]) * mass_flow_rate_cold
    saturation_enthalpy_g = (saturation_enthalpy_g - data_h_cold[0]) * mass_flow_rate_cold

    # Liquid part
    m_l = (saturation_temperature - T3) / (saturation_enthalpy_l - data_Q_cold[0])
    c_l = T3 - m_l * data_Q_cold[0]

    # Gas part
    m_g = (T4 - saturation_temperature) / (data_Q_cold[-1] - saturation_enthalpy_g)
    c_g = T4 - m_g * data_Q_cold[-1]

    # Calculation of the temperature vectors for the step by step cycle
    cold_liquid_Q = np.linspace(data_Q_cold[0], saturation_enthalpy_l, Model)
    cold_liquid_temperatures = m_l * cold_liquid_Q + c_l
    cold_gas_Q = np.linspace(saturation_enthalpy_g, data_Q_cold[-1], Model)
    cold_gas_temperatures = m_g * cold_gas_Q + c_g
    cold_saturation_temperatures = np.full(Model, saturation_temperature)
    cold_saturation_Q = np.linspace(saturation_enthalpy_l, saturation_enthalpy_g, Model)

    # Initialize the driver
    driver = webdriver.Chrome(service=service, options=options)

    # Open the NIST Fluid Properties web page for water
    driver.get(
        "https://webbook.nist.gov/cgi/fluid.cgi?ID=C7732185&TUnit=C&PUnit=bar&DUnit=mol%2Fl&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm&Type=IsoBar&RefState=DEF&Action=Page")

    # Wait for the page to load fully
    WebDriverWait(driver, 60).until(
        EC.presence_of_element_located((By.NAME, "P")))

    # Locate the temperature and pressure input fields
    pressure_input = driver.find_element(By.NAME, "P")
    # Locate the temperature input fields by their 'name' attributes
    TLow_input = driver.find_element(By.NAME, "TLow")
    THigh_input = driver.find_element(By.NAME, "THigh")

    # Clear and input values for temperatures and pressure
    TLow_input.clear()
    TLow_input.send_keys(str(T1st))

    THigh_input.clear()
    THigh_input.send_keys(str(T2st))

    pressure_input.clear()
    pressure_input.send_keys(str(120))

    # Submit the form by pressing RETURN in one of the input fields
    pressure_input.send_keys(Keys.RETURN)

    # Wait for the page to load after form submission
    time.sleep(5)  # Adjust sleep time as needed or use WebDriverWait

    # Get the page source after form submission
    page_source = driver.page_source

    # Use BeautifulSoup to parse the HTML page source
    soup = BeautifulSoup(page_source, 'html.parser')

    # Find the table using BeautifulSoup
    table = soup.find('table', {'class': 'small', 'border': '1'})

    # Extract the columns of the highest isobar
    data_s_hot = []
    data_h_hot = []
    data_T_hot = []
    data_mu_hot = []
    data_k_hot = []
    data_cp_hot = []
    data_V_hot = []

    for row in table.find_all('tr')[1:]:  # Skip the header row
        cols = row.find_all('td')
        data_T_hot.append(float(cols[0].text.strip()))
        data_V_hot.append(float(cols[3].text.strip()))
        data_h_hot.append(float(cols[5].text.strip()))
        data_s_hot.append(float(cols[6].text.strip()))
        data_cp_hot.append(float(cols[8].text.strip()))
        data_mu_hot.append(float(cols[11].text.strip()))
        data_k_hot.append(float(cols[12].text.strip()))

    driver.close()  # Closes the current active tab

    data_h_hot = np.array(data_h_hot)
    data_h_hot = data_h_hot * 1000
    mass_flow_rate_hot = (mass_flow_rate_cold * (data_h_cold[-1] - data_h_cold[0])) / (data_h_hot[-1] - data_h_hot[0])
    data_Q_hot = (data_h_hot - data_h_hot[0]) * mass_flow_rate_hot

    # Hot fluid side
    m_h = (T2st - T1st) / (data_Q_hot[-1] - data_Q_hot[0])
    c_h = T2st - m_h * data_Q_hot[-1]

    # Calculation of the temperature vectors for the step by step cycle
    hot_liquid_Q = np.linspace(data_Q_hot[0], saturation_enthalpy_l, Model)
    hot_gas_Q = np.linspace(saturation_enthalpy_g, data_Q_hot[-1], Model)
    hot_saturation_Q = np.linspace(saturation_enthalpy_l, saturation_enthalpy_g, Model)
    hot_liquid_temperatures = m_h * hot_liquid_Q + c_h
    hot_gas_temperatures = m_h * hot_gas_Q + c_h
    hot_saturation_temperatures = m_h * hot_saturation_Q + c_h

    # Calculation of the pinch point
    # Calculate differences for each type
    liquid_difference = hot_liquid_temperatures - cold_liquid_temperatures
    saturation_difference = hot_saturation_temperatures - cold_saturation_temperatures
    gas_difference = hot_gas_temperatures - cold_gas_temperatures

    # Combine all differences into one vector
    T_difference = np.concatenate([liquid_difference, saturation_difference, gas_difference])

    # Calculate the pinch point as the minimum difference
    pinch_point = min(T_difference)

    """
    # Check plots
    # Plot t2_vector against h2_vector
    plt.plot(hot_liquid_Q, hot_liquid_temperatures, linestyle='-', label='T vs Q')
    plt.plot(hot_gas_Q, hot_gas_temperatures, linestyle='-', label='T vs Q')
    plt.plot(hot_saturation_Q, hot_saturation_temperatures, linestyle='-', label='T vs Q')
    plt.plot(cold_liquid_Q, cold_liquid_temperatures, linestyle='-', label='T vs Q')
    plt.plot(cold_gas_Q, cold_gas_temperatures, linestyle='-', label='T vs Q')
    plt.plot(cold_saturation_Q, cold_saturation_temperatures, linestyle='-', label='T vs Q')

    # Add labels and title
    plt.xlabel('Energy exchanged')
    plt.ylabel('Temperature')
    plt.title('Temperature vs. Energy Exchanged')

    # Add legend and grid
    plt.legend()
    plt.grid(True)

    # Display the plot
    plt.show()
    """

    # Plot creation
    hot_Q = np.concatenate([hot_liquid_Q, hot_saturation_Q, hot_gas_Q])
    hot_T = np.concatenate([hot_liquid_temperatures, hot_saturation_temperatures, hot_gas_temperatures])
    cold_Q = np.concatenate([cold_liquid_Q, cold_saturation_Q, cold_gas_Q])
    cold_T = np.concatenate([cold_liquid_temperatures, cold_saturation_temperatures, cold_gas_temperatures])

    # Calculation of relevant parameters for the for cycle
    p = Lc / N_t
    b = p - t
    De = 2 * b / Phi
    Lw = Lh + Dp
    Ach = b * Lw
    N_cp = (N_t - 1) / 2 * N_p

    # Cycle over the liquid sides on the left
    Areas_liquid = []
    i = None
    for i in range(0, len(hot_liquid_Q) - 1):
        Delta_T1 = hot_liquid_temperatures[i] - cold_liquid_temperatures[i]
        Delta_T2 = hot_liquid_temperatures[i + 1] - cold_liquid_temperatures[i + 1]
        Delta_Tlog = (Delta_T1 - Delta_T2) / math.log(Delta_T1 / Delta_T2)
        middle_h = (hot_liquid_Q[i] + hot_liquid_Q[i + 1]) / 2
        j = None
        for j in range(0, len(data_Q_cold) - 1):
            if (data_Q_cold[j] - middle_h) * (data_Q_cold[j + 1] - middle_h) < 0:
                k_cold = (middle_h - data_Q_cold[j]) / (data_Q_cold[j + 1] - data_Q_cold[j]) * (
                            data_k_cold[j + 1] - data_k_cold[j]) + data_k_cold[j]
                cp_cold = (middle_h - data_Q_cold[j]) / (data_Q_cold[j + 1] - data_Q_cold[j]) * (
                            data_cp_cold[j + 1] - data_cp_cold[j]) + data_cp_cold[j]
                mu_cold = (middle_h - data_Q_cold[j]) / (data_Q_cold[j + 1] - data_Q_cold[j]) * (
                            data_mu_cold[j + 1] - data_mu_cold[j]) + data_mu_cold[j]
                break
        j = None
        for j in range(0, len(data_Q_hot) - 1):
            if (data_Q_hot[j] - middle_h) * (data_Q_hot[j + 1] - middle_h) < 0:
                k_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (
                            data_k_hot[j + 1] - data_k_hot[j]) + data_k_hot[j]
                cp_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (
                            data_cp_hot[j + 1] - data_cp_hot[j]) + data_cp_hot[j]
                mu_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (
                            data_mu_hot[j + 1] - data_mu_hot[j]) + data_mu_hot[j]
                break
        Gc_hot = mass_flow_rate_hot / (Ach * N_cp)
        Gc_cold = mass_flow_rate_cold / (Ach * N_cp)
        Re_hot = Gc_hot * De / mu_hot
        Re_cold = Gc_cold * De / mu_cold
        Pr_hot = cp_hot * 1000 * mu_hot / k_hot
        Pr_cold = cp_cold * 1000 * mu_cold / k_cold
        C_h, n = get_heat_transfer_coefficients(beta, Re_hot)
        h_hot = k_hot * C_h * Re_hot ** n * Pr_hot ** (1 / 3) / De
        C_h, n = get_heat_transfer_coefficients(beta, Re_cold)
        h_cold = k_cold * C_h * Re_cold ** n * Pr_cold ** (1 / 3) / De
        U = (1 / h_hot + t / k_SS + 1 / h_cold) ** (-1)
        Areas_liquid.append((hot_liquid_Q[i + 1] - hot_liquid_Q[i]) / U / Delta_Tlog)

    # Cycle over the gas sides on the right
    Areas_gas = []
    i = None
    for i in range(0, len(hot_gas_Q) - 1):
        Delta_T1 = hot_gas_temperatures[i] - cold_gas_temperatures[i]
        Delta_T2 = hot_gas_temperatures[i + 1] - cold_gas_temperatures[i + 1]
        Delta_Tlog = (Delta_T1 - Delta_T2) / math.log(Delta_T1 / Delta_T2)
        middle_h = (hot_gas_Q[i] + hot_gas_Q[i + 1]) / 2
        j = None
        for j in range(0, len(data_Q_cold) - 1):
            if (data_Q_cold[j] - middle_h) * (data_Q_cold[j + 1] - middle_h) < 0:
                k_cold = (middle_h - data_Q_cold[j]) / (data_Q_cold[j + 1] - data_Q_cold[j]) * (
                            data_k_cold[j + 1] - data_k_cold[j]) + data_k_cold[j]
                cp_cold = (middle_h - data_Q_cold[j]) / (data_Q_cold[j + 1] - data_Q_cold[j]) * (
                            data_cp_cold[j + 1] - data_cp_cold[j]) + data_cp_cold[j]
                mu_cold = (middle_h - data_Q_cold[j]) / (data_Q_cold[j + 1] - data_Q_cold[j]) * (
                            data_mu_cold[j + 1] - data_mu_cold[j]) + data_mu_cold[j]
        j = None
        for j in range(0, len(data_Q_hot) - 1):
            if (data_Q_hot[j] - middle_h) * (data_Q_hot[j + 1] - middle_h) < 0:
                k_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (
                            data_k_hot[j + 1] - data_k_hot[j]) + data_k_hot[j]
                cp_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (
                            data_cp_hot[j + 1] - data_cp_hot[j]) + data_cp_hot[j]
                mu_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (
                            data_mu_hot[j + 1] - data_mu_hot[j]) + data_mu_hot[j]
        Gc_hot = mass_flow_rate_hot / (Ach * N_cp)
        Gc_cold = mass_flow_rate_cold / (Ach * N_cp)
        Re_hot = Gc_hot * De / mu_hot
        Re_cold = Gc_cold * De / mu_cold
        Pr_hot = cp_hot * 1000 * mu_hot / k_hot
        Pr_cold = cp_cold * 1000 * mu_cold / k_cold
        C_h, n = get_heat_transfer_coefficients(beta, Re_hot)
        h_hot = k_hot * C_h * Re_hot ** n * Pr_hot ** (1 / 3) / De
        C_h, n = get_heat_transfer_coefficients(beta, Re_cold)
        h_cold = k_cold * C_h * Re_cold ** n * Pr_cold ** (1 / 3) / De
        U = (1 / h_hot + t / k_SS + 1 / h_cold) ** (-1)
        Areas_gas.append((hot_gas_Q[i + 1] - hot_gas_Q[i]) / U / Delta_Tlog)

    # Cycle over the saturation sides in the middle
    Areas_saturation = []
    i = None
    for i in range(0, len(hot_saturation_Q) - 1):
        Delta_T1 = hot_saturation_temperatures[i] - cold_saturation_temperatures[i]
        Delta_T2 = hot_saturation_temperatures[i + 1] - cold_saturation_temperatures[i + 1]
        Delta_Tlog = (Delta_T1 - Delta_T2) / math.log(Delta_T1 / Delta_T2)
        middle_h = (hot_saturation_Q[i] + hot_saturation_Q[i + 1]) / 2
        k_cold = (middle_h - saturation_enthalpy_l) / (saturation_enthalpy_g - saturation_enthalpy_l) * (
                    gas_k - liquid_k) + liquid_k
        cp_cold = (middle_h - saturation_enthalpy_l) / (saturation_enthalpy_g - saturation_enthalpy_l) * (
                    gas_cp - liquid_cp) + liquid_cp
        mu_cold = (middle_h - saturation_enthalpy_l) / (saturation_enthalpy_g - saturation_enthalpy_l) * (
                    gas_mu - liquid_mu) + liquid_mu
        j = None
        for j in range(0, len(data_Q_hot) - 1):
            if (data_Q_hot[j] - middle_h) * (data_Q_hot[j + 1] - middle_h) < 0:
                k_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (
                            data_k_hot[j + 1] - data_k_hot[j]) + data_k_hot[j]
                cp_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (
                            data_cp_hot[j + 1] - data_cp_hot[j]) + data_cp_hot[j]
                mu_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (
                            data_mu_hot[j + 1] - data_mu_hot[j]) + data_mu_hot[j]
        Gc_hot = mass_flow_rate_hot / (Ach * N_cp)
        Gc_cold = mass_flow_rate_cold / (Ach * N_cp)
        Re_hot = Gc_hot * De / mu_hot
        Re_cold = Gc_cold * De / mu_cold
        Pr_hot = cp_hot * 1000 * mu_hot / k_hot
        Pr_cold = cp_cold * 1000 * mu_cold / k_cold
        C_h, n = get_heat_transfer_coefficients(beta, Re_hot)
        h_hot = k_hot * C_h * Re_hot ** n * Pr_hot ** (1 / 3) / De
        h_cold = k_cold * 0.26 * Re_cold ** 0.65 * Pr_cold ** 0.4 / De
        U = (1 / h_hot + t / k_SS + 1 / h_cold) ** (-1)
        Areas_saturation.append((hot_saturation_Q[i + 1] - hot_saturation_Q[i]) / U / Delta_Tlog)

    Area = np.sum(Areas_gas) + np.sum(Areas_liquid) + np.sum(Areas_saturation)
    return Area, pinch_point, cold_T, cold_Q, hot_T, hot_Q
