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
import shutil
import os

# Detect paths dynamically
chromedriver_path = os.getenv("CHROMEDRIVER_PATH", shutil.which("chromedriver"))
chromium_path = os.getenv("CHROMIUM_PATH", shutil.which("chromium-browser"))

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

def calculate_area_s(D_t, D_s, N_t, N_p, P_T, C, B, k_SS, Model,
                   saturation_temperature, p2, T3, T4, T1st, T2st,
                   mass_flow_rate_cold, Q_th, saturation_enthalpy_l, saturation_enthalpy_g, Link):
    """Calculate the total heat exchange area."""

    saturation_enthalpy_l = saturation_enthalpy_l * 1000
    saturation_enthalpy_g = saturation_enthalpy_g * 1000

    # Set up the Selenium WebDriver in headless mode
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
    pressure_input.send_keys(str(p2))

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
    data_d_cold = []

    for row in table.find_all('tr')[1:]:  # Skip the header row
        cols = row.find_all('td')
        data_T_cold.append(float(cols[0].text.strip()))
        data_d_cold.append(float(cols[2].text.strip()))
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
            liquid_d= data_d_cold[i]
            gas_d = data_d_cold[i + 1]
            break

    #Save the Heat values
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
    cold_gas_Q = np.linspace(saturation_enthalpy_g , data_Q_cold[-1], Model)
    cold_gas_temperatures = m_g * cold_gas_Q + c_g
    cold_saturation_temperatures = np.full(Model, saturation_temperature)
    cold_saturation_Q = np.linspace(saturation_enthalpy_l, saturation_enthalpy_g, Model)

    # Initialize the driver
    driver = webdriver.Chrome(service=service, options=options)

    # Open the NIST Fluid Properties web page for water
    driver.get("https://webbook.nist.gov/cgi/fluid.cgi?ID=C7732185&TUnit=C&PUnit=bar&DUnit=mol%2Fl&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm&Type=IsoBar&RefState=DEF&Action=Page")

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
    data_d_hot = []

    for row in table.find_all('tr')[1:]:  # Skip the header row
        cols = row.find_all('td')
        data_T_hot.append(float(cols[0].text.strip()))
        data_d_hot.append(float(cols[2].text.strip()))
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

    # Calculation of relevant parameters for the for cycle
    #U_t = (4 * mass_flow_rate_hot * N_p / N_t) / (math.pi * D_t ** 2 * )
    D_o = D_t + C
    a_s = C * B * D_s / P_T
    D_e = 4 * (P_T ** 2 - math.pi / 4 * D_o ** 2) / (math.pi * D_o)
    # D_e = 1.10 / D_o * (P_T ** 2 - 0.917 * D_o ** 2)
    D_m = (D_t - D_o) / math.log(D_t / D_o)
    t = D_o - D_t
    R_t = 0
    R_o = 0

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
                k_cold = (middle_h - data_Q_cold[j]) / (data_Q_cold[j + 1] - data_Q_cold[j]) * (data_k_cold[j + 1] - data_k_cold[j]) + data_k_cold[j]
                cp_cold = (middle_h - data_Q_cold[j]) / (data_Q_cold[j + 1] - data_Q_cold[j]) * (data_cp_cold[j + 1] - data_cp_cold[j]) + data_cp_cold[j]
                mu_cold = (middle_h - data_Q_cold[j]) / (data_Q_cold[j + 1] - data_Q_cold[j]) * (data_mu_cold[j + 1] - data_mu_cold[j]) + data_mu_cold[j]
                d_cold = (middle_h - data_Q_cold[j]) / (data_Q_cold[j + 1] - data_Q_cold[j]) * (data_d_cold[j + 1] - data_d_cold[j]) + data_d_cold[j]
                break
        j = None
        for j in range(0, len(data_Q_hot) - 1):
            if (data_Q_hot[j] - middle_h) * (data_Q_hot[j + 1] - middle_h) < 0:
                k_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (data_k_hot[j + 1] - data_k_hot[j]) + data_k_hot[j]
                cp_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (data_cp_hot[j + 1] - data_cp_hot[j]) + data_cp_hot[j]
                mu_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (data_mu_hot[j + 1] - data_mu_hot[j]) + data_mu_hot[j]
                d_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (data_d_hot[j + 1] - data_d_hot[j]) + data_d_hot[j]
                break
        U_t = (4 * mass_flow_rate_hot * N_p / N_t) / (math.pi * d_hot * D_t ** 2)
        Re_hot = d_hot * U_t * D_t / mu_hot
        Gc_cold = mass_flow_rate_cold / a_s
        Re_cold = D_e * Gc_cold / mu_cold
        Pr_hot = cp_hot * 1000 * mu_hot / k_hot
        Pr_cold = cp_cold * 1000 * mu_cold / k_cold
        #f_cold = (0.790 * math.log(Re_cold) - 1.64) ** (-2)
        f_hot = (0.790 * math.log(Re_hot) - 1.64) ** (-2)
        #Nu_cold = ((f_cold / 8) * (Re_cold - 1000) * Pr_cold) / (1 + 12.7 * (f_cold / 8)**0.5 * (Pr_cold**(2/3) - 1))
        Nu_hot = ((f_hot / 8) * (Re_hot - 1000) * Pr_hot) / (1 + 12.7 * (f_hot / 8) ** 0.5 * (Pr_hot ** (2 / 3) - 1))
        #h_cold = k_cold * Nu_cold / D_e
        h_hot = k_hot * Nu_hot / D_t
        #h_hot = k_hot * 0.023 * Re_hot ** 0.8 * Pr_hot ** 0.4 / D_t
        h_cold = k_cold * 0.36 * Re_cold ** 0.55 * Pr_cold ** (1 / 3) / D_e
        U = (D_o / (h_hot * D_t) + 1 / h_cold + t * D_o / (k_SS * D_m) + R_o + R_t * D_o / D_t)**-1
        Area = (hot_liquid_Q[i + 1] - hot_liquid_Q[i]) / U / Delta_Tlog
        Areas_liquid.append(Area)

    # Cycle over the saturation sides in the middle
    Areas_saturation = []
    i = None
    for i in range(0, len(hot_saturation_Q) - 1):
        Delta_T1 = hot_saturation_temperatures[i] - cold_saturation_temperatures[i]
        Delta_T2 = hot_saturation_temperatures[i + 1] - cold_saturation_temperatures[i + 1]
        Delta_Tlog = (Delta_T1 - Delta_T2) / math.log(Delta_T1 / Delta_T2)
        middle_h = (hot_saturation_Q[i] + hot_saturation_Q[i + 1]) / 2
        k_cold = (middle_h - saturation_enthalpy_l) / (saturation_enthalpy_g - saturation_enthalpy_l) * (gas_k - liquid_k) + liquid_k
        cp_cold = (middle_h - saturation_enthalpy_l) / (saturation_enthalpy_g - saturation_enthalpy_l) * (gas_cp - liquid_cp) + liquid_cp
        mu_cold = (middle_h - saturation_enthalpy_l) / (saturation_enthalpy_g - saturation_enthalpy_l) * (gas_mu - liquid_mu) + liquid_mu
        quality = 1 - (middle_h - saturation_enthalpy_l) / (saturation_enthalpy_g - saturation_enthalpy_l)
        j = None
        for j in range(0, len(data_Q_hot) - 1):
            if (data_Q_hot[j] - middle_h) * (data_Q_hot[j + 1] - middle_h) < 0:
                k_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (data_k_hot[j + 1] - data_k_hot[j]) + data_k_hot[j]
                cp_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (data_cp_hot[j + 1] - data_cp_hot[j]) + data_cp_hot[j]
                mu_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (data_mu_hot[j + 1] - data_mu_hot[j]) + data_mu_hot[j]
        U_t = (4 * mass_flow_rate_hot * N_p / N_t) / (math.pi * d_hot * D_t ** 2)
        Re_hot = d_hot * U_t * D_t / mu_hot
        Gc_cold = mass_flow_rate_cold / a_s
        Re_cold = D_e * Gc_cold / mu_cold
        Pr_hot = cp_hot * 1000 * mu_hot / k_hot
        Pr_cold = cp_cold * 1000 * mu_cold / k_cold
        f_hot = (0.790 * math.log(Re_hot) - 1.64) ** (-2)
        Nu_hot = ((f_hot / 8) * (Re_hot - 1000) * Pr_hot) / (1 + 12.7 * (f_hot / 8) ** 0.5 * (Pr_hot ** (2 / 3) - 1))
        h_hot = k_hot * Nu_hot / D_t
        # h_hot = k_hot * 0.023 * Re_hot ** 0.8 * Pr_hot ** 0.4 / D_t
        h_f = 0.021 * Re_cold ** 0.8 * Pr_cold ** 0.43 * k_cold / D_e
        h_cold = h_f * (1 + quality * (liquid_d / gas_d - 1)) ** 0.5
        U = (D_o / (h_hot * D_t) + 1 / h_cold + t * D_o / (k_SS * D_m) + R_o + R_t * D_o / D_t)**-1
        Area = (hot_saturation_Q[i + 1] - hot_saturation_Q[i]) / U / Delta_Tlog
        Areas_saturation.append(Area)

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
                k_cold = (middle_h - data_Q_cold[j]) / (data_Q_cold[j + 1] - data_Q_cold[j]) * (data_k_cold[j + 1] - data_k_cold[j]) + data_k_cold[j]
                cp_cold = (middle_h - data_Q_cold[j]) / (data_Q_cold[j + 1] - data_Q_cold[j]) * (data_cp_cold[j + 1] - data_cp_cold[j]) + data_cp_cold[j]
                mu_cold = (middle_h - data_Q_cold[j]) / (data_Q_cold[j + 1] - data_Q_cold[j]) * (data_mu_cold[j + 1] - data_mu_cold[j]) + data_mu_cold[j]
                d_cold = (middle_h - data_Q_cold[j]) / (data_Q_cold[j + 1] - data_Q_cold[j]) * (data_d_cold[j + 1] - data_d_cold[j]) + data_d_cold[j]
                break
        j = None
        for j in range(0, len(data_Q_hot) - 1):
            if (data_Q_hot[j] - middle_h) * (data_Q_hot[j + 1] - middle_h) < 0:
                k_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (data_k_hot[j + 1] - data_k_hot[j]) + data_k_hot[j]
                cp_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (data_cp_hot[j + 1] - data_cp_hot[j]) + data_cp_hot[j]
                mu_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (data_mu_hot[j + 1] - data_mu_hot[j]) + data_mu_hot[j]
                d_hot = (middle_h - data_Q_hot[j]) / (data_Q_hot[j + 1] - data_Q_hot[j]) * (data_d_hot[j + 1] - data_d_hot[j]) + data_d_hot[j]
                break
        U_t = (4 * mass_flow_rate_hot * N_p / N_t) / (math.pi * d_hot * D_t ** 2)
        Re_hot = d_hot * U_t * D_t / mu_hot
        Gc_cold = mass_flow_rate_cold / a_s
        Re_cold = D_e * Gc_cold / mu_cold
        Pr_hot = cp_hot * 1000 * mu_hot / k_hot
        Pr_cold = cp_cold * 1000 * mu_cold / k_cold
        f_hot = (0.790 * math.log(Re_hot) - 1.64) ** (-2)
        Nu_hot = ((f_hot / 8) * (Re_hot - 1000) * Pr_hot) / (1 + 12.7 * (f_hot / 8) ** 0.5 * (Pr_hot ** (2 / 3) - 1))
        h_hot = k_hot * Nu_hot / D_t
        # h_hot = k_hot * 0.023 * Re_hot ** 0.8 * Pr_hot ** 0.4 / D_t
        h_cold = k_cold * 0.36 * Re_cold ** 0.55 * Pr_cold ** (1 / 3) / D_e
        U = (D_o / (h_hot * D_t) + 1 / h_cold + t * D_o / (k_SS * D_m) + R_o + R_t * D_o / D_t)**-1
        Area = (hot_gas_Q[i + 1] - hot_gas_Q[i]) / U / Delta_Tlog
        Areas_gas.append(Area)

    Area = np.sum(Areas_gas) + np.sum(Areas_liquid) + np.sum(Areas_saturation)
    return Area

