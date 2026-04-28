# written by Giorgia Trabucco
#---------------------------------------------

import os
import json
import math
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from matplotlib.ticker import MultipleLocator

#---------------------------------------------

# Define constant function
def constant(x, c):
    return c

#---------------------------------------------

# temperature analysis
runs = [10586, 10587]

err_temp = []

for run in runs:
    path = f".../run_{run}"
    if run == runs[0]: state = "off"
    else: state = "on"

    indices = range(0, 100001, 1000)
    all_dfs = []  # to store each DataFrame
    rows = []

    for i in indices:
        file = f"{path}/temperatures_{i}.json"
        if not os.path.exists(file):
            print(f"Skipping missing file: {file}")
            continue

        with open(file) as f:
            data = json.load(f)
        
        dm08_chip1 = {'DM08': data['DM08']['1'], 'chip': 1}
        # Convert to DataFrame
        df = pd.DataFrame([dm08_chip1])

        """# remove unwanted keys
        for unwanted in ['CCBOARD', 'PCCA', 'PCCB']:
            data.pop(unwanted, None)

        # create DataFrame from this file
        #df = pd.DataFrame(data)
        df = pd.DataFrame(data).reset_index().rename(columns={'index': 'chip'})"""
        
        # add time/index info if you want to track which file it came from
        df["step"] = i  

        # move the step column to front
        df = df[["step"] + [col for col in df.columns if col != "step"]]

        all_dfs.append(df)

        for dm, values in data.items():

            if dm in ["CCBOARD", "PCCA", "PCCB"]:
                continue

            for chip, temp in values.items():
                rows.append({
                    "run": run,
                    "step": i,
                    "DM": dm,
                    "chip": int(chip),
                    "temperature": temp
                })
    df_tot = pd.DataFrame(rows)
    #print(df_tot)

    # Combine all into a single DataFrame
    combined_df = pd.concat(all_dfs, ignore_index=True)

    #-----------------------------------

    # Convert from wide to long format
    df_long = combined_df.melt(
        id_vars=['step', 'chip'],
        var_name='DM',
        value_name='temperature'
    )

    #mean and std dev
    mean_temp = df_long['temperature'].mean()
    std_temp = df_long['temperature'].std()
    """mean_temp = df_tot['temperature'].mean()
    std_temp = df_tot['temperature'].std()"""
    err_temp.append(std_temp)
    #err_temp.append(std_temp / math.sqrt(1000.))

    print("Mean temperature:", mean_temp)
    print("Standard deviation:", std_temp)
    print("temperature error:", err_temp[-1])

    # Compute min and max for each DM + chip combination
    stats = (
        df_long
        .groupby(['DM', 'chip'])['temperature']
        .agg(['min', 'max'])
        .reset_index()
    )
    #print(stats)

    # Plot
    plt.figure(figsize=(10,6))

    # Create one line per (DM, chip) combination
    for (dm, chip), group in df_long.groupby(['DM', 'chip']):
        plt.plot(group['step'], group['temperature'], label=f"{dm} - chip{chip}")

        # Get min and max for this DM + chip
        row = stats[(stats['DM'] == dm) & (stats['chip'] == chip)].iloc[0]
        min_val, max_val = row['min'], row['max']


    plt.xlabel("Cycles", fontsize=22)
    plt.ylabel("Temperature (°C)", fontsize=22)
    plt.title(f"{dm} - chip{chip} - TEC{state} - max={row['max']:.2f}°C, min={row['min']:.2f}°C", fontsize=20)
    
    """plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    # Create text block (min/max summary)
    for _, row in stats.iterrows():
        text_summary = f"max={row['max']:.2f}°C\nmin={row['min']:.2f}°C" """

    # Add text box under the legend
    """plt.gcf().text(
        1.05, 0.5, text_summary, 
        fontsize=10, 
        va='top', 
        ha='left', 
        transform=plt.gca().transAxes,
        bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.4')
    )"""
    plt.tight_layout()
    plt.savefig(f".../temps_{state}.png", dpi=300)

    data_dm08 = {'DM08': combined_df['DM08']}
    df = pd.DataFrame(data_dm08).reset_index().rename(columns={'index': 'chip'})
    #print(df)

#---------------------------------------------
#---------------------------------------------

# this script takes data from one sm with TECs on and off,
# the point is to prove that the time resolution isn't affected by some find of noise due to TECs.    

file_off = ".../plots_run_10586/peaks.pq"
file_on = ".../plots_run_10587/peaks.pq"

df_off = pd.read_parquet(file_off)
df_on = pd.read_parquet(file_on)

#------------------------------------

# sigma difference: I want it to be compatible with 0

sigma_off = df_off[df_off["chip"] == 1][["bar", "num_events", "sigma_dt3", "sigma_err_dt3"]]
sigma_on = df_on[df_on["chip"] == 1][["bar", "num_events", "sigma_dt3", "sigma_err_dt3"]]

merged = pd.merge(sigma_off, sigma_on, on="bar", suffixes=("_off", "_on"))

# Compute difference
merged["delta_sigma"] = merged["sigma_dt3_on"] - merged["sigma_dt3_off"]

# Propagate the error
merged["delta_sigma_err"] = np.sqrt(merged["sigma_err_dt3_on"]**2 + merged["sigma_err_dt3_off"]**2)
#merged["delta_sigma_err"] = np.sqrt((merged["sigma_dt3_on"]**2 / merged["num_events_on"]) + (merged["sigma_dt3_off"]**2 / merged["num_events_off"]))

# Final result
result = merged[["bar", "delta_sigma", "delta_sigma_err"]]

print(result)

# Extract x, y, and error
x = result["bar"]
y = result["delta_sigma"]
yerr = result["delta_sigma_err"]

# Perform the fit with error bars as weights
popt, pcov = curve_fit(constant, x, y, sigma=yerr, absolute_sigma=True)
c_fit = popt[0]
c_err = np.sqrt(pcov[0][0])

# Plot data with error bars
plt.figure(figsize=(10, 6))
plt.errorbar(x, y, yerr=yerr, fmt='o', capsize=4, label="Data")

# Plot fitted constant line
plt.axhline(c_fit, color='red', linestyle='--', label=f"Fit: c = {c_fit:.1e} ± {c_err:.1e}")

# Customize the plot
plt.gca().xaxis.set_major_locator(MultipleLocator(1))

plt.xlabel("Bar", fontsize=22)
plt.ylabel("Time resolution difference (ns)", fontsize=22)
plt.title(r"$(\sigma_t^{on} - \sigma_t^{off})$ vs Bar", fontsize=20)
#plt.title("Time resolution difference (on-off) vs Bar", fontsize=20)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()

# Save the plot
plt.savefig(".../delta_sigma_vs_bar.png", dpi=300)

#------------------------------------
# Peak energy ratio: I want it to be one (if so we chose the right overvoltage and breakdownvoltage)

peak_off = df_off[df_off["chip"] == 1][["bar", "num_events", "peak", "peak_err", "sigma", "sigma_err"]]
peak_on = df_on[df_on["chip"] == 1][["bar", "num_events", "peak", "peak_err", "sigma", "sigma_err"]]

# uncertainty due to temperature variation
volt_var = 37.5   # mV/°C
OV = [3000, 2741] # mV
dOV = [e * volt_var  for e in err_temp]
dov = [x / y for x, y in zip(dOV, OV)]

# ecceptance band
dR = math.sqrt(dov[0]**2 + dov[1]**2)
print("Acceptance band:", dR)

peak_off["peak_err"] = peak_off["peak_err"] + (dov[0] * peak_off["peak"])
peak_on["peak_err"] = peak_on["peak_err"] + (dov[1] * peak_on["peak"])

merged = pd.merge(peak_off, peak_on, on="bar", suffixes=("_off", "_on"))

# Compute ratio
merged["ratio_peak"] = merged["peak_on"] / merged["peak_off"]

# Propagate the error
merged["ratio_peak_err"] = merged["ratio_peak"] * np.sqrt(
    (merged["peak_err_on"] / merged["peak_on"])**2 +
    (merged["peak_err_off"] / merged["peak_off"])**2
)

# Final result
result = merged[["bar", "ratio_peak", "ratio_peak_err"]]

print(result)

# Extract x, y, and error
x = result["bar"]
y = result["ratio_peak"]
yerr = result["ratio_peak_err"]

# Perform the fit with error bars as weights
popt, pcov = curve_fit(constant, x, y, sigma=yerr, absolute_sigma=True)
c_fit = popt[0]
c_err = np.sqrt(pcov[0][0])

# Plot data with error bars
plt.figure(figsize=(10, 6))
plt.errorbar(x, y, yerr=yerr, fmt='o', capsize=4, label="Data")

# Plot fitted constant line
plt.axhline(c_fit, color='red', linestyle='--', label=f"Fit: c = {c_fit:.4f} ± {c_err:.4f}")

# Acceptance band around the fit
"""xbar = [-0.5, 15.5]
plt.fill_between(
    xbar,
    c_fit - dR * c_fit,
    c_fit + dR * c_fit,
    color='red',
    alpha=0.2,
    label=r"Acceptance band (temperature)"
)"""

# Customize the plot
"""plt.xlim(-0.5, 15.5)
plt.ylim(0.985, 1.005)"""
plt.gca().xaxis.set_major_locator(MultipleLocator(1))

plt.xlabel("Bar", fontsize=22)
plt.ylabel("Energy ratio", fontsize=22)
plt.title(r"$(\mu_E^{on}/\mu_E^{off})$ vs Bar", fontsize=20)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()



# Save the plot
plt.savefig("/home/cmsdaq/DAQ/mtd_daq/data/outputs/plots/tofhir/plots_run_10586/ratio_peak_vs_bar.png", dpi=300)
