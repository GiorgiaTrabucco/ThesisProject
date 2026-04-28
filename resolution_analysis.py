# written by Giorgia Trabucco
#--------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import ROOT as R
import math

#--------------------------------------------------

# function that plots the bar map of the bad energy resolution bars
def plot_bar_map(name, title, data, mu, sigma, sigma_type, output_path, ymin=None, ymax=None):
    [...]
    
#--------------------------------------------------

# this script takes data of time resolution, it computes a mu and a sigma of the distribution
# the point is to set boundaries on the accepted time resolution
# I want to use it for the study of bad energy resolution arrays

#--------------------------------------------------

# Save data in a list of dataframes 
chips = [0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
runs = [x for x in range(10540, 10565) if x not in (10553,)]

dfs = []
for run, chip in zip(runs, chips):
    file = f".../peaks.pq"
    df = pd.read_parquet(file)
    dfs.append(df[df["chip"] == chip])

#------------------------------------------

# histograms of all the bars: relative energy resolution
color0 = R.kOrange
canvas0 = R.TCanvas("canvas0", "histogram", 1300, 1000)
canvas0.SetBottomMargin(0.15)
canvas0.SetLeftMargin(0.15)

canvas0.cd()

bin_edges0 = np.array([1.52+i*0.1 for i in range(24)])
hist0 = R.TH1F("hist0", "", len(bin_edges0)-1, bin_edges0)

min_hist = 10
max_hist = 0
index = 0
for df in dfs:
    df['en_res'] = 100. * df['sigma'].values / df['peak'].values
    df['en_res_err'] = df['en_res'] * np.sqrt((df["sigma_err"] / df["sigma"])**2 + (df["peak_err"] / df["peak"])**2)

    values = df['en_res']
    index2 = 0
    for value in values:
        hist0.Fill(value)
        if value < min_hist: min_hist = value
        if value > max_hist: max_hist = value

        if value < 2: print(f"run{runs[index]}_chip{chips[index]}_bar{index2}: {value}")
        elif value > 3.5: print(f"BAD: run{runs[index]}_chip{chips[index]}_bar{index2}")
        index2 = index2 + 1
    index = index + 1
print(f"\n energy res hist: min={min_hist}, max={max_hist} \n")

#hist0.SetStats(0)
hist0.SetFillStyle(3001)
hist0.SetLineColor(color0)
hist0.SetFillColor(color0)
hist0.GetXaxis().SetTitle("Energy Resolution (%)")
hist0.GetYaxis().SetTitle("Number of crystal bars")
hist0.GetXaxis().SetTitleSize(0.05)
hist0.GetYaxis().SetTitleSize(0.05)
hist0.GetXaxis().SetNdivisions(5)
hist0.Draw()

mu0 = hist0.GetMean()
sigma0 = hist0.GetStdDev()

canvas0.SaveAs(f".../energy_res.png")

#------------------------------------------

# data outside mu +- 3 sigma: compute again mu and sigma
df_out0 = pd.DataFrame() 
df_filtered0 = pd.DataFrame() 

for df in dfs:
    mask = (100. * df['sigma'] / df['peak'] < mu0 - 3 * sigma0) | (100. * df['sigma'] / df['peak'] > mu0 + 3 * sigma0)
    df_out0 = pd.concat([df_out0, df[mask]], ignore_index=True)
    df_filtered0 = pd.concat([df_filtered0, df[~mask]], ignore_index=True)

# redo the histograms
bin_edges0 = np.array([2.+i*0.1 for i in range(20)])
hist01 = R.TH1F("hist01", "Energy Resolution", len(bin_edges0)-1, bin_edges0)
canvas0.cd()

for value in 100. * df_filtered0["sigma"].values / df_filtered0["peak"].values:
        hist01.Fill(value)

hist01.Draw("hist")
hist01.SetFillStyle(3001)
hist01.SetLineColor(color0)
hist01.SetFillColor(color0)
hist01.GetXaxis().SetTitle("Energy Resolution (%)")
hist01.GetYaxis().SetTitle("Number of crystal bars")
hist01.GetXaxis().SetTitleSize(0.05)
hist01.GetYaxis().SetTitleSize(0.05)

"""mu0 = hist01.GetMean()
sigma0 = hist01.GetStdDev()"""

canvas0.SaveAs(f".../energy_res_cleaned.png")

#------------------------------------------

# histograms of all the bars: time resolution
color = R.kViolet
#canvas = R.TCanvas("canvas", "histogram", 3000, 1000)
canvas = R.TCanvas("canvas", "histogram", 1300, 1000)
#canvas.Divide(2, 1)
canvas.SetLeftMargin(0.15)

# sigma_corr1
#canvas.cd(1)
canvas.cd()

bin_edges = np.array([0.052+i*0.0015 for i in range(23)])
hist1 = R.TH1F("hist1", "", len(bin_edges)-1, bin_edges)
min_hist = 1
max_hist = 0
for i in range(len(dfs)):
    values = dfs[i]['sigma_dt3'].values

    for value in values:
        hist1.Fill(value)
        if value > max_hist: max_hist = value
        if value < min_hist: min_hist = value

#print(f"\n time res hist: min={min_hist}, max={max_hist} \n")

#hist1.SetStats(0)
hist1.SetFillStyle(3001)
hist1.SetLineColor(color)
hist1.SetFillColor(color)
hist1.GetXaxis().SetTitle("Time Resolution (ns)")
hist1.GetYaxis().SetTitle("Number of crystal bars")
hist1.GetXaxis().SetTitleSize(0.05)
hist1.GetYaxis().SetTitleSize(0.05)
hist1.GetXaxis().SetNdivisions(5)
hist1.Draw()
#hist1.Draw("same")

mu1 = hist1.GetMean()
sigma1 = hist1.GetStdDev()

# sigma_corr2
"""canvas.cd(2)

hist2 = R.TH1D("hist2", "Time Resolution phase corrected", len(bin_edges)-1, bin_edges)

for i in range(len(dfs)):
    values = dfs[i]['sigma_dt4'].values

    for value in values:
        #print(value)
        hist2.Fill(value)
        #hist2.SetStats(0)
        hist2.SetFillStyle(3001)
        hist2.SetLineColor(color)
        hist2.SetFillColor(color)

hist2.Draw()
hist2.GetXaxis().SetTitle("Time Resolution (ns)")
hist2.GetYaxis().SetTitle("Number of crystal bars")
hist2.GetXaxis().SetTitleSize(0.05)
hist2.GetYaxis().SetTitleSize(0.05)

mu2 = hist2.GetMean()
sigma2 = hist2.GetStdDev()"""

canvas.SaveAs(f".../time_res.png")

#------------------------------------------

# data outside mu +- 3 sigma: compute again mu and sigma
df_out1 = pd.DataFrame() 
#df_out2 = pd.DataFrame()
df_filtered1 = pd.DataFrame() 
#df_filtered2 = pd.DataFrame() 

for df in dfs:
    values1 = df['sigma_dt3'].values
    #values2 = df['sigma_dt4'].values

    # Build masks for outliers
    mask1 = (values1 < mu1 - 3 * sigma1) | (values1 > mu1 + 3 * sigma1)
    #mask2 = (values2 < mu2 - 3 * sigma2) | (values2 > mu2 + 3 * sigma2)

    # Append rows where condition is True
    df_out1 = pd.concat([df_out1, df[mask1]], ignore_index=True)
    #df_out2 = pd.concat([df_out2, df[mask2]], ignore_index=True)

    df_filtered1 = pd.concat([df_filtered1, df[~mask1]], ignore_index=True)
    #df_filtered2 = pd.concat([df_filtered2, df[~mask2]], ignore_index=True)
    
#print(df_out1)
#print(df_out2)

# redo the histograms
bin_edges = np.array([0.05+i*0.002 for i in range(11)])
hist3 = R.TH1F("hist3", "Time Resolution energy corrected", len(bin_edges)-1, bin_edges)
canvas.cd(1)

for value in df_filtered1["sigma_dt3"].values:
        hist3.Fill(value)

hist3.Draw("hist")
hist3.SetFillStyle(3001)
hist3.SetLineColor(color)
hist3.SetFillColor(color)
hist3.GetXaxis().SetTitle("Time Resolution (ns)")
hist3.GetYaxis().SetTitle("Number of crystal bars")
hist3.GetXaxis().SetTitleSize(0.05)
hist3.GetYaxis().SetTitleSize(0.05)

"""mu1 = hist3.GetMean()
sigma1 = hist3.GetStdDev()
"""

"""hist4 = R.TH1F("hist4", "Time Resolution energy corrected", len(bin_edges)-1, bin_edges)

canvas.cd(2)

for value in df_filtered2["sigma_dt4"].values:
        hist4.Fill(value)

hist4.Draw("hist")
hist4.SetFillStyle(3001)
hist4.SetLineColor(color)
hist4.SetFillColor(color)
hist4.GetXaxis().SetTitle("Time Resolution (ns)")
hist4.GetYaxis().SetTitle("Number of crystal bars")
hist4.GetXaxis().SetTitleSize(0.05)
hist4.GetYaxis().SetTitleSize(0.05)

mu2 = hist4.GetMean()
sigma2 = hist4.GetStdDev()"""

canvas.SaveAs(f".../time_res_cleaned.png")

#------------------------------------------

# plot the map of the tray board (RU0, TRAY6)
df_out_en0 = pd.DataFrame()
df_out_t10 = pd.DataFrame()
df_tot0 = pd.DataFrame()
dm_bad0 = set()
dm_bad10 = set()

for df in dfs:
    df_tot0 = pd.concat([df_tot0, df], ignore_index=True)

    # --- en ---
    mumu = df['en_res']
    sigmama = mumu * df['en_res_err']
    df['num_sigma_en'] = (mumu - mu0) / np.sqrt(sigma0**2 + sigmama**2)
    mask_en = df['num_sigma_en'] > 3.
    selected_en0 = df[mask_en]
    df_out_en0 = pd.concat([df_out_en0, selected_en0], ignore_index=True)
    dm_bad0.update(selected_en0['DM'].astype(int).unique())

    # --- t1 ---
    df['num_sigma_dt'] = (df['sigma_dt3'] - mu1) / np.sqrt(sigma1**2 + df['sigma_err_dt3']**2)
    mask_t1 = df['num_sigma_dt'] > 3.
    selected_t1 = df[mask_t1]
    df_out_t10 = pd.concat([df_out_t10, selected_t1], ignore_index=True)
    dm_bad10.update(selected_t1['DM'].astype(int).unique())

output_path = f".../plots_run_{runs[-1]}"

plot_bar_map("bad_en_res_map_tray", "Bad Energy Resolution map", df_out_en0, mu0, sigma0, "sigma", output_path)
plot_bar_map("bad_time_res1_map_tray", "Bad Time Resolution map", df_out_t10, mu1, sigma1, "sigma_dt3", output_path)

"""plot_bar_map("bad_en_res_map_tray", "Bad Energy Resolution map", df_tot0, mu0, sigma0, "sigma", output_path)
plot_bar_map("bad_time_res1_map_tray", "Bad Time Resolution map", df_tot0, mu1, sigma1, "sigma_dt3", output_path)"""

#------------------------------------------

# study of bad energy resolution arrays
# board with 10 dms: 6 of them are bad (we want to identify them without knowing a priori)
# What I expect is to find more than 2 affected channels for a bad sm
#check the enenrgy resolution
#check the time resolution

#--------------------------------------------------

# Save data in a list of dataframes (bad board)
chips = [1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1]
runs = [x for x in range(10593, 10616) if x not in (10594,10597,10598)]

dfs2 = []
for run, chip in zip(runs, chips):
    file = f".../peaks.pq"
    df = pd.read_parquet(file)
    dfs2.append(df[df["chip"] == chip])
    #if run == 10615: print(dfs2[-1])

#------------------------------------------

# histograms of all the bars (bad board): relative energy resolution
color0 = R.kOrange
canvas0 = R.TCanvas("canvas5", "histogram", 1300, 1000)
canvas0.SetBottomMargin(0.15)
canvas0.SetLeftMargin(0.15)

# sigma_corr1
canvas0.cd()

bin_edges0 = np.array([1.9+i*0.25 for i in range(26)])
hist5 = R.TH1F("hist5", "", len(bin_edges0)-1, bin_edges0)

min_hist = 100
max_hist = 0
for df in dfs2:
    df['en_res'] = 100. * df['sigma'].values / df['peak'].values
    df['en_res_err'] = df['en_res'] * np.sqrt((df["sigma_err"] / df["sigma"])**2 + (df["peak_err"] / df["peak"])**2)
    values = df['en_res']

    for value in values:
        hist5.Fill(value)
        if value < 0: print(runs[i])
        if value < min_hist: min_hist = value 
        if value > max_hist: max_hist = value

print(f"\n BAD energy res hist: min={min_hist}, max= {max_hist} \n")

#hist5.SetStats(0)
hist5.SetFillStyle(3001)
hist5.SetLineColor(color0)
hist5.SetFillColor(color0)
hist5.GetXaxis().SetTitle("Energy Resolution (%)")
hist5.GetYaxis().SetTitle("Number of crystal bars")
hist5.GetXaxis().SetTitleSize(0.05)
hist5.GetYaxis().SetTitleSize(0.05)
hist5.GetXaxis().SetNdivisions(5)
hist5.Draw()

mu_bad0 = hist5.GetMean()
sigma_bad0 = hist5.GetStdDev()

canvas0.SaveAs(f".../energy_res_bad_board.png")

#------------------------------------------

# histograms of all the bars (bad board): time resolution
color = R.kViolet
canvas = R.TCanvas("canvas6", "histogram", 1300, 1000)
canvas0.SetBottomMargin(0.15)
canvas.SetLeftMargin(0.15)

# sigma_corr1
canvas.cd()

bin_edges = np.array([0.053+i*0.0012 for i in range(22)])
hist6 = R.TH1F("hist6", "", len(bin_edges)-1, bin_edges)

min_hist = 100
max_hist = 0
for i in range(len(dfs2)):
    values = dfs2[i]['sigma_dt3'].values

    for value in values:
        hist6.Fill(value)
        if value < 0: print(runs[i])
        if value < min_hist: min_hist = value 
        if value > max_hist: max_hist = value

print(f"\n BAD time res hist: min={min_hist}, max= {max_hist} \n")

#hist6.SetStats(0)
hist6.SetFillStyle(3001)
hist6.SetLineColor(color)
hist6.SetFillColor(color)   
hist6.GetXaxis().SetTitle("Time Resolution (ns)")
hist6.GetYaxis().SetTitle("Number of crystal bars")
hist6.GetXaxis().SetTitleSize(0.05)
hist6.GetYaxis().SetTitleSize(0.05)
hist6.GetXaxis().SetNdivisions(5)
hist6.Draw()

mu_bad1 = hist6.GetMean()
sigma_bad1 = hist6.GetStdDev()

canvas.SaveAs(f".../time_res_bad_board.png")

#------------------------------------------

# data bad energy res
df_out_en = pd.DataFrame()
df_out_t1 = pd.DataFrame()
df_out_t2 = pd.DataFrame()
df_tot = pd.DataFrame()

dm_bad = set()
dm_bad1 = set()
dm_bad2 = set()
for df in dfs2:
    df_tot = pd.concat([df_tot, df], ignore_index=True)
    # --- en ---
    mumu = df['en_res']
    sigmama = mumu * df['en_res_err']
    df['num_sigma_en'] = (mumu - mu0) / np.sqrt(sigma0**2 + sigmama**2)
    mask_en = df['num_sigma_en'] > 3.
    selected_en = df[mask_en]
    df_out_en = pd.concat([df_out_en, selected_en], ignore_index=True)
    dm_bad.update(selected_en['DM'].astype(int).unique())

    # --- t1 ---
    df['num_sigma_dt'] = (df['sigma_dt3'] - mu1) / np.sqrt(sigma1**2 + df['sigma_err_dt3']**2)
    mask_t1 = df['num_sigma_dt'] > 3.
    selected_t1 = df[mask_t1]
    df_out_t1 = pd.concat([df_out_t1, selected_t1], ignore_index=True)
    dm_bad1.update(selected_t1['DM'].astype(int).unique())

    # --- t2 ---
    #mask_t2 = df['sigma_dt4'] > mu2 + 3 * sigma2
    """mask_t2 = (df['sigma_dt4'] - mu2) / np.sqrt(sigma2**2 + df['sigma_err_dt4']**2) > 3.
    selected_t2 = df[mask_t2]
    df_out_t2 = pd.concat([df_out_t2, selected_t2], ignore_index=True)
    dm_bad2.update(selected_t2['DM'].astype(int).unique())"""


print("Bad energy resolution")
counts = df_out_en.groupby(['DM', 'chip']).size().reset_index(name='#channels')
#counts = counts[counts['#channels'] > 2]
print(counts)

"""print("Bad time resolution (energy corrected) channels")
print(df_out_t1)
print("Bad time resolution (phase corrected) channels")
print(df_out_t2)"""

print("Bad time resolution (energy corrected)")
counts1 = df_out_t1.groupby(['DM', 'chip']).size().reset_index(name='#channels')
#counts1 = counts1[counts1['#channels'] > 2]
print(counts1)

"""print("Bad time resolution (phase corrected)")
counts2 = df_out_t2.groupby(['DM', 'chip']).size().reset_index(name='#channels')
#counts2 = counts2[counts2['#channels'] > 2]
print(counts2)"""

"""# Find DM values that have both chip 0 and chip 1
valid_dms = counts2.groupby('DM')['chip'].nunique()
valid_dms = valid_dms[valid_dms > 1].index
counts2 = counts2[counts2['DM'].isin(valid_dms)]
print(counts2)"""

dm_bad = list(map(int, dm_bad))
dm_bad1 = list(map(int, dm_bad1))
#dm_bad2 = list(map(int, dm_bad2))
print(f"Bad energy resolution DMs: {dm_bad}")
print(f"Bad time resolution (energy corrected) DMs: {dm_bad1}")
#print(f"Bad time resolution (phase corrected) DMs: {dm_bad2}")

#--------------------------------------------------
# plot the bar map of the bad energy resolution bars
output_path = f".../plots_run_{runs[-1]}"

# maps: colored bars are bad ones
plot_bar_map("bad_en_res_map", "Bad Energy Resolution map", df_out_en, mu0, sigma0, "sigma", output_path)
plot_bar_map("bad_time_res1_map", "Bad Time Resolution map", df_out_t1, mu1, sigma1, "sigma_dt3", output_path)

"""plot_bar_map("bad_en_res_map", "Bad Energy Resolution map", df_tot, mu0, sigma0, "sigma", output_path)
plot_bar_map("bad_time_res1_map", "Bad Time Resolution map", df_tot, mu1, sigma1, "sigma_dt3", output_path)"""

