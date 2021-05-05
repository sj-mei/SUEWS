import numpy as np
import pandas as pd
import supy as sp
import matplotlib.pyplot as plt
from scipy import stats

plt.close("all")

"""
# old SUEWS output files.
# suews output
df1_suews = pd.read_csv("C:/Users/nx902220/Documents/GitHub_SUEWS/SUEWS/Test/BaseRun/2020b/Output/test1_2004_SUEWS_60.txt",header=0,delimiter=r"\s+")
df2_suews = pd.read_csv("C:/Users/nx902220/Documents/GitHub_SUEWS/SUEWS/Test/BaseRun/2020b/Output/test1_2005_SUEWS_60.txt",header=0,delimiter=r"\s+")
df_suews = df1_suews.append(df2_suews)
print("df_suews:\n",df_suews)

# debug output
df1_debug = pd.read_csv("C:/Users/nx902220/Documents/GitHub_SUEWS/SUEWS/Test/BaseRun/2020b/Output/test1_2004_debug_60.txt",header=0,delimiter=r"\s+")
df2_debug = pd.read_csv("C:/Users/nx902220/Documents/GitHub_SUEWS/SUEWS/Test/BaseRun/2020b/Output/test1_2005_debug_60.txt",header=0,delimiter=r"\s+")
df_debug = df1_debug.append(df2_debug)
print("df_debug:\n",df_debug)
"""

### read in

# suews output
df_suews = pd.read_csv("C:/Users/nx902220/Documents/GitHub_SUEWS/SUEWS/Test/BaseRun/2020b/Output/Kc1_2011_SUEWS_60.txt",header=0,delimiter=r"\s+")
print("df_suews:\n",df_suews)

# debug output
df_debug = pd.read_csv("C:/Users/nx902220/Documents/GitHub_SUEWS/SUEWS/Test/BaseRun/2020b/Output/Kc1_2011_debug_60.txt",header=0,delimiter=r"\s+")
print("df_debug:\n",df_debug)

# KCL observations
df_obs = pd.read_csv("C:/Users/nx902220/Documents/GitHub_SUEWS/SUEWS/python_analysis/Obs_London_KCL_1h.txt",header=0,delimiter=r"\s+")
print("df_obs:\n",df_obs)


def compose_date(years, months=1, days=1, weeks=None, hours=None, minutes=None,
                 seconds=None, milliseconds=None, microseconds=None, nanoseconds=None):
    years = np.asarray(years) - 1970
    months = np.asarray(months) - 1
    days = np.asarray(days) - 1
    types = ('<M8[Y]', '<m8[M]', '<m8[D]', '<m8[W]', '<m8[h]',
             '<m8[m]', '<m8[s]', '<m8[ms]', '<m8[us]', '<m8[ns]')
    vals = (years, months, days, weeks, hours, minutes, seconds,
            milliseconds, microseconds, nanoseconds)
    return sum(np.asarray(v, dtype=t) for t, v in zip(types, vals)
               if v is not None)

### create datetime arrays 

# suews
df_suews['datetime'] = compose_date(df_suews['Year'], days=df_suews['DOY'], hours=df_suews['Hour'], minutes=df_suews['Min'])
print("df_suews:\n",df_suews.loc[:,'datetime'])
# debug
df_debug['datetime'] = compose_date(df_debug['Year'], days=df_debug['DOY'], hours=df_debug['Hour'], minutes=df_debug['Min'])
print("df_debug:\n",df_debug.loc[:,'datetime'])
# obs
df_obs['datetime'] = compose_date(df_obs['Year'], days=df_obs['DOY'], hours=df_obs['Hour'], minutes=df_obs['Min'])
print("df_obs:\n",df_obs.loc[:,'datetime'])

### time slice

start = '2011-01-01'
end = '2011-12-01'

time_slice = np.logical_and(df_debug['datetime'] >= start, df_debug['datetime'] <= end)
df_debug = df_debug[time_slice]
df_suews = df_suews[time_slice]
time_slice = np.logical_and(df_obs['datetime'] >= start, df_obs['datetime'] <= end)
df_obs = df_obs[time_slice]

### plot

linestyle="-"
marker=""

linestyle_obs = "--"
marker_obs = ""

## radiation timeseries
plt.figure(figsize=(12,8))
plt.plot(df_suews['datetime'],df_suews['QN'],linestyle=linestyle,marker=marker,color="black",label=r'$Q_N$')
plt.plot(df_obs['datetime'],df_obs['QN'],linestyle=linestyle_obs,marker=marker_obs,color="black")
plt.plot(df_suews['datetime'],df_suews['Lup'],linestyle=linestyle,marker=marker,color="blue",label=r'$L\uparrow$')
plt.plot(df_obs['datetime'],df_obs['Lup'],linestyle=linestyle_obs,marker=marker_obs,color="blue")
plt.plot(df_suews['datetime'],df_suews['Ldown'],linestyle=linestyle,marker=marker,color="red",label=r'$L\downarrow$')
plt.plot(df_obs['datetime'],df_obs['Ldown'],linestyle=linestyle_obs,marker=marker_obs,color="red")
plt.plot(df_suews['datetime'],df_suews['Kup'],linestyle=linestyle,marker=marker,color="purple",label=r'$K\uparrow$')
plt.plot(df_obs['datetime'],df_obs['Kup'],linestyle=linestyle_obs,marker=marker_obs,color="purple")
plt.plot(df_suews['datetime'],df_suews['Kdown'],linestyle=linestyle,marker=marker,color="orange",label=r'$K\downarrow$')
plt.plot(df_obs['datetime'],df_obs['Kdown'],linestyle=linestyle_obs,marker=marker_obs,color="orange")
plt.xlabel("Time")
plt.ylabel(r"Energy Flux Density (Wm$^{-2}$)")
plt.xticks(rotation='vertical')
#plt.ylim(-100,900)
plt.legend(loc=2)
plt.tight_layout()

## temperature timeseries
plt.figure(figsize=(12,8))
plt.plot(df_suews['datetime'],df_suews['T2'],linestyle=linestyle,marker=marker,label=r'$T_2$')
plt.plot(df_suews['datetime'],df_suews['Ts'],linestyle=linestyle,marker=marker,label=r'$T_s$')
plt.plot(df_suews['datetime'],df_suews['Tsurf'],linestyle=linestyle,marker=marker,label=r'$T_{surf}$')
plt.plot(df_obs['datetime'],df_obs['Tair'],linestyle=linestyle,marker=marker,label=r'Observed $T_{air}$')
plt.xlabel("Time")
plt.ylabel(r"T ($^{\circ}C$)")
plt.xticks(rotation='vertical')
#plt.ylim(-5,9)
plt.legend(loc=2)
plt.tight_layout()

## albedo and emissivity timeseries
plt.figure(figsize=(12,8))
plt.plot(df_debug['datetime'],df_debug['alb_spc'],linestyle=linestyle,marker=marker,label=r'$\alpha$')
plt.plot(df_obs['datetime'],df_obs['Kup']/df_obs['Kdown'],linestyle=linestyle,marker=marker,label=r'Observed $\alpha$')
plt.plot(df_debug['datetime'],1-df_debug['emiss_spc'],linestyle=linestyle,marker=marker,label=r'$1-\epsilon$')
plt.xlabel("Time")
plt.ylabel(r"Albedo and 1-Emissivity")
plt.xticks(rotation='vertical')
plt.ylim(0,0.3)
plt.legend(loc=2)
plt.tight_layout()

## Create an Zenith column in obs interpolated from suews
df_obs['Zenith'] = np.interp(df_obs["datetime"],df_suews["datetime"],df_suews["Zenith"])

## Zenith timeseries
plt.figure(figsize=(8,6))
plt.plot(df_suews['datetime'],df_suews['Zenith'],linestyle='solid',marker='',label=r'SUEWS')
plt.xlabel("Time")
plt.ylabel(r"Zenith")
plt.legend(loc=2)
plt.tight_layout()

## observed and suews panel albedo Kup:Kdown
fig, (ax1, ax2) = plt.subplots(1, 2)
# obs
# calculate linear regression: observed
df_obs_daytime = df_obs[(df_obs['Zenith'] <= 80)]
slope, intercept, r_value, p_value, std_err = stats.linregress(df_obs_daytime['Kdown'].dropna(),df_obs_daytime['Kup'].dropna())
line_x = np.arange(df_obs_daytime['Kdown'].min(), df_obs_daytime['Kdown'].max())
line_y = slope*line_x + intercept
# plot
ax1.plot(df_obs_daytime['Kdown'],df_obs_daytime['Kup'],linestyle="none",marker="o")
ax1.plot(line_x, line_y,label='$%.3fx + %.3f$, $R^2=%.3f$' % (slope, intercept, r_value**2))
ax1.set_xlabel(r'$K\downarrow$ (Wm$^{-2}$)')
ax1.set_ylabel(r'$K\uparrow$ (Wm$^{-2}$)')
ax1.set_xlim(0,)
ax1.set_ylim(0,)
ax1.legend(loc=2)
ax1.set_title('Observed')
# spartacus
# calculate linear regression: spartacus
df_suews_daytime = df_suews[(df_suews['Zenith'] <= 80)]
slope, intercept, r_value, p_value, std_err = stats.linregress(df_suews_daytime['Kdown'].dropna(),df_suews_daytime['Kup'].dropna())
line_x = np.arange(df_suews_daytime['Kdown'].min(), df_suews_daytime['Kdown'].max())
line_y = slope*line_x + intercept
ax2.plot(df_suews_daytime['Kdown'],df_suews_daytime['Kup'],linestyle="none",marker="o")
ax2.plot(line_x, line_y,label='$%.3fx + %.3f$, $R^2=%.3f$' % (slope, intercept, r_value**2))
ax2.set_xlabel(r'$K\downarrow$ (Wm$^{-2}$)')
ax2.set_ylabel(r'$K\uparrow$ (Wm$^{-2}$)')
ax2.set_xlim(0,)
ax2.set_ylim(0,)
ax2.legend(loc=2)
ax2.set_title('Spartacus')
# overall formatting
#fig.suptitle('Horizontally stacked subplots')
fig.set_figheight(4)
fig.set_figwidth(8)
plt.tight_layout()
