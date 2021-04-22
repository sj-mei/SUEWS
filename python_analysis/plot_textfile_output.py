import numpy as np
import pandas as pd
import supy as sp
import matplotlib.pyplot as plt

plt.close("all")

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

# suews
df_suews['datetime'] = compose_date(df_suews['Year'], days=df_suews['DOY'], hours=df_suews['Hour'], minutes=df_suews['Min'])
print("df_suews:\n",df_suews.loc[:,'datetime'])
# debug
df_debug['datetime'] = compose_date(df_debug['Year'], days=df_debug['DOY'], hours=df_debug['Hour'], minutes=df_debug['Min'])
print("df_debug:\n",df_debug.loc[:,'datetime'])

# radiation
plt.figure(figsize=(12,8))
plt.plot(df_suews['datetime'], df_suews['QN'],label=r'$Q_N$')
plt.plot(df_suews['datetime'], df_suews['Lup'],label=r'$L\uparrow$')
plt.plot(df_suews['datetime'], df_suews['Ldown'],label=r'$L\downarrow$')
plt.plot(df_suews['datetime'], df_suews['Kup'],label=r'$K\uparrow$')
plt.plot(df_suews['datetime'], df_suews['Kdown'],label=r'$K\downarrow$')
plt.xlabel("Time")
plt.ylabel(r"Energy Flux Density (Wm$^{-2}$)")
plt.xticks(rotation='vertical')
plt.ylim(-80,380)
plt.legend(loc=2)
plt.tight_layout()

# temperature
plt.figure(figsize=(12,8))
plt.plot(df_suews['datetime'], df_suews['T2'],label=r'$T_2$')
plt.plot(df_suews['datetime'], df_suews['Ts'],label=r'$T_s$')
plt.plot(df_suews['datetime'], df_suews['Tsurf'],label=r'$T_{surf}$')
plt.xlabel("Time")
plt.ylabel(r"T ($^{\circ}C$)")
plt.xticks(rotation='vertical')
plt.ylim(-5,9)
plt.legend(loc=2)
plt.tight_layout()

# albedo and emissivity
plt.figure(figsize=(12,8))
plt.plot(df_debug['datetime'], df_debug['alb_spc'],label=r'$\alpha$')
plt.plot(df_debug['datetime'], 1-df_debug['emiss_spc'],label=r'$1-\epsilon$')
plt.xlabel("Time")
plt.ylabel(r"Albedo and 1-Emissivity")
plt.xticks(rotation='vertical')
plt.ylim(0,)
plt.legend(loc=2)
plt.tight_layout()