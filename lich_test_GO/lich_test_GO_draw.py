import pandas as pd, matplotlib.pyplot as plt

import matplotlib as mpl

mpl.rcParams.update({
    "figure.dpi": 300,
    "figure.figsize": (3.2, 2.4),    
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.linewidth": 0.8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.minor.size": 1.5,
    "ytick.minor.size": 1.5,
    "axes.spines.top": True,
    "axes.spines.right": True,
    "axes.spines.left": True,
    "axes.spines.bottom": True,
    "legend.frameon": False,
})


df = pd.read_csv("./frameCount.csv")
#frame,ncubic,nhex,nmix,nci,nhi,nch,nint,nliq
for col in ["ncubic","nhex","nmix","nci","nhi","nch","nint","nliq"]:
    plt.plot(df["frame"], df[col], label=col)
plt.legend()
plt.xlabel("ps")
plt.ylabel("count")
plt.show()
plt.savefig("frameCount.png")

