import matplotlib.pyplot as plt
import pandas as pd

csv_path = 'figures_data/figure8_ablation_case3_power.csv'

df = pd.read_csv(csv_path)
df.plot.bar(rot=0, title="Peak Theta Power by Ablation Experiment", ylabel='[V^2/Hz]', xlabel='Experiment')
plt.show()
