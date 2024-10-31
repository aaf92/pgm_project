import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

# Edge data
data = pd.read_csv("/ix/djishnu/Aaron_F/PGM_project/edge_data.csv")
edge_values = data["posterior_edge"].tolist()

# Plot Distribution
plt.figure(figsize=(10, 6))
sns.histplot(edge_values, kde=True, bins=30, color='blue')  # kde=True adds a density curve
plt.title('Posterior Edge Values', fontsize=16)
plt.xlabel('edge_probability', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.grid(visible=True, linestyle='--', alpha=0.5)
plt.savefig('/ix/djishnu/Aaron_F/PGM_project/20241030_Resutls/edge_distribution_il6.pdf',format='pdf')
