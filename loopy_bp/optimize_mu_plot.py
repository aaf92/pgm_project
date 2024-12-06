import pandas as pd
import matplotlib.pyplot as plt

# Read the data from the CSV file
df = pd.read_csv("/ix/djishnu/Aaron_F/PGM_project/loopy_bp/mu_score.csv")

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(df['mu'], df['score'], marker='o', label='Score vs. Mu', color='b')

# Add titles and labels
plt.title("Score vs Mu", fontsize=14)
plt.xlabel("Mu Parameter", fontsize=12)
plt.ylabel("Score", fontsize=12)
plt.tight_layout()

# save
plt.savefig("/ix/djishnu/Aaron_F/PGM_project/loopy_bp/mu_score_plot.png")
