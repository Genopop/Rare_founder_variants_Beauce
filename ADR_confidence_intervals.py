import pandas as pd
import numpy as np
import os
from scipy.stats import t

# === Load ancestor information ===
print("Get ancestor info")
sum_ancetres = pd.read_csv(
    os.path.join(data_dir, "dates_locations.txt"), sep=" "
)
ancestors_infos = pd.DataFrame({
    'ind': sum_ancetres['IndID'],
    'region': sum_ancetres['RegionIDMariage'],
    'date': sum_ancetres['DateMariage']
})

# === Load regional ancestor files (Outputs from ancestors_diversity_ratio.py) ===
print("Get region files")
total_beauce = pd.read_csv(os.path.join(data_dir, "Beauce_ADR.csv"))
total_beauce['pro'] = 'Beauce'
total_sag = pd.read_csv(os.path.join(data_dir, "SLSJ_ADR.csv"))
total_sag['pro'] = 'SLSJ'
total_montreal = pd.read_csv(os.path.join(data_dir, "Mtl_ADR.csv"))
total_montreal['pro'] = 'Montreal'

# Combine all regions
total_df = pd.concat([total_beauce, total_montreal, total_sag], ignore_index=True)
total_df = total_df[total_df['ancestor_id'] != 0]

# Merge with ancestor metadata
total_infos = pd.merge(total_df, ancestors_infos, left_on='ancestor_id', right_on='ind')
total_infos = total_infos[(total_infos['date'] >= 1630) & (total_infos['date'] <= 1950)]

# === Process decades and filter iterations ===
print("Group by decades")
total_infos['decade'] = (np.ceil(total_infos['date'] / 10) * 10).astype(int)
iteration_cols = [col for col in total_infos.columns if col.startswith('iteration_')]
total_infos['pro'] = total_infos['pro'].astype('category')
total_infos = total_infos[total_infos[iteration_cols].sum(axis=1) > 0]

# === Reshape data to long format for ancestor counts ===
print("Get unique ancestors")
long_df = total_infos.melt(
    id_vars=['ancestor_id', 'decade', 'pro'],
    value_vars=iteration_cols,
    var_name='iteration',
    value_name='ancestor_value'
)
long_df = long_df[long_df['ancestor_value'] != 0]
# Count unique ancestors per iteration
unique = (
    long_df
    .groupby(['decade', 'pro', 'iteration'], as_index=False)
    .agg(n_unique=('ancestor_id', 'nunique'))
)

# === Compute total number of ancestor occurrences per iteration ===
print("Get total number of apparitions")
df_year_grouped = (
    total_infos
    .groupby(['decade', 'pro'], as_index=False)
    .agg({col: 'sum' for col in iteration_cols})
)
totals_long = df_year_grouped.melt(
    id_vars=['decade', 'pro'],
    value_vars=iteration_cols,
    var_name='iteration',
    value_name='total_count'
)

# === Compute ADR (ancestor diversity ratio) ===
print("Calculate ADR")
joined = pd.merge(unique, totals_long, on=['decade', 'pro', 'iteration'], how='left')
joined['ratio'] = joined['n_unique'] / joined['total_count']

# Pivot back to wide format to calculate statistics
ratios_wide = joined.pivot_table(
    index=['decade', 'pro'],
    columns='iteration',
    values='ratio'
).reset_index()

# === Compute confidence intervals for ADR per decade/probands' region ===
print("Get confidence intervals")
def compute_stats(row):
    values = row.filter(like='iteration_').dropna()
    n = len(values)
    mean = values.mean()
    sd = values.std(ddof=1)
    if n > 1:
        ci = t.interval(0.95, df=n-1, loc=mean, scale=sd / np.sqrt(n))
    else:
        ci = (np.nan, np.nan)
    return pd.Series({'mean': mean, 'ci_lower': ci[0], 'ci_upper': ci[1]})

stats_df = ratios_wide.copy()
stats = stats_df.apply(compute_stats, axis=1)
df_stats = pd.concat([stats_df[['decade', 'pro']], stats], axis=1)
df_stats = df_stats[(df_stats['decade'] >= 1640) & (df_stats['decade'] <= 1940)]

# === Export final result ===
output_path = os.path.join(data_dir, "ADR_beauce-slsj-mtl.csv")
df_stats.to_csv(output_path, index=False)
print("Output saved")
