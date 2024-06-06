import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.interpolate import interp1d
from project_filters import project_filters


moh_label = ', '.join(project_filters['Serotype'])
skygrid_label = ', '.join(project_filters['Serotype']) + ' (Skygrid)'

# Read the first CSV file into a DataFrame
moh_data = pd.read_csv('moh_months.csv')

# Ensure that the 'Year' and 'Month' columns are available for processing
if 'Year' in moh_data.columns and 'Month' in moh_data.columns:
    # Create a new DataFrame with the desired columns and calculations
    # Prepare the list of column names based on the serotypes
    serotype_columns = [f'{s}.cases' for s in project_filters['Serotype']]

    # Calculate the sum of the selected serotype columns
    moh_data['Total.Cases'] = moh_data[serotype_columns].sum(axis=1)

    filtered_moh_data = pd.DataFrame({
        'Date': (moh_data['Year'] + moh_data['Month'] / 12),
        'Type': ', '.join(project_filters['Serotype']),
        'log.Pop': moh_data['Total.Cases']
    })

else:
    print("Error: 'Year' and/or 'Month' columns are missing from the dataset.")
    
# # Read the second tab-separated file into a DataFrame
skygrid_data = pd.read_csv('skygrid_data.txt', sep='\t')
skygrid_data = skygrid_data[['time', 'mean']]
skygrid_data = pd.DataFrame({
    'Date': skygrid_data['time'],
    'Type': ', '.join(project_filters['Serotype']) + ' (Skygrid)',
    'log.Pop': skygrid_data['mean']
})

# Merge the information from the second DataFrame into the first one
appended_df = pd.concat([filtered_moh_data, skygrid_data], ignore_index=True)

# Find the largest Date value for D3_SkyGrid
max_date_d3_skygrid = appended_df[appended_df['Type'] == skygrid_label]['Date'].max()

# Create a mask to filter D3 rows
mask_d3 = (appended_df['Type'] == moh_label) & (appended_df['Date'] <= max_date_d3_skygrid)

# Filter D3 rows based on the mask
appended_df = appended_df[mask_d3 | (appended_df['Type'] == skygrid_label)]

# Specify x-values for interpolation based on 'D3' type
specified_x_values = appended_df[appended_df['Type'] == moh_label]['Date'].values

# Extract 'D3' values for specified x-values
d3_values = appended_df[appended_df['Type'] == moh_label].set_index('Date')['log.Pop'].reindex(specified_x_values)

# Extract 'D3_SkyGrid' values for existing x and y values
existing_x_values = appended_df[appended_df['Type'] == skygrid_label]['Date'].values
existing_y_values = appended_df[appended_df['Type'] == skygrid_label]['log.Pop'].values

# Perform linear interpolation for 'D3_SkyGrid' values
interpolator = interp1d(existing_x_values, existing_y_values, kind='linear', fill_value="extrapolate")
interpolated_values = interpolator(specified_x_values)

# Calculate the Pearson correlation coefficient and p-value
correlation_coefficient, p_value = pearsonr(interpolated_values, d3_values.dropna())

# Print the results
print(f"Pearson Correlation Coefficient: {correlation_coefficient}")
print(f"P-value: {p_value}")

# Set Seaborn theme to 'viridis'
sns.set_theme(style="whitegrid", palette="viridis")

# Create a figure and axis
fig, ax1 = plt.subplots()

# Plot the 'D3' values on the left axis
sns.lineplot(data=appended_df[appended_df['Type'] == moh_label], x='Date', 
             y='log.Pop', ax=ax1, label=moh_label, legend=False)

# Create a twin axis on the right
ax2 = ax1.twinx()

# Plot the 'D3_SkyGrid' values on the right axis
sns.lineplot(data=appended_df[appended_df['Type'] == skygrid_label], x='Date', 
             y='log.Pop', ax=ax2, label=skygrid_label, color='orange')

# Set labels and legend
ax1.set_ylabel('Cases: ' + moh_label)
ax2.set_ylabel('log.Pop: ' + skygrid_label)
# Adjust the legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2)

# Add Pearson correlation coefficient and P-value to the plot
text_x = ax1.get_xlim()[0]  # Get the left x-axis limit for placing text
text_y = ax1.get_ylim()[1]  # Get the upper y-axis limit for placing text
ax1.text(text_x, text_y, f'Pearson Coefficient: {correlation_coefficient:.2f}\nP-value: {p_value:.2e}', verticalalignment='top')



# Set the lower x-limit to 2012
plt.xlim(left=project_filters['Year'][0], right=project_filters['Year'][-1])

# Show the plot
plt.show() 