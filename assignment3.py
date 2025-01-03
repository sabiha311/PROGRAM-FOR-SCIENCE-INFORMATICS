# 1 A .import pandas as pd

# Load the Excel file
gene_expression_data = pd.read_excel('Gene_Expression_Data.xlsx')

# Load the CSV file
gene_information = pd.read_csv('Gene_Information.csv')

# Load the TSV file
sample_information = pd.read_csv('Sample_Information.tsv', sep='\t')

# Display the first few rows of each dataframe
print("Gene Expression Data:")
print(gene_expression_data.head())

print("\nGene Information:")
print(gene_information.head())

print("\nSample Information:")
print(sample_information.head())
B.import pandas as pd

 1 b .import pandas as pd

# Load the data
gene_expression_data = pd.read_excel('Gene_Expression_Data.xlsx')
sample_information = pd.read_csv('Sample_Information.tsv', sep='\t')

# Display the first few rows of each dataframe for reference
print("Original Gene Expression Data:")
print(gene_expression_data.head())

print("\nSample Information:")
print(sample_information.head())

# Assuming the first column in sample_information contains the original sample names
# and a column named 'Phenotype' contains the corresponding phenotypes.

# Create a mapping from original sample names to new names based on phenotype
sample_mapping = dict(zip(sample_information.iloc[:, 0], sample_information['group']))

# Rename columns in gene expression data
# Exclude the first column if it's an identifier (like gene names)
original_columns = gene_expression_data.columns[1:]  # Adjust as needed
new_columns = [sample_mapping.get(col, col) for col in original_columns]  # Map to new names

# Adding suffixes to ensure unique names
new_columns_with_suffixes = []
seen = {}
for col in new_columns:
    if col in seen:
        seen[col] += 1
        new_columns_with_suffixes.append(f"{col}_{seen[col]}")
    else:
        seen[col] = 1
        new_columns_with_suffixes.append(col)

# Update the column names in the gene expression data
gene_expression_data.columns = [gene_expression_data.columns[0]] + new_columns_with_suffixes

# Display the updated gene expression data
print("\nUpdated Gene Expression Data:")
print(gene_expression_data.head())
1.c .import pandas as pd

# Load the data again if needed
gene_expression_data = pd.read_excel('Gene_Expression_Data.xlsx')
sample_information = pd.read_csv('Sample_Information.tsv', sep='\t')

# Create a mapping from original sample names to new names based on phenotype
# The column name is changed from 'Phenotype' to 'group'
sample_mapping = dict(zip(sample_information.iloc[:, 0], sample_information['group']))

# Rename columns in gene expression data
original_columns = gene_expression_data.columns[1:]  # Adjust if necessary
new_columns = [sample_mapping.get(col, col) for col in original_columns]

# Adding suffixes to ensure unique names
new_columns_with_suffixes = []
seen = {}
for col in new_columns:
    if col in seen:
        seen[col] += 1
        new_columns_with_suffixes.append(f"{col}_{seen[col]}")
    else:
        seen[col] = 1
        new_columns_with_suffixes.append(col)

# Update the column names in the gene expression data
gene_expression_data.columns = [gene_expression_data.columns[0]] + new_columns_with_suffixes

# Split the data based on phenotype
tumor_data = gene_expression_data.loc[:, gene_expression_data.columns.str.contains('tumor', case=False)]
normal_data = gene_expression_data.loc[:, gene_expression_data.columns.str.contains('normal', case=False)]

# Optionally, include the gene names (assuming they are in the first column)
tumor_data = pd.concat([gene_expression_data.iloc[:, 0], tumor_data], axis=1)
normal_data = pd.concat([gene_expression_data.iloc[:, 0], normal_data], axis=1)

# Display the split data
print("\nTumor Data:")
print(tumor_data.head())

print("\nNormal Data:")
print(normal_data.head())

1d..import pandas as pd

# Load the data again if needed (assuming previous steps are done)
gene_expression_data = pd.read_excel('Gene_Expression_Data.xlsx')
sample_information = pd.read_csv('Sample_Information.tsv', sep='\t')

# Create a mapping from original sample names to new names based on phenotype
# Change 'Phenotype' to 'group' to match the actual column name
sample_mapping = dict(zip(sample_information.iloc[:, 0], sample_information['group']))

# Rename columns in gene expression data
original_columns = gene_expression_data.columns[1:]  # Adjust if necessary
new_columns = [sample_mapping.get(col, col) for col in original_columns]

# Adding suffixes to ensure unique names
new_columns_with_suffixes = []
seen = {}
for col in new_columns:
    if col in seen:
        seen[col] += 1
        new_columns_with_suffixes.append(f"{col}_{seen[col]}")
    else:
        seen[col] = 1
        new_columns_with_suffixes.append(col)

# Update the column names in the gene expression data
gene_expression_data.columns = [gene_expression_data.columns[0]] + new_columns_with_suffixes

# Split the data based on phenotype
tumor_data = gene_expression_data.loc[:, gene_expression_data.columns.str.contains('tumor', case=False)]
normal_data = gene_expression_data.loc[:, gene_expression_data.columns.str.contains('normal', case=False)]

# Optionally, include the gene names (assuming they are in the first column)
tumor_data = pd.concat([gene_expression_data.iloc[:, 0], tumor_data], axis=1)
normal_data = pd.concat([gene_expression_data.iloc[:, 0], normal_data], axis=1)

# Compute the average expression for each probe
# Exclude the first column (Probe_ID) when calculating the mean

# Change here: Select only numeric columns for calculating mean
tumor_average = tumor_data.select_dtypes(include='number').mean(axis=1)
normal_average = normal_data.select_dtypes(include='number').mean(axis=1)

# Create DataFrames for average expressions
tumor_average_df = pd.DataFrame({
    'Probe': tumor_data.iloc[:, 0],  # First column is probes
    'Average Expression (Tumor)': tumor_average
})

normal_average_df = pd.DataFrame({
    'Probe': normal_data.iloc[:, 0],  # First column is probes
    'Average Expression (Normal)': normal_average
1 e.import pandas as pd

# Load the data again if needed (assuming previous steps are done)
gene_expression_data = pd.read_excel('Gene_Expression_Data.xlsx')
sample_information = pd.read_csv('Sample_Information.tsv', sep='\t')

# Create a mapping from original sample names to new names based on phenotype
sample_mapping = dict(zip(sample_information.iloc[:, 0], sample_information['group']))

# Rename columns in gene expression data
original_columns = gene_expression_data.columns[1:]  # Adjust if necessary
new_columns = [sample_mapping.get(col, col) for col in original_columns]

# Adding suffixes to ensure unique names
new_columns_with_suffixes = []
seen = {}
for col in new_columns:
    if col in seen:
        seen[col] += 1
        new_columns_with_suffixes.append(f"{col}_{seen[col]}")
    else:
        seen[col] = 1
        new_columns_with_suffixes.append(col)

# Update the column names in the gene expression data
gene_expression_data.columns = [gene_expression_data.columns[0]] + new_columns_with_suffixes

# Split the data based on phenotype
tumor_data = gene_expression_data.loc[:, gene_expression_data.columns.str.contains('tumor', case=False)]
normal_data = gene_expression_data.loc[:, gene_expression_data.columns.str.contains('normal', case=False)]

# Optionally, include the gene names (assuming they are in the first column)
tumor_data = pd.concat([gene_expression_data.iloc[:, 0], tumor_data], axis=1)
normal_data = pd.concat([gene_expression_data.iloc[:, 0], normal_data], axis=1)

# Compute the average expression for each probe
tumor_average = tumor_data.drop(columns=[tumor_data.columns[0]]).mean(axis=1)
normal_average = normal_data.drop(columns=[normal_data.columns[0]]).mean(axis=1)
# Create DataFrames for average expressions
tumor_average_df = pd.DataFrame({
    'Probe': tumor_data.iloc[:, 0],  # First column is probes
    'Average Expression (Tumor)': tumor_average
})

normal_average_df = pd.DataFrame({
    'Probe': normal_data.iloc[:, 0],  # First column is probes
    'Average Expression (Normal)': normal_average
})

# Merge the two average DataFrames on the Probe column
merged_averages = pd.merge(tumor_average_df, normal_average_df, on='Probe')

# Compute fold change: (Tumor - Control) / Control
merged_averages['Fold Change'] = (merged_averages['Average Expression (Tumor)'] - merged_averages['Average Expression (Normal)']) / merged_averages['Average Expression (Normal)']

# Display the results
print("\nFold Change for each Probe:")
print(merged_averages[['Probe', 'Fold Change']].head())
1 f.Display the first few rows of merged_averages to check fold change values
print("Merged Averages DataFrame:")
print(merged_averages.head())

# Display the first few rows of gene_information to check identifiers
print("\nGene Information DataFrame:")
print(gene_information.head())

# Filter for absolute fold change greater than 5
significant_genes = merged_averages[merged_averages['Fold Change'].abs() > 5]

# Check if significant_genes is empty
if significant_genes.empty:
    print("\nNo genes with absolute fold change greater than 5 found.")
else:
    # Merge with gene information to get gene details
    result = pd.merge(significant_genes, gene_information, left_on='Probe', right_on='Gene_ID', how='left')

    # Display the significant genes
    print("\nSignificant Genes with Fold Change Magnitude Greater than 5:")
    print(result[['Probe', 'Fold Change', 'Gene_Name']].head())  # Adjust column names as needed
1 g .import pandas as pd

# Load the gene information
gene_information = pd.read_csv('Gene_Information.csv')

# Assuming merged_averages has been created with average expressions and fold changes
# Filter for absolute fold change greater than 5
significant_genes = merged_averages[merged_averages['Fold Change'].abs() > 5]

# Check if there are any significant genes
if not significant_genes.empty:
    # Merge with gene information to get gene details
    result = pd.merge(significant_genes, gene_information, left_on='Probe', right_on='Gene_ID', how='left')

    # Add a new column for expression status
    result['Expression Status'] = result.apply(
        lambda row: 'Higher in Tumor' if row['Average Expression (Tumor)'] > row['Average Expression (Normal)']
        else 'Higher in Normal' if row['Average Expression (Tumor)'] < row['Average Expression (Normal)']
        else 'Equal', axis=1
    )

    # Display the significant genes with expression status
    print("\nSignificant Genes with Fold Change Magnitude Greater than 5 and Expression Status:")
    print(result[['Probe', 'Fold Change', 'Gene_Name', 'Expression Status']])
else:
    print("No significant genes found with fold change magnitude greater than 5.")

2 a 2A.import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your significant genes DataFrame (result from part 1g)
# result = pd.read_csv('your_significant_genes_data.csv')  # Uncomment if loading from CSV

# 1. Overview of the Data
print(result.info())
print(result.describe())

# 2. Visualize Distributions

# Histogram of Fold Change
plt.figure(figsize=(10, 6))
sns.histplot(result['Fold Change'], bins=30, kde=True)
plt.title('Distribution of Fold Change')
plt.xlabel('Fold Change')
plt.ylabel('Frequency')
plt.show()

# Box Plot of Fold Change
plt.figure(figsize=(10, 6))
sns.boxplot(x=result['Fold Change'])
plt.title('Box Plot of Fold Change')
plt.xlabel('Fold Change')
plt.show()

# 3. Average Expression Distributions
plt.figure(figsize=(12, 6))
sns.histplot(result['Average Expression (Tumor)'], bins=30, kde=True, label='Tumor', color='red', alpha=0.5)
sns.histplot(result['Average Expression (Normal)'], bins=30, kde=True, label='Normal', color='blue', alpha=0.5)
plt.title('Distribution of Average Expressions')
plt.xlabel('Expression Level')
plt.ylabel('Frequency')
plt.legend()
plt.show()

# Box Plot for Average Expressions
plt.figure(figsize=(12, 6))
sns.boxplot(data=result[['Average Expression (Tumor)', 'Average Expression (Normal)']])
plt.title('Box Plot of Average Expressions')
plt.ylabel('Expression Level')
plt.xticks([0, 1], ['Tumor', 'Normal'])
plt.show()

# 4. Correlation Analysis (if applicable)
# Only applicable if there are more numeric variables to analyze
correlation_matrix = result[['Fold Change', 'Average Expression (Tumor)', 'Average Expression (Normal)']].corr()
plt.figure(figsize=(8, 6))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt='.2f')
plt.title('Correlation Matrix')
plt.show()

# 5. Categorical Analysis: Expression Status Distribution
plt.figure(figsize=(8, 6))
sns.countplot(x='Expression Status', data=result)
plt.title('Expression Status Distribution')
plt.xlabel('Expression Status')
plt.ylabel('Count')
plt.xticks(rotation=45)
plt.show()

2 b.import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your significant genes DataFrame (result from part 1g)
# result = pd.read_csv('your_significant_genes_data.csv')  # Uncomment if loading from CSV

# 1. Count the number of DEGs by chromosome
deg_counts = result['Chromosome'].value_counts()

# 2. Create a histogram (bar plot) of DEGs by chromosome
plt.figure(figsize=(10, 6))
sns.barplot(x=deg_counts.index, y=deg_counts.values, palette='viridis')
plt.title('Distribution of Differentially Expressed Genes (DEGs) by Chromosome')
plt.xlabel('Chromosome')
plt.ylabel('Number of DEGs')
plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
plt.show()
2.c .C.import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your significant genes DataFrame (result from part 1g)
# result = pd.read_csv('your_significant_genes_data.csv')  # Uncomment if loading from CSV

# 1. Ensure there is a 'Sample Type' column indicating Normal or Tumor
# If 'Sample Type' is not present, create it based on your data structure
# For example, if you have separate Tumor and Normal DataFrames used to create 'result':
result['Sample Type'] = result.apply(lambda row: 'Tumor' if row['Average Expression (Tumor)'] > row['Average Expression (Normal)'] else 'Normal', axis=1) 
# This line adds the missing 'Sample Type' column based on the expression values.

# Alternatively, if 'Sample Type' is in a different column, rename it:
# result = result.rename(columns={'YourColumnName': 'Sample Type'})

# 2. Count the number of DEGs by chromosome and sample type
deg_counts = result.groupby(['Chromosome', 'Sample Type']).size().reset_index(name='Count')

# 3. Create a bar plot of DEGs by chromosome segregated by sample type
plt.figure(figsize=(12, 6))
sns.barplot(x='Chromosome', y='Count', hue='Sample Type', data=deg_counts, palette='Set2')
plt.title('Distribution of Differentially Expressed Genes (DEGs) by Chromosome Segregated by Sample Type')
plt.xlabel('Chromosome')
plt.ylabel('Number of DEGs')
plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
plt.legend(title='Sample Type')
plt.show()
2.d.#D.import pandas as pd
import matplotlib.pyplot as plt

# Load your significant genes DataFrame (result from part 1g)
# result = pd.read_csv('your_significant_genes_data.csv')  # Uncomment if loading from CSV

# 1. Count the number of upregulated and downregulated genes
upregulated_count = result[result['Fold Change'] > 1].shape[0]  # Adjust the threshold if needed
downregulated_count = result[result['Fold Change'] < -1].shape[0]  # Adjust the threshold if needed

# 2. Calculate total DEGs
total_degs = upregulated_count + downregulated_count

# 3. Calculate percentages
# Check if total_degs is zero to avoid ZeroDivisionError
if total_degs != 0:
    upregulated_percentage = (upregulated_count / total_degs) * 100
    downregulated_percentage = (downregulated_count / total_degs) * 100
else:
    upregulated_percentage = 0  # or any other suitable value
    downregulated_percentage = 0  # or any other suitable value

# 4. Prepare data for plotting
percentages = [upregulated_percentage, downregulated_percentage]
labels = ['Upregulated in Tumor', 'Downregulated in Tumor']

# 5. Create bar chart
plt.figure(figsize=(8, 5))
plt.bar(labels, percentages, color=['blue', 'orange'])
plt.title('Percentage of DEGs Upregulated and Downregulated in Tumor Samples')
plt.ylabel('Percentage (%)')
plt.ylim(0, 100)  # Set y-axis limits
plt.show()

2.e .import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the raw gene expression data
# gene_expression_data = pd.read_excel('Gene_Expression_Data.xlsx')  # Uncomment if loading from Excel
# gene_expression_data = pd.read_csv('Gene_Expression_Data.csv')  # Uncomment if loading from CSV

# Assume the first column is gene names and set it as index
gene_expression_data.set_index(gene_expression_data.columns[0], inplace=True)

# Create a heatmap
plt.figure(figsize=(12, 8))
sns.heatmap(gene_expression_data, cmap='viridis')
plt.title('Heatmap of Gene Expression by Sample')
plt.xlabel('Samples')
plt.ylabel('Genes')
plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
plt.show()
2.f import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the raw gene expression data
gene_expression_data = pd.read_excel('Gene_Expression_Data.xlsx')  # Uncomment and adjust the file path if needed
# gene_expression_data = pd.read_excel('Gene_Expression_Data.xlsx')  # Uncomment if loading from Excel
# gene_expression_data = pd.read_csv('Gene_Expression_Data.csv')  # Uncomment if loading from CSV

# Set the first column (gene names) as index
gene_expression_data.set_index(gene_expression_data.columns[0], inplace=True)

# Create a clustermap
sns.clustermap(gene_expression_data, cmap='viridis', figsize=(12, 8))
plt.title('Clustermap of Gene Expression by Sample')
plt.show()
2.g from the visualisations several key findings are obtained to make a gene expression data .The heatmap hihlighted revealing cluster of genes awith high expression profiles.
the analysis of DEGs showed a notable percentage of genes that were upregulated in tumor samples compared to normal samples, as illustrated in the bar chart. This suggests potential pathways that may be activated in tumorigenesis. .
The clustermap further reinforced these insights by organizing genes and samples based on their expression similarities, allowing us to identify groups of differentially expressed genes (DEGs) that are significantly altered between tumor and normal samples..
Overall, these findings provide a comprehensive overview of gene expression alterations associated with tumor samples, highlighting potential targets for further investigation and therapeutic intervention..



