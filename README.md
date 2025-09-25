# Gene Expression Explorer

An interactive web app for gene expression data analysis and visualization.
Upload your dataset, define sample groups easily, view heatmaps, volcano plots, gene boxplots, and download results—all powered by Streamlit and Python.

## Features

- **Universal File Upload:** Supports CSV, TSV, Excel (.xlsx), and compressed (.gz) files.
- **Flexible Group Assignment:** Automatically detects sample groups from column names or allows manual selection/upload of metadata.
- **Advanced Visualizations:** Clustered heatmaps, volcano plots, and gene expression boxplots.
- **Customizable Analysis:** Adjust p-value and log2 fold change thresholds, filter by minimum gene expression.
- **Downloadable Results:** Export analysis summaries and significance flags as CSV.
- **User-Friendly Interface:** Simple data upload, immediate data preview, and easy navigation.

## Example Data Format

| Gene_ID   | Sample1 | Sample2 | ... |
|-----------|---------|---------|-----|
| GENE001   | 12.5    | 8.3     | ... |
| GENE002   | 5.7     | 4.2     | ... |
| ...       | ...     | ...     | ... |

- Genes as rows (first column = gene names/IDs)
- Samples as columns (any naming pattern)
- Expression values should be numeric
- At least two samples per group for statistical analysis

## Installation

`git clone https://github.com/aahmeddev/universal-gene-expression-explorer.git`

`cd universal-gene-expression-explorer`

`pip install -r requirements.txt`


# Usage
`streamlit run geo_app.py`

- Upload your gene expression data file using the sidebar.
- Define your sample groups (auto-detection, manual, or metadata file).
- Adjust analysis parameters (optional) in the sidebar.
- View clustered heatmaps, volcano plots, and gene boxplots.
- Download your results as CSV for further analysis.

# Technologies Used
- Streamlit
- Pandas
- Numpy
- Scipy
- Seaborn
- Matplotlib
- Plotly
