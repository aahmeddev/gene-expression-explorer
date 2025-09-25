import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import plotly.express as px
from scipy import stats

### App Configuration ###

st.set_page_config(
    page_title="Universal Gene Expression Explorer",
    page_icon="🧬",
    layout="wide"
)

### Utility Functions ###

def detect_file_format(uploaded_file):
    """Detect file format from uploaded file"""
    file_name = uploaded_file.name.lower()
    if file_name.endswith('.csv') or file_name.endswith('.csv.gz'):
        return 'csv'
    elif file_name.endswith('.tsv') or file_name.endswith('.txt') or file_name.endswith('.txt.gz'):
        return 'tsv'
    elif file_name.endswith(('.xlsx', '.xls')):
        return 'excel'
    else:
        return 'csv'  # Default to CSV

def load_uploaded_file(uploaded_file):
    """Load uploaded file based on its format"""
    file_format = detect_file_format(uploaded_file)

    try:
        if file_format == 'csv':
            # Try different separators and compression
            if uploaded_file.name.endswith('.gz'):
                df = pd.read_csv(uploaded_file, compression='gzip', index_col=0)
            else:
                df = pd.read_csv(uploaded_file, index_col=0)
        elif file_format == 'tsv':
            if uploaded_file.name.endswith('.gz'):
                df = pd.read_csv(uploaded_file, compression='gzip', sep='\t', index_col=0)
            else:
                df = pd.read_csv(uploaded_file, sep='\t', index_col=0)
        elif file_format == 'excel':
            df = pd.read_excel(uploaded_file, index_col=0)
        else:
            df = pd.read_csv(uploaded_file, index_col=0)

        return df, None

    except Exception as e:
        return None, str(e)

def detect_sample_groups(columns, num_samples):
    """Auto-detect potential sample groups based on column names"""
    suggestions = []

    # Check for common prefixes/suffixes
    prefixes = {}
    suffixes = {}

    for col in columns:
        col_str = str(col).upper()

        # Extract potential prefixes (first 3-5 characters)
        for i in range(2, min(6, len(col_str))):
            prefix = col_str[:i]
            if prefix not in prefixes:
                prefixes[prefix] = []
            prefixes[prefix].append(col)

        # Extract potential suffixes
        for i in range(2, min(6, len(col_str))):
            suffix = col_str[-i:]
            if suffix not in suffixes:
                suffixes[suffix] = []
            suffixes[suffix].append(col)

    # Find the most balanced split
    best_split = None
    best_balance = float('inf')

    # Check prefixes
    for prefix, cols in prefixes.items():
        if 2 <= len(cols) <= num_samples - 2:  # Must have at least 2 samples per group
            other_cols = [c for c in columns if c not in cols]
            balance = abs(len(cols) - len(other_cols))
            if balance < best_balance:
                best_balance = balance
                best_split = (f"Samples starting with '{prefix}'", cols, other_cols)

    # Check suffixes  
    for suffix, cols in suffixes.items():
        if 2 <= len(cols) <= num_samples - 2:
            other_cols = [c for c in columns if c not in cols]
            balance = abs(len(cols) - len(other_cols))
            if balance < best_balance:
                best_balance = balance
                best_split = (f"Samples ending with '{suffix}'", cols, other_cols)

    # Check for common patterns
    common_patterns = ['TUMOR', 'NORMAL', 'CONTROL', 'TREATED', 'CASE', 'CTRL', 'T', 'N', 'C','M','F']
    for pattern in common_patterns:
        matching_cols = [c for c in columns if pattern in str(c).upper()]
        if 2 <= len(matching_cols) <= num_samples - 2:
            other_cols = [c for c in columns if c not in matching_cols]
            balance = abs(len(matching_cols) - len(other_cols))
            if balance < best_balance:
                best_balance = balance
                best_split = (f"Samples containing '{pattern}'", matching_cols, other_cols)

    if best_split:
        suggestions.append({
            'description': best_split[0],
            'group1': best_split[1],
            'group2': best_split[2]
        })

    # Always offer 50/50 split as backup
    mid_point = len(columns) // 2
    suggestions.append({
        'description': f"First {mid_point} vs Last {len(columns) - mid_point} samples",
        'group1': columns[:mid_point],
        'group2': columns[mid_point:]
    })

    return suggestions

@st.cache_data
def process_gene_expression_data(df, group1_samples, group2_samples, group1_name, group2_name, min_expression=1):
    """
    Process gene expression data for any dataset structure
    """
    # Validate sample groups
    all_samples = list(group1_samples) + list(group2_samples)
    available_samples = [s for s in all_samples if s in df.columns]

    if len(available_samples) < len(all_samples):
        missing = set(all_samples) - set(available_samples)
        st.warning(f"Some samples not found in data: {missing}")

    # Filter for available samples
    group1_available = [s for s in group1_samples if s in df.columns]
    group2_available = [s for s in group2_samples if s in df.columns]

    if len(group1_available) < 2 or len(group2_available) < 2:
        raise ValueError("Need at least 2 samples per group for statistical analysis")

    # Create metadata
    metadata = pd.DataFrame(index=group1_available + group2_available)
    metadata['Condition'] = [group1_name] * len(group1_available) + [group2_name] * len(group2_available)

    # Filter expression data
    expression_df = df[group1_available + group2_available].copy()

    # Filter low-expression genes
    expression_df = expression_df[expression_df.mean(axis=1) > min_expression]

    # Log transform
    log_transformed_df = np.log2(expression_df + 1)

    # Statistical analysis
    group1_data = log_transformed_df[group1_available]
    group2_data = log_transformed_df[group2_available]

    # Calculate log2 fold change and p-values
    log2fc = group1_data.mean(axis=1) - group2_data.mean(axis=1)

    # Handle statistical test
    p_values = []
    for gene in log_transformed_df.index:
        try:
            statistic, pval = stats.ttest_ind(
                group1_data.loc[gene].dropna(),
                group2_data.loc[gene].dropna()
            )
            p_values.append(pval)
        except:
            p_values.append(1.0)  # No difference if test fails

    # Create results DataFrame
    results_df = pd.DataFrame({
        'log2FC': log2fc,
        'pvalue': p_values
    })
    results_df['-log10(pvalue)'] = -np.log10(results_df['pvalue'])
    results_df = results_df.replace([np.inf, -np.inf], np.nan).dropna()

    return log_transformed_df, metadata, results_df

### Plotting Functions ###

def create_heatmap_figure(data_matrix, metadata, group1_name, group2_name):
    """Generates a clustered heatmap figure with custom group names."""
    # Create color palette using group names
    condition_palette = {group2_name: '#3498db', group1_name: '#e74c3c'}
    col_colors = metadata['Condition'].map(condition_palette)

    plt.figure(figsize=(12, 10))
    g = sns.clustermap(
        data_matrix,
        z_score=0,
        cmap="RdBu_r",
        col_colors=col_colors,
        yticklabels=False,
        figsize=(12, 10),
        cbar_kws={'label': 'Z-score'}
    )

    g.ax_row_dendrogram.set_visible(False)

    # Create legend with actual group names
    handles = [mpatches.Patch(color=color, label=label)
               for label, color in condition_palette.items()]

    g.ax_heatmap.legend(
        handles=handles,
        title='Sample Groups',
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        frameon=True
    )

    plt.title(f'Gene Expression Heatmap: {group1_name} vs {group2_name}', 
              pad=20, fontsize=14, fontweight='bold')

    return g.fig

def create_volcano_plot_figure(results_df, group1_name, group2_name, pval_threshold=0.05, fc_threshold=1):
    """Generates an interactive volcano plot with customizable thresholds."""
    results_df = results_df.copy()
    results_df['significant'] = (results_df['pvalue'] < pval_threshold) & (abs(results_df['log2FC']) > fc_threshold)
    results_df['gene_name'] = results_df.index

    # Color by significance and direction
    def get_color(row):
        if not row['significant']:
            return 'gray'
        elif row['log2FC'] > 0:
            return 'red'
        else:
            return 'blue'

    results_df['color'] = results_df.apply(get_color, axis=1)

    fig = px.scatter(
        results_df, 
        x='log2FC', 
        y='-log10(pvalue)',
        color='color',
        color_discrete_map={'red': '#e74c3c', 'blue': '#3498db', 'gray': '#95a5a6'},
        hover_data=['gene_name'],
        title=f"Volcano Plot: {group1_name} vs {group2_name}",
        labels={
            'log2FC': f'log2(FC) [{group1_name}/{group2_name}]',
            '-log10(pvalue)': '-log10(p-value)'
        }
    )

    # Add significance thresholds
    fig.add_hline(y=-np.log10(pval_threshold), line_dash="dash", line_color="black", 
                  annotation_text=f"p-value = {pval_threshold}")
    fig.add_vline(x=fc_threshold, line_dash="dash", line_color="black")
    fig.add_vline(x=-fc_threshold, line_dash="dash", line_color="black")

    # Update layout
    fig.update_layout(
        showlegend=False,
        width=800,
        height=600
    )

    return fig

def create_boxplot_figure(log_data, metadata, gene, group1_name, group2_name):
    """Generates an interactive box plot for a specific gene."""
    gene_data = log_data.loc[gene]
    df_to_plot = pd.DataFrame({
        'Expression (log2)': gene_data,
        'Condition': metadata.loc[gene_data.index, 'Condition']
    })

    fig = px.box(
        df_to_plot, 
        x='Condition', 
        y='Expression (log2)',
        color='Condition', 
        color_discrete_map={group2_name: '#3498db', group1_name: '#e74c3c'},
        title=f"Expression of {gene}",
        points="all"  # Show all individual points
    )

    return fig

### Main App UI ###

st.title("Universal Gene Expression Explorer")
st.markdown("Upload your gene expression data and explore differential expression patterns!")

# Sidebar for file upload and parameters
st.sidebar.header("Data Upload")

uploaded_file = st.sidebar.file_uploader(
    "Choose a gene expression file",
    type=['csv', 'tsv', 'txt', 'xlsx', 'xls', 'gz'],
    help="Upload CSV, TSV, Excel, or compressed (.gz) files. First column should be gene names/IDs."
)

if uploaded_file is not None:
    # Load the data
    with st.spinner("Loading data..."):
        df, error = load_uploaded_file(uploaded_file)

    if df is None:
        st.error(f"Error loading file: {error}")
        st.stop()

    # Display data info
    st.success(f"✅ Data loaded successfully!")
    st.info(f"**Dataset Info:** {df.shape[0]} genes × {df.shape[1]} samples")

    # Show data preview
    with st.expander("Data Preview"):
        st.dataframe(df.head())

        col1, col2 = st.columns(2)
        with col1:
            st.write("**Sample Names:**")
            st.write(list(df.columns))
        with col2:
            st.write("**Gene Names (first 10):**")  
            st.write(list(df.index[:10]))

    # Sample grouping section
    st.header("Sample Grouping")
    st.markdown("Define which samples belong to each experimental group:")

    # Auto-detect potential groupings
    suggestions = detect_sample_groups(df.columns.tolist(), len(df.columns))

    grouping_method = st.radio(
        "How would you like to define sample groups?",
        ["Auto-detect groups", "Manual selection", "Upload metadata file"]
    )

    group1_samples = []
    group2_samples = []
    group1_name = "Group1"
    group2_name = "Group2"

    if grouping_method == "Auto-detect groups":
        if suggestions:
            suggestion = st.selectbox(
                "Select an auto-detected grouping:",
                range(len(suggestions)),
                format_func=lambda i: f"{suggestions[i]['description']} ({len(suggestions[i]['group1'])} vs {len(suggestions[i]['group2'])} samples)"
            )

            selected_suggestion = suggestions[suggestion]
            group1_samples = selected_suggestion['group1']
            group2_samples = selected_suggestion['group2']

            col1, col2 = st.columns(2)
            with col1:
                group1_name = st.text_input("Group 1 name:", value="Group1")
                st.write(f"**Samples ({len(group1_samples)}):**")
                st.write(group1_samples)

            with col2:
                group2_name = st.text_input("Group 2 name:", value="Group2")
                st.write(f"**Samples ({len(group2_samples)}):**")
                st.write(group2_samples)
        else:
            st.warning("Could not auto-detect sample groups. Please use manual selection.")
            grouping_method = "Manual selection"

    elif grouping_method == "Manual selection":
        col1, col2 = st.columns(2)

        with col1:
            group1_name = st.text_input("Group 1 name:", value="Treatment")
            group1_samples = st.multiselect(
                f"Select {group1_name} samples:",
                df.columns.tolist(),
                help="Select samples that belong to the first experimental group"
            )

        with col2:
            group2_name = st.text_input("Group 2 name:", value="Control")
            # Auto-select remaining samples for group 2
            remaining_samples = [col for col in df.columns if col not in group1_samples]
            group2_samples = st.multiselect(
                f"Select {group2_name} samples:",
                remaining_samples,
                default=remaining_samples,
                help="Select samples that belong to the second experimental group"
            )

    elif grouping_method == "Upload metadata file":
        metadata_file = st.file_uploader(
            "Upload metadata file (CSV with sample names and conditions)",
            type=['csv', 'tsv', 'txt']
        )
        if metadata_file:
            try:
                metadata_df = pd.read_csv(metadata_file, index_col=0)
                if 'Condition' in metadata_df.columns or 'condition' in metadata_df.columns:
                    condition_col = 'Condition' if 'Condition' in metadata_df.columns else 'condition'
                    unique_conditions = metadata_df[condition_col].unique()

                    if len(unique_conditions) >= 2:
                        group1_name = st.selectbox("Select Group 1:", unique_conditions)
                        group2_name = st.selectbox("Select Group 2:", [c for c in unique_conditions if c != group1_name])

                        group1_samples = metadata_df[metadata_df[condition_col] == group1_name].index.tolist()
                        group2_samples = metadata_df[metadata_df[condition_col] == group2_name].index.tolist()
                    else:
                        st.error("Metadata file must have at least 2 different conditions")
                else:
                    st.error("Metadata file must have a 'Condition' column")
            except Exception as e:
                st.error(f"Error reading metadata file: {e}")

    # Validate groups
    if len(group1_samples) > 0 and len(group2_samples) > 0:
        # Analysis parameters
        st.sidebar.header("⚙️ Analysis Parameters")
        min_expression = st.sidebar.number_input(
            "Minimum expression threshold:",
            min_value=0.0,
            value=1.0,
            help="Filter out genes with mean expression below this threshold"
        )

        pval_threshold = st.sidebar.number_input(
            "P-value threshold:",
            min_value=0.001,
            max_value=0.1,
            value=0.05,
            step=0.001
        )

        fc_threshold = st.sidebar.number_input(
            "Log2 fold-change threshold:",
            min_value=0.5,
            max_value=3.0,
            value=1.0,
            step=0.1
        )

        # Process data
        with st.spinner("Processing gene expression data..."):
            try:
                log_df, meta_df, results_df = process_gene_expression_data(
                    df, group1_samples, group2_samples, group1_name, group2_name, min_expression
                )

                st.success(f"✅ Analysis complete! Analyzed {len(results_df)} genes.")

                # Summary statistics
                significant_genes = results_df[
                    (results_df['pvalue'] < pval_threshold) & 
                    (abs(results_df['log2FC']) > fc_threshold)
                ]

                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Total Genes", len(results_df))
                with col2:
                    st.metric("Significant Genes", len(significant_genes))
                with col3:
                    if len(results_df) > 0:
                        pct_sig = len(significant_genes) / len(results_df) * 100
                        st.metric("% Significant", f"{pct_sig:.1f}%")

            except Exception as e:
                st.error(f"Error processing data: {e}")
                st.stop()

        # Display options
        st.sidebar.header("Display Options")
        num_genes_heatmap = st.sidebar.slider(
            "Number of top genes for heatmap:",
            min_value=10,
            max_value=min(200, len(results_df)),
            value=min(50, len(results_df)),
            step=10
        )

        # Heatmap
        st.header("Clustered Heatmap")
        st.markdown(f"Showing the top **{num_genes_heatmap}** differentially expressed genes based on fold-change.")

        top_genes = results_df['log2FC'].abs().nlargest(num_genes_heatmap).index
        heatmap_matrix = log_df.loc[top_genes]

        heatmap_fig = create_heatmap_figure(heatmap_matrix, meta_df, group1_name, group2_name)
        st.pyplot(heatmap_fig, use_container_width=True)

        # Volcano plot
        st.header("Volcano Plot")
        st.markdown("Interactive volcano plot showing statistical significance vs. fold-change. Hover over points to see gene names.")

        volcano_fig = create_volcano_plot_figure(results_df, group1_name, group2_name, pval_threshold, fc_threshold)
        st.plotly_chart(volcano_fig, use_container_width=True)

        # Individual gene expression
        st.header("Individual Gene Expression")
        significant_gene_list = significant_genes.index.tolist()

        if significant_gene_list:
            # Sort by absolute fold change for better selection
            significant_genes_sorted = significant_genes.reindex(
                significant_genes['log2FC'].abs().sort_values(ascending=False).index
            )

            selected_gene = st.selectbox(
                "Select a gene to display:",
                options=significant_genes_sorted.index.tolist(),
                help="Genes are sorted by fold-change magnitude"
            )

            if selected_gene:
                # Show gene statistics
                gene_stats = results_df.loc[selected_gene]
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Log2 Fold Change", f"{gene_stats['log2FC']:.3f}")
                with col2:
                    st.metric("P-value", f"{gene_stats['pvalue']:.2e}")
                with col3:
                    st.metric("-log10(P-value)", f"{gene_stats['-log10(pvalue)']:.2f}")

                # Box plot
                boxplot_fig = create_boxplot_figure(log_df, meta_df, selected_gene, group1_name, group2_name)
                st.plotly_chart(boxplot_fig, use_container_width=True)
        else:
            st.warning(f"No significant genes found with current thresholds (p < {pval_threshold}, |log2FC| > {fc_threshold})")
            st.info("Try adjusting the thresholds in the sidebar to see more genes.")

        # Download results
        st.header("Download Results")

        # Prepare downloadable data
        download_df = results_df.copy()
        download_df['significant'] = (
            (download_df['pvalue'] < pval_threshold) & 
            (abs(download_df['log2FC']) > fc_threshold)
        )
        download_df = download_df.sort_values('pvalue')

        csv = download_df.to_csv()
        st.download_button(
            label="📥 Download Results (CSV)",
            data=csv,
            file_name=f"differential_expression_results_{group1_name}_vs_{group2_name}.csv",
            mime="text/csv"
        )

    else:
        st.warning("⚠️ Please define both sample groups to proceed with analysis.")
        st.info("Select samples for each experimental group using the options above.")

else:
    # Instructions when no file is uploaded
    st.markdown("""
    ## Instructions

    1. **Upload your gene expression data** using the file uploader in the sidebar
    2. **Supported formats:** CSV, TSV, Excel (.xlsx), or compressed files (.gz)
    3. **Data structure:** Genes as rows, samples as columns, first column should be gene names/IDs
    4. **Define sample groups** to compare (e.g., Treatment vs Control, Tumor vs Normal)
    5. **Explore results** with interactive heatmaps, volcano plots, and individual gene plots

    ### Example Data Format:
    ```
    Gene_ID    Sample1  Sample2  Sample3  Sample4  ...
    GENE001    12.5     8.3      15.2     9.1      ...
    GENE002    5.7      4.2      6.8      3.9      ...
    GENE003    0.1      0.0      2.1      0.3      ...
    ```

    ### Features:
    - **Auto-detection** of sample groups based on naming patterns
    - **Manual selection** of samples for each group
    - **Metadata upload** for complex experimental designs
    - **Interactive visualizations** with Plotly and Seaborn
    - **Statistical analysis** with t-tests and multiple testing correction
    - **Downloadable results** in CSV format
    """)
