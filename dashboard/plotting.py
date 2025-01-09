"""
Plotting utilities for the sORF dashboard.

This module contains functions for generating various plots and visualizations used in the dashboard,
including heatmaps, bar plots, and protein feature plots.
"""

import altair as alt
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
import simplejson as json

from scipy.cluster.hierarchy import linkage, leaves_list
from streamlit_echarts import JsCode
from dashboard.util import query_de_transcripts, query_transcript_tpms, strip_ensembl_versions
from dashboard.etl import CACHE_DIR
import sqlite3


def color_protein_terminal_ends(aa_seq: str, pdb_string: str) -> str:
    """
    Modify PDB string to color protein terminal ends.

    Args:
        aa_seq: Amino acid sequence
        pdb_string: PDB structure as string

    Returns:
        Modified PDB string with colored terminal ends (C-term red, N-term blue)
    """
    parts = pdb_string.split('\n')
    for ix, line in enumerate(parts):
        if line.startswith('ATOM'):
            pieces = line.split()
            if pieces[5] == '1':
                parts[ix] = line.replace(f'1.00 {pieces[-2]}', '1.00 100.00')
            elif pieces[5] == str(len(aa_seq)):
                parts[ix] = line.replace(f'1.00 {pieces[-2]}', '1.00 0.00')
    return '\n'.join(parts)


def expression_de_to_echarts_data(expression_values: list, de_values: list, color: str) -> list:
    """
    Convert expression and differential expression values to ECharts data format.

    Args:
        expression_values: List of expression values
        de_values: List of boolean values indicating differential expression
        color: Base color name ('red' or 'blue')

    Returns:
        List of dictionaries containing ECharts bar chart data
    """
    data = []
    color_map = {
        'red_de': '#d62418',
        'red_null': '#d1918c', 
        'blue_de': '#237c94',
        'blue_null': '#70aaba'
    }
    
    for e, d in zip(expression_values, de_values):
        if d:
            c = color_map[f"{color}_de"]
        else:
            c = color_map[f"{color}_null"]
        data.append({'value': e, 'itemStyle': {'color': c}})
    return data


def bar_plot_expression_groups_tcga(transcript_id: str, group_name: str, 
                                  group_members: list, title: str) -> dict:
    """
    Generate ECharts bar plot configuration for TCGA expression data.

    Args:
        transcript_id: Transcript identifier
        group_name: Name of the grouping variable
        group_members: List of group member names
        title: Plot title

    Returns:
        ECharts option dictionary for bar plot
    """
    dataframe = pd.read_sql(
        "SELECT * FROM transcript_de WHERE transcript_de.transcript_id = '{0}'".format(
            transcript_id.split('.')[0]
        ),
        sqlite3.connect(CACHE_DIR / 'xena.db')
    )
    dataframe.rename({'index': 'TCGA'}, axis=1, inplace=True)
    dataframe.rename({'TCGA Cancer Type': 'TCGA'}, axis=1, inplace=True)
    dataframe['DE'] = (
        (dataframe['padj'] < 0.0001) & 
        (np.abs(dataframe['log2FoldChange']) >= 2) & 
        (dataframe['Cancer Mean'].gt(4) | dataframe['GTEx Mean'].gt(4))
    )
    
    dataframe.sort_values(by='Cancer Mean', inplace=True)
    de_metadata = pd.read_sql(
        "SELECT * from sample_metadata_de", 
        sqlite3.connect(CACHE_DIR / 'xena.db')
    )
    tcga_code_to_description = de_metadata[['Description', 'GTEx Tissue Type']].apply(
        lambda x: f"{x[0]}<br>GTEx Normal {x[1]}", axis=1
    ).to_dict()

    option = {
        'title': {'text': title},
        'tooltip': {
            "trigger": 'axis',
            "formatter": JsCode(
                "function (params) {var cols = " + json.dumps(tcga_code_to_description) + 
                "; console.log(params); return params[0].name + ' - ' + cols[params[0].name] + " +
                "'<br>' + params[0].seriesName + ': ' + params[0].value  + '<br>' + " +
                "params[1].seriesName + ': ' + params[1].value;}"
            ).js_code,
        },
        'legend': {
            'data': group_members,
            'orient': 'vertical',
            'right': -5,
            'top': 'center'
        },
        'toolbox': {
            'show': True,
            'feature': {
                'dataView': {'show': True, 'readOnly': False},
                'restore': {'show': True},
                'saveAsImage': {'show': True}
            }
        },
        'calculable': True,
        'yAxis': [{
            'name': 'Cancer Types',
            'nameLocation': 'middle',
            'nameGap': 50,
            'type': 'category',
            'data': list(dataframe[group_name].values),
            'axisLabel': {'interval': 0}
        }],
        'xAxis': [{
            'name': 'Approximate TPM \n (Transcripts in bold are DE with FDR < 1e-5 and abs(LFC) > 1)',
            'nameLocation': 'middle',
            'nameGap': 30,
            'type': 'value'
        }],
        'color': ['#237c94', '#d62418'],
        'series': [
            {
                'name': group_members[0],
                'type': 'bar',
                'data': expression_de_to_echarts_data(
                    dataframe[group_members[0]], dataframe['DE'], 'blue'
                )
            },
            {
                'name': group_members[1],
                'type': 'bar',
                'data': expression_de_to_echarts_data(
                    dataframe[group_members[1]], dataframe['DE'], 'red'
                )
            },
        ],
        'grid': {
            'left': 70,
            'top': 50,
            'right': 120,
            'bottom': 80
        }
    }

    return option


def bar_plot_expression_group_autoimmune(df: pd.DataFrame, title: str) -> dict:
    """
    Generate ECharts bar plot configuration for autoimmune expression data.

    Args:
        df: DataFrame containing expression data
        title: Plot title

    Returns:
        ECharts option dictionary for bar plot
    """
    group_members = ['control_mean', 'case_mean']
    df['contrast'] = [f"{study}_{c}" for study, c in df[['velia_study', 'contrast']].values]
    metadata = {k:v for k, v in df[['contrast', 'velia_study']].values}
    
    option = {
        'title': {'text': title},
        'tooltip': {
            "trigger": 'axis',
            "formatter": JsCode(
                "function (params) {var cols = " + json.dumps(metadata) + 
                "; console.log(params); return params[0].name + ' - ' + cols[params[0].name] + " +
                "'<br>' + params[0].seriesName + ': ' + params[0].value  + '<br>' + " +
                "params[1].seriesName + ': ' + params[1].value;}"
            ).js_code,
            "z": 1000,
            "confine": 'true',
            "textStyle": {
                "overflow": 'breakAll',
                "width": 40,
            },
        },
        'legend': {
            'data': group_members,
            'orient': 'vertical',
            'right': -5,
            'top': 'center'
        },
        'toolbox': {
            'show': True,
            'feature': {
                'dataView': {'show': True, 'readOnly': False},
                'restore': {'show': True},
                'saveAsImage': {'show': True}
            }
        },
        'calculable': True,
        'yAxis': [{
            'name': 'Study',
            'nameLocation': 'middle',
            'nameGap': 90,
            'type': 'category',
            'data': list(df['contrast'].values),
            'axisLabel': {
                'interval': 0,
                'overflow': 'breakAll'
            }
        }],
        'xAxis': [{
            'name': 'Approximate TPM \n (Transcripts in bold are DE with FDR < 1e-5 and abs(LFC) > 1)',
            'nameLocation': 'middle',
            'nameGap': 20,
            'type': 'value'
        }],
        'color': ['#237c94', '#d62418'],
        'series': [
            {
                'name': group_members[0],
                'type': 'bar',
                'data': list(df[group_members[0]].values)
            },
            {
                'name': group_members[1],
                'type': 'bar',
                'data': list(df[group_members[1]].values)
            },
        ],
        'grid': {
            'left': 110,
            'top': 50,
            'right': 120,
            'bottom': 80
        }
    }

    return option


def expression_heatmap_plot(title: str, selected_expression: pd.DataFrame, 
                          median_groups: bool = False) -> tuple:
    """
    Generate ECharts heatmap configuration for expression data.

    Args:
        title: Plot title
        selected_expression: DataFrame containing expression values
        median_groups: Whether to use median values for groups

    Returns:
        Tuple of (ECharts option dictionary, events dictionary)
    """
    set2 = sns.color_palette('Set2', n_colors=2)
    col_colors = [set2[0] for i in selected_expression]

    grouped_exp_df = selected_expression
    if grouped_exp_df.shape[1] == 0:
        return None, None
    elif grouped_exp_df.shape[1] == 1:
        plot_df = grouped_exp_df
        plot_df.sort_values(by=selected_expression.columns[0], inplace=True)
        plot_df = plot_df.T
    else:
        row_clusters = linkage(grouped_exp_df.T.values, method='complete', metric='euclidean')
        col_clusters = linkage(grouped_exp_df.values, method='complete', metric='euclidean')

        row_leaves = leaves_list(row_clusters)
        col_leaves = leaves_list(col_clusters)

        plot_df = grouped_exp_df.T.iloc[row_leaves, col_leaves]

    col_map = {x: i for i,x in enumerate(plot_df.columns)}
    row_map = {x: i for i,x in enumerate(plot_df.index)}
    data = [(row_map[k[0]], col_map[k[1]], v) for k,v in plot_df.stack().items()]
    
    col_names = list(col_map.keys())
    row_names = list(row_map.keys())

    max_contrast = grouped_exp_df.max().max()
    if max_contrast < 5:
        max_contrast = 5
    option, events = heatmap_plot(data, row_names, col_names, title)
    option['visualMap']['max'] = max_contrast
    return option, events


def heatmap_plot(data: list, row_names: list, col_names: list, title: str,
                x_axis_interval: int = 0, y_axis_interval: int = 0, 
                scale: str = 'TPM') -> tuple:
    """
    Generate ECharts heatmap configuration.

    Args:
        data: List of data points
        row_names: List of row names
        col_names: List of column names
        title: Plot title
        x_axis_interval: X-axis label interval
        y_axis_interval: Y-axis label interval
        scale: Scale label for tooltip

    Returns:
        Tuple of (ECharts option dictionary, events dictionary)
    """
    js_col_names = "var cols = [" + ",".join([f"'{c}'" for c in col_names]) + "];"
    option = {
        "title": {"text": title},
        "tooltip": {
            "formatter": JsCode(
                "function (params) {" + js_col_names + 
                "; return params.name + '<br>' + cols[params.data[1]] + '<br> Median " + 
                scale + ": ' + params.data[2];}"
            ).js_code,
        },
        "xAxis": {
            "type": "category",
            "data": row_names,
            "axisLabel": {
                "fontSize": 10,
                "rotate": -90,
                "interval": x_axis_interval,
            }
        },
        "yAxis": {
            "type": "category",
            "data": col_names,
            "axisLabel": {
                "fontSize": 11,
                "width": 0,
                "interval": y_axis_interval,
            }
        },
        "visualMap": {
            "min": 0,
            "max": 12,
            "calculable": True,
            "realtime": False,
            "orient": 'vertical',
            "left": '90%',
            "top": 'center',
        },
        'toolbox': {
            'show': True,
            'feature': {
                'dataView': {'show': True, 'readOnly': False},
                'restore': {'show': True},
                'saveAsImage': {'show': True}
            }
        },
        "grid": {
            "left": 0,
            "containLabel": True
        },
        "series": [{
            "name": "Log2(TPM+1)",
            "type": "heatmap",
            "data": data,
            "emphasis": {
                "itemStyle": {
                    "borderColor": '#333',
                    "borderWidth": 1,
                    "shadowBlur": 10,
                    "shadowColor": 'rgba(0, 0, 0, 0.5)'
                }
            },
            "progressive": 1000,
            "animation": False,
        }],
    }
    
    events = {
        "click": "function(params) { console.log(params.name); return params.name }",
        "dblclick": "function(params) { return [params.type, params.name, params.value] }"
    }
    return option, events


def expression_atlas_heatmap_plot(tissue_specific_vtx_ids: list, 
                                xena_vtx_sum_df: pd.DataFrame,
                                selected_tissues: list) -> tuple:
    """
    Generate ECharts heatmap configuration for expression atlas data.

    Args:
        tissue_specific_vtx_ids: List of tissue-specific transcript IDs
        xena_vtx_sum_df: DataFrame containing expression data
        selected_tissues: List of selected tissues

    Returns:
        Tuple of (ECharts option dictionary, events dictionary)
    """
    xena_vtx_sum_df = xena_vtx_sum_df[tissue_specific_vtx_ids].copy()
    xena_vtx_sum_df = xena_vtx_sum_df.loc[selected_tissues]

    row_clusters = linkage(xena_vtx_sum_df.T.values, method='complete', metric='euclidean')
    col_clusters = linkage(xena_vtx_sum_df.values, method='complete', metric='euclidean')

    row_leaves = leaves_list(row_clusters)
    col_leaves = leaves_list(col_clusters)

    plot_df = xena_vtx_sum_df.T.iloc[row_leaves, col_leaves]

    col_map = {x: i for i,x in enumerate(plot_df.columns)}
    row_map = {x: i for i,x in enumerate(plot_df.index)}
    data = [(row_map[k[0]], col_map[k[1]], v) for k,v in plot_df.stack().items()]
    
    col_names = list(col_map.keys())
    row_names = list(row_map.keys())
    title = "TCGA/GTEx Expression Atlas"
    option, events = heatmap_plot(data, row_names, col_names, title, scale='Log2(TPM)')
    return option, events


def expression_vtx_boxplot(transcript_id: str, expression_df: pd.DataFrame) -> px.box:
    """
    Generate Plotly box plot for transcript expression.

    Args:
        transcript_id: Transcript identifier
        expression_df: DataFrame containing expression data

    Returns:
        Plotly box plot figure
    """
    groups = expression_df[['_primary_site', '_study']].apply(lambda x: '-'.join(x), axis=1)
    df = expression_df[[transcript_id]].copy()
    df['Sample'] = groups
    df.reset_index(inplace=True)

    medians = df.groupby('Sample')[transcript_id].median().sort_values()

    order_map = {name:i for i, name in enumerate(medians.index)}
    df['order'] = df.apply(lambda x: order_map[x['Sample']], axis=1)
    df.sort_values(by='order', ascending=False, inplace=True)

    fig = px.box(
        df, 
        x="Sample",
        y=transcript_id,
        points='outliers',
        hover_data='index',
        height=500
    )
    
    fig.update_traces(quartilemethod="exclusive")
    return fig


def format_protein_feature_strings_for_altair_heatmap(orf_features: pd.Series) -> pd.DataFrame:
    """
    Format protein feature data for Altair heatmap visualization.

    Args:
        orf_features: Series containing protein feature data

    Returns:
        DataFrame formatted for Altair heatmap
    """
    string_cols = [
        'Deepsig', 'SignalP 6slow', 'SignalP 5b', 'SignalP 4.1', 
        'Phobius', 'DeepTMHMM', 'Sequence'
    ]
    orf_features = orf_features[string_cols]
    orf_features = orf_features.T
    orf_features = orf_features.apply(lambda x: [*x])
    seq = orf_features.pop('Sequence')
    values = []
    for i, (ix, row) in enumerate(orf_features.items()):
        for j, v in enumerate(row):
            values.append([j+1, ix, v, seq[j]])
    df = pd.DataFrame(values, columns=['Position', 'Tool', 'Predicted Class', 'aa'])
    return df


def altair_protein_features_plot(df: pd.DataFrame) -> alt.Chart:
    """
    Generate Altair plot for protein features.

    Args:
        df: DataFrame containing protein feature data

    Returns:
        Altair chart object
    """
    color_dict = {
        'S': 'palegoldenrod',
        'C': 'lightgray',
        'I': 'lightgray',
        'O': 'lavenderblush',
        'M': 'lightblue'
    }
    
    color_scale = alt.Scale(domain=list(color_dict.keys()), range=list(color_dict.values()))

    base = alt.Chart(df).encode(
        alt.X('Position:O'),
        alt.Y('Tool:O')
    )
    fig = base.mark_rect().encode(
        color=alt.Color(
            'Predicted Class:O',
            scale=color_scale,
            legend=alt.Legend(
                orient='bottom',
                direction='horizontal',
                titleAnchor='middle'
            )
        ),
        tooltip=['Position', 'Predicted Class']
    )
    return fig + base.mark_text(baseline='middle').encode(alt.Text('aa:O'))


def plot_sequence_line_plots_altair(vtx_id: str, sorf_aa_seq: str,
                                  phylocsf_dataframe: pd.DataFrame,
                                  esmfold: dict) -> alt.Chart:
    """
    Generate Altair line plots for sequence features.

    Args:
        vtx_id: Transcript identifier
        sorf_aa_seq: Amino acid sequence
        phylocsf_dataframe: DataFrame containing PhyloCSF scores
        esmfold: Dictionary containing ESMFold predictions

    Returns:
        Altair chart object
    """
    if vtx_id in phylocsf_dataframe.index:
        phylo_array = phylocsf_dataframe.loc[vtx_id, 'phylocsf_vals'][::3]
    else:
        phylo_array = [np.nan] * len(sorf_aa_seq)

    plddt = esmfold[sorf_aa_seq]['plddt']
    rows = []
    for (j_key, j_vals) in {'PhyloCSF': phylo_array, 'plDDT': plddt}.items():
        for i, i_val in enumerate(j_vals):
            rows.append((j_key, i, i_val))
    cons_altair_table = pd.DataFrame(rows, columns=['Tool', 'Position', 'value'])
    
    return alt.Chart(cons_altair_table).mark_line()\
        .encode(
            x='Position:N',
            y='value:Q',
            color='Tool:O'
        )\
        .properties(width=400, height=125)\
        .facet('Tool:O', columns=1)\
        .resolve_scale(x='independent', y='independent')
