import pandas as pd
from streamlit_echarts import st_echarts, JsCode
import streamlit as st
import altair as alt
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
from plotly import graph_objects as go


def plot_transcripts_differential_expression_barplot(transcript_ids, de_tables_dict, title):
    """
    Parameters
    ----------
    transcript_ids : List[string]
        The transcript IDs to plot differential expression for. (Assumes they exist in data).
    de_tables_dict : Dict
        Dictionary of transcript_id keys and dataframe values. Dataframes contain tables of DE 
          results for different cancer types.
          
    """
    sum_expression_cancer = pd.DataFrame(pd.DataFrame([de_tables_dict[tid]['Cancer Average'] for tid in transcript_ids if len(de_tables_dict[tid])>0]).sum(axis=0), columns = ['Sum'])
    sum_expression_cancer['condition'] = 'Cancer'
    sum_expression_normal = pd.DataFrame(pd.DataFrame([de_tables_dict[tid]['GTEx Average'] for tid in  transcript_ids if len(de_tables_dict[tid])>0]).sum(axis=0), columns = ['Sum'])
    sum_expression_normal['condition'] = 'GTEx'
    de = pd.DataFrame([de_tables_dict[tid]['padj']<0.00001 for tid in transcript_ids if len(de_tables_dict[tid])>0]).sum(axis=0)
    result = pd.concat([sum_expression_cancer, sum_expression_normal])#, '# DE Transcripts':de})
    result['TCGA'] = result.index
    result['# DE Transcripts'] = [de.loc[i] for i in result.index]
    # Define the bar plot using plotly express
    return bar_plot_expression_groups(result, 'TCGA', ['GTEx', 'Cancer'], title)


def plot_vtx_transcripts_heatmap(vtx_id, vtx_id_to_transcripts, xena_expression, xena_metadata):
    selected_transcripts_exact = vtx_id_to_transcripts.loc[vtx_id, 'transcripts_exact']
    selected_transcripts_overlapping = vtx_id_to_transcripts.loc[vtx_id, 'transcripts_overlapping']
    selected_transcripts = np.concatenate([selected_transcripts_exact, selected_transcripts_overlapping])        
    xena_overlap = xena_expression.columns.intersection(selected_transcripts)
    set2 = sns.color_palette('Set2', n_colors=2)
    if len(xena_overlap)>1:
        col_colors = [set2[0] if i in selected_transcripts_exact else set2[1] for i in xena_overlap]
        selected_expression = xena_expression[xena_overlap]
        groups = list(map(lambda x: '-'.join(map(str, x)), xena_metadata[['_primary_site', '_study']].values))
        heatmap_figure = sns.clustermap(selected_expression.groupby(groups).median(), col_colors=col_colors,
                                cmap='coolwarm', cbar_kws={'label': 'Log(TPM+0.001)'}, center=1, vmin=-3, vmax=6)
    elif len(xena_overlap) == 1:
        # st.write('Only 1 transcript')
        col_colors = [set2[0] if i in selected_transcripts_exact else set2[1] for i in xena_overlap]
        selected_expression = xena_expression[list(xena_overlap)]
        # selected_expression['Null'] = -1
        # col_colors.append(set2[1])
        print(selected_expression.shape, col_colors)
        groups = list(map(lambda x: '-'.join(map(str, x)), xena_metadata[['_primary_site', '_study']].values))
        heatmap_figure, ax = plt.subplots()
        sns.heatmap(selected_expression.groupby(groups).median(), ax=ax,
                                cmap='coolwarm', cbar_kws={'label': 'Log(TPM+0.001)'}, center=1, vmin=-3, vmax=6)
    else:
        return 'No transcripts in TCGA/GTEx/TARGET found containing this sORF'
    return heatmap_figure


def plot_structure_plddt(plddt, ax):
    """
    """
    fs=8
    n_aa = len(plddt)
    ax.plot(range(n_aa), plddt)
    ax.set_title("ESMFold Prediction Confidence", fontsize=fs)
    ax.set_xlabel("Amino Acid", fontsize=fs)
    ax.set_ylabel("plDDT", fontsize=fs)
    ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), fontsize=fs)
    ax.set_yticks(ax.get_yticks(), ax.get_yticklabels(), fontsize=fs)
    ax.set_xlim(0, n_aa)
    ax.set_ylim([0, 100])
    ax.hlines(50, 0, n_aa, linestyle='--', color='r', alpha=0.5, label='Low Confidence')
    ax.hlines(90, 0, n_aa, linestyle='--', color='g', alpha=0.5, label='High Confidence')
    return ax


def expression_de_to_echarts_data(expression_values, de_values, color):
    """
    """
    data = []
    color_map = {'red_de': '#d62418',
                 'red_null': '#d1918c',
                 'blue_de': '#237c94',
                 'blue_null': '#70aaba'
                }
    # colors = ['black', 'red']
    for e, d in zip(expression_values, de_values):
        if d:
            c = color_map[f"{color}_de"]
        else:
            c = color_map[f"{color}_null"]
        data.append({'value': e, 'itemStyle': {'color': c}})
    return data


def bar_plot_expression_groups(dataframe, group_name, group_members, title):
    if dataframe.shape[0]>0:
        dataframe = dataframe.copy()
        new_df = pd.DataFrame()
        for g, subdf in dataframe.groupby('condition'):
            new_df[g] = subdf['Sum']
        new_df[group_name] = subdf[group_name]
        new_df['DE'] = subdf['# DE Transcripts']
        dataframe = new_df.copy()
    else:
        return None
    option = {
      'title': {'text': title},
      'tooltip': {'trigger': 'axis'},
      'legend': {'data': group_members},
      'toolbox': {
        'show': True,
        'feature': {
          'dataView': { 'show': True, 'readOnly': False },
          # 'magicType': { show: true, type: ['line', 'bar'] },
          'restore': { 'show': True },
          'saveAsImage': { 'show': True }
        }
      },
      'calculable': True,
      'xAxis': [
        {
          'type': 'category',
          'data': list(dataframe[group_name].values),
          'axisLabel': { 'interval': 0, 'rotate': 90}
        }
      ],
      'yAxis': [
        {'type': 'value'}
      ],
        'color': ['#237c94', '#d62418'],
      'series': [
        {
          'name': group_members[0],
          'type': 'bar',
          'data': expression_de_to_echarts_data(dataframe[group_members[0]], dataframe['DE'], 'blue')
          # markPoint: {
          #   data: [
          #     { type: 'max', name: 'Max' },
          #     { type: 'min', name: 'Min' }
          #   ]
          # },
        },
        {
          'name': group_members[1],
          'type': 'bar',
          'data': expression_de_to_echarts_data(dataframe[group_members[1]], dataframe['DE'], 'red')
          # markPoint: {
          #   data: [
          #     { type: 'max', name: 'Max' },
          #     { type: 'min', name: 'Min' }
          #   ]
          # },
        },
      ]
    }

    return option


def expression_heatmap_plot(vtx_id, vtx_id_to_transcripts, xena_expression, xena_metadata):
    """
    """
    # Plot transcript expression levels
    selected_transcripts_exact = vtx_id_to_transcripts.loc[vtx_id, 'transcripts_exact']
    selected_transcripts_overlapping = vtx_id_to_transcripts.loc[vtx_id, 'transcripts_overlapping']
    selected_transcripts = np.concatenate([selected_transcripts_exact, selected_transcripts_overlapping])        
    xena_overlap = xena_expression.columns.intersection(selected_transcripts)
    set2 = sns.color_palette('Set2', n_colors=2)
    if len(xena_overlap)>1:
        col_colors = [set2[0] if i in selected_transcripts_exact else set2[1] for i in xena_overlap]
        selected_expression = xena_expression[xena_overlap]
        groups = list(map(lambda x: '-'.join(map(str, x)), xena_metadata[['_primary_site', '_study']].values))
        fig = sns.clustermap(selected_expression.groupby(groups).median(), col_colors=col_colors,
                                  cmap='coolwarm', cbar_kws={'label': 'Log(TPM+0.001)'}, center=1, vmin=-3, vmax=6)
        
    elif len(xena_overlap) == 1:
        st.write('Only 1 transcript')
        col_colors = [set2[0] if i in selected_transcripts_exact else set2[1] for i in xena_overlap]
        selected_expression = xena_expression[list(xena_overlap)]
        print(selected_expression.shape, col_colors)
        groups = list(map(lambda x: '-'.join(map(str, x)), xena_metadata[['_primary_site', '_study']].values))
        fig, ax = plt.subplots()
        sns.heatmap(selected_expression.groupby(groups).median(), ax=ax,
                                  cmap='coolwarm', cbar_kws={'label': 'Log(TPM+0.001)'}, center=1, vmin=-3, vmax=6)
    
    else:
        fig = None

    return fig


def expression_de_plot(vtx_id, vtx_id_to_transcripts, de_tables_dict):
    """
    """
    tids = vtx_id_to_transcripts.loc[vtx_id, 'transcripts_exact']
    sum_expression_cancer = pd.DataFrame(pd.DataFrame([de_tables_dict[tid]['Cancer Average'] for tid in tids if len(de_tables_dict[tid])>0]).sum(axis=0), columns = ['Sum'])
    sum_expression_cancer['condition'] = 'Cancer'
    sum_expression_normal = pd.DataFrame(pd.DataFrame([de_tables_dict[tid]['GTEx Average'] for tid in  tids if len(de_tables_dict[tid])>0]).sum(axis=0), columns = ['Sum'])
    sum_expression_normal['condition'] = 'GTEx'
    de = pd.DataFrame([de_tables_dict[tid]['padj']<0.001 for tid in tids if len(de_tables_dict[tid])>0]).sum(axis=0)
    result = pd.concat([sum_expression_cancer, sum_expression_normal])#, '# DE Transcripts':de})
    result['TCGA'] = result.index
    result['# DE Transcripts'] = [de.loc[i] for i in result.index]

    # Define the bar plot using Altair
    fig = px.bar(result, y='TCGA', x='Sum', color='condition', barmode='group')
    fig.add_trace(go.Scatter(
        y=result['TCGA'],
        x=result['Sum'],
        mode="text",
        name="# DE Transcripts",
        text=result['# DE Transcripts'],
        textposition="top left",
        showlegend=False
    ))

    return fig, result


def expression_atlas_heatmap_plot(xena_tau_df, data, col_names, row_names, values):
    """
    """
    tissue_specific_vtx_ids = list(xena_tau_df[xena_tau_df['tau'].between(*values)].index)
    row_df = pd.DataFrame(row_names)

    specific_df = row_df[row_df.apply(lambda x: True if x[0] in tissue_specific_vtx_ids else False, axis=1)]
    specific_df.reset_index(inplace=True)
    subset_data = []

    for d in data:
        if d[0] in specific_df['index'].values:
            new_idx = int(specific_df[specific_df['index'] == d[0]].index[0])
            subset_data.append((new_idx, d[1], d[2]))
    
    row_names = list(specific_df[0])

    js_col_names = "var cols = [" + ",".join([f"'{c}'" for c in col_names]) + "];"
    #st.write(js_col_names)

    option = {
        "tooltip": {
            "formatter": JsCode("function (params) {" + js_col_names + "; return params.name + '<br>' + cols[params.data[1]] + '<br> Log2(TPM+1): ' + params.data[2];}").js_code,
        },
        "xAxis": {
            "type": "category", 
            "data": row_names, 
            "axisLabel": {
                "fontSize": 10,
                "rotate": -90,
                "interval": 1,
            }
            },
        "yAxis": {
            "type": "category", 
            "data": col_names,
            "axisLabel": {
                "fontSize": 11,
                "width": 0,
                "interval": 0,
            } 
            },
        "visualMap": {
            "min": 0,
            "max": 12,
            "calculable": True,
            "realtime": False,
            "orient": 'vertical',
            "left": '95%',
            "top": 'center'
        },
        "grid": {
            "left": '20%',
        },
        "series": [
            {
                "name": "Log2(TPM+1)",
                "type": "heatmap",
                "data": subset_data,
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
            }
        ],
    }
    
    events = {
        "click": "function(params) { console.log(params.name); return params.name }",
        "dblclick": "function(params) { return [params.type, params.name, params.value] }"
    }

    return option, events, tissue_specific_vtx_ids


def expression_vtx_boxplot(vtx_id, expression_df):
    """
    """

    df = expression_df[['primary disease or tissue', vtx_id]].copy()
    df.reset_index(inplace=True)

    medians = df.groupby('primary disease or tissue')[vtx_id].median().sort_values()

    order_map = {name:i for i, name in enumerate(medians.index)}
    df['order'] = df.apply(lambda x: order_map[x['primary disease or tissue']], axis=1)
    df.sort_values(by='order', ascending=False, inplace=True)

    fig = px.box(df, x="primary disease or tissue", 
                 y=vtx_id,  
                 points='outliers', hover_data='index',
                 category_orders={''}, height=500)
    
    fig.update_traces(quartilemethod="exclusive") # or "inclusive", or "linear" by default

    return fig


def format_protein_feature_strings_for_altair_heatmap(orf_features):
    # orf_features = orf_string_table.T[vtx_id]
    orf_features = orf_features.T
    orf_features = orf_features.apply(lambda x: [*x])
    seq = orf_features.pop('Sequence')
    values = []
    for i, (ix, row) in enumerate(orf_features.items()):
        for j, v in enumerate(row):
            values.append([j, ix, v, seq[j]])
    df = pd.DataFrame(values, columns = ['Position', 'Tool', 'Predicted Class', 'aa'])
    return df


def altair_protein_features_plot(df):
    color_dict = {
        'S': 'palegoldenrod',
        'C': 'lightgray',
        'I': 'lightgray',
        'O': 'lightblue'
    }
    
    color_scale = alt.Scale(domain=list(color_dict.keys()), range=list(color_dict.values()))

    base = alt.Chart(df).encode(
        alt.X('Position:O'),
        alt.Y('Tool:O')
    )
    fig = base.mark_rect().encode(
        color=alt.Color('Predicted Class:O', scale=color_scale, legend=alt.Legend(
        orient='none',
        legendX=130, legendY=-60,
        direction='horizontal',
        titleAnchor='middle')),
        tooltip = ['Position'])
    return fig+base.mark_text(baseline='middle').encode(alt.Text('aa:O'))

def plot_sequence_line_plots_altair(vtx_id, sorf_aa_seq, phylocsf_dataframe, kibby, esmfold):
    phylo_array = phylocsf_dataframe.loc[vtx_id, 'phylocsf_vals'][::3]
    kib = [i*100 for i in kibby.loc[vtx_id, 'conservation']]
    plddt = esmfold[sorf_aa_seq]['plddt']
    rows = []
    for (j_key, j_vals) in {'PhyloCSF': phylo_array, 'plDDT': plddt, 'Kibby Conservation': kib}.items():
        for i, i_val in enumerate(j_vals):
            rows.append((j_key, i, i_val))
    cons_altair_table = pd.DataFrame(rows, columns = ['Tool', 'Position', 'value'])
    return alt.Chart(cons_altair_table).mark_line().encode(x='Position:N',
                                                           y='value:Q', color='Tool:O').properties( width=400, height=125
                                                                                                            ).facet('Tool:O', columns=1)

