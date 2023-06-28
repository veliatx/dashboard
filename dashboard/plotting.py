import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from streamlit_echarts import st_echarts
import streamlit as st
import plotly.express as px
from plotly import graph_objects as go


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


def bar_plot_expression_groups(dataframe, group_name, group_members):
    """
    """
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
      'title': {'text': 'Expression Data'},
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
