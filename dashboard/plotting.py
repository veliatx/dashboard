import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from streamlit_echarts import st_echarts, JsCode
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
            #"containLabel": True,
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
    st.write(expression_df.head())
    option = {
      "title": [
          {"text": f"GTEx Expression {vtx_id} ", "left": "center"},
          {
              "text": "",
              "borderColor": "#999",
              "borderWidth": 1,
              "textStyle": {"fontWeight": "normal", "fontSize": 14, "lineHeight": 20},
              "left": "10%",
              "top": "90%",
          },
      ],
      "dataset": [
          {
              "source": [
                  [
                      850,
                      740,
                      900,
                      1070,
                      930,
                      850,
                      950,
                      980,
                      980,
                      880,
                      1000,
                      980,
                      930,
                      650,
                      760,
                      810,
                      1000,
                      1000,
                      960,
                      960,
                  ],
                  [
                      960,
                      940,
                      960,
                      940,
                      880,
                      800,
                      850,
                      880,
                      900,
                      840,
                      830,
                      790,
                      810,
                      880,
                      880,
                      830,
                      800,
                      790,
                      760,
                      800,
                  ],
                  [
                      880,
                      880,
                      880,
                      860,
                      720,
                      720,
                      620,
                      860,
                      970,
                      950,
                      880,
                      910,
                      850,
                      870,
                      840,
                      840,
                      850,
                      840,
                      840,
                      840,
                  ],
                  [
                      890,
                      810,
                      810,
                      820,
                      800,
                      770,
                      760,
                      740,
                      750,
                      760,
                      910,
                      920,
                      890,
                      860,
                      880,
                      720,
                      840,
                      850,
                      850,
                      780,
                  ],
                  [
                      890,
                      840,
                      780,
                      810,
                      760,
                      810,
                      790,
                      810,
                      820,
                      850,
                      870,
                      870,
                      810,
                      740,
                      810,
                      940,
                      950,
                      800,
                      810,
                      870,
                  ],
              ]
          },
          {
              "transform": {
                  "type": "boxplot",
                  "config": {"itemNameFormatter": "expr {value}"},
              }
          },
          {"fromDatasetIndex": 1, "fromTransformResult": 1},
      ],
      "tooltip": {"trigger": "item", "axisPointer": {"type": "shadow"}},
      "grid": {"left": "10%", "right": "10%", "bottom": "15%"},
      "xAxis": {
          "type": "category",
          "boundaryGap": True,
          "nameGap": 30,
          "splitArea": {"show": False},
          "splitLine": {"show": False},
      },
      "yAxis": {
          "type": "value",
          "name": "TPM",
          "splitArea": {"show": True},
      },
      "series": [
          {"name": "boxplot", "type": "boxplot", "datasetIndex": 1},
          {"name": "outlier", "type": "scatter", "datasetIndex": 2},
      ],
    }
    return option
    st_echarts(option, height="400px")


def expression_vtx_boxplot2(vtx_id, expression_df):
    """
    """

    df = expression_df[['primary disease or tissue', vtx_id]]
    df.reset_index(inplace=True)
    fig = px.box(df, x="primary disease or tissue", y=vtx_id, labels={'primary disease or tissue': 'GTEx Primary Tissue'}, points='all', hover_data='index')
    fig.update_traces(quartilemethod="exclusive") # or "inclusive", or "linear" by default

    return fig


def features_to_int_arrays(row):
    char_to_int = {'S': 1, 'C': 0, 'O': 2, 'M': 3, 'I': 0}
    convert_to_integer_array = lambda x: np.array([char_to_int[i] for i in x])
    textrow = row.apply(lambda x: [i for i in x]).copy()
    row = row.apply(convert_to_integer_array)
    imdf = pd.DataFrame(index = [i for i in row.name])
    textdf = pd.DataFrame(index = [i for i in row.name])
    for ix, row in row.items():
        imdf[ix] = row
        textdf[ix] = textrow.loc[ix]
    imdf = imdf.T
    textdf = textdf.T
    imdf.loc['aa'] = imdf.columns
    imdf.columns = np.arange(imdf.shape[1])
    return imdf, textdf
