import altair as alt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
import simplejson as json
import streamlit as st
import streamlit.components.v1 as components

from plotly import graph_objects as go
from scipy.cluster.hierarchy import linkage, leaves_list
from streamlit_echarts import JsCode
import streamlit.components.v1 as components

def color_protein_terminal_ends(aa_seq, pdb_string):
    """
    Modify pdb string to color terminal ends.
        C-term is red, and N-term is blue.
    aa_seq : str
        Amino acid sequence
    pdb_string : str
        PDB Structure as string.
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

def plot_transcripts_differential_expression_barplot(transcript_ids, de_tables_dict, de_metadata, title):
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
    
    de = pd.DataFrame([(de_tables_dict[tid]['padj']<0.00001) & (np.abs(de_tables_dict[tid]['log2FC'])>1) & ((de_tables_dict[tid]['Cancer Average'] > 1) | (de_tables_dict[tid]['GTEx Average'] > 1)) for tid in transcript_ids if len(de_tables_dict[tid])>0]).sum(axis=0)

    result = pd.concat([sum_expression_cancer, sum_expression_normal])#, '# DE Transcripts':de})
    result['TCGA'] = result.index
    result['# DE Transcripts'] = [de.loc[i] for i in result.index]
    result['GTEx Normal Tissue'] = de_metadata['GTEx']
    # Define the bar plot using plotly express
    return bar_plot_expression_groups(result, 'TCGA', ['GTEx', 'Cancer'], de_metadata, title)


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


def bar_plot_expression_groups(dataframe, group_name, group_members, de_metadata, title):
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
    
    dataframe.sort_values(by='Cancer', inplace=True)
    tcga_code_to_description = de_metadata[['Description', 'GTEx']].apply(lambda x: f"{x[0]}<br>GTEx Normal {x[1]}", axis=1).to_dict()
    option = {
      'title': {'text': title},
      'tooltip': {
          "trigger": 'axis',
          #"formatter": JsCode("function (params) {console.log(params)}").js_code,
          "formatter": JsCode("function (params) {var cols = " + json.dumps(tcga_code_to_description) + "; console.log(params); return params[0].name + ' - ' + cols[params[0].name] + '<br>' + params[0].seriesName + ': ' + params[0].value  + '<br>' + params[1].seriesName + ': ' + params[1].value;}").js_code,
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
          'dataView': { 'show': True, 'readOnly': False },
          'restore': { 'show': True },
          'saveAsImage': { 'show': True }
        }
      },
      'calculable': True,
      'yAxis': [
        {
          'name': 'Cancer Types',
          'nameLocation': 'middle',
          'nameGap': 50,
          'type': 'category',
          'data': list(dataframe[group_name].values),
          'axisLabel': { 'interval': 0}#, 'rotate': 90}
        }
      ],
      'xAxis': [
        {'name': 'Approximate TPM \n (Transcripts in bold are DE with FDR < 1e-5 and abs(LFC) > 1)',
         'nameLocation': 'middle',
         'nameGap': 30,
         'type': 'value'}
      ],
        'color': ['#237c94', '#d62418'],
      'series': [
        {
          'name': group_members[0],
          'type': 'bar',
          'data': expression_de_to_echarts_data(dataframe[group_members[0]], dataframe['DE'], 'blue')
        },
        {
          'name': group_members[1],
          'type': 'bar',
          'data': expression_de_to_echarts_data(dataframe[group_members[1]], dataframe['DE'], 'red')
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


def expression_heatmap_plot(vtx_id, vtx_id_to_transcripts, xena_expression, xena_metadata, title):
    """
    """
    # Plot transcript expression levels
    selected_transcripts_exact = vtx_id_to_transcripts.loc[vtx_id, 'transcripts_exact']
    selected_transcripts_overlapping = vtx_id_to_transcripts.loc[vtx_id, 'transcripts_overlapping']
    selected_transcripts = np.concatenate([selected_transcripts_exact, selected_transcripts_overlapping])        
    xena_overlap = xena_expression.columns.intersection(selected_transcripts)
    set2 = sns.color_palette('Set2', n_colors=2)
    
    col_colors = [set2[0] if i in selected_transcripts_exact else set2[1] for i in xena_overlap]
    
    selected_expression = xena_expression[xena_overlap]
    groups = list(map(lambda x: '-'.join(map(str, x)), xena_metadata[['_primary_site', '_study']].values))
    grouped_exp_df = selected_expression.groupby(groups).median()
    
    if grouped_exp_df.shape[1] == 0:
        return None, None
    elif grouped_exp_df.shape[1] == 1:
        plot_df = grouped_exp_df
        plot_df.sort_values(by=xena_overlap[0], inplace=True)
        plot_df = plot_df.T
    else:
        row_clusters = linkage(grouped_exp_df.T.values, method='complete', metric='euclidean')
        col_clusters = linkage(grouped_exp_df.values, method='complete', metric='euclidean')

        # compute the leaves order
        row_leaves = leaves_list(row_clusters)
        col_leaves = leaves_list(col_clusters)

        # reorder the DataFrame according to the clusters
        plot_df = grouped_exp_df.T.iloc[row_leaves, col_leaves]

    col_map = {x: i for i,x in enumerate(plot_df.columns)}
    row_map = {x: i for i,x in enumerate(plot_df.index)}
    data = [(row_map[k[0]], col_map[k[1]], v) for k,v in plot_df.stack().items()]
    
    col_names = list(col_map.keys())
    row_names = list(row_map.keys())
    
    js_col_names = "var cols = [" + ",".join([f"'{c}'" for c in col_names]) + "];"
    
    updated_row_names = []
    for r in row_names:
        if r in selected_transcripts_exact:
            updated_row_names.append(f'**{r}')
        else:
            updated_row_names.append(r)
    max_contrast = grouped_exp_df.max().max()
    if max_contrast < 5:
        max_contrast = 5
    option = {
        "title": {"text": title},
        "tooltip": {
            "formatter": JsCode("function (params) {" + js_col_names + "; return params.name + '<br>' + cols[params.data[1]] + '<br> Median TPM: ' + params.data[2];}").js_code,
        },
        "xAxis": {
            "type": "category", 
            "data": updated_row_names, 
            "axisLabel": {
                "fontSize": 10,
                "rotate": -90,
                "interval": 0,
            }
            },
        "yAxis": {
            "type": "category", 
            "data": col_names,
            "axisLabel": {
                "fontSize": 10,
                "width": 0,
                "interval": 0,
            } 
            },
        "visualMap": {
            "min": 0.25,
            "max": max_contrast,
            "calculable": True,
            "realtime": False,
            #"inRange": {
            #    "color": [
            #        '#5782bc', 
            #        '#7e9ac2', 
            #        '#a3b4cd', 
            #        '#cad0dd', 
            #        '#efeef1', 
            #        '#f7eae8', 
            #        '#e6c5c3', 
            #        '#d7a09d', 
            #        '#c87e7b', 
            #        '#b95b5a'
            #    ]
            #},
            "orient": 'vertical',
            "left": '90%',
            "top": 'center'
        },
        'toolbox': {
            'show': True,
            'feature': {
                'dataView': { 'show': True, 'readOnly': False },
                'restore': { 'show': True },
                'saveAsImage': { 'show': True }
            }
        },
        "grid": {
            "left": '30%',
            "bottom": '15%'
        },
        "series": [
            {
                "name": "Log2(TPM+1)",
                "type": "heatmap",
                "data": data,
                "borderColor": '#333',
                "borderWidth": 1,
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

    return option, events


def expression_atlas_heatmap_plot(xena_tau_df, values, xena_vtx_sum_df):

    """
    """
    tissue_specific_vtx_ids = list(xena_tau_df[xena_tau_df['tau'].between(*values)].index)
    #row_df = pd.DataFrame(row_names)

    #specific_df = row_df[row_df.apply(lambda x: True if x[0] in tissue_specific_vtx_ids else False, axis=1)]
    #specific_df.reset_index(inplace=True)
    xena_vtx_sum_df = xena_vtx_sum_df[tissue_specific_vtx_ids].copy()

    # compute the clusters
    row_clusters = linkage(xena_vtx_sum_df.T.values, method='complete', metric='euclidean')
    col_clusters = linkage(xena_vtx_sum_df.values, method='complete', metric='euclidean')

    # compute the leaves order
    row_leaves = leaves_list(row_clusters)
    col_leaves = leaves_list(col_clusters)

    # reorder the DataFrame according to the clusters
    plot_df = xena_vtx_sum_df.T.iloc[row_leaves, col_leaves]

    col_map = {x: i for i,x in enumerate(plot_df.columns)}
    row_map = {x: i for i,x in enumerate(plot_df.index)}
    data = [(row_map[k[0]], col_map[k[1]], v) for k,v in plot_df.stack().items()]
    
    col_names = list(col_map.keys())
    row_names = list(row_map.keys())

    js_col_names = "var cols = [" + ",".join([f"'{c}'" for c in col_names]) + "];"

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
            "top": 'center',
            #"inRange": {
            #    "color": [
            #        '#5782bc', 
            #        '#7e9ac2', 
            #        '#a3b4cd', 
            #        '#cad0dd', 
            #        '#efeef1', 
            #        '#f7eae8', 
            #        '#e6c5c3', 
            #        '#d7a09d', 
            #        '#c87e7b', 
            #        '#b95b5a'
            #    ]
            #},
        },
        "grid": {
            "left": '20%',
        },
        "series": [
            {
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
            }
        ],
    }
    
    events = {
        "click": "function(params) { console.log(params.name); return params.name }",
        "dblclick": "function(params) { return [params.type, params.name, params.value] }"
    }

    return option, events, tissue_specific_vtx_ids


def expression_vtx_boxplot(transcript_id, expression_df):
    """
    """
    if '**' in transcript_id:
        transcript_id = transcript_id[2:]
    groups = expression_df[['_primary_site', '_study']].apply(lambda x: '-'.join(x), axis = 1)
    df = expression_df[[transcript_id]].copy()
    df['Sample'] = groups
    df.reset_index(inplace=True)

    medians = df.groupby('Sample')[transcript_id].median().sort_values()

    order_map = {name:i for i, name in enumerate(medians.index)}
    df['order'] = df.apply(lambda x: order_map[x['Sample']], axis=1)
    df.sort_values(by='order', ascending=False, inplace=True)

    fig = px.box(df, x="Sample", 
                 y=transcript_id,  
                 points='outliers', hover_data='index',
                 category_orders={''}, height=500)
    
    fig.update_traces(quartilemethod="exclusive") # or "inclusive", or "linear" by default

    return fig


def format_protein_feature_strings_for_altair_heatmap(orf_features):
    """
    """
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
    """
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
        color=alt.Color('Predicted Class:O', scale=color_scale, legend=alt.Legend(
        orient='none',
        legendX=45, legendY=150,
        direction='horizontal',
        titleAnchor='middle')),
        tooltip = ['Position', 'Predicted Class'])
    return fig+base.mark_text(baseline='middle').encode(alt.Text('aa:O'))


def plot_sequence_line_plots_altair(vtx_id, sorf_aa_seq, phylocsf_dataframe, kibby, esmfold):
    # try:
    if vtx_id in phylocsf_dataframe.index:
        phylo_array = phylocsf_dataframe.loc[vtx_id, 'phylocsf_vals'][::3]
    else:
        phylo_array = [np.nan]*len(sorf_aa_seq)
    kib = [i*100 for i in kibby.loc[vtx_id, 'conservation']]
    plddt = esmfold[sorf_aa_seq]['plddt']
    rows = []
    for (j_key, j_vals) in {'PhyloCSF': phylo_array, 'plDDT': plddt, 'Kibby Conservation': kib}.items():
        for i, i_val in enumerate(j_vals):
            rows.append((j_key, i, i_val))
    cons_altair_table = pd.DataFrame(rows, columns = ['Tool', 'Position', 'value'])
    return alt.Chart(cons_altair_table).mark_line()\
                                       .encode(x='Position:N', y='value:Q', color='Tool:O')\
                                       .properties(width=400, height=125)\
                                       .facet('Tool:O', columns=1)\
                                       .resolve_scale(x='independent', y='independent')


