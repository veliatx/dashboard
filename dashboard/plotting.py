import matplotlib.pyplot as plt
import pandas as pd
from streamlit_echarts import st_echarts
import streamlit as st

def plot_structure_plddt(plddt, ax):
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
    st_echarts(option, height="400px")