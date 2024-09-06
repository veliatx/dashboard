from typing import List, Any, Dict, Tuple
from PIL import Image
import random

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from streamlit_plotly_events import plotly_events
from streamlit_echarts import st_echarts
from scipy.cluster.hierarchy import linkage, leaves_list
import plotly.graph_objects as go

from dashboard.etl import DB_CONNECTION_STRING
from dashboard.util import filter_dataframe_dynamic
from dashboard import plotting
from dashboard.data_load import *

from expression_atlas_db import queries, base


blood_tissues = [
    # Blood tissues don't have a great location in body. Display on bottom left.
    'NEUTROPHIL', 'BLOOD', 'PBMC', 'DENDRITIC_CELL', 'MONOCYTES', 'MONOCYTE',
    'HEMATOPOIETIC_STEM_CELL', 'GRANULOCYTIC_PRECURSOR', 'MONOCYTIC_PRECURSOR',
    'ERYTHROID_PRECURSOR', 'DENDRITIC_CELLS', 
]
tissue_positions = {
    # Each tissue needs to be manually filled in. If not, it will be displayed
    # on bottom right. 
    'OLIGODENDROCYTE': (0.5, 0.86),
    'BRAIN': (0.5, 0.85),
    'LUNG': (0.48, 0.7),
    'LIVER': (0.5, 0.62),
    'SIGMOID_COLON': (0.52, 0.52),
    'ASCENDING_COLON': (0.52, 0.55),
    'THYROID': (0.5, 0.79),
    'SKELETAL_MUSCLE': (0.48, 0.25),
    'VISCERAL_ADIPOSE_TISSUE': (0.48, 0.58),
    'HEART': (0.52, 0.73),
    'HEART_MUSCLE': (0.52, 0.72),
    'DUODENUM': (0.5, 0.63),
    'KIDNEY': (0.48, 0.59),
    'JEJUNUM': (0.5, 0.64),
    'STOMACH_MUSCULARIS': (0.5, 0.65),
    'RECTUM': (0.5, 0.51),
    'ILEUM': (0.5, 0.61),
    'BONE': (0.475, 0.42),
    'SYNOVIAL_TISSUE': (0.475, 0.32),
    'CARTILAGE': (0.475, 0.33),
    'SALIVARY_GLAND': (0.49, 0.82),
    'ADIPOCYTE': (0.46, 0.54),
    'SKIN': (0.565, 0.49),
    'PANCREATIC_ISLET': (0.52, 0.6),
    'COLON': (0.52, 0.53),
    'ATRIAL_FIBROBLAST': (0.52, 0.73),
}


def update_queue(
    session: base._Session,
    update_df: pd.DataFrame,
    update_columns: List[str],
) -> None:
    """Callback update queue."""
    edit_df = queries.update_studyqueue(
        session,
        update_df,
        update_columns=update_columns,
    )
    st.session_state.pop('editor_queue_df')
    st.session_state['edited_result_df'] = edit_df

def submit_queue(
    Session: Any,
    placeholder: st.delta_generator.DeltaGenerator,
) -> None:
    """Callback submit queue."""
    srp_id = st.session_state['srp_id']
    geo_id = st.session_state['geo_id']
    category = st.session_state['category']
    disease = st.session_state['disease']
    tissue = st.session_state['tissue']
    contrast = st.session_state['contrast']
    requestor = st.session_state['requestor']
    technology = st.session_state['technology']
    comments = st.session_state['comments']
    priority = st.session_state['priority']

    if not srp_id.startswith(('SRP', 'ERP', 'DRP')) or len(srp_id) != 9:
        st.session_state['submit_results'] = (
            False,
            'SRP ID must start with (SRP, ERP, or DRP) and be 9 characters.',
        )
    elif geo_id != '' and not geo_id.startswith('GSE'):
        st.session_state['submit_results'] = (
            False,
            'GEO ID misformatted, must start with GSE.',
        )
    elif len(category) == 0 or category.startswith('ex:'):
        st.session_state['submit_results'] = (
            False,
            'Must submit valid request.',
        )
    elif len(disease) == 0 or disease.startswith('ex:'):
        st.session_state['submit_results'] = (
            False,
            'Must submit valid disease.',
        )
    elif len(tissue) == 0 or tissue.startswith('ex:'):
        st.session_state['submit_results'] = (
            False,
            'Must submit valid tissue.',
        )
    elif len(contrast) == 0 or contrast.startswith('ex:'):
        st.session_state['submit_results'] = (
            False,
            'Must submit valid contrast.',
        )
    elif len(requestor) == 0 or requestor.startswith('ex:'):
        st.session_state['submit_results'] = (
            False,
            'Must submit valid requestor.',
        )
    elif not technology:
        st.session_state['submit_results'] = (
            False,
            'Must submit valid technology.',
        )
    elif comments.startswith('ex:'):
        comments = ''
    else:
        with placeholder, Session.begin() as session:
            with st.spinner('Thinking...'):
                results = queries.submit_studyqueue(
                    session,
                    srp_id,
                    category,
                    technology.upper(),
                    disease,
                    tissue,
                    contrast,
                    requestor,
                    priority,
                    comments,
                    geo_id=geo_id if len(geo_id) > 0 else None,
                )
        st.session_state['submit_results'] = results

def set_random_tissue_positions(
    tissue_positions: Dict[str, Tuple[float, float]] = tissue_positions,
) -> Dict[str, Tuple[float, float]]:
    """Randomly jitter tissues positions around x position."""
    random.seed(42)
    set_tissue_positions = tissue_positions.copy()
    for t in tissue_positions.keys():
        x, y = tissue_positions[t]
        set_tissue_positions[t] = (x+random.uniform(-0.01,0.01), y)
    return set_tissue_positions

def plot_human_body_summary(
    summary_level: str,
    sizes: pd.Series,
    summary_text: pd.Series = pd.Series(),
    h: float = 800,
    w: float = 1000,
    tissue_positions: Dict[str, Tuple[float, float]] = set_random_tissue_positions(),
    blood_tissues: List[str] = blood_tissues,
    scale_size: float = 5, 
    background_source: str = 'human_body.jpg',
) -> go.Figure:
    """
    Plot human body with tissues, diseases, contrasts overlayed on figure. 
    """
    extra_tissues = [
        t for t in sizes.index \
        if t not in blood_tissues and t not in tissue_positions.keys()
    ]
    sorted_tissues_positions = dict(
        sorted(
            tissue_positions.items(), 
            key=lambda x: x[1][1], 
            reverse=False,
        )
    )

    # Split tissues by left or right positioning, tissues without a label outside
    # from blood positioned on bottom left. 
    left_tissues = {k:v for k,v in sorted_tissues_positions.items() if v[0] < 0.5}
    right_tissues = {k:v for k,v in sorted_tissues_positions.items() if v[0] >= 0.5}

    left_y = np.linspace(0.3, 0.95, len(left_tissues))
    right_y = np.linspace(0.3, 0.95, len(right_tissues))
    blood_y1 = np.linspace(0.05, 0.1, len(blood_tissues))
    blood_y2 = np.linspace(0.05, 0.275, len(blood_tissues))
    extra_y1 = np.linspace(0.05, 0.1, len(extra_tissues))
    extra_y2 = np.linspace(0.05, 0.25, len(extra_tissues))

    fig = go.Figure()
    fig.add_layout_image(
        {
            'source': Image.open(background_source),
            'x': 0,
            'y': 1,
            'sizex':1,
            'sizey':1,
            'sizing':'stretch',
            'layer': 'below',
            'opacity': 0.5,
        }
    )
    fig.update_layout(
        {
            'title': f'Expression Atlas {summary_level.capitalize()}',
            'xaxis_range': [0,1],
            'yaxis_range': [0,1],
            'xaxis_showgrid': False,
            'yaxis_showgrid': False,
            'xaxis_fixedrange': True,
            'yaxis_fixedrange': True,
            'xaxis_visible': False,
            'yaxis_visible': False,
            'width':w,
            'height':h,
            'template': 'plotly_white',
        },
    )

    for y, (t, (p_x, p_y)) in zip(left_y, left_tissues.items()):
        circle_size = sizes.loc[t]*scale_size if t in sizes.index else 1
        fig.add_trace(
            go.Scatter(
                x=[p_x, 0.25],
                y=[p_y, y],
                hovertext=[
                    '', 
                    f'{sizes.get(t, 0)} {summary_level}<br />'+
                    '<br />'.join(summary_text.get(t, [])),
                ],
                mode="lines+text+markers",
                text=['', t],
                textposition='middle left',
                line={'color': 'grey'},
                marker={
                    'color': ['white', 'grey'],
                    'size': [3, circle_size],

                },
                showlegend=False,
            )
        )
    for y, (t, (p_x, p_y)) in zip(right_y, right_tissues.items()):
        circle_size = sizes.loc[t]*scale_size if t in sizes.index else 1
        fig.add_trace(
            go.Scatter(
                x=[p_x, 0.75],
                y=[p_y, y],
                hovertext=[
                    '', 
                    f'{sizes.get(t, 0)} {summary_level}<br />'+
                    '<br />'.join(summary_text.get(t, [])),
                ],
                mode="lines+text+markers",
                text=['', t],
                textposition='middle right',
                line={'color': 'grey'},
                marker={
                    'color': ['white', 'grey'],
                    'size': [3, circle_size],

                },
                showlegend=False,
            )
        )
    for y1, y2, t in zip(blood_y1, blood_y2, blood_tissues):
        circle_size = sizes.loc[t]*scale_size if t in sizes.index else 1
        fig.add_trace(
            go.Scatter(
                x=[0.48, 0.25],
                y=[y1, y2],
                hovertext=[
                    '', 
                    f'{sizes.get(t, 0)} {summary_level}<br />'+
                    '<br />'.join(summary_text.get(t, [])),
                ],
                mode="lines+text+markers",
                text=['', t],
                textposition='middle left',
                line={'color': 'grey'},
                marker={
                    'color': ['white', 'grey'],
                    'size': [3, circle_size],

                },
                showlegend=False,
            )
        )
    for y1, y2, t in zip(extra_y1, extra_y2, extra_tissues):
        circle_size = sizes.loc[t]*scale_size if t in sizes.index else 1
        fig.add_trace(
            go.Scatter(
                x=[0.52, 0.75],
                y=[y1, y2],
                hovertext=[
                    '', 
                    f'{sizes.get(t, 0)} {summary_level}<br />'+
                    '<br />'.join(summary_text.get(t, [])),
                ],
                mode="lines+text+markers",
                text=['', t],
                textposition='middle right',
                line={'color': 'grey'},
                marker={
                    'color': ['white', 'grey'],
                    'size': [3, circle_size],

                },
                showlegend=False,
            )
        )
    return fig

def plot_bar_summary(
    counts: pd.Series,
    summary_level: str,
    counts_information: pd.Series = pd.Series(),
    counts_information_level: str = 'Meta',
    w: int = 1000,
    h: int = 1000,
) -> go.Figure:
    """Plot bar plot with information on number of diseases/contrasts etc. 
    """
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            y=counts.index,
            x=counts.values,
            orientation='h',
            hovertext=[
                f'{counts_information_level.capitalize()}:<br />'+'<br / >'.join(
                    counts_information.get(s, [])
                ) for s in counts.index
        ],

        )
    )
    fig.update_layout(
        {
            'yaxis_autorange': 'reversed',
            'title': f'Expression Atlas {summary_level.capitalize()}',
            'title_x': 0.5,
            'xaxis_showgrid': False,
            'yaxis_showgrid': False,
            'xaxis_fixedrange': True,
            'yaxis_fixedrange': True,
            'yaxis_tickmode': 'linear',
            'width':w,
            'height':h,
            'template': 'plotly_white',
        }
    )
    return fig    


def expression_atlas_summary() -> None:
    """Landing page for expression atlas.
    """
    Session = base.configure(DB_CONNECTION_STRING)
    queue_display_cols = [
            "id", "created_at", "requestor", "priority", "srp_id", "geo_id", "pmid", "tissue", 
            "disease", "category", "contrast", "status", "quality", "technology", "comments",
            "data_released", "title", "description",
        ]
    tab1, tab2, tab3 = st.tabs(['Expression Atlas Summary', 'Queue', 'Request Dataset'])
    with tab1:
        with Session.begin() as session:
            if 'summary_dfs' not in st.session_state.keys():
                study_df, contrast_df, sample_df = queries.build_expression_atlas_summary_dfs(session)
                st.session_state['summary_dfs'] = (study_df, contrast_df, sample_df,)
            else:
                study_df, contrast_df, sample_df = st.session_state['summary_dfs']
        st.subheader(f"Studies: {study_df['velia_id'].nunique()}")
        st.divider()
        human_study_col, bar_study_col = st.columns(2)
        study_tissue_vals = study_df.loc[
            :,['velia_id','primary_tissue']
        ].drop_duplicates()['primary_tissue'].value_counts()
        study_tissue_diseases = study_df.loc[
            :,['velia_id','primary_tissue', 'primary_disease']
        ].drop_duplicates().groupby('primary_tissue')['primary_disease'].agg(list)
        
        study_disease_counts = study_df.loc[
            :,['velia_id', 'primary_disease']
        ].drop_duplicates()['primary_disease'].value_counts()
        study_disease_names = study_df.loc[
            :,['velia_id', 'primary_disease']
        ].drop_duplicates().groupby('primary_disease')['velia_id'].agg(list)
        
        with human_study_col:
            fig_study_body = plot_human_body_summary(
                'studies', 
                study_tissue_vals, 
                summary_text=study_tissue_diseases,
                background_source='assets/human_body.jpg',
                h=700,
            )
            st.plotly_chart(fig_study_body, use_container_width=True)
            st.dataframe(study_df, hide_index=True)
        with bar_study_col:
            fig_study_bar = plot_bar_summary(
                study_disease_counts,
                'studies',
                counts_information=study_disease_names,
                counts_information_level='studies',
                h=700
            )
            st.plotly_chart(fig_study_bar, use_container_width=True)
            st.dataframe(study_df.loc[:,['primary_tissue','primary_disease']].value_counts())
        
        st.subheader(f"Contrasts: {contrast_df['contrast_name'].size}")
        st.divider()

        human_contrast_col, bar_contrast_col = st.columns(2)
        contrast_tissue_vals = contrast_df.loc[
            :,['velia_id','primary_tissue', 'contrast_name']
        ].drop_duplicates()['primary_tissue'].value_counts()
        contrast_tissue_diseases = contrast_df.loc[
            :,['velia_id','primary_tissue', 'contrast_name']
        ].drop_duplicates().groupby('primary_tissue')['contrast_name'].agg(list)
        
        contrast_study_counts = contrast_df.loc[
            :,['velia_id', 'contrast_name']
        ].drop_duplicates()['contrast_name'].value_counts()
        contrast_study_names = contrast_df.loc[
            :,['velia_id', 'contrast_name']
        ].drop_duplicates().groupby('contrast_name')['velia_id'].agg(list)

        with human_contrast_col:
            fig_contrast_body = plot_human_body_summary(
                'contrasts', 
                contrast_tissue_vals, 
                summary_text=contrast_tissue_diseases,
                background_source='assets/human_body.jpg',
                h=700,
                scale_size=2.5,
            )
            st.plotly_chart(fig_contrast_body, use_container_width=True)
            st.dataframe(contrast_df, hide_index=True)
        with bar_contrast_col:
            fig_contrast_bar = plot_bar_summary(
                contrast_study_counts.head(50),
                'contrasts top 50',
                counts_information=contrast_study_names,
                counts_information_level='studies',
                h=700
            )
            st.plotly_chart(fig_contrast_bar, use_container_width=True)

        st.subheader(f"Samples: {sample_df['srx_id'].nunique()}")
        st.divider()

        human_sample_col, df_sample_col = st.columns(2)
        sample_tissue_vals = sample_df.loc[
            :,['srx_id','primary_tissue']
        ].drop_duplicates()['primary_tissue'].value_counts()
        
        with human_sample_col:
            fig_sample_body = plot_human_body_summary(
                'samples', 
                sample_tissue_vals, 
                background_source='assets/human_body.jpg',
                h=700,
                scale_size=0.1,
            )
            st.plotly_chart(fig_sample_body, use_container_width=True)
        with df_sample_col:
            st.dataframe(sample_df, hide_index=True)

    with tab2:
        not_editable_columns = [
            "id",
            "created_at",
            "requestor",
            "srp_id",
            "geo_id",
        ]
        editable_columns = [
            "tissue",
            "disease",
            "contrast",
            "category",
            "status",
            "technology",
            "quality",
            "comments",
            "priority",
        ]
        with Session.begin() as session:
            queue_df = queries.fetch_studyqueue(session)
            queue_df.loc[:,editable_columns] = queue_df.loc[:,editable_columns].astype(str)
            queue_df = filter_dataframe_dynamic(queue_df.loc[:,queue_display_cols], f'queue')
            edit_queue = st.toggle("Edit queue", value=False)
            if not edit_queue:
                st.dataframe(queue_df, hide_index=True)
            else:
                with st.form('edit_dataset'):
                    st.write('Click on cell to edit:')
                    edited_queue_df = st.data_editor(
                        queue_df.loc[:, not_editable_columns + editable_columns].assign(
                            delete = False,
                        ),
                        disabled=not_editable_columns,
                        hide_index=True,
                        key='editor_queue_df',
                    )
                    stage_edits = st.form_submit_button("Stage Edits")
                    if stage_edits:
                        st.session_state['edited_queue_df'] = edited_queue_df
                if 'edited_queue_df' in st.session_state:
                    if queue_df.index.equals(st.session_state['edited_queue_df'].index):
                        delete_ids = queue_df.loc[:, not_editable_columns + editable_columns].assign(
                            delete = False,
                        ).compare(
                            st.session_state['edited_queue_df'], result_names=("left", "right",), keep_equal=False
                        ).index
                        st.write('Staging queue for editing:')
                        st.dataframe(
                            st.session_state['edited_queue_df'].loc[
                                delete_ids
                            ], 
                            hide_index=True,
                        )
                        submit = st.button(
                            "Edit Queue",
                            on_click=update_queue,
                            args=(
                                session, 
                                st.session_state['edited_queue_df'].loc[
                                    delete_ids
                                ],
                                editable_columns,
                            ),
                            disabled=len(delete_ids) == 0,
                        )
                        if submit and 'edited_result_df' in st.session_state:
                            st.write('Edited:')
                            st.dataframe(
                                st.session_state['edited_result_df'].loc[
                                    :,not_editable_columns + editable_columns
                                ]
                            )

    with tab3:
        with st.form("request_dataset"):
            srp_id = st.text_input(
                "SRP ID:", "ex: SRP,ERP,DRP...",
                key="srp_id",
            )
            geo_id = st.text_input(
                "GEO ID:", "ex: GSE...",
                key="geo_id",
            )
            technology = st.selectbox(
                "Technology:", ("Bulk", "10X-5", "10X-3", "SMART-SEQ", "Other"), 
                index=None,
                key="technology",
            )
            category = st.text_input(
                "Study Category:", 
                "ex: uP discovery efforts in functional dyspesia disease vs control.",
                key="category",
            )
            disease = st.text_input(
                "Study Disease(s):", 
                "ex: Diabetic nephropathy, functional dyspesia, etc.",
                key="disease",
            )
            tissue = st.text_input(
                "Study Tissue(s):", 
                "ex: liver, duodenum, etc.",
                key="tissue",
            )
            contrast = st.text_input(
                "Desired Contrast(s):", 
                "ex: FUNCTIONAL_DYSPESIA_VS_CONTROL",
                key="contrast",
            )
            requestor = st.text_input(
                "Requestor:", 
                "ex: Long, Steve, etc.",
                key="requestor",
            )
            comments = st.text_input(
                "Comments:",
                "ex: Looks good!",
                key="comments",
            )
            priority = st.selectbox(
                "Priority:", ("HIGH", "MEDIUM", "LOW", "HOLD"),
                index=None,
                key="priority",
            )

            placeholder = st.empty()

            submit = st.form_submit_button(
                "Request",
                on_click=submit_queue,
                args=(Session, placeholder),
            )
        if submit and 'submit_results' in st.session_state:
            results = st.session_state['submit_results']
            if results[0] and isinstance(results[1], pd.DataFrame):
                st.write(f'{srp_id} already exists in queue:')
                st.dataframe(results[1].loc[:, results[1].columns.isin(queue_display_cols)])
            elif not results[0] and isinstance(results[1], pd.DataFrame):
                st.write(f'{srp_id} successfully submitted to queue:')
                st.dataframe(results[1].loc[:, results[1].columns.isin(queue_display_cols)])
            else:
                st.error(f'{srp_id} submission to queue failed, inquire about {srp_id}: {results[1]}')
            st.session_state.pop('submit_results')
