import base64
import json
import pathlib
import os
import time

from os import listdir
from os import path
from os.path import isfile, join
import numpy as np
import dash
import dash_core_components as dcc
import dash_html_components as html
from io import BytesIO
from Bio import SeqIO
import pandas as pd
import plotly.graph_objs as go

from keras.models import load_model
from PIL import Image
from dash.dependencies import Input, Output, State
import pickle
import plotly.express as px
import dash_bootstrap_components as dbc
import dash_table
from plotly.colors import n_colors
from Bio import SeqIO
import metagenomics_processor
from MetaG_backend import backend
import metagenomics_counter
import pyVectorizer
# from helpers import load_mnist, parse_image, numpy_to_b64, create_img, label_mapping
import run_tsne

# Setting up upload path
UPLOAD_DIRECTORY = "fasta_files"
if not os.path.exists(UPLOAD_DIRECTORY):
    os.makedirs(UPLOAD_DIRECTORY)


# get relative data folder
PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("data").resolve()



# Import settings data from csv
SETTINGS_PATH = PATH.joinpath("settings").resolve()
settings_df = pd.read_csv(SETTINGS_PATH.joinpath("sample_settings.csv"))
settings_df['Index'] = range(1, len(settings_df) + 1)
settings_df = settings_df[["Index","FASTA File ", "Iterations","k in kmer","Autoencoder","2D","3D","Activation"]]
#print(settings_df)

PAGE_SIZE = 5

with open(PATH.joinpath("demo_intro.md"), "r") as file:
    demo_intro_md = file.read()

with open(PATH.joinpath("PCATSNE_intro.md"), "r") as file:
    PCATSNE_intro_md = file.read()

with open(PATH.joinpath("Taxa_intro.md"), "r") as file:
    Taxa_intro_md = file.read()


def numpy_to_b64(array, scalar=True):
    # Convert from 0-1 to 0-255
    if scalar:
        array = np.uint8(255 * array)

    im_pil = Image.fromarray(array)
    buff = BytesIO()
    im_pil.save(buff, format="png")
    im_b64 = base64.b64encode(buff.getvalue()).decode("utf-8")

    return im_b64


# Methods for creating components in the layout code
def Card(children, **kwargs):
    return html.Section(children, className="card-style")


def NamedSlider(name, short, min, max, step, val, marks=None):
    if marks:
        step = None
    else:
        marks = {i: i for i in range(min, max + 1, step)}

    return html.Div(
        style={"margin": "25px 5px 30px 0px"},
        children=[
            f"{name}",
            html.Div(
                style={"margin-left": "5px"},
                children=[
                    dcc.Slider(
                        id=f"slider-{short}",
                        min=min,
                        max=max,
                        marks=marks,
                        step=step,
                        value=val,
                    )
                ],
            ),
        ],
    )


def NamedInlineRadioItems(name, short, options, val, **kwargs):
    return html.Div(
        id=f"div-{short}",
        style={"display": "inline-block"},
        children=[
            f"{name}:",
            dcc.RadioItems(
                id=f"radio-{short}",
                options=options,
                value=val,
                labelStyle={"display": "inline-block", "margin-right": "7px"},
                style={"display": "inline-block", "margin-left": "7px"},
            ),
        ],
    )


AE_Dict = {'3': ['[32,x,32]', '[32,4,x,4,32]', '[32,16,x,16,32]'], '4': ['[136,64,x,64,136]', '[136,32,x,32,136]', '[136,16,x,16,136]']}
kmer = list(AE_Dict.keys())
AE_Options = AE_Dict[kmer[0]]
uploaded_fasta_files = [f.split('.')[0] for f in listdir("fasta_files") if isfile(join("fasta_files", f))]

def create_layout(app):
    # Actual layout of the app
    return html.Div(
        className="row",
        style={"max-width": "100%", "font-size": "1.5rem", "padding": "0px 0px"},
        children=[
            # Header
            html.Div(
                className="row header",
                id="app-header",
                style={"background-image": 'url("/assets/white_bg.jpg")', 'background-size': '1600px 1000px'},
                children=[
                    html.Div(
                        [
                            html.Img(
                                src=app.get_asset_url("logo.png"),
                                className="logo",
                                id="plotly-image",
                            ),
                            html.H3(
                                "MetaG",
                                #className="header_title",
                                id="meta-app-title",
                                style={'font-family': 'cursive', 'color':'Black', 'padding-left':"200px !important", 'text-shadow': '2px 2px #ff0000'},
                            ),
                            
                        ],
                        className="three columns header_img",
                    ),
                    html.Div(
                        [
                            html.H5(
                                "Comprehensive Metagenomics Analysis Tool",
                                className="header_title",
                                id="app-title",
                                style={'font-family': 'arial', 'color':'Black'},
                            )
                        ],
                        className="nine columns header_title_container",
                    ),
                ],
            ),
            html.Div(
                className="row background",
                children=[
                dcc.Loading(
            id="loading-1",
            type="default",
            children=html.Div(id="loading-output-1")
        ),
                html.Div(
                    className="twelve columns",
                    id="tab-div",
                    
                        children=[
                            Card(dcc.Tabs(id="main-tabs", value='tab1', children=[dcc.Tab(label='Basic Autoencoder Analysis', value='tab1',
                                children=[
                                    # Demo Description
            html.Div(
                className="row background",
                id="pcatsne-explanation",
                style={"padding": "50px 45px"},
                children=[
                    
                    html.Div(
                        id="description-text", children=dcc.Markdown(demo_intro_md)
                    ),
                ],
            ),
            # Body
            html.Div(
                className="row background",
                style={"padding": "10px"},
                children=[
                    html.Div(
                        style={"border": "1px solid", "box-shadow": "5px 5px 5px #888888", "border-radius":"5px", "background-color":"#E6E5E6"},
                        className="three columns",
                        children=[
                            Card(dcc.Tabs(id="tabs", value='tab1', children=[dcc.Tab(label='Settings', value='tab1', id ="basic-analysis",
                                children=

                                [
                                    html.Hr(),
                                    html.P("Metagenomic sample"),
                                    dcc.Dropdown(
                                        id="fasta-dropdown",
                                        style={"margin-top": "50px"},
                                        searchable=False,
                                        clearable=False,
                                        options=[],
                                        placeholder="Select a dataset",
                                        value="Synthetic1",
                                    ),
                                    html.Br(),
                                    html.Hr(),
                                    html.P("Choose Training Sample Percentage"),
                                    dcc.Input(id="training_downsample", type="number", value=60, min=10, max=100, step=5),
                                    html.Br(),
                                    html.Hr(),
                                    
                                     html.Div([
                                        html.P("Choose Iterations"),
                                        dcc.Slider(
                                            id = "iterations-slider",
                                            min=200,
                                            max=7000,
                                            step=50,
                                            value=2950,
                                        ),
                                        html.Div(id='slider-output-container')
                                    ]),
                                     html.Hr(),
                                     html.P("Choose k in k-mers"),
                                    dcc.Dropdown(
                                        id = "kmer",
                                        options=[{'label':k, 'value':k} for k in kmer],
                                        value='4'
                                    )  ,
                                    html.Hr(),
                                    html.P("Choose Autoencoder"),
                                    dcc.Dropdown(
                                        id='ae-dropdown',
                                        value= "[136,64,x,64,136]"
                                    ),
                                     html.Hr(),
                                    html.P("Choose Dimensions"),
                                    dcc.RadioItems(
                                        id = "dimension",
                                        inputStyle ={'margin-right':'30px'},
                                        options=[
                                            {'label': '2D    ', 'value': '2d'},
                                            {'label': '3D    ', 'value': '3d'},
                                           
                                        ],
                                        value='2d'
                                    )  ,
                                     html.Hr(),
                                    html.P("Activation Function"),
                                    dcc.RadioItems(
                                        id = "activation",
                                        inputStyle ={'margin-right':'30px'},
                                        options=[
                                            {'label': 'Sigmoid', 'value': 'sigmoid'},
                                            {'label': 'Tanh', 'value': 'tanh'},
                                        ],
                                        value='sigmoid'
                                    )  ,
                                    html.Hr(),
                                    dbc.Button("Submit", id="submit-btn", color="primary", className="mr-1", style={"color":"White", "background-color":"#32CD32", "margin": "0 auto", "display": "block"}),
                                    #dbc.Alert(id='submit-alert'),
                                    html.Hr(),
                                   html.Div( id='submit-alert'),
                                    
                                ]),
                                dcc.Tab(label='Data', value='tab2', children=[
                        
                        html.Br(),
                        html.P(
                            "FASTA file containing the sequence reads of the metagenomic sample is required to do the analysis and prompt the visualizations."
                        ),
                         html.Br(),
                        html.Img(src='assets/reads.PNG', style={"width":"100%", "height":"50%", "padding-bottom":"20px"}),
                        html.Br(),
                        html.P("Please upload your FASTA file containing metagenomic shotgun sequence reads in order to get the visualizations."),
                        html.Br(),
                        html.P("After uploading, please go to the settings tab and select the uploaded FASTA file from the dropdown menu and add the settings."),
                        
                        
                        html.Span(
                            className="control-label", children=[html.H5("Upload the sequence file")]
                        ),
                        dcc.Upload(
                            id="upload-data",
                            className="upload-component",
                            children=html.Div(
                                ["Drag and Drop or ", html.A("Select Files")]
                            ),
                            style={"padding-bottom":"400px"},
                        ),
                        html.Div(id='upload-output'),
                        ]),
                                ])
                            )
                        ],
                    ),
                     html.Div(
                        className="nine columns",
                        
                        children=[
                            html.Div(
                            style={"border": "1px solid", "box-shadow": "5px 5px 5px #888888", "border-radius":"5px", "padding":"10px 10px 70px 10px", "background-color":"#E6E5E6"},
                            children=[
                            html.H3(
                                        className="graph-title",
                                        children="Current Results",
                                    ),
                            dash_table.DataTable(
                                id='table-paging-and-sorting',
                                columns=[
                                    {'name': i, 'id': i, 'deletable': True} for i in settings_df.columns
                                ],
                                page_current=0,
                                page_size=PAGE_SIZE,
                                page_action='custom',
                                row_selectable="single",
                                selected_rows=[0],
                                sort_action='custom',
                                sort_mode='single',
                                sort_by=[],
                                style_as_list_view=True,
                                style_header={'backgroundColor': 'rgb(30, 30, 30)', 'fontWeight': 'bold'},
                                
                                style_cell={
                                    'backgroundColor': 'rgb(50, 50, 50)',
                                    'color': 'white'
                                },
                            ),
                            ],
                            ),
                        ],
                    ),
                     html.Br(),
                    
                    html.Br(),
                     html.Div(
                         className="nine columns",
                         style = {"padding-top":"20px"},
                         children=[
                    html.Div(
                        
                        style={"border": "1px solid", "box-shadow": "5px 5px 5px #888888", "border-radius":"5px", "padding":"10px 10px 10px 10px", "background-color":"#E6E5E6"},
                        children=[
                            html.H3(    id = "dr-graph-title",
                                        className="graph-title",
                                        children="Sequence data reduced Representation ",
                                    ),
                            dcc.Graph(id="graph-3d-plot-tsne", style={"height": "98vh", "margin-bottom":"50px"})
                        ],
                    ),
                ],
            ),

                    html.Div(
                         className="nine columns",
                         style = {"padding-top":"20px"},
                         children=[
                    html.Div(
                        style={"border": "1px solid", "box-shadow": "5px 5px 5px #888888", "border-radius":"5px", "padding":"10px 10px 10px 10px", "background-color":"#E6E5E6"},
                        children=[
                            html.H3(
                                        className="graph-title",
                                        children="Metagenomic Nucleotide Frequencies",
                                        
                                    ),
                            dcc.Graph(id='example-graph-2', style={"height": "98vh","margin-bottom":"50px"})
                        ],
                    ),
                ],
            ),
                    html.Div(
                         className="nine columns",
                         style = {"padding-top":"20px"},
                         children=[
                     html.Div(
                        
                        style={"border": "1px solid", "box-shadow": "5px 5px 5px #888888", "border-radius":"5px", "padding":"10px 10px 10px 10px", "background-color":"#E6E5E6", "margin-left":"320px","margin-right":"-320px"},
                        children=[
                            html.H3(
                                        className="graph-title",
                                        children="Sequence Content",
                                        
                                    ),
                            dcc.Graph(id='ACGT_donut_chart', style={"height": "98vh", "margin-bottom":"50px"})
                        ],
                    ),
                ],
            ),
                    html.Div(
                        className="twelve columns",
                        id="footer",
                        children=[
                            html.Div(style={"height":"100px"}),
                        ],
                    ),
                    html.Div(
                        className="three columns",
                        id="euclidean-distance",
                        children=[
                            Card(
                                style={"padding": "5px"},
                                children=[
                                    
                                    html.Div(id="div-plot-click-image"),
                                    html.Div(id="div-plot-click-wordemb"),
                                ],
                            )
                        ],
                    ),
                    
                ],
            ),

                                ],
                            ),
                            dcc.Tab(label='Taxanomic Analysis', value='tab2', children=[
                                html.Div(
                                style={"padding": "10px"},
                                children=[
                                html.Div(
                                    className="row background",
                                    id="taxa-explanation",
                                    style={"padding": "50px 45px"},
                                    children=[
                                        
                                        html.Div(
                                            id="taxa", children=dcc.Markdown(Taxa_intro_md)
                                        ),
                                    ],
                                ),
                                html.Div(
                         className="twelve columns",
                         style = {"padding-top":"20px"},
                         children=[
                            
                            html.Br(),
                                    html.Div(
                                        
                                        children=[
                                    html.Div(
                                        
                                        style={"border": "1px solid", "box-shadow": "5px 5px 5px #888888", "border-radius":"5px", "padding":"10px 10px 10px 10px", "background-color":"#E6E5E6"},
                                        children=[
                                            html.H3(    id = "hybrid-graph-title",
                                                        className="graph-title",
                                                        children="Reduced Representation with Taxonomy Details ",
                                                    ),
                                            dcc.Graph(id="hybrid-plot", style={"height": "98vh", "margin-bottom":"50px"})
                                        ],
                                    ),
                                ],
                            ),
                            
                            html.Div(
                                style={"border": "1px solid", "box-shadow": "5px 5px 5px #888888", "border-radius":"5px", "padding":"10px 10px 10px 10px", "background-color":"#E6E5E6"},
                                children=[
                                dcc.Graph(id='domain_donut_chart', style={"height": "49vh","width": "49%",  "display":"inline-block", "padding-right":"2%"}),
                                dcc.Graph(id='phylum_donut_chart', style={"height": "49vh","width": "49%", "display":"inline-block"})
                            ]),
                            html.Div(
                                style={"border": "1px solid", "box-shadow": "5px 5px 5px #888888", "border-radius":"5px", "padding":"10px 10px 10px 10px", "background-color":"#E6E5E6"},
                                children=[
                                dcc.Graph(id='class_donut_chart', style={"height": "49vh","width": "49%",  "display":"inline-block", "padding-right":"2%"}),
                                dcc.Graph(id='order_donut_chart', style={"height": "49vh", "width": "49%","display":"inline-block"})
                            ]),
                            html.Div(
                                style={"border": "1px solid", "box-shadow": "5px 5px 5px #888888", "border-radius":"5px", "padding":"10px 10px 10px 10px", "background-color":"#E6E5E6"},
                                children=[
                                dcc.Graph(id='family_donut_chart', style={"height": "49vh","width": "49%",  "display":"inline-block", "padding-right":"2%"}),
                                dcc.Graph(id='genus_donut_chart', style={"height": "49vh","width": "49%", "display":"inline-block"})
                            ]),
                            html.Div(
                                style={"border": "1px solid", "box-shadow": "5px 5px 5px #888888", "border-radius":"5px", "padding":"10px 10px 10px 10px", "background-color":"#E6E5E6"},
                                children=[
                                dcc.Graph(id='species_donut_chart', style={"height": "75vh",  "align":"center"}),
                
                            ]),
                            
                            
                ],
            ),
                                    html.Div(
                                                className="twelve columns",
                                                style = {"padding-top":"20px"},
                                                children=[
                                                    html.Div(
                                                        
                                                        style={"border": "1px solid", "box-shadow": "5px 5px 5px #888888", "border-radius":"5px", "padding":"10px 10px 10px 10px", "background-color":"#E6E5E6"},
                                                        children=[
                                                            html.H3(
                                                                        className="graph-title",
                                                                        children="GC content Distribution",
                                                                        
                                                                    ),
                                                            dcc.Graph(id='violin_plot', style={"height": "98vh", "margin-bottom":"50px"})
                                                        ],
                                                    ),
                                        ],
                                    ),
                                ]
                                ),
                            ],
                            ),
                            dcc.Tab(label='t-SNE / PCA Analysis', value='tab3', children=[
                                html.Div(
                                style={"padding": "10px"},
                                children=[
                               
                                html.Div(
                                    className="row background",
                                    id="demo-explanation",
                                    style={"padding": "50px 45px"},
                                    children=[
                                        
                                        html.Div(
                                            id="pcatsne-intro", children=dcc.Markdown(PCATSNE_intro_md)
                                        ),
                                    ],
                                ),
                                # Body
                                html.Div(
                                    className="row background",
                                    style={"padding": "10px"},
                                    children=[
                                        html.Div(
                                            style={"border": "1px solid", "box-shadow": "5px 5px 5px #888888", "border-radius":"5px", "background-color":"#E6E5E6"},
                                            className="three columns",
                                            children=[
                                                Card(id="pcatsne-tabs", 
                                                    children=

                                                    [   html.H5("t-SNE Settings"),
                                                      
                                                        html.Hr(),
                                                        html.P("Metagenomic sample"),
                                                        dcc.Dropdown(
                                                            id="pcatsne-fasta-dropdown",
                                                            searchable=False,
                                                            clearable=False,
                                                            options=[],
                                                            placeholder="Select a dataset",
                                                            value="Synthetic1",
                                                        ),
                                                        html.Hr(),
                                                        html.P("Kmer"),
                                                        dcc.Dropdown(
                                                            id = "pcatsne-kmer",
                                                            options=[{'label':k, 'value':k} for k in kmer],
                                                            value='4'
                                                        )  ,
                                                        html.Hr(),
                                                        html.P("Plot Dimension"),
                                                        dcc.Dropdown(
                                                            id = "pcatsne-dim",
                                                            options=[{'label': '2D    ', 'value': '2d'},{'label': '3D    ', 'value': '3d'}],
                                                            value='2d'
                                                        )  ,
                                                        html.Hr(),
                                                        NamedSlider(
                                                            name="Number Of Iterations",
                                                            short="iterations",
                                                            
                                                            min=250,
                                                            max=1000,
                                                            step=None,
                                                            val=500,
                                                            marks={
                                                                i: str(i) for i in [250, 500, 750, 1000]
                                                            },
                                                        ),
                                                        NamedSlider(
                                                            name="Perplexity",
                                                            short="perplexity",
                                                            
                                                            min=3,
                                                            max=100,
                                                            step=None,
                                                            val=30,
                                                            marks={i: str(i) for i in [3, 10, 30, 50, 100]},
                                                        ),
                                                        NamedSlider(
                                                            name="Initial PCA Dimensions",
                                                            short="pca-dimension",
                                                            
                                                            min=20,
                                                            max=100,
                                                            step=10,
                                                            val=50,
                                                            marks={i: str(i) for i in [20,30,40,50,60,70,80,90,100]},
                                                        ),
                                                        NamedSlider(
                                                            name="Learning Rate",
                                                            short="learning-rate",
                                                           
                                                            min=20,
                                                            max=200,
                                                            step=None,
                                                            val=100,
                                                            marks={i: str(i) for i in [20,40,60,80,100,120, 140,160, 180, 200]},
                                                        ),
                                    
                                                        
                                                        dbc.Button("Submit", id="pcatsne-submit-btn", color="primary", className="mr-1", style={"color":"White", "background-color":"#32CD32", "margin": "0 auto", "display": "block"}),
                                                        #dbc.Alert(id='submit-alert'),
                                                        html.Hr(),
                                                        html.Div( id='pcatsne-submit-alert'),
                                                        
                                                    ],
                                                    
                                                    )
                                                
                                            ],
                                        ),
                                        
                                        html.Div(
                                            className="nine columns",
                                            children=[
                                        html.Div(
                                            
                                            style={"border": "1px solid", "box-shadow": "5px 5px 5px #888888", "border-radius":"5px", "padding":"10px 10px 10px 10px", "background-color":"#E6E5E6"},
                                            children=[
                                                html.H3(    id = "tsne-graph-title",
                                                            className="graph-title",
                                                            children="t-SNE Representation ",
                                                        ),
                                                dcc.Graph(id="tsne-plot", style={"height": "98vh", "margin-bottom":"50px"})
                                            ],
                                        ),
                                    ],
                                ),
                                
                                       
                                        
                                        
                                        
                                    ],
                                ),
                                ])])
                            ]
                            )
                            )
                        ]
                
                )
                ]

            ),
            
            
        ],
    )


def demo_callbacks(app):
    @app.callback(
    dash.dependencies.Output('slider-output-container', 'children'),
    [dash.dependencies.Input('iterations-slider', 'value')])
    def update_output(value):
        return 'You have selected {} iterations'.format(value)


    def preprocessing():
        print("preprocessing started ...")



    @app.callback(
    dash.dependencies.Output('ae-dropdown', 'options'),
    [dash.dependencies.Input('kmer', 'value')])
    def update_ae_dropdown(k):
        return [{'label': i, 'value': i} for i in AE_Dict[k]]

    @app.callback(
    dash.dependencies.Output('fasta-dropdown', 'options'),
    [Input("upload-data", "filename"),Input("upload-data", "contents")])
    def update_fasta_dropdown(uploaded_filenames,uploaded_file_contents):
        uploaded_fasta_files = [f.split('.')[0] for f in listdir("fasta_files") if isfile(join("fasta_files", f))]
        if uploaded_filenames is not None:
            file_name = '_'.join(uploaded_filenames.split(".")[:-1])
            mod_uploaded_filenames = file_name.replace(" ","_") + ".fasta"
            uploaded_fasta_files.append(mod_uploaded_filenames.split('.')[0])
        
        return [{'label': i, 'value': i} for i in uploaded_fasta_files]

    @app.callback(
    dash.dependencies.Output('pcatsne-fasta-dropdown', 'options'),
    [Input("upload-data", "filename"),Input("upload-data", "contents")])
    def update_fasta_dropdown(uploaded_filenames,uploaded_file_contents):
        uploaded_fasta_files = [f.split('.')[0] for f in listdir("fasta_files") if isfile(join("fasta_files", f))]
        if uploaded_filenames is not None:
            file_name = '_'.join(uploaded_filenames.split(".")[:-1])
            mod_uploaded_filenames = file_name.replace(" ","_") + ".fasta"
            uploaded_fasta_files.append(mod_uploaded_filenames.split('.')[0])
        
        return [{'label': i, 'value': i} for i in uploaded_fasta_files]

    def save_file(name, content):
        """Decode and store a file uploaded with Plotly Dash."""
        #print(content)
        #print(name)
        with open(os.path.join(UPLOAD_DIRECTORY, name), "w") as fp:
            fp.write(content)
    
    
    @app.callback(
    Output("upload-output", "children"),
    [Input("upload-data", "filename"), Input("upload-data", "contents")])
    def update_output(uploaded_filenames, uploaded_file_contents):
        """Save uploaded files and regenerate the file list."""
        
        if uploaded_filenames is not None and uploaded_file_contents is not None:
            decoded_content = base64.b64decode(uploaded_file_contents.split(',')[-1]).decode('utf-8')
            file_extention = uploaded_filenames.split(".")[-1]
            file_name = '_'.join(uploaded_filenames.split(".")[:-1])
            uploaded_filenames = file_name.replace(" ","_") + ".fasta"
            save_file(uploaded_filenames, decoded_content)
        return html.Div()

    @app.callback(Output("loading-output-1", "children"), [Input("tab-div", "children")])
    def input_triggers_spinner(children):
        time.sleep(2)
        return None

        


    @app.callback(
    Output('table-paging-and-sorting', 'data'),
    [
     Input('table-paging-and-sorting', "page_size"),
     Input('table-paging-and-sorting', 'sort_by'),
     Input('table-paging-and-sorting', 'selected_rows'),
     Input('table-paging-and-sorting', "page_current"),]
     )
    def update_table(page_size, sort_by, selected_rows,page_current ):
        if len(sort_by):
            # Import settings data from csv
            SETTINGS_PATH = PATH.joinpath("settings").resolve()
            settings_df = pd.read_csv(SETTINGS_PATH.joinpath("sample_settings.csv"))
            settings_df['Index'] = range(1, len(settings_df) + 1)
            settings_df = settings_df[["Index","FASTA File ", "Iterations","k in kmer","Autoencoder","2D","3D","Activation"]]
            dff = settings_df.sort_values(
                sort_by[0]['column_id'],
                ascending=sort_by[0]['direction'] == 'asc',
                inplace=False
            )
        else:

            SETTINGS_PATH = PATH.joinpath("settings").resolve()
            settings_df = pd.read_csv(SETTINGS_PATH.joinpath("sample_settings.csv"))
            settings_df['Index'] = range(1, len(settings_df) + 1)
            settings_df = settings_df[["Index","FASTA File ", "Iterations","k in kmer","Autoencoder","2D","3D","Activation"]]
            dff = settings_df

        return dff.iloc[
            page_current*page_size:(page_current+ 1)*page_size
        ].to_dict('records')


    def generate_figure_image(embedding_df):

        df = pd.DataFrame(dict(x=embedding_df["a0"].tolist(),
                               y=embedding_df["a1"].tolist(),
                               z=embedding_df["a2"].tolist()))

        figure = px.scatter_3d(df, x="x", y="y", z="z")
        figure.update_traces(marker_line=dict(width=1, color='DarkSlateGray'))

        return figure


    def generate_figure_2d_image(embedding_df):

        df = pd.DataFrame(dict(x=embedding_df["a0"].tolist(),
                               y=embedding_df["a1"].tolist(),))

        figure = px.scatter(df, x="x", y="y",render_mode='webgl')
        figure.update_traces(marker_line=dict(width=1, color='DarkSlateGray'))

        return figure

    def generate_hybrid_fig(embedding_df):

        df = pd.DataFrame(dict(x=embedding_df["a0"].tolist(),
                               y=embedding_df["a1"].tolist(),
                               z=embedding_df["a2"].tolist(),
                               color = embedding_df["lin_0"].tolist()))

        figure = px.scatter_3d(df, x="x", y="y", z="z", color = "color")
        figure.update_traces(marker_line=dict(width=1, color='DarkSlateGray'))

        return figure
    
    def generate_hybrid_2d_fig(embedding_df):

        df = pd.DataFrame(dict(x=embedding_df["a0"].tolist(),
                               y=embedding_df["a1"].tolist(),
                               color = embedding_df["lin_0"].tolist()))

        figure = px.scatter(df, x="x", y="y",color = "color", render_mode='webgl')
        figure.update_traces(marker_line=dict(width=1, color='DarkSlateGray'))

        return figure
    
     

    @app.callback(
        [Output("graph-3d-plot-tsne", "figure"),
         Output("ACGT_donut_chart", "figure"),
         Output("example-graph-2", "figure"),
         Output("dr-graph-title", "children"),
         Output("hybrid-plot", "figure"),
         Output("domain_donut_chart", "figure"),
         Output("phylum_donut_chart", "figure"),
         Output("class_donut_chart", "figure"),
         Output("order_donut_chart", "figure"),
         Output("family_donut_chart", "figure"),
         Output("genus_donut_chart", "figure"),
         Output("species_donut_chart", "figure"),
         Output("violin_plot", "figure"),
         Output("submit-alert", "children")],
        [Input("submit-btn", "n_clicks"),
         Input('table-paging-and-sorting', 'sort_by'),
         Input('table-paging-and-sorting', 'selected_rows'),
         Input('table-paging-and-sorting', "page_size"),
         Input('table-paging-and-sorting', "page_current")],
        [
            
            State("fasta-dropdown", "value"),
            State("iterations-slider", "value"),
            State("kmer", "value"),
            State("ae-dropdown", "value"),
            State("dimension", "value"),
            State("activation", "value"),
            State("training_downsample", "value"),
            
        ],
    )
    def display_2d_scatter_plot(
        btn1,
        sort_by,
        selected_rows,
        page_size,
        current_page,
        fasta_name,
        iterations,
        kmer,
        ae,
        dim,
        act,
        training_percentage,
        
    ):
            
            print(current_page)
            if selected_rows and not btn1:
                return load_from_table(sort_by, selected_rows, current_page, page_size)
            input_url = [
                "fasta_files",
                str(fasta_name) + ".fasta",
            ]
            input_file_path = PATH.joinpath(*input_url)

            data_url = [
                "results",
                str(fasta_name),
                str(dim)+"_$$_"+str(kmer)+"_$$_"+str(iterations)+"_$$_"+str(ae)+"_$$_"+str(act)+".csv",
            ]
            
            results_path = PATH.joinpath(*data_url)
            folder_created_now = False
            if not os.path.exists("results/"+str(fasta_name)):
                folder_created_now = True
                os.makedirs("results/"+str(fasta_name))

            freq_counter_url = [
                "results",
                str(fasta_name),
            ]
            freq_counter_path = PATH.joinpath(*freq_counter_url)
            donut_fig = calc_acgt_count(input_file_path, int(kmer), freq_counter_path)
            freq_fig = freq_count(input_file_path, int(kmer), freq_counter_path)

            try:
                print("available")
                embedding_df = pd.read_csv(
                    results_path, index_col=0, encoding="ISO-8859-1"
                )

            except FileNotFoundError as error:
                print(
                    error,
                    "\nThe dataset was not found. Please generate it using generate_demo_embeddings.py",
                )

                act = str(act)
                ui_ae = ae
                ae = ae[1:-1].split(",")
                
                ind = ae.index('x')
                ae[ind] = int(dim[0])
                for i in range(len(ae)):
                    ae[i] = int(ae[i])
                
                # l = metagenomics_processor.process(input_file_path, int(kmer), int(iterations), act, ae, results_path)
                

                

                training_percentage = training_percentage/100
                if folder_created_now:
                    l = backend.process(input_file_path, int(kmer), int(iterations), act, ae,training_percentage)
                else:
                    import re
                    for f in os.listdir('./results/'+str(fasta_name)):
                        print("inside 1")
                        if re.match('3d', f):
                            print("inside 3d")
                            embedding_df = pd.read_csv(
                                "./results/"+str(fasta_name)+"/"+str(f), index_col=0, encoding="ISO-8859-1"
                            )
                            l = backend.process_ae(input_file_path, int(kmer), int(iterations), act, ae,
                                                        training_percentage, embedding_df)
                            break
                        elif re.match('2d', f):
                            print("inside 2d")
                            embedding_df = pd.read_csv(
                                "./results/"+str(fasta_name)+"/"+str(f), index_col=0, encoding="ISO-8859-1"
                            )
                            l = backend.process_ae(input_file_path, int(kmer), int(iterations), act, ae,
                                                        training_percentage, embedding_df)
                            break

                        else:
                            l = backend.process(input_file_path, int(kmer), int(iterations), act, ae,
                                                           training_percentage)

                print("out of process")

                l.to_csv(results_path)
                embedding_df = pd.read_csv(
                    results_path, index_col=0, encoding="ISO-8859-1"
                )
                print("read csv")
                settings_path = 'settings/sample_settings.csv'
               
               
                settings = {'FASTA File ': [str(fasta_name) + ".fasta"],
                            'Iterations': [iterations],
                            'k in kmer': [kmer],
                            'Autoencoder': [str(ui_ae).replace(" ", "")],
                            '2D': 'Yes' if dim == "2d" else "",
                            '3D': 'Yes' if dim == "3d" else "",
                            'Activation': [act]
                            }
                print(settings)
                print(pd.DataFrame(settings))
                pd.DataFrame(settings).to_csv(settings_path, mode='a', header=not (path.exists(settings_path)), index=False)


            
            sampled_embedding_df = embedding_df[embedding_df['sampled']]
 
            sampled_embedding_df.loc[:,sampled_embedding_df.columns.str.contains('lin')]=sampled_embedding_df.loc[:,sampled_embedding_df.columns.str.contains('lin')].apply(lambda x : x.fillna(value="NA"))
            
            dim = int(dim[0])
            if dim==2:
                figure = generate_figure_2d_image(sampled_embedding_df)
                graph_title = "Sequence data reduced to 2D Representation "
                hybrid_plot = generate_hybrid_2d_fig(sampled_embedding_df)
            else:
                figure = generate_figure_image(sampled_embedding_df)
                graph_title = "Sequence data reduced to 3D Representation "
                hybrid_plot = generate_hybrid_fig(sampled_embedding_df)
            msg = "Sucsessful !"
            st = {'width':'100%','color'
            :'Green','background-color':'#D3FFC5 ',' border-style':'solid','border-radius': '2px', 'padding':'3px 3px 2px 2px', 'text-align': 'center',"border": "1px solid", "box-shadow": "5px 5px 5px #888888", "border-radius":"5px",}
            domain_donut_figure = generate_domain_donut(embedding_df)
            phylum_donut_figure = generate_phylum_donut(embedding_df)
            class_donut_figure = generate_class_donut(embedding_df)
            order_donut_figure = generate_order_donut(embedding_df)
            family_donut_figure = generate_family_donut(embedding_df)
            genus_donut_figure = generate_genus_donut(embedding_df)
            species_donut_figure = generate_species_donut(embedding_df)
            violin_plot = generate_violin_plot(embedding_df)
            return figure, donut_fig, freq_fig, graph_title, hybrid_plot,domain_donut_figure,phylum_donut_figure,class_donut_figure,order_donut_figure,family_donut_figure,genus_donut_figure,species_donut_figure, violin_plot,html.Div(msg, style = st)

    def calc_acgt_count(input_path, k, output_path):
        try:
            with open(output_path/'_acgt.json') as json_file:
                data = json.load(json_file)

        except FileNotFoundError:
            metagenomics_counter.count_process(input_path, k, output_path)
            with open(output_path/'_acgt.json') as json_file:
                data = json.load(json_file)

        donut_labels = ['A', 'C', 'G', 'T']
        donut_values = data
        donut_fig = go.Figure(data=[go.Pie(labels=donut_labels, values=donut_values, hole=.3)])
        GC_value = (donut_values[1]+donut_values[2])/sum(donut_values)
        GC_value = round(GC_value*100, 3)
        GC_content = "GC-content = " + str(GC_value) + "%"
        donut_fig.update_layout(title=GC_content)
        return donut_fig


    def freq_count(input_path, k, output_path):
        file_name = str(k)+"_freqs.json"
        try:
            with open(output_path/file_name) as json_file:
                data = json.load(json_file)

        except FileNotFoundError:
            metagenomics_counter.freq_process(input_path, k, output_path)
            with open(output_path/file_name) as json_file:
                data = json.load(json_file)

        k_mer_list = pyVectorizer.all_kmers(k)
        freq_fig = go.Figure(data=[go.Scatter(
            # x=[i for i in range(len(data[1]))],
            x=k_mer_list,
            y=data[1],
            mode='lines+markers')])
        return freq_fig


    def load_from_table(sort_by,selected_rows, current_page, page_size):
        print("inside load from table")
        # Import settings data from csv
        SETTINGS_PATH = PATH.joinpath("settings").resolve()
        settings_df = pd.read_csv(SETTINGS_PATH.joinpath("sample_settings.csv"))
        settings_df['Index'] = range(1, len(settings_df) + 1)
        settings_df = settings_df[["Index","FASTA File ", "Iterations","k in kmer","Autoencoder","2D","3D","Activation"]]
            
        if len(sort_by):
            dff = settings_df.sort_values(
                sort_by[0]['column_id'],
                ascending=sort_by[0]['direction'] == 'asc',
                inplace=False
            )
        else:
            dff = settings_df
        # if selected_rows is not None and len(selected_rows) :
            # print(dff.iloc[[selected_rows[0]]]['FASTA File '].tolist())
        print(dff.iloc[[(current_page*page_size)+selected_rows[0]]])
        fasta_name = dff.iloc[[(current_page*page_size)+selected_rows[0]]]['FASTA File '].tolist()[0][:-6]
        iterations = dff.iloc[[(current_page*page_size)+selected_rows[0]]]['Iterations'].tolist()[0]
        kmer = dff.iloc[[(current_page*page_size)+selected_rows[0]]]['k in kmer'].tolist()[0]
        ae = dff.iloc[[(current_page*page_size)+selected_rows[0]]]['Autoencoder'].tolist()[0]
        dim2 = dff.iloc[[(current_page*page_size)+selected_rows[0]]]['2D'].tolist()[0]
        dim3 = dff.iloc[[(current_page*page_size)+selected_rows[0]]]['3D'].tolist()[0]
        act = dff.iloc[[(current_page*page_size)+selected_rows[0]]]['Activation'].tolist()[0]

        print(fasta_name)

        dim=""
        if  isinstance(dim2, str):
            dim = "2d"
        else:
            dim = "3d"

        freq_counter_url = [
            "results",
            str(fasta_name),
        ]
        freq_counter_path = PATH.joinpath(*freq_counter_url)

        freq_fig = freq_count("",int(kmer),freq_counter_path)
        donut_fig = calc_acgt_count("", int(kmer), freq_counter_path)

        scatter_url = [
            "results",
            str(fasta_name),
            str(dim) + "_$$_" + str(kmer) + "_$$_" + str(iterations) + "_$$_" + str(ae) + "_$$_" + str(act) + ".csv",
        ]
        results_path = PATH.joinpath(*scatter_url)
        embedding_df = pd.read_csv(
            results_path, index_col=0, encoding="ISO-8859-1"
        )
       
        sampled_embedding_df = embedding_df[embedding_df['sampled']]
        sampled_embedding_df.loc[:,sampled_embedding_df.columns.str.contains('lin')]=sampled_embedding_df.loc[:,sampled_embedding_df.columns.str.contains('lin')].apply(lambda x : x.fillna(value="NA"))

        

        if int(dim[0]) == 2:
            scatter_figure = generate_figure_2d_image(sampled_embedding_df)
            graph_title = "Sequence data reduced to 2D Representation "
            hybrid_plot = generate_hybrid_2d_fig(sampled_embedding_df)
        else:
            scatter_figure = generate_figure_image(sampled_embedding_df)
            graph_title = "Sequence data reduced to 3D Representation "
            hybrid_plot = generate_hybrid_fig(sampled_embedding_df)
        domain_donut_figure = generate_domain_donut(embedding_df)
        phylum_donut_figure = generate_phylum_donut(embedding_df)
        class_donut_figure = generate_class_donut(embedding_df)
        order_donut_figure = generate_order_donut(embedding_df)
        family_donut_figure = generate_family_donut(embedding_df)
        genus_donut_figure = generate_genus_donut(embedding_df)
        species_donut_figure = generate_species_donut(embedding_df)
        violin_plot = generate_violin_plot(embedding_df)
        

        return scatter_figure, donut_fig, freq_fig, graph_title, hybrid_plot,domain_donut_figure,phylum_donut_figure,class_donut_figure,order_donut_figure,family_donut_figure,genus_donut_figure,species_donut_figure, violin_plot, html.Div()

    def generate_taxa_donut():
        
        labels = ['mycobacterium avium','mycobacterium bovis','mycobacterium tuberculosis','Clostridium botulinum','Clostridium perfringens', 'Clostridium tetani']
        values = [4500, 2500, 1053, 500, 20,80]
        fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=.3)])
        fig.update_traces(marker=dict(line=dict(color='#000000', width=2)))
        return fig
    def generate_domain_donut(df):
        df = df.groupby('lin_7').size().reset_index(name='counts')
        labels = list(df['lin_7'])
        values = list(df['counts'])
        fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=.3)])
        fig.update_layout(title="Domain Summary")
        fig.update_traces(marker=dict(line=dict(color='#000000', width=2)))
        return fig
    def generate_phylum_donut(df):
        df = df.groupby('lin_5').size().reset_index(name='counts')
        labels = list(df['lin_5'])
        values = list(df['counts'])
        fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=.3)])
        fig.update_layout(title="Phylum Summary")
        fig.update_traces(marker=dict(line=dict(color='#000000', width=2)))
        return fig
    def generate_class_donut(df):
        df = df.groupby('lin_4').size().reset_index(name='counts')
        labels = list(df['lin_4'])
        values = list(df['counts'])
        fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=.3)])
        fig.update_layout(title="Class Summary")
        fig.update_traces(marker=dict(line=dict(color='#000000', width=2)))
        return fig
    def generate_order_donut(df):
        df = df.groupby('lin_3').size().reset_index(name='counts')
        labels = list(df['lin_3'])
        values = list(df['counts'])
        fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=.3)])
        fig.update_layout(title="Order Summary")
        fig.update_traces(marker=dict(line=dict(color='#000000', width=2)))
        return fig
    def generate_family_donut(df):
        df = df.groupby('lin_2').size().reset_index(name='counts')
        labels = list(df['lin_2'])
        values = list(df['counts'])
        fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=.3)])
        fig.update_layout(title="Family Summary")
        fig.update_traces(marker=dict(line=dict(color='#000000', width=2)))
        return fig
    def generate_genus_donut(df):
        df = df.groupby('lin_1').size().reset_index(name='counts')
        labels = list(df['lin_1'])
        values = list(df['counts'])
        fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=.3)])
        fig.update_layout(title="Genus Summary")
        fig.update_traces(marker=dict(line=dict(color='#000000', width=2)))
        return fig
    def generate_species_donut(df):
        df = df.groupby('lin_0').size().reset_index(name='counts')
        labels = list(df['lin_0'])
        values = list(df['counts'])
        fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=.3)])
        fig.update_layout(title="Species Summary")
        fig.update_traces(marker=dict(line=dict(color='#000000', width=2)))
        return fig

    def generate_violin_plot(df):
        species = df.lin_0.unique().tolist()

        #gc_content=(g+c)*100/(a+t+g+c)
        #df = pd.read_csv("https://raw.githubusercontent.com/plotly/datasets/master/violin_data.csv")

        fig = go.Figure()


        for s in species:
            fig.add_trace(go.Violin(x=df['lin_0'][df['lin_0'] == s],
                                    y=((df['g']+df['c'])*100/(df['a']+df['c']+df['g']+df['t']))[df['lin_0'] == s],
                                    name=s,
                                    box_visible=True,
                                    meanline_visible=True))


        return fig


    @app.callback(
         [Output("tsne-plot", "figure"),
         Output("tsne-graph-title", "children"),],
        [Input("pcatsne-submit-btn", "n_clicks"),],
        [
            State("pcatsne-fasta-dropdown", "value"),
            State("slider-iterations", "value"),
            State("slider-perplexity", "value"),
            State("slider-pca-dimension", "value"),
            State("slider-learning-rate", "value"),
            State("pcatsne-kmer", "value"),
            State("pcatsne-dim", "value"),
        ],
    )
    def display_tsne_plot(
        pca_tsne_btn,
        fasta_name,
        iterations,
        perplexity,
        pca_dim,
        lr,
        kmer,
        dim,
    ):

            input_url = [
                "fasta_files",
                str(fasta_name) + ".fasta",
            ]
            input_file_path = PATH.joinpath(*input_url)

            data_url = [
                "results",
                str(fasta_name),
                "tsne",
                str(dim)+"_$$_"+str(kmer)+"_$$_"+"pca"+"_$$_"+str(pca_dim)+"_$$_"+"tsne"+"_$$_"+str(iterations)+"_$$_"+str(perplexity)+"_$$_"+str(lr)+".csv",
            ]
            results_path = PATH.joinpath(*data_url)

            if not os.path.exists("results/"+str(fasta_name)+"/tsne"):
                os.makedirs("results/"+str(fasta_name)+"/tsne")


            try:
                embedding_df = pd.read_csv(
                    results_path, index_col=0, encoding="ISO-8859-1"
                )

            except FileNotFoundError as error:
                print(
                    error,
                    "\nThe dataset was not found. Please generate it using generate_demo_embeddings.py",
                )

                run_tsne.process(input_file_path, int(kmer), int(pca_dim), int(dim[0]), int(iterations), int(lr), int(perplexity), results_path)

                print("out of process")
                embedding_df = pd.read_csv(
                    results_path, index_col=0, encoding="ISO-8859-1"
                )
                print("read csv")



            dim = int(dim[0])
            if dim==2:
                figure = generate_figure_2d_image(embedding_df)
                graph_title = "Sequence data reduced to 2D Representation "

            else:
                figure = generate_figure_image(embedding_df)
                graph_title = "Sequence data reduced to 3D Representation "


            return figure, graph_title



