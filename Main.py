# coding=utf-8
'''This is the main file for our Pyruvate_model simulation. It contains the GUI, user inputs, the import of
Model_mergeable and the interactive dashboard. Please make sure that Model_mergeable.py is located in the same
folder as this file to run.'''


'''List of symbols:
cA      acetate concentration (g L)1)
cA,0    acetate concentration in the feed (g L)1)
cG      glucose concentration (g L)1)
cG,0    glucose concentration in the feed (g L)1)
cP      pyruvate concentration (g L)1)
cP,max  critical pyruvate concentration above which
        reaction cannot proceed (g L)1)
cX      biomass concentration (g L)1)
KI      inhibition constant for pyruvate production
        (g L)1)
K_I A   inhibition constant for biomass growth on
        acetate (g L)1)
KP      saturation constant for pyruvate production
        (g L)1)
KP      inhibition constant of Jerusalimsky (g L)1)
K_S A   Monod growth constant for acetate (g L)1)
K_S G   Monod growth constant for glucose (g L)1)
mA      maintenance coefficient for growth on acetate
        (g g)1 h)1)
mG      maintenance coefficient for growth on glucose
        (g g)1 h)1)
n       constant of extended Monod kinetics (Levenspiel) (–)
q_V     volumetric flow rate (L h)1)
q_VA    volumetric flow rate of acetate (L h)1)
q_VG    volumetric flow rate of glucose (L h)1)
r_A     specific rate of acetate consumption
        (g g)1 h)1)
r_G     specific rate of glucose consumption
        (g g)1 h)1)
r_P     specific rate of pyruvate production
        (g g)1 h)1)
r_P,max maximum specific rate of pyruvate production
        (g g)1 h)1)
t       time (h)
V       reaction (broth) volume (L)
Y_PG    yield coefficient pyruvate from glucose
        (g g)1)
Y_XA    yield coefficient biomass from acetate (g g)1)
Y_XA,max maximum yield coefficient biomass from acetate
        (g g)1)
Y_XG    yield coefficient biomass from glucose (g g)1)
Y_XG,max maximum yield coefficient biomass from glucose
        (g g)1)'''

# libraries and imports
import numpy as np
from scipy.integrate import odeint
import plotly.graph_objects as go
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from tkinter import *
import Model_mergeable
import webbrowser

#GUI integration
root = Tk ()
root.title ( "pyruvate model" )
theLabel = Label ( root, text="please create your process" )
root.geometry ( "800x300" )

# Labels for different user inputs
labelBiomass = Label ( master=root, bg='#F0E68C', text='Starting Biomass in g/L:' )
labelBiomass.place ( x=54, y=30, width=140, height=27 )
labelPyruvate = Label ( master=root, bg="light blue", text="Starting Pyruvate:" )
labelPyruvate.place ( x=54, y=60, width=140, height=27 )
labelGlucose = Label ( master=root, bg="light green", text="Starting Glucose:" )
labelGlucose.place ( x=54, y=90, width=140, height=27 )
labelAcetate = Label ( master=root, bg="light yellow", text="Starting Acetate:" )
labelAcetate.place ( x=54, y=120, width=140, height=27 )
labelVolume = Label ( master=root, bg="light pink", text="Starting Volume:" )
labelVolume.place ( x=54, y=150, width=140, height=27 )
labelVG = Label ( master=root, bg="#008B8B", text="Volumetric Glucose flowrate:" )
labelVG.place ( x=450, y=30, width=190, height=27 )
labelVA = Label ( master=root, bg="#7FFFd4", text="Volumetric Acetate flowrate:" )
labelVA.place ( x=450, y=60, width=190, height=27 )
labelFG = Label ( master=root, bg="#87CeFa", text="Feed Glucose Concentration:" )
labelFG.place ( x=450, y=90, width=190, height=27 )
labelFA = Label ( master=root, bg="#b0C4de", text="Feed Acetate Concentration:" )
labelFA.place ( x=450, y=120, width=190, height=27 )

# Entries for user inputs
entryBiomass = Entry ( master=root, bg='white' )
entryBiomass.place ( x=230, y=30, width=40, height=27 )
entryPyruvate = Entry ( master=root, bg="white" )
entryPyruvate.place ( x=230, y=60, width=40, height=27 )
entryGlucose = Entry ( master=root, bg="white" )
entryGlucose.place ( x=230, y=90, width=40, height=27 )
entryAcetate = Entry ( master=root, bg="white" )
entryAcetate.place ( x=230, y=120, width=40, height=27 )
entryVolume = Entry ( master=root, bg="white" )
entryVolume.place ( x=230, y=150, width=40, height=27 )
entryVG = Entry ( master=root, bg="white" )
entryVG.place ( x=700, y=30, width=40, height=27 )
entryVA = Entry ( master=root, bg="white" )
entryVA.place ( x=700, y=60, width=40, height=27 )
entryFG = Entry ( master=root, bg="white" )
entryFG.place ( x=700, y=90, width=40, height=27 )
entryFA = Entry ( master=root, bg="white" )
entryFA.place ( x=700, y=120, width=40, height=27 )

# Button for calculation
buttonCalculate = Button ( master=root, bg='#ADFF2F', text='lets calculate', command=root.quit )
buttonCalculate.place ( x=320, y=230, width=100, height=27 )
# heading to the window
theLabel.pack ()
# activation of the window
root.mainloop ()

# Model inputs that are used for setting timeframe and time-intervall for Model
t = np.arange ( 0, 50, 0.1 )

# converting user inputs into variables that are used in the model
data = dict ()
data["q_vg"] = 0  # Batch Feed Flow = 0
data["q_va"] = 0  # Batch Feed Flow = 0
data["c_gnull"] = float ( entryFG.get () )  # concentration
data["c_anull"] = float ( entryFA.get () )  # concentration
q_vgfed = float ( entryVG.get () )  # Fedbatch GlucoseFeed Flow is user definable
q_vafed = float ( entryVA.get () )  # Fedbatch AcetateFeed Flow is user definable

# User input for Starting-values of differential equations of our model
Biomass_init = float ( entryBiomass.get () )  # User input for Biomass
Pyruvate_init = float ( entryPyruvate.get () )  # User input for Pyruvate
Glucose_init = float ( entryGlucose.get () )  # User input for Glucose
Acetate_init = float ( entryAcetate.get () )  # User input for Acetate
Volume_init = float ( entryVolume.get () )  # User input for Volume

# import of our Model which is accessed from file: Model_mergeable.py
result = Model_mergeable.calculate ( data, Biomass_init, Pyruvate_init, Glucose_init, Acetate_init, Volume_init,
                                     q_vafed, q_vgfed, t )

# Display the result of the simulation graphically
fig = go.Figure ()

# Add Traces
fig.add_trace ( go.Scatter ( x=t, y=result[:, 0], name='Biomass Concentration', line=dict ( color='black' ) ) )
fig.add_trace ( go.Scatter ( x=t, y=result[:, 2], name='Glucose Concentration', line=dict ( color='red' ) ) )
fig.add_trace ( go.Scatter ( x=t, y=result[:, 1], name='Pyruvate Concentration', line=dict ( color='green' ) ) )
fig.add_trace ( go.Scatter ( x=t, y=result[:, 3], name='Acetate Concentration', line=dict ( color='blue' ) ) )

# Add figure title, x- and y-axis title, and legend title
fig.update_layout(xaxis_title='Time [h]',
                  legend_title='Legend Title',
                  yaxis_title='Concentrations [gL^-1]')

webbrowser.open('http://127.0.0.1:8050/')  # to open dash page automatically


# Making an interactive Dashboard
app = dash.Dash ( __name__ )

# Creating the layout for dashboard
app.layout = html.Div(children=[html.H1(children='Modeling of the pyruvate production with Escherichia coli '
                                                     'in a fed-batch bioreactor', style={'textAlign': 'center'}),
                                  html.Div(children='Simulation of Pyruvate Model 9', style={'textAlign': 'center'}),
                                  html.Div(children=' Özge Kilic ',
                                             style={'textAlign': 'center'}),
                                  dcc.Graph(id='Simulation of Pyruvate Model 9', figure=fig),
                                  html.Div(children='Change Parameters for Fed-Batch Process:',
                                             style={'textAlign': 'left', 'font-weight':'Bold', 'font-size': 25}),
                                  html.Div([
                                      html.Span(children='Initial Biomass Concentration:',
                                                  style={'width': '15%', 'display': 'inline-block'}),
                                      dcc.Input(id="InitialBiomassConcentration", type="number",
                                                  placeholder="Initial Biomass Concentration",
                                                  style={'margin-right': '10px'}),
                                      html.Span(children='Volumetric Glucose Flowrate:',
                                                  style={'width': '15%', 'display': 'inline-block'}),
                                      dcc.Input(id="VolumetricGlucoseFlowrate", type="number",
                                                  placeholder="Volumetric Glucose Flowrate",
                                                  style={'margin-right': '10px'}),
                                      html.Span(children='Initial Volume:',
                                                  style={'width': '15%', 'display': 'inline-block'}),
                                      dcc.Input(id="InitialVolume", type="number",
                                                  placeholder="Initial Volume",
                                                  style={'margin-right': '10px'})]),
                                  html.Div([
                                      html.Span(children='Initial Glucose Concentration:',
                                                  style={'width': '15%', 'display': 'inline-block'}),
                                      dcc.Input(id="InitialGlucoseConcentration", type="number",
                                                  placeholder="Initial Glucose Concentration",
                                                  style={'margin-right': '10px'}),
                                      html.Span(children='Volumetric Acetate Flowrate:',
                                                  style={'width': '15%', 'display': 'inline-block'}),
                                      dcc.Input(id="VolumetricAcetateFlowrate", type="number",
                                                  placeholder="Volumetric Acetate Flowrate",
                                                  style={'margin-right': '10px'}),
                                      html.Span(children='Feed Glucose Concentration:',
                                                  style={'width': '15%', 'display': 'inline-block'}),
                                      dcc.Input(id="FeedGlucoseConcentration", type="number",
                                                  placeholder="Feed Glucose Concentration",
                                                  style={'margin-right': '10px'})]),
                                  html.Div([
                                      html.Span(children='Initial Pyruvate Concentration:',
                                                  style={'width': '15%', 'display': 'inline-block'}),
                                      dcc.Input(id="InitialPyruvateConcentration", type="number",
                                                  placeholder="Initial Pyruvate Concentration",
                                                  style={'margin-right': '10px'}),
                                      html.Span(children='Feed Acetate Concentration:',
                                                  style={'width': '15%', 'display': 'inline-block'}),
                                      dcc.Input(id="FeedAcetateConcentration", type="number",
                                                  placeholder="Feed Acetate Concentration",
                                                  style={'margin-right': '10px'}),
                                      html.Span(children='Initial Acetate Concentration:',
                                                  style={'width': '15%', 'display': 'inline-block'}),
                                      dcc.Input(id="InitialAcetateConcentration", type="number",
                                                  placeholder="Initial Acetate Concentration",
                                                  style={'margin-right': '10px'})]),

                                  html.Div([
                                      html.Span(style={'width': '40%', 'display': 'inline-block'}),
                                      html.Button('Recalculate', id='recalculate-val')
                                  ], style={'padding': '10px'})])

# To recalculate the model when user input new value
@app.callback (
    Output ( "Simulation of Pyruvate Model 9", "figure" ),
    Input ( "recalculate-val", "n_clicks" ),
    state=[State ( "InitialBiomassConcentration", "value" ),
           State ( "InitialGlucoseConcentration", "value" ),
           State ( "InitialPyruvateConcentration", "value" ),
           State ( "InitialAcetateConcentration", "value" ),
           State ( "VolumetricGlucoseFlowrate", "value" ),
           State ( "VolumetricAcetateFlowrate", "value" ),
           State ( "FeedGlucoseConcentration", "value" ),
           State ( "FeedAcetateConcentration", "value" ),
           State ( "InitialVolume", "value" ),
     ])

# define the function to update the graph with the new value
def update_graph(recalValue,IBC, IGC, IPC, IAC, VGF, VAF, FGC, FAC, IV):
    data = dict ()
    data["q_vg"] = 0  # Batch Feed Flow = 0
    data["q_va"] = 0  # Batch Feed Flow = 0
    data["c_gnull"] = FGC  # concentration
    data["c_anull"] = FAC  # concentration
    result = Model_mergeable.calculate (data, IBC, IPC, IGC, IAC, IV,
                                         VAF, VGF, t )
    fig = go.Figure ()

    # Add Traces
    fig.add_trace ( go.Scatter ( x=t, y=result[:, 0], name='Biomass Concentration', line=dict ( color='black' ) ) )
    fig.add_trace ( go.Scatter ( x=t, y=result[:, 2], name='Glucose Concentration', line=dict ( color='red' ) ) )
    fig.add_trace ( go.Scatter ( x=t, y=result[:, 1], name='Pyruvate Concentration', line=dict ( color='green' ) ) )
    fig.add_trace ( go.Scatter ( x=t, y=result[:, 3], name='Acetate Concentration', line=dict ( color='blue' ) ) )

    # Add figure title, x- and y-axis title, and legend title
    fig.update_layout ( xaxis_title='Time [h]',
                        legend_title='Legend Title',
                        yaxis_title='Concentrations [gL^-1]' )

    return fig

# to run dash server
if __name__ == '__main__':
    app.run_server ( debug=False )
