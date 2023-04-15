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

# libraries
import numpy as np
from scipy.integrate import odeint
import plotly.graph_objects as go
import dash
import dash_core_components as dcc
import dash_html_components as html
from tkinter import *

root = Tk()
root.title("pyruvate model")
theLabel = Label(root, text="please create your process")
root.geometry("800x300")

# Label Biomass
labelBiomass = Label(master=root, bg='#F0E68C', text='Starting Biomass in g/L:')
labelBiomass.place(x=54, y=30, width=140, height=27)
labelPyruvate = Label(master=root, bg="light blue", text="Starting Pyruvate:")
labelPyruvate.place(x=54, y=60, width=140, height=27)
labelGlucose = Label(master=root, bg="light green", text="Starting Glucose:")
labelGlucose.place(x=54, y=90, width=140, height=27)
labelAcetate = Label(master=root, bg="light yellow", text="Starting Acetate:")
labelAcetate.place(x=54, y=120, width=140, height=27)
labelVolume = Label(master=root, bg="green", text="Starting Volume:")
labelVolume.place(x=54, y=150, width=140, height=27)
labelVG = Label(master=root, bg="#008B8B",text= "Volumetric Glucose flowrate:" )
labelVG.place(x=450, y=30, width=190, height=27)
labelVA = Label(master=root, bg="#7FFFd4",text= "Volumetric Acetate flowrate:" )
labelVA.place(x=450, y=60, width=190, height=27)
labelFG = Label(master=root, bg="#87CeFa",text= "Feed Glucose Concentration:" )
labelFG.place(x=450, y=90, width=190, height=27)
labelFA = Label(master=root, bg="#b0C4de",text= "Feed Acetate Concentration:" )
labelFA.place(x=450, y=120, width=190, height=27)

# Entry for the Biomass
entryBiomass = Entry(master=root, bg='white')
entryBiomass.place(x=230, y=30, width=40, height=27)
entryPyruvate = Entry(master=root, bg="white")
entryPyruvate.place(x=230, y=60, width=40, height=27)
entryGlucose = Entry(master=root, bg="white")
entryGlucose.place(x=230, y=90, width=40, height=27)
entryAcetate = Entry(master=root, bg="white")
entryAcetate.place(x=230, y=120, width=40, height=27)
entryVolume = Entry(master=root, bg="white")
entryVolume.place(x=230, y=150, width=40, height=27)
entryVG = Entry(master=root, bg="white")
entryVG.place(x=700, y=30, width=40, height=27)
entryVA = Entry(master=root, bg="white")
entryVA.place(x=700, y=60, width=40, height=27)
entryFG = Entry(master=root, bg="white")
entryFG.place(x=700, y=90, width=40, height=27)
entryFA = Entry(master=root, bg="white")
entryFA.place(x=700, y=120, width=40, height=27)

# Button for calculate
buttonCalculate = Button(master=root, bg='#ADFF2F', text='lets calculate', command=root.quit)
buttonCalculate.place(x=320, y=230, width=100, height=27)
#heading in the window
theLabel.pack()
# activation the window
root.mainloop()

# Model for the fed-batch simulation
def pyruvate_fedbatch_model(x0, time, data):
    #start_values
    c_x = x0[0]
    c_p = x0[1]
    c_g = x0[2]
    c_a = x0[3]
    V = x0[4]

    #constants
    #unchangable constants
    c_pmax = 63.6
    Y_xgmax = 0.918 #+/- 0.101
    Y_xamax = 1.52  #+/- 0.752
    Y_pg = 0.932    #+/- 0.380
    mu_max = 1.69   #+/- 0.287
    K_gs = 0.110    #+/- 0.105
    K_as = 0.322    #+/- 0.194
    m_g = 0.0490    #+/- 0.0177
    m_a = 0.147     #+/- 0.0660
    alpha = 0.610   #+/- 0.334
    beta = 0.214    #+/- 0.0289
    K_p = 3.93      #+/- 1.32

    #algebraic_equations
    mu = mu_max * c_g / (K_gs + c_g) * c_a / (K_as + c_a) * K_p / (c_p + K_p)
    '''r_p = (alpha * dcxdt + beta * c_x) * (1 - c_p / c_pmax)'''
    Y_xg = (Y_xgmax * mu) / (Y_xgmax * m_g + mu)
    Y_xa = (Y_xamax * mu) / (Y_xamax * m_a + mu)
    r_g = mu / Y_xg
    r_a = mu / Y_xa

    #differential_equations
    dVdt = data["q_vg"] + data["q_va"]
    dcxdt = - dVdt / V * c_x + mu * c_x
    r_p = (alpha * dcxdt + beta * c_x) * (1 - c_p / c_pmax)
    dcgdt = - dVdt / V * c_g + data["q_vg"] / V * data["c_gnull"] - r_g * c_x - r_p * c_x
    dcadt = - dVdt / V * c_a + data["q_va"] / V * data["c_anull"] - r_a * c_x
    dcpdt = - dVdt / V * c_p + r_p * c_x * Y_pg

    return [dcxdt, dcpdt, dcgdt, dcadt, dVdt]

#setting timeframe and time-intervall
t = np.arange(0, 50, 0.1)

data = dict()
# user defined constants are 0 for batch phase
data["q_vg"] = 0
data["q_va"] = 0
data["c_gnull"] = float(entryFG.get())   #concentration
data["c_anull"] = float(entryFA.get())      #concentration

#User input for Starting-values of differential equations
Biomass_init = float(entryBiomass.get())
Pyruvate_init = float(entryPyruvate.get())
Glucose_init = float(entryGlucose.get())
Acetate_init = float(entryAcetate.get())
Volume_init = float(entryVolume.get())


#adding initial values to group start_values
start_values = [Biomass_init, Pyruvate_init, Glucose_init, Acetate_init, Volume_init]

#simulation of the pyruvate_model
simdata = odeint(pyruvate_fedbatch_model, start_values, t, (data,))

# creating index to find out at what point glucose is dropping to 0
index = 0
while(True):
    if(index >= simdata.shape[0] -1 or simdata[index][2] <= 0):
        break
    index += 1

data["q_vg"] = float(entryVG.get())
data["q_va"] = float(entryVA.get())

#taking batch data just until glucose=0
simdatanew = simdata[:index, :]

#starting fed-batch-model with ending data of batch-model
simdatafed = odeint(pyruvate_fedbatch_model, simdata[index], t, (data,))

#appending batch data and fed_batch data to a single result
result = np.append(simdatanew, simdatafed, axis=0)

print(result.shape)
print(result)

# Display the result of the simulation graphically
fig = go.Figure()

# Add Traces
fig.add_trace(go.Scatter(x=t, y=result[:, 0], name='Biomass Concentration', line =dict( color='black')))
fig.add_trace(go.Scatter(x=t, y=result[:, 2], name='Glucose Concentration', line=dict(color='red')))
fig.add_trace(go.Scatter(x=t, y=result[:, 1], name='Pyruvate Concentration', line=dict(color='green')))
fig.add_trace(go.Scatter(x=t, y=result[:, 3], name='Acetate Concentration', line=dict(color='blue')))

# Add figure title, x- and y-axis title, and legend title
fig.update_layout ( xaxis_title='Time [h]',
                    legend_title='Legend Title',
                    yaxis_title='Concentrations [gL^-1]' )

app = dash.Dash ( __name__ )

app.layout = html.Div ( children=[html.H1 ( children='Modeling of the pyruvate production with Escherichia coli '
                                                     'in a fed-batch bioreactor', style={'textAlign': 'center'} ),
                                  html.Div ( children='''
                                Simulation of Pyruvate Model 9''', style={'textAlign': 'center'} ),
                                  html.Div ( children= '''  ''',
                                             style={'textAlign':'center'}),
                                  dcc.Graph ( id='Simulation of Pyruvate Model 9', figure=fig )
                                  ])

if __name__ == '__main__':
    app.run_server ( debug=False )

