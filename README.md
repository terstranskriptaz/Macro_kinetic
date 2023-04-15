# Macro_kinetic
Macro Kinetic Growth Model for E.coli 
Model simulated in Python 3.8 and the following libraries of Python have been used:
    numpy,
    odeint,
    tkinter,
    plotly.graph_objects,
    webbrowser,
    dash
    

The publication
    "Modeling of the pyruvate production with Escherichia coli in a
    fed-batch bioreactor" by Zelić, B., Vasić-Rački, Đ., Wandrey, C., Takors, R.
    in 2004.
    PMID: 15085423  DOI: 10.1007/s00449-004-0358-0
    
    https://pubmed.ncbi.nlm.nih.gov/15085423/
  
Model abstract:
   In this study 10 different unstructered model has been developed to model cell growth, substrate cG,0 consumption,
    and product formation of the pyruvate producing strain Escherichia coli YYC202 ldhA::Kan strain used in 
    fed-batch processes. With the estimations of the parameters from the fed-batch fermentation with a qVG=10 mL/h 
    glucose feed rate developed feeding strategies with different glucose feed rates (qVG = 10 and 30  mL/h ). 
    
   The main reason to choose model 9 in 10 different model is; the model represents an approximate analogy 
    to the non-competitive substrate inhibition, well reflective experimental data of the prediction biomass 
    and pyruvate curves, lowest residual sum of squares, mechanistically correction, acceptable confidence in interval 
    and it contains all experimentally observed effects like growth inhibition by pyruvate and 
    pyruvate inhibition product formation. 
    
    
To run the code two files are necessary: 
    Main.py and Model_mergaeble.py 

The code enables the user to estimate pyruvate production by E.coli based on a batch-fedbatch simulation. 
We integrated a GUI where the user can give inputs for the simulation.
Furthermore we developed an interactive dashboard to change inputs and simulate the graph with optimised values. 
