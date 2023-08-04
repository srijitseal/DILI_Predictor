import streamlit as st

st.set_page_config(
    
    page_title="DILI Predictor Documentation",
    page_icon="🧊",
    layout="wide",
    initial_sidebar_state="expanded"
    )

st.markdown(
"""
    ### DILI Predictor  Technology Overview 

    DILI Predictor is an open-source app framework built specifically for
    predicting human drug induced liver injury.
    
    DILI Predictor employs eleven proxy-DILI labels from in vitro (e.g., mitochondrial toxicity, bile 
    salt export pump inhibition) and in vivo (e.g., preclinical rat hepatotoxicity studies) datasets 
    along with structural fingerprints and physicochemical parameters as features.
    
    All our data and models are publicly available at GitHub and further details are published at bioaRxiv.
"""
)


