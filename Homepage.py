import streamlit as st

st.set_page_config(
    page_title="DILI Predictor",
    page_icon="Logo.png",
    layout="wide",
    initial_sidebar_state="expanded"
)

left_col, right_col = st.columns(2)

right_col.write("# Welcome to DILI Predictor")
right_col.write("v2.1.0")
right_col.write("Created by Srijit Seal, Dominic Williams, Ola Spjuth and Andreas Bender")
left_col.image("Logo.png")



st.sidebar.success("")

st.markdown(
"""
    DILI Predictor is an open-source app framework built specifically for
    human drug-induced liver injury (DILI)  

    DILI Predictor employs eleven proxy-DILI labels from in vitro 
    (e.g., mitochondrial toxicity, bile salt export pump inhibition) 
    and in vivo (e.g., preclinical rat hepatotoxicity studies) 
    datasets along with pharmacokinetic parameters, structural fingerprints and physicochemical 
    parameters as features.
    
    Select from the sidebar to predict DILI for a single molecule!
    For bulk jobs, or local use: use code from Github page: https://github.com/srijitseal/DILI_Predictor
    
    ### Want to learn more?
    - Check out our paper at [bioRxiv](Appearing here soon!)
    """
)
st.markdown("---")

left_col, right_col = st.columns(2)

right_col.image("Logo.png")

left_col.markdown(
        """
        ### Usage
        On the left pane is the main menu for navigating to 
        the following pages in the PK Predictor application:
        - **Home Page:** You are here!
        - **Documentation:** Explains how the algorithm works and the overview of the models
        - **Submit single SMILES:** You can enter the smiles of the query compound here to obtain a detailed analysis of the predicted DILI and predicitons for proxy-DILI labels and PK parameters used by the model.
        """
    )
st.markdown("---")


left_info_col, right_info_col = st.columns(2)

left_info_col.markdown(

        f"""
        ### Authors
        
        ##### Srijit Seal 
        - Email:  <seal@broadinstitute.org>
        - GitHub: https://github.com/srijitseal
        ##### Andreas Bender
        - Email: <ab454@cam.ac.uk>
        """,
        unsafe_allow_html=True,
    )

right_info_col.markdown(
        """
        ### Funding
        - Cambridge Centre for Data-Driven Discovery and Accelerate Programme for Scientific Discovery
         """
    )

right_info_col.markdown(
        """
        ### License
        Apache License 2.0
        """
    )
