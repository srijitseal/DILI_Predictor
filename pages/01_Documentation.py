import streamlit as st

st.set_page_config(
    
    page_title="PKSmart Predictor",
    page_icon="🧊",
    layout="wide",
    initial_sidebar_state="expanded"
    )

st.image("pages/DrugWiseEndpoints.png")



st.markdown(
"""
    ### PKSmart PK Prediction Technology Overview 

    PKSmart is an open-source app framework built specifically for
    predicting human and animal PK parameters.
    
    PKSmart PK Predictor technology uses chemical structure, physicochemical properties
    and predicted animal properties to predict human PK parameters for small molecule drug
    discovery. 

    Our models are trained on chemicall diverse datasets and employ an applicability domain analysis
    and fold error estimation to give a predicted value along with a interval within which the model
    predicts the error to lie.

    We also provide you with the overlayed chemical space showing where your compound lies with respect
    to our training data. 

    All our data and models are publicly available at GitHub and further details are published at bioaRxiv.
"""
)

st.sidebar.success("")

st.markdown("---")

col_1, col_2, col_3 = st.columns([1,6,1])
col_2.image("logo_bottom.png")

col_1, col_3, col_4, col_2 = st.columns([1, 4, 4,1])
col_3.markdown(
        """
        ### Human PK Parameters
        PKSmart's PK Predictor can be used to predict Human PK parameters of :
         
        - steady-state volume of distribution VDss(L/kg)
        - clearance CL (mL/min/kg)
        - half-Life t½ (h) 
        - fraction unbound in plasma fu (dimensionless)  
        - mean residence time MRT (h)"""
)
col_4.markdown(  """          
        ###  Animal PK parameters 
        PK Predictor can be used to predict animal PK parameters for dog, monkey and rat for:

        - steady-state volume of distribution VDss(L/kg)
        - clearance CL (mL/min/kg)
        - fraction unbound in plasma fu (dimensionless) 
                
        """
    )
st.markdown("---")

col_1, col_2, col_3 = st.columns([1,8,1])

col_2.image("./pages/Modelling_workflow.png")
col_2.markdown(
        """
        ### How PKSmart works:
        
        Chemical Structure, Physicochemical property and predicted Animal PK parameters
         
        - The PK Predictor is based on the compound library of the Lombardo dataset. This involves 1,283 compounds with annotated human PK parameters and 372 unique compounds in the animal dataset

        - PKSmart generates Morgan fingerprints for structural features and Mordred descriptors as physicpchemical features. 
        - Feature selection methods filter and search the most distinct fingerprint bits and descriptors relevant to the compounds in the dataset.
        - Random Forest models build on the 372 compounds in the animal dataset allow using these selected chemical structure and physicochemical properties to predict animal PK parameters of 1,283 compounds in the human dataset.
        - Based on the combinaton of structural fingerprints, physicochemical properties and predicted animal PK parameters, PKSmart trains a model on the human dataset.
        - The final human PK predictor predicts the desired PK parameters. Extensive knowledge of out-of-fold nested cross validation predictions allow us to estimate fold error based on structural similarity of compound to the training data.
        - PKSmart also enables overlaying a query compound on the physicochemical space of the training data thus enabling applicability domain analysis.

"""
)
col_2.image("./pages/ChemicalSpace.png")

st.markdown("---")

col_1, col_2, col_3 = st.columns([1,8,1])
col_2.markdown(
        """
        ### How well does PKSmart work?
        
        Our models in the nested cross validation outperform many state-of-art publicly available tools for PK prediction. 
        - Steady-state volume of distribution VDss(L/kg) achieved a **R2 = 0.55**, with **76%** compounds within in 3-fold error.
        - clearance CL (mL/min/kg) achieved a **R2 = 0.31**, with **71%** compounds within in 3-fold error.
        - half-Life t½ (h) achieved a **R2 = 0.31**, with **76%** compounds within in 3-fold error.
        - fraction unbound in plasma fu (dimensionless) achieved a **R2 = 0.56**, with **71%** compounds within in 3-fold error.
        - mean residence time MRT (h) achieved a **R2 = 0.29**, with **70%** compounds within in 3-fold error.

"""
)
col_2.image("./pages/NCV_mean_results.png")

st.markdown("---")

col_1, col_2, col_3 = st.columns([1,8,1])
col_2.markdown(
        """
        ### How do we estimate fold error?
        
        Our results in the nested cross validation allow us to analyse the relation between fold error and structural similarity of a compound. 
        - We then determine the similarity of the query compound to the entire human training data.
        - Next, we use a Kernel ridge regression fit to estimate the fold error based on this structural similarity.
        - This process is repeated for each endpoint of human PK parameters.

"""
)
col_2.image("./pages/Applicability Domain_how.png")
col_2.image("./pages/Applicability Domain.png")
st.markdown("---")

left_info_col, right_info_col = st.columns(2)

left_info_col.markdown(

        f"""
        ### Authors
        
        ##### Srijit Seal [![Follow](https://img.shields.io/twitter/follow/srijitseal?style=social)](https://www.twitter.com/srijitseal)
        - Email:  <ss2686@cam.ac.uk>
        - GitHub: https://github.com/srijitseal
        ##### Andreas Bender [![Follow](https://img.shields.io/twitter/follow/AndreasBenderUK?style=social)](https://www.twitter.com/AndreasBenderUK)
        - Email: <ab454@cam.ac.uk>
        """,
        unsafe_allow_html=True,
    )

right_info_col.markdown(
        """
        ### Funding
        - Cambridge Centre for Data Driven Discovery and Accelerate Programme for Scientific Discovery under the project title “Theoretical, Scientific, and Philosophical Perspectives on Biological Understanding in the Age of Artificial Intelligence”, made possible by a donation from Schmidt Futures
        - Cambridge Commonwealth, European and International Trust
        - Boak Student Support Fund (Clare Hall)
        - Jawaharlal Nehru Memorial Fund
        - Allen, Meek and Read Fund
        - Trinity Henry Barlow (Trinity College)
         """
    )

right_info_col.markdown(
        """
        ### License
        Apache License 2.0
        """
    )