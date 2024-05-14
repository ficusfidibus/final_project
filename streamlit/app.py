
import streamlit as st
import pandas as pd
from padelpy import from_sdf
import pubchempy as pcp
import pickle
import glob
from io import StringIO
import os

# set pathes for outputs
sdf_output_path = "./data/sdf/output"
descriptors_output_path = "./data/descriptors.csv"

# creates the folder if not exists
if not os.path.exists(sdf_output_path):
    os.makedirs(sdf_output_path)

for file in glob.glob(f"{sdf_output_path}/*.sdf"):
    os.remove(file)

def split_into_files(input, output_path):
    # import libraries
    import os
    
    # create output dict
    chem_info = {"pubchemid": [],
                #"iupac_name": [],
                "isomeric_smiles": []
                }
    
    # get the file content and translate it from binary to normal
    content = StringIO(input.getvalue().decode("utf-8"))

    blocks = content.read().split("$$$$")

    for i, block in enumerate(blocks, start=1):
        # check if block is empty
        if not block.strip():
            continue
        # writes the blocks into a new file + "$$$$", which is needed, that padelpy recognize the end
        with open(f"{output_path}/output_{i:03}.sdf", "w") as writer:
            writer.write(block.strip())
            writer.write("\n\n$$$$")

        with open(f"{output_path}/output_{i:03}.sdf", "r") as reader:
            query = reader.read()

        try:
            result = pcp.get_compounds(query, "sdf")
            # get the compound ID from the result
            compound_id = result[0].cid
            # get the name
            #compound_name = result[0].iupac_name
            try:
                # get the SMILES
                compound_smile = result[0].isomeric_smiles
            except:
                compound_smile = None
            # create df
            chem_info["pubchemid"] += [compound_id]
            #chem_info["iupac_name"] += [compound_name]
            chem_info["isomeric_smiles"] += [compound_smile]

        except:
            result = "PubChem information couldn't be loaded"
            chem_info["pubchemid"] += [result]
            #chem_info["iupac_name"] += [None]
            chem_info["isomeric_smiles"] += [None]

        

    return chem_info


def get_sdf(input_type, input):

    result = pcp.get_compounds(input, input_type)
    
    # get the compound ID from the result
    compound_id = result[0].cid
    # get the name
    compound_name = result[0].iupac_name
    # get the SMILES
    compound_smile = result[0].isomeric_smiles

    # create df
    chem_info = {"pubchemid": [compound_id],
                #"iupac_name": [compound_name],
                "isomeric_smiles": [compound_smile]
                }

    # get the SDF from from PubCHem based on cid
    sdf = pcp.get_sdf(compound_id)

    # maybe not useful
    if sdf:
        # write the response into a SDF file
        with open(f"{sdf_output_path}/test_molecule.sdf", "w") as writer:
            writer.write(sdf)
        st.write("SDF_saved")
    else:
        print("no SDF found")
    
    return sdf, chem_info

def calc_descriptors(input_path, output_path):
    import glob
    from padelpy import from_sdf

    sdf_paths = glob.glob(f"{input_path}/*.sdf")
    sdf_paths.sort()

    not_calculated = []
    descriptors = []
    indices = []

    for i, path in enumerate(sdf_paths):
        st.write(f"Calculation of descriptors, molecule {i+1}")
        try:
            descriptors += from_sdf(path,
                                fingerprints=False,
                                descriptors=True,
                                timeout=60,
                                threads=-1)
        except:
            st.write(f"Descriptors for molecule {i+1} coudnt be calculated")
            not_calculated += [i]
            continue

    descriptors_df = pd.DataFrame(descriptors)
    descriptors_df.to_csv(output_path, index=False)
    # need to load the df from csv to get dtypes as numeric. Because of None values the whole df is object...
    descriptors_df = pd.read_csv(output_path)
    return descriptors_df, not_calculated

def load_model(model_nr="model_1"):
    model_file_path = {"model_1": "./model1_desc_nusvc_01.pkl"}

    # load model from pickle file
    with open(model_file_path[model_nr], 'rb') as file:  
        model = pickle.load(file)
    return model

########################################################################################################################            

allowed_dtypes = ["PubChemID", "SMILES notation", "InChI-Key", "Molecule Name", "SDF File"]
dtypes_dict = {"SMILES notation": "smiles", "InChI-Key": "inchikey", "Molecule Name": "name", "PubChemID": "cid"}
chem_info_df = pd.DataFrame()

st.title("Anti-Biofilm activity predictor")

# Select input datatype
dtype_selection = st.radio("Select input dataype", allowed_dtypes)
st.write(f"You selected: {dtype_selection}")

# Text area
input = st.text_area("Enter SMILES, InChI-Key, Molecule Name, PubChemID, ChemblID")
# safe outcome if searching for more
safe = st.toggle("Safe prediction for iterative queries")

# File upload
file = st.file_uploader("""Upload a SDF file (More than one entry per file is possible: delimiter '\$\$\$\$')""")

# Button
if st.button("Predict anti Biofilm activity"):

    with st.spinner("Prediction is processing..."):

        # if file or input
        if file:  
            # check what data is uploaded and start if .sdf
            if file.name[-3:] == "sdf":
                # get the file content and translate it from binary to normal
                chem_info = split_into_files(file, sdf_output_path)

        elif input:
            st.write(f"You entered: {input}")
            sdf, chem_info = get_sdf(input=input, input_type=dtypes_dict[dtype_selection])
             
        # Calulate descriptors
        if safe:
            try: 
                former_desc_df = pd.read_csv(descriptors_output_path)
                descriptors_df, not_calculated = calc_descriptors(sdf_output_path, descriptors_output_path)
                desc_safe_df = pd.concat([former_desc_df, descriptors_df])
                desc_safe_df.to_csv(descriptors_output_path, index=False)

            # safe if first run
            except:
                descriptors_df, not_calculated = calc_descriptors(sdf_output_path, descriptors_output_path)
                descriptors_df.to_csv(descriptors_output_path, index=False)
        else:
            descriptors_df, not_calculated = calc_descriptors(sdf_output_path, descriptors_output_path)

        # load the model
        model = load_model()
        # prediction
        prediction = model.predict(descriptors_df)

        chem_info_df = pd.DataFrame(chem_info)

        chem_info_df.loc[~chem_info_df.index.isin(not_calculated), "Prediction"] = prediction
        chem_info_df.loc[chem_info_df["Prediction"] == 1, "Prediction"] = "active"
        chem_info_df.loc[chem_info_df["Prediction"] == 0, "Prediction"] = "inactive"
        # if safe file toggled safe the result of prediction into an csv file. If former search is done, 

        if safe:
            # Concat with former search
            try: 
                former_info_df = pd.read_csv("./prediction.csv")
                chem_info_df = pd.concat([former_info_df, chem_info_df])
                chem_info_df.to_csv("./prediction.csv", index=False)
            # safe if first run
            except:
                chem_info_df.to_csv("./prediction.csv", index=False)

        st.markdown(f"""#### Prediction:""")
        st.dataframe(chem_info_df)

expand = st.expander("Download CSV's")

with expand:
    try:
        with open(descriptors_output_path, 'rb') as f:
                st.download_button('Download descriptor CSV', f, file_name='descriptors.csv')

        st.download_button('Download prediction CSV', chem_info_df.to_csv(index=False), file_name="prediction.csv")
    except:
        print("no data yet")
    
if st.button("Reset queries"):
    for file in glob.glob(f"./*.csv"):
        os.remove(file)
    for file in glob.glob(f"./data/*.csv"):
        os.remove(file)



