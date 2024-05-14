def split_into_files(input_path, output_path):
    # importing libraries needed
    import os
    import pubchempy as pcp
    # create output dict
    chem_info = {"pubchemid": [],
                "iupac_name": [],
                "isomeric_smiles": []
                }

    with open(input_path, "r") as input:
        content = input.read()

        blocks = content.split("$$$$")
        # creates the folder if not exists
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        for i, block in enumerate(blocks, start=1):
            # check if block is empty
            if not block.strip():
                continue
            try:
                result = pcp.get_compounds(block, "sdf")
                # get the compound ID from the result
                compound_id = result[0].cid
                # get the name
                compound_name = result[0].iupac_name
                # get the SMILES
                compound_smile = result[0].isomeric_smiles
                # create df
                chem_info["pubchemid"] += [compound_id]
                chem_info["iupac_name"] += [compound_name]
                chem_info["isomeric_smiles"] += [compound_smile]

            except:
                result = "PubChem information couldn't be loaded"
                chem_info["pubchemid"] += [result]
                chem_info["iupac_name"] += [result]
                chem_info["isomeric_smiles"] += [result]

           
            # create the name for the file to be used
            file_name = input_path.split("/")[-1].split(".")[0]

            # writes the blocks into a new file + "$$$$", which is needed, that padelpy recognize the end
            with open(f"{output_path}/{file_name}{i:03}.sdf", "w") as writer:
                writer.write(block.strip())
                writer.write("\n\n$$$$")
    return chem_info


