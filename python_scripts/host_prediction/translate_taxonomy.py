import os # Importing built-in module for handling files and directories 
import sys # Built in python library to access system specific parameters such as home directory, path of current script etc...  
         # This is used when we need a file or an object from the OS (Operating System) level. For example os.path here returns absolute paths for different files and directories 
import pandas as pd    # A fast data structure in python to handle large scale real-time analytics which uses open source library   Pandas, we can use this when dealing with table like file formats (CSV). It is used by thousands of developers around the world.    
                       
# Function definition for extracting genus part from NCBI taxonomy using RaFAH's output  # Define a function to get ncbi_genus and return it if found in s else 'none'.  
def get_ncbi_genus(s):     
    '''This method is used as helper functions for getting the genus part from NCBI taxonomy.'''     # Adding comments inside this line, explaining what function does 
                                                            assertion are done here to ensure that file extensions and names of input files match our expectations  
                                                                                  os module in python allows us interaction with operating system level functionalities such as creating absolute path name etc...    sys library is a built-in Python Module which provides access to some variables used by the interpreter, like LANG, PATH ,SYSTEM32,...,etc.  # using these modules we can handle command line arguments passed while running our script
def get_gtdb_upto_genus(s):   '''This method is defined as helper function for extracting complete taxonomic classification up to genus from GTDB'''    assertion are done here, similar purpose like above about file extensions and names of input files.  This will be used when reading the translation table csv
     # Main Function: Execution starts at this point if script is executed as a main program otherwise it won’t start executing further   def in Python which declares block (i.e., code within {}) are called blocks, here we have defined our functionality to be performed when calling the python file ‘py'
    for ss in s.split(";"):
        if ss.startswith("g__"): return ss[3:]
    return "none"

# Extract the complete taxonomic classification up to genus from GTDB
def get_gtdb_upto_genus(s):
    return s.split(";s__")[0]

def main():

    # set the filename of the rafah output file containing the best prediction per contig
    rafah_filename = os.path.abspath(sys.argv[1])
    assert rafah_filename.endswith(".tsv")

    # set the filename of the translation table with NCBI and GTDB taxonomic IDs
    translation_filename = os.path.abspath(sys.argv[2])
    assert translation_filename.endswith(".csv")
    
    # set the output filename
    out_filename = os.path.abspath(sys.argv[3])
    assert out_filename.endswith(".csv")

    # read the input tables using pandas
    translate_df = pd.read_csv(translation_filename)
    rafah_df = pd.read_csv(rafah_filename, sep='\t')

    # we are only interested in the first 3 columns
    rafah_df = rafah_df.iloc[:, :3]

    # create a dictionary for efficient matching of NCBI to GTDB
    translation_dict = {get_ncbi_genus(r["ncbi_taxonomy"]):get_gtdb_upto_genus(r["gtdb_taxonomy"]) for i, r in translate_df.iterrows()}

    # translate all NCBI genus predictions to GTDB.
    # We don't manipulate the table while iterating over it because this could lead to problems.
    translated_hosts = []
    for i, row in rafah_df.iterrows():
        translated_hosts.append(translation_dict[row["Predicted_Host"]] if row["Predicted_Host"] in translation_dict else "no_translation") 

    # add the translated predictions to the dataframe and save the table.
    rafah_df["Predicted_Host"] = translated_hosts
    rafah_df.to_csv(out_filename,index=False)

if __name__ == "__main__":
    main()

