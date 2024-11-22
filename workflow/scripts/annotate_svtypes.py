import pandas as pd
import numpy as np

# paths
input_svtable = snakemake.input["all_SVs"]
caller_tables = snakemake.input["caller_tables"]
out_svtable = snakemake.output["all_SVs"]
# test paths
# input_svtable = ("/Users/asherpreskasteinberg/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/"
#                  "lab_notebook/APS017.1_3x3_SV_analysis/SVChordinator_troubleshoot/SVChordinator_test/somatic_SVs"
#                  "/sample.filtered_ensemble.annotated.tsv")
# caller_tables = [
#     "/Users/asherpreskasteinberg/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/lab_notebook/APS017.1_3x3_SV_analysis/SVChordinator_troubleshoot/SVChordinator_test/raw_SVs/sample/sample.nanomonsv.tsv",
#     "/Users/asherpreskasteinberg/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/lab_notebook/APS017.1_3x3_SV_analysis/SVChordinator_troubleshoot/SVChordinator_test/raw_SVs/sample/sample.SAVANA.tsv",
#     "/Users/asherpreskasteinberg/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/lab_notebook/APS017.1_3x3_SV_analysis/SVChordinator_troubleshoot/SVChordinator_test/raw_SVs/sample/sample.Severus.tsv"
# ]
# out_svtable = ("/Users/asherpreskasteinberg/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/"
#                  "lab_notebook/APS017.1_3x3_SV_analysis/SVChordinator_troubleshoot/SVChordinator_test/somatic_SVs"
#                  "/test.reannotated.tsv")

# load in the input from the different individual callers as one unified pandas dataframe
all_callers = pd.DataFrame()
for i in np.arange(0, len(caller_tables)):
    table_path = caller_tables[i]
    caller_table = pd.read_csv(table_path, sep="\t")
    caller_table["caller"] = "caller_{i}".format(i=i)
    all_callers = pd.concat([all_callers, caller_table])
# rename columns
all_callers.columns = ["chrom", "pos", "ID", "SV_Type", "Strands", "BP_notation", "caller"]
# load input table of consensus SV calls
all_svtable = pd.read_csv(input_svtable, sep="\t")
# search for BNDs (the unresolved breakpoint)
bnd_table = all_svtable[all_svtable["SV_Type"] == "BND"]
# reannotate translocations first because they're easy
tra_table = bnd_table[bnd_table["chrom1"] != bnd_table["chrom2"]]
tra_table = tra_table[tra_table["SV_Type"] == "TRA"]
# now let's re-annotate unresolved breakpoints to their closest match
# based on strands and annotation from individual callers
unresolved_table = bnd_table[bnd_table["chrom1"] == bnd_table["chrom2"]]
sv_types = []
for i, row in unresolved_table.iterrows():
    call_ids = row["SV_callers"].split(",")
    # get breakpoint ids
    bp_ids = []
    for call in call_ids:
        terms = call.split("_")
        bp_id = "_".join(terms[1:])
        bp_ids.append(bp_id)
    # search for these ids in the all_callers table
    bp_df = all_callers[all_callers["ID"].isin(bp_ids)]
    # iterate through until we can annotate
    sv_type = "BND"
    j = 0
    while sv_type == "BND" and j < len(bp_df):
        bp_row = bp_df.iloc[j]
        strands = bp_row["Strands"]
        if bp_row["SV_Type"] != "BND":
            sv_type = bp_row["SV_Type"]
        elif strands == "++" or strands == "--":
            sv_type = "INV"
        elif strands == "+-":
            sv_type = "DEL"
        elif strands == "-+":
            sv_type = "DUP"
        j = j + 1
    # append sv type
    sv_types.append(sv_type)

# make new table
reannotated_table = unresolved_table
reannotated_table["SV_Type"] = sv_types

# stitch everything together into final output
resolved_table = all_svtable[all_svtable["SV_Type"] != "BND"]
# add translocations
final_table = pd.concat([resolved_table, tra_table, reannotated_table])
# write output
final_table.to_csv(out_svtable, sep="\t", index=False)





