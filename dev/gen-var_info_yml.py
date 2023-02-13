# this script is to generate the var_info.yml file
# import common libs
# %%
import json
import yaml
import pandas as pd

from pathlib import Path

# %%
# load var info csv files
p_dir_csv = Path("../docs/source/proc_var_info/")
dict_df_var_info = {}
for p_fn_csv in p_dir_csv.glob("*.csv"):
    print(p_fn_csv)
    df_var_info = pd.read_csv(p_fn_csv)
    # print(df_var_info)
    dict_df_var_info[p_fn_csv.stem] = df_var_info

# merge all dataframes
df_var_info = (
    pd.concat(
        dict_df_var_info,
        names=["var_type"],
        # ignore_index=True,
    )
    .reset_index()
    .set_index(["variable", "var_type"])
    .drop(columns="level_1")
    .rename(index=lambda x: x.replace("df_", ""))
    .reset_index("var_type")
)
# p_fn_csv = Path("../docs/source/proc_var_info/df_state.csv")
# df_var_info = pd.read_csv(p_fn_csv)
for var in df_var_info.index.unique():
    len_var = df_var_info.loc[var].shape
    if len(len_var) > 1:
        if len_var[0] > 1:
            df_var = df_var_info.loc[var].dropna(how="all", axis=1).drop_duplicates()
            len_var = df_var.shape
            print(var, len_var[0])
            print(df_var.loc[var].apply(lambda x: x.dropna().unique()))
            print()

df_var_info.loc["DensSnow_Water"].dropna(how="all", axis=1).drop_duplicates()
# %%
# load "loading path" json file
p_fn_json = Path("../src/supy/var2siteselect.json")
with open(p_fn_json, "r") as json_file:
    json_data = json.load(json_file)

# %%
# Create a dictionary with the JSON data stored under new sub-keys
yaml_data = {}
for key in df_var_info.index.unique():
    print(key)
    try:
        df_var = (
            df_var_info.loc[key]
            .dropna(how="all", axis=1)
            .drop_duplicates()
        )
        if df_var['var_type'][0] == 'output':
            df_var=df_var.apply(lambda x: x.dropna().unique())
            # if length of unique values is 1, then reduce to a scalar
            df_var=df_var.apply(lambda x: x[0] if len(x)==1 else x.tolist())

        print(df_var)
        dict_var = df_var.to_dict()
    except ValueError:
        dict_var = df_var_info.loc[key].dropna().to_dict()
    yaml_data[key] = dict_var
    if key in json_data:
        yaml_data[key]["path_loading"] = json_data[key]


# Write the YAML data to a file
with open("var_info.yml", "w") as yaml_file:
    yaml.dump(
        yaml_data,
        yaml_file,
        default_flow_style=False,
        sort_keys=False,
    )
