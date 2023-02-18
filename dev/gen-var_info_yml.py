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
        names=["type"],
        # ignore_index=True,
    )
    .reset_index()
    .set_index(["variable", "type"])
    .drop(columns="level_1")
    .rename(index=lambda x: x.replace("df_", ""))
    .reset_index("type")
)

# debugging info: show all keys with multiple values
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
        df_var = df_var_info.loc[key].dropna(how="all", axis=1).drop_duplicates()
        if df_var["type"][0] == "output":
            df_var = df_var.apply(lambda x: x.dropna().unique())
            # if length of unique values is 1, then reduce to a scalar
            df_var = df_var.apply(lambda x: x[0] if len(x) == 1 else x.tolist())

        print(df_var)
        dict_var = df_var.to_dict()
    except ValueError:
        dict_var = df_var_info.loc[key].dropna().to_dict()
    # replace all keys to lower case
    yaml_data[key] = {k.lower(): v for k, v in dict_var.items()}
    if key in json_data:
        yaml_data[key]["loading path"] = json_data[key]

# re-organise the yaml data
# merge groups under the key 'type'
# remove the key 'group' from the sub-keys
for k, v in yaml_data.items():
    dict_var = v
    if "group" in dict_var:
        dict_var["type"] = {"output": dict_var["group"]}
        del dict_var["group"]
    else:
        dict_var["type"] = {"input": dict_var["type"]}

    # rename the key 'dimensionality' to 'data dimensions':
    if "dimensionality" in dict_var:
        dict_var["data dimensions"] = [
            dict_var["dimensionality"],
            {"remarks": dict_var["dimensionality remarks"]},
        ]
        del dict_var["dimensionality"]
        del dict_var["dimensionality remarks"]

    # add a key 'scheme' to setting-related variables
    if "physics scheme" not in dict_var and "input" in dict_var["type"]:
        if "state" in dict_var["type"]["input"]:
            dict_var["physics scheme"] = {
                "scheme to add": [
                    "code 1",
                    "code 2",
                ],
            }

    # write the data back to the yaml_data
    yaml_data[k] = dict_var


# Write the YAML data to a file
with open("var_info.yml", "w") as yaml_file:
    yaml.dump(
        yaml_data,
        yaml_file,
        default_flow_style=False,
        sort_keys=False,
    )
# %%
yaml_data[key]
# %%
sorted(yaml_data)
# %%
import supy as sp
df_state,df_forcing=sp.load_SampleData()
# %%
df_state.shape
# %%
