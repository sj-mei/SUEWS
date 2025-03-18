import numpy as np
import pandas as pd
import copy
from .supy_driver import suews_driver as sd


##############################################################################
# post-processing part
# get variable information from Fortran
def get_output_info_df():
    from packaging.version import parse as LooseVersion

    size_var_list = sd.output_size()
    var_list_x = [np.array(sd.output_name_n(i)) for i in np.arange(size_var_list) + 1]

    df_var_list = pd.DataFrame(var_list_x, columns=["var", "group", "aggm", "outlevel"])

    # strip leading and trailing spaces
    fun_strip = lambda x: x.decode().strip()
    if LooseVersion(pd.__version__) >= LooseVersion("2.1.0"):
        # if pandas version is 2.1.0 or above, we can use `df.map`
        df_var_list = df_var_list.map(fun_strip)
    else:
        # otherwise, we need to use `df.applymap`
        df_var_list = df_var_list.applymap(fun_strip)

    df_var_list_x = df_var_list.replace(r"^\s*$", np.nan, regex=True).dropna()
    var_dfm = df_var_list_x.set_index(["group", "var"])
    return var_dfm


# get variable info as a DataFrame
# save `var_df` for later use
var_df = get_output_info_df()

# dict as var_df but keys in lowercase
var_df_lower = {group.lower(): group for group in var_df.index.levels[0].str.strip()}

#  generate dict of functions to apply for each variable
dict_func_aggm = {
    "T": "last",
    "A": "mean",
    "S": "sum",
    "L": "last",
}
var_df["func"] = var_df.aggm.apply(lambda x: dict_func_aggm[x])

# dict of resampling ruls:
#  {group: {var: agg_method}}
dict_var_aggm = {
    group: var_df.loc[group, "func"].to_dict() for group in var_df.index.levels[0]
}


# generate index for variables in different model groups
def gen_group_cols(group_x):
    # get correct group name by cleaning and swapping case
    group = group_x.replace("dataoutline", "").replace("line", "")
    # print group
    group = var_df_lower[group]
    header_group = np.apply_along_axis(
        list, 0, var_df.loc[["datetime", group]].index.values
    )[:, 1]

    # generate MultiIndex if not `datetimeline`
    if not group_x == "datetimeline":
        index_group = pd.MultiIndex.from_product(
            [[group], header_group], names=["group", "var"], sortorder=None
        )
    else:
        index_group = header_group

    return index_group


# merge_grid: useful for both `dict_output` and `dict_state`
def pack_df_grid(dict_output):
    # pack all grid and times into index/columns
    df_xx = pd.DataFrame.from_dict(dict_output, orient="index")
    # pack
    df_xx0 = df_xx.map(pd.Series)
    df_xx1 = df_xx0.map(pd.DataFrame.from_dict)
    df_xx2 = pd.concat(
        {
            grid: pd.concat(df_xx1[grid].to_dict()).unstack().dropna(axis=1)
            for grid in df_xx1.columns
        }
    )
    # drop redundant levels
    df_xx2.columns = df_xx2.columns.droplevel(0)
    # regroup by `grid`
    df_xx2.index.names = ["grid", "time"]
    gb_xx2 = df_xx2.groupby(level="grid")
    # merge results of each grid
    xx3 = gb_xx2.agg(lambda x: tuple(x.values)).map(np.array)

    return xx3


# generate MultiIndex for variable groups
def gen_index(varline_x):
    var_x = varline_x.replace("dataout", "").replace("block", "").replace("line", "")
    group = var_df_lower[var_x]
    var = var_df.loc[group].index.tolist()
    mindex = pd.MultiIndex.from_product([[group], var], names=["group", "var"])
    return mindex


# generate one MultiIndex from a whole dict
def gen_MultiIndex(dict_x):
    x_keys = dict_x.keys()
    mindex = pd.concat([gen_index(k).to_frame() for k in x_keys]).index
    return mindex


# generate one Series from a dict entry
def gen_Series(dict_x, varline_x):
    m_index = gen_index(varline_x)
    res_Series = pd.Series(dict_x[varline_x], index=m_index)
    return res_Series


# merge a whole dict into one Series
def comb_gen_Series(dict_x):
    x_keys = dict_x.keys()
    res_Series = pd.concat([gen_Series(dict_x, k) for k in x_keys])
    return res_Series


# pack up output of `run_suews`
def pack_df_output_line(dict_output):
    import pickle

    pickle.dump(dict_output, open("dict_output.pkl", "wb"))
    print("dict_output saved to dict_output.pkl")
    # TODO: add output levels as in the Fortran version
    df_output = pd.DataFrame(dict_output).T
    # df_output = pd.concat(dict_output).to_frame().unstack()
    # set index level names
    index = df_output.index.set_names(["datetime", "grid"])
    # clean columns
    columns = gen_MultiIndex(df_output.iloc[0])
    values = np.apply_along_axis(np.hstack, 1, df_output.values)
    df_output = pd.DataFrame(values, index=index, columns=columns)
    return df_output


def pack_df_state(dict_state):
    df_state = pd.DataFrame(dict_state).T
    # df_state = pd.concat(dict_state).to_frame().unstack()
    # set index level names
    df_state.index = df_state.index.set_names(["datetime", "grid"])

    return df_state


def pack_df_output_array(dict_output_array, df_forcing):
    grid_list = list(dict_output_array.keys())
    grid_start = grid_list[0]
    col_df = gen_MultiIndex(dict_output_array[grid_start])
    dict_df = {}
    for grid in grid_list:
        array_grid = np.hstack([v[:, 5:] for v in dict_output_array[grid].values()])
        df_grid = pd.DataFrame(array_grid, columns=col_df, index=df_forcing.index)

        dict_df.update({grid: df_grid})

    # join results of all grids
    df_grid_res = pd.concat(dict_df, keys=dict_df.keys())

    # set index level names
    df_grid_res.index.set_names(["grid", "datetime"], inplace=True)

    return df_grid_res


def pack_df_output_block(dict_output_block, df_forcing_block):
    col_df = gen_MultiIndex(dict_output_block)
    val_df = np.hstack([ar[:, 5:] for ar in dict_output_block.values()])
    ind_df = df_forcing_block.index.rename("datetime")
    df_out = pd.DataFrame(val_df, columns=col_df, index=ind_df)

    return df_out


# resample supy output
def resample_output(df_output, freq="60T", dict_aggm=dict_var_aggm):
    # get grid and group names
    list_grid = df_output.index.get_level_values("grid").unique()
    list_group = df_output.columns.get_level_values("group").unique()

    # resampling output according to different rules defined in dict_aggm
    # note the setting in .resample: (closed='right',label='right')
    # which is to conform to SUEWS convention
    # that timestamp refer to the ending of previous period
    df_rsmp = pd.concat(
        {
            grid: pd.concat(
                {
                    group: df_output.loc[grid, group]
                    .resample(
                        freq,
                        closed="right",
                        label="right",
                    )
                    .agg(dict_aggm[group])
                    for group in list_group
                },
                axis=1,
                names=["group", "var"],
            )
            for grid in list_grid
        },
        names=["grid"],
    )

    # clean results
    df_rsmp = df_rsmp.dropna(how="all", axis=0)

    return df_rsmp


def proc_df_rsl(df_output, debug=False):
    try:
        # if we work on the whole output with multi-index columns
        df_rsl_raw = df_output["RSL"].copy()
    except:
        # if we directly work on the RSL output
        df_rsl_raw = df_output.copy()

    try:
        # drop unnecessary columns if existing
        df_rsl_data = df_rsl_raw.drop(["Year", "DOY", "Hour", "Min", "Dectime"], axis=1)
    except:
        df_rsl_data = df_rsl_raw

    # retrieve data for plotting
    df_rsl = df_rsl_data.iloc[:, : 30 * 4]
    df_rsl.columns = (
        df_rsl.columns.str.split("_")
        .map(lambda l: tuple([l[0], int(l[1])]))
        .rename(["var", "level"])
    )
    df_rsl_proc = df_rsl.stack()
    if debug:
        # retrieve debug variables
        df_rsl_debug = df_rsl_data.iloc[:, 120:]
        return df_rsl_proc, df_rsl_debug
    else:
        return df_rsl_proc


def is_numeric(obj):
    """
    Check if an object is numeric.

    Parameters:
    obj (object): The object to be checked.

    Returns:
    bool: True if the object is numeric, False otherwise.
    """
    if isinstance(obj, (int, float, complex)):
        return True
    if isinstance(obj, np.ndarray):
        return np.issubdtype(obj.dtype, np.number)
    return False


def inspect_dts_structure(dts_obj):
    """Inspect debug object structure once and cache the property paths"""
    props = {}
    for attr in dir(dts_obj):
        if not attr.startswith("_") and not callable(getattr(dts_obj, attr)):
            val = getattr(dts_obj, attr)
            # If it's a nested object, get its properties too
            if hasattr(val, "__dict__"):
                sub_props = [
                    sub_attr
                    for sub_attr in dir(val)
                    if not sub_attr.startswith("_")
                    and not callable(getattr(val, sub_attr))
                ]
                props[attr] = sub_props
            else:
                props[attr] = None
    return props


def fast_pack_dts(dts_obj, structure=None):
    """Pack debug object using cached structure or inspect if needed"""
    if structure is None:
        structure = inspect_dts_structure(dts_obj)

    result = {}
    for attr, sub_props in structure.items():
        val = getattr(dts_obj, attr)
        if sub_props:
            # Nested object
            result[attr] = {sub_prop: getattr(val, sub_prop) for sub_prop in sub_props}
        else:
            # Direct value
            result[attr] = val

    return result


def pack_dts_batch(dts_objs, structure):
    """Process multiple derived type objects at once"""
    return [fast_pack_dts(obj, structure) for obj in dts_objs]


def pack_dts_selective(dts_obj, needed_vars):
    """Pack only specified variables from derived type object

    Args:
        dts_obj: Fortran derived type object
        needed_vars: dict of {state_name: [var_names]} to extract
    """
    result = {}
    get = getattr

    for state, vars in needed_vars.items():
        state_obj = get(dts_obj, state)
        if vars:  # If specific variables are requested
            result[state] = {var: get(state_obj, var) for var in vars}
        else:  # If None, get all non-private attributes
            result[state] = get(state_obj)

    return result


def pack_dict_dts(dts):
    """Fast version that assumes well-formed DTS objects"""
    dict_dts = {}
    for attr in dir(dts):
        # Skip magic methods with one check
        if not attr.startswith("_"):
            val = getattr(dts, attr)
            # Most values should be numeric, so check that first
            if isinstance(val, (int, float, bool, np.ndarray)):
                dict_dts[attr] = val
            # Only recurse if not numeric and not callable
            elif not callable(val):
                dict_dts[attr] = pack_dict_dts(val)
    return dict_dts


def has_dict(d):
    """
    Check if a dictionary contains any values that are dictionaries.

    Parameters:
    d (dict): The dictionary to check.

    Returns:
    bool: True if the dictionary contains any values that are dictionaries, False otherwise.
    """
    return any(isinstance(v, dict) for v in d.values())


def pack_df_dts_raw(dict_dts):
    """
    Packs a dictionary of DTS-derived information into a pandas DataFrame.

    Args:
        dict_debug (dict): A dictionary containing debug information.

    Returns:
        pandas.DataFrame: A DataFrame containing the packed debug information.
    """

    # import pdb; pdb.set_trace()
    dict_df_dts = {}
    for k, v in dict_dts.items():
        if has_dict(v):
            dict_df_dts[k] = pack_df_dts_raw(v)
        else:
            dict_df_dts[k] = pd.Series(v)

    df_dts = pd.concat(dict_df_dts, axis=0)

    return df_dts


def pack_df_dts(dict_dts):
    """
    Packs the DTS-derived dictionary into a DataFrame.

    Args:
        dict_dts (dict): The DTS-derived dictionary containing the DTS-derived information.

    Returns:
        pandas.DataFrame: The packed DataFrame with DTS-derived information.

    """
    df_dts_raw = pack_df_dts_raw(dict_dts)
    df_dts_raw.index = df_dts_raw.index.rename(
        [
            "datetime",
            "grid",
            "step",
            "group",
            "var",
        ]
    )
    df_dts = (
        df_dts_raw.unstack(level=["group", "var"])
        .sort_index(level=0, axis=1)
        .dropna(axis=1, how="all")
    )
    return df_dts
