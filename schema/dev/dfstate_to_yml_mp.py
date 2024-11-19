import pandas as pd
import yaml
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap

IGNORECOLSDICT = {
    "Model Use": [
        'Driver',
        'NARP',
        'BEERS',
        'STEBBS',
        'waterdist',
        'anemsn',
        'RSL',
        'ESTM',
        'resist',
        'SPARTACUS-S',
        'atmmoisstab',
        'biogenCO2',
        'OHM',
        'LUMPS',
        'dailystate',
        'ANOHM',
        'EHC',
        'SNOW'
    ],
    "ModelTranslation": [
        'toSUEWS',
        'toSTEBBS',
        'toSPARTACUS-S'
    ],
    "Methods": [
        'CBLUse',
        'SnowUse',
        'NetRadiationMethod',
        'BaseTMethod',
        'EmissionsMethod',
        'StorageHeatMethod',
        'OHMIncQF',
        'StabilityMethod',
        'RoughLenHeatMethod',
        'RoughLenMomMethod',
        'SMDMethod',
        'WaterUseMethod',
        'DiagMethod',
        'stebbsmethod'
    ]
}

class ModelParameters:
    def __init__(self, csv_file, encoding='utf-8'):
        self.df = pd.read_csv(csv_file, encoding=encoding)
        self.df = self.remove_ignored_columns(ignoreColsdict=IGNORECOLSDICT)
        self.df = self.convert_to_multiindex()

    def remove_ignored_columns(self, ignoreColsdict):
        columns_to_ignore = [col for cols in ignoreColsdict.values() for col in cols]
        df_dropped = self.df.copy()
        df_dropped.drop(columns=columns_to_ignore, inplace=True, errors='ignore')
        return df_dropped

    def convert_to_multiindex(self):
        df_mi = self.df.copy()
        df_mi.set_index(['Model', 'Input Type', 'Category', 'Scale', 'SuPy Input'], inplace=True)
        return df_mi


    def filter_parameters(self, model=None, input_type=None, category=None, scale=None, contains=None):
        params = self.df.copy()
        if model is not None:
            params = params[self.df.index.get_level_values('Model') == model]
        if input_type is not None:
            params = params[params.index.get_level_values('Input Type') == input_type]
        if category is not None:
            params = params[params.index.get_level_values('Category') == category]
        if scale is not None:
            params = params[params.index.get_level_values('Scale') == scale]
        if contains is not None:
            params = params[params.index.get_level_values('SuPy Input').str.contains(contains)]
        
        self.df.drop(index=params.index, inplace=True)
        
        params = params.index.get_level_values('SuPy Input').tolist()
        params = {'stebbs_'+param: None for param in params}

        return params

class YamlEditor:
    def __init__(self, yaml_file_path):
        self.yaml_file_path = yaml_file_path
        self.config_data = None
        self.read_yaml_file()

    def read_yaml_file(self):
        with open(self.yaml_file_path, 'r') as file:
            self.config_data = yaml.safe_load(file)[0]

    def add_method(self, method, method_type, value=0):
        self.config_data['model'][method_type][method] = value

    def add_siteInfo(self, siteInfo:dict):
        def update_dict(original, update):
            if isinstance(original, dict):
                for key, value in update.items():
                    if isinstance(value, dict) and key in original:
                        update_dict(original[key], value)
                    elif isinstance(value, list):
                        for item in value:
                            update_dict(original[key], item)
                    else:
                        if key in original:
                            original[key].update(value)
                        else:
                            original[key] = value

        if 'site' not in self.config_data:
            self.config_data['site'] = {}
        if 'properties' not in self.config_data['site'][0]:
            self.config_data['site'][0]['properties'] = {}

        update_dict(self.config_data['site'][0]['properties'], siteInfo)

    def add_all_parameters(self, model_name, parameters):
        if 'site' not in self.config_data:
            self.config_data['site'] = {}
        if 'properties' not in self.config_data['site'][0]:
            self.config_data['site'][0]['properties'] = {}
        if model_name not in self.config_data['site'][0]['properties']:
            self.config_data['site'][0]['properties'][model_name] = {}
        self.config_data['site'][0]['properties'][model_name] = parameters
    
    def convert_to_commented_map(self):
        def convert_dict_to_commented_map(d):
            if isinstance(d, dict):
                commented_map = CommentedMap()
                for k, v in d.items():
                    commented_map[k] = convert_dict_to_commented_map(v)
                return commented_map
            elif isinstance(d, list):
                return [convert_dict_to_commented_map(item) for item in d]
            else:
                return d

        self.config_data = convert_dict_to_commented_map(self.config_data)
        self.config_data = CommentedMap(self.config_data)

    def write_yaml_file(self):
        with open('./schema/dev/config-suews_update.yml', 'w') as file:
            yaml.dump(self.config_data, file)
    
    def search_dict(self, dct, key, path=None):
        """
        Recursively search for a key in a dictionary that may contain sub-dictionaries.
        
        :param d: The dictionary to search.
        :param key: The key to search for.
        :param path: The current path to the key.
        :return: The path to the key as a list of keys, or None if the key is not found.
        """
        if path is None:
            path = []
        if key in dct:
            return path + [key]
        for k, v in dct.items():
            if isinstance(v, (dict, CommentedMap)):
                result = self.search_dict(v, key, path + [k])
                if result is not None:
                    return result
            elif isinstance(v, list):
                for i, item in enumerate(v):
                    if isinstance(item, (dict, CommentedMap)):
                        result = self.search_dict(item, key, path + [k, i])
                        if result is not None:
                            return result
        return None
    
    def make_comment(self, path, value):
        """
        Update the value of a key using the path to the key.
        
        :param path: The path to the key as a list of keys.
        :param value: The new value to set.
        """
        d = self.config_data
        for key in path[:-1]:
            d = d[key]
        d.yaml_add_eol_comment(value, path[-1])

    def write_yaml_file_with_comments(self, comments=None):
        yaml = YAML()
        yaml.indent(mapping=2, sequence=4, offset=2)
        self.convert_to_commented_map()
        if comments:
            for key, value in comments.items():
                key_path = self.search_dict(self.config_data, key)
                if key_path is not None:
                    self.make_comment(key_path, value)
                    
        with open('./schema/dev/config-suews_update.yml', 'w') as file:
            yaml.dump(self.config_data, file)


if __name__ == '__main__':
    # Read in the data in multi-index format
    supyParameters = ModelParameters(csv_file='./schema/dev/2024-11-suews_stebbs_translation(Inputs).csv', encoding='latin1')

    # Collect some parameters from the STEBBS model
    stebbs_methods = supyParameters.filter_parameters(
        model='STEBBS',
        input_type='Method'
    )
    stebbs_building_params = supyParameters.filter_parameters(
        model='STEBBS',
        input_type='Parameter',
        scale='Building'
    )
    stebbs_wallroof_params = supyParameters.filter_parameters(
        model='STEBBS',
        input_type='Parameter',
        scale='Facet',
        contains='Wall'
    )
    stebbs_window_params = supyParameters.filter_parameters(
        model='STEBBS',
        input_type='Parameter',
        scale='Facet',
        contains='Window'
    )
    stebbs_floor_params = supyParameters.filter_parameters(
        model='STEBBS',
        input_type='Parameter',
        scale='Facet',
        contains='Floor'
    )
    stebbs_qf_params = supyParameters.filter_parameters(
        model='STEBBS',
        category='QF'
    )
    stebbs_dhw_params = supyParameters.filter_parameters(
        model='STEBBS',
        category='QF DHW'
    )
    stebbs_temp_params = supyParameters.filter_parameters(
        model='STEBBS',
        category='QF Temp'
    )
    stebbs_population_params = supyParameters.filter_parameters(
        model='STEBBS',
        category='Population'
    )
    stebbs_appliance_params = supyParameters.filter_parameters(
        model='STEBBS',
        category='Appliance'
    )
    stebbs_dynamic_params = supyParameters.filter_parameters(
        model='STEBBS',
        category='Dynamic'
    )
    siteInfo_update = {'land_cover': {'bldgs':{
        'WallRoof': stebbs_wallroof_params,
        'Window': stebbs_window_params,
        'Floor': stebbs_floor_params,
        'QF': stebbs_qf_params,
        'QF DHW': stebbs_dhw_params,
        'QF Temp': stebbs_temp_params,
        'Population': stebbs_population_params,
        'Appliance': stebbs_appliance_params,
        'Dynamic': stebbs_dynamic_params
        }
    }}
    siteInfo_update['land_cover']['bldgs'].update(stebbs_building_params)

    # Edit the yaml file
    yamlEditor = YamlEditor(yaml_file_path='./schema/dev/config-suews.yml')
    for method in stebbs_methods:
        yamlEditor.add_method(method=method, method_type='physics')
    
    yamlEditor.add_siteInfo(siteInfo=siteInfo_update)
    yamlEditor.add_all_parameters(model_name='stebbs', parameters=supyParameters.filter_parameters(model='STEBBS'))
    comments = {
        'lat': 'Latitude',
        'lng': 'Longitude',
        'alt': 'Altitude',
        'site': 'Site information',
    }


    # Add missing columns:
    missing_columns = pd.read_csv('./df_state_test(in).csv')
    columns = missing_columns.columns.tolist()

    # Check if the columns are in the dataframe
    for col in columns:
        path = yamlEditor.search_dict(yamlEditor.config_data, col)
        if path is None:
            print(f'{col} not in the yaml file')


    yamlEditor.write_yaml_file_with_comments(comments=comments)