from schema.dev.def_config_suews import *

config = SUEWSConfig()
# Convert model dump to YAML format
with open('./schema/dev/config-suews-default.yml', 'w') as file:
    yaml.dump(config.model_dump(exclude_none=True), file, sort_keys=False, allow_unicode=True)