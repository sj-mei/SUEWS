from pathlib import Path
from supy.data_model.core import SUEWSConfig
import json

def main():
    # Generate JSON Schema from the canonical SUEWSConfig model
    schema = SUEWSConfig.model_json_schema()
    path_schema = Path(__file__) / "schema.json"
    path_schema.write_text(json.dumps(schema, indent=2, ensure_ascii=False))
    print(f"Schema written to {path_schema}")

if __name__ == "__main__":
    main()