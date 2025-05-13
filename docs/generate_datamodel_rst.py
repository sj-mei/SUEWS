import inspect
import importlib
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Type, Union, get_args, get_origin, Literal

from pydantic import BaseModel
from pydantic.fields import FieldInfo

# Add the project root to sys.path to allow importing supy
# Assuming this script is in 'docs/generate_datamodel_rst.py'
# And supy is in 'src/supy'
PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC_PATH = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC_PATH))

# Try to import the data models from supy
# This relies on supy.data_model being an importable path
# and __init__.py files being set up correctly.
try:
    import supy.data_model
except ImportError as e:
    print(f"Error importing supy.data_model: {e}")
    print("Please ensure that 'src' is in the Python path and supy.data_model is accessible.")
    sys.exit(1)


# --- Helper function to get user-friendly type names ---
def get_user_friendly_type_name(type_hint: Any) -> str:
    """Generates a user-friendly string representation of a type hint."""
    origin = get_origin(type_hint)
    args = get_args(type_hint)

    if origin is Union:
        # Handle Optional[X] as X (Optional)
        if len(args) == 2 and type(None) in args:
            non_none_arg = next(arg for arg in args if arg is not type(None))
            return f"{get_user_friendly_type_name(non_none_arg)} (Optional)"
        return " | ".join(get_user_friendly_type_name(arg) for arg in args)
    if origin is list or origin is List:
        return f"List of {get_user_friendly_type_name(args[0])}" if args else "List"
    if origin is dict or origin is Dict:
        if args and len(args) == 2:
            return f"Mapping from {get_user_friendly_type_name(args[0])} to {get_user_friendly_type_name(args[1])}"
        return "Mapping"
    if hasattr(type_hint, "__name__"):
        # Check if it's a ValueWithDOI wrapper
        if type_hint.__name__ == "ValueWithDOI" and args:
             return f"Value (type: {get_user_friendly_type_name(args[0])}) with DOI/Reference"
        return type_hint.__name__
    return str(type_hint)


# --- Main RST Generation Logic ---
def generate_rst_for_model(
    model_class: Type[BaseModel],
    output_dir: Path,
    processed_models: set[Type[BaseModel]],
    all_supy_models: Dict[str, Type[BaseModel]]
) -> None:
    """
    Generates an .rst file for a given Pydantic model, focusing on user configuration.
    """
    if model_class in processed_models:
        return
    processed_models.add(model_class)

    model_name = model_class.__name__
    rst_content = []

    # Title
    rst_content.append(model_name.replace("_", " ").title())
    rst_content.append("=" * len(rst_content[-1]))
    rst_content.append("")

    # Model Docstring (if any)
    if model_class.__doc__:
        rst_content.append(inspect.cleandoc(model_class.__doc__))
        rst_content.append("")

    rst_content.append("**Parameters:**")
    rst_content.append("")

    for field_name, field_info in model_class.model_fields.items():
        field_type_hint = field_info.annotation
        user_type_name = get_user_friendly_type_name(field_type_hint)

        # Start option block
        rst_content.append(f".. option:: {field_name} <{user_type_name}>")
        rst_content.append("")

        # Description
        description = getattr(field_info, 'description', None)
        if description:
            rst_content.append(f"   {description.strip()}")
            rst_content.append("")

        # Unit (Placeholder - requires units in Pydantic model or parsing from description)
        # Example: unit = field_info.json_schema_extra.get('unit') if field_info.json_schema_extra else None
        unit = None # Placeholder
        if hasattr(field_info, 'json_schema_extra') and isinstance(field_info.json_schema_extra, dict):
            unit = field_info.json_schema_extra.get('unit')
        if unit:
            rst_content.append(f"   :Unit: {unit}")
        else:
            # Try to parse from description, e.g., "Some value [unit]"
            if description and '[' in description and ']' in description:
                try:
                    parsed_unit = description[description.rfind('[')+1:description.rfind(']')]
                    if len(parsed_unit) < 10 and not ' ' in parsed_unit: # Basic sanity check
                         rst_content.append(f"   :Unit: {parsed_unit}")
                except:
                    pass # Ignore parsing errors

        # Default value
        default_value = "Not specified"
        if field_info.default is not None and field_info.default != inspect.Parameter.empty:
            default_value = f"``{field_info.default!r}``"
        elif field_info.default_factory is not None:
            try:
                default_value = f"``{field_info.default_factory()!r}`` (generated)"
            except Exception:
                default_value = "Dynamically generated"

        rst_content.append(f"   :Default: {default_value}")

        # Constraints
        constraints_desc = []
        # Standard Pydantic v2 constraint attributes
        constraint_attrs_map = {
            'gt': '>',
            'ge': '>=',
            'lt': '<',
            'le': '<=',
            'min_length': 'Minimum length',
            'max_length': 'Maximum length',
            'multiple_of': 'Must be a multiple of',
            'pattern': 'Must match regex pattern',
        }

        for attr, desc_prefix in constraint_attrs_map.items():
            if hasattr(field_info, attr):
                value = getattr(field_info, attr)
                if value is not None: # Ensure the constraint is actually set
                    constraints_desc.append(f"{desc_prefix}: ``{value!r}``")

        # Check for enum/Literal constraints from the type hint itself
        origin_type = get_origin(field_type_hint)
        args = get_args(field_type_hint)
        if origin_type is Literal:
            constraints_desc.append(f"Allowed values: {', '.join(f'``{arg!r}``' for arg in args)}")
        elif (origin_type is Union and
              any(get_origin(arg) is Literal for arg in args)):
            # Handle Union of Literals, e.g. Optional[Literal['a', 'b']]
            literal_args_combined = []
            for union_arg in args:
                if get_origin(union_arg) is Literal:
                    literal_args_combined.extend(get_args(union_arg))
            if literal_args_combined:
                 constraints_desc.append(f"Allowed values: {', '.join(f'``{arg!r}``' for arg in set(literal_args_combined))}")

        if constraints_desc:
            rst_content.append(f"   :Constraints: {'; '.join(constraints_desc)}")


        # Link to nested models
        # Check if the raw type or any type argument is a Pydantic model we know
        possible_model_types = [field_type_hint] + list(get_args(field_type_hint))
        nested_model_to_document = None
        for pt in possible_model_types:
            origin_pt = get_origin(pt) or pt # get actual type if it's a generic alias
            if hasattr(origin_pt, '__name__') and origin_pt.__name__ in all_supy_models and issubclass(origin_pt, BaseModel) and origin_pt != model_class:
                nested_model_to_document = origin_pt
                break

        if nested_model_to_document:
            nested_model_name = nested_model_to_document.__name__
            rst_content.append("")
            rst_content.append(f"   Details for this parameter group can be found in the :doc:`{nested_model_name.lower()}` section.")
            # Recursively generate RST for the nested model
            generate_rst_for_model(nested_model_to_document, output_dir, processed_models, all_supy_models)

        rst_content.append("") # Blank line after each option

    # Write to file
    rst_file_path = output_dir / f"{model_name.lower()}.rst"
    with open(rst_file_path, "w") as f:
        f.write("\n".join(rst_content))
    print(f"Generated: {rst_file_path}")


def get_all_models_in_module(module) -> Dict[str, Type[BaseModel]]:
    """
    Inspects a module and returns a dictionary of Pydantic models.
    """
    models = {}
    for name, obj in inspect.getmembers(module):
        if inspect.isclass(obj) and issubclass(obj, BaseModel) and obj.__module__ == module.__name__:
            models[name] = obj
    return models

def main():
    output_dir_name = "yaml_input"
    docs_source_path = PROJECT_ROOT / "docs" / "source"
    output_dir = docs_source_path / output_dir_name
    output_dir.mkdir(exist_ok=True)

    # Clean up previously generated RST files in the output directory
    print(f"Cleaning up previously generated .rst files in {output_dir}...")
    for rst_file in output_dir.glob("*.rst"):
        try:
            rst_file.unlink()
            print(f"  Deleted {rst_file}")
        except OSError as e:
            print(f"  Error deleting {rst_file}: {e}")
    print("Cleanup complete.")
    print("") # Add a blank line for better log readability

    processed_models = set()
    all_supy_data_models = {}

    # Discover models in supy.data_model and its submodules
    data_model_module_root = Path(supy.data_model.__file__).parent
    for py_file in data_model_module_root.glob("*.py"):
        if py_file.name == "__init__.py":
            module_name_to_import = "supy.data_model"
        else:
            module_name_to_import = f"supy.data_model.{py_file.stem}"

        try:
            module = importlib.import_module(module_name_to_import)
            all_supy_data_models.update(get_all_models_in_module(module))
        except ImportError as e:
            print(f"Could not import {module_name_to_import}: {e}")
            continue

    # Define top-level models for documentation (these will be the starting points)
    # Users should provide a list of "main" configuration models
    # For now, let's try with a few common ones if they exist
    top_level_model_names = ["Site", "SiteProperties", "LandCover", "SnowParams", "InitialStates", "Model"] # Add more as needed

    models_to_process = []
    for name in top_level_model_names:
        if name in all_supy_data_models:
            models_to_process.append(all_supy_data_models[name])
        else:
            print(f"Warning: Top-level model '{name}' not found in supy.data_model.")

    if not models_to_process:
        print("No top-level models found to process. Exiting.")
        print(f"Available models: {list(all_supy_data_models.keys())}")
        return

    for model_class in models_to_process:
        generate_rst_for_model(model_class, output_dir, processed_models, all_supy_data_models)

    # Create an index.rst for the generated files
    index_rst_path = output_dir / "yaml_input.rst"
    with open(index_rst_path, "w") as f:
        f.write(".. _yaml_input:\n\n")
        f.write("User Input Configuration Reference\n")
        f.write("===================================\n\n")
        f.write("This section provides a detailed reference for all parameters required to configure the SUEWS model.\n\n")
        f.write(".. toctree::\n")
        f.write("   :maxdepth: 2\n")
        f.write("   :caption: Configuration Sections:\n\n")
        # Sort files for consistent toctree
        generated_files = sorted([p.stem for p in output_dir.glob("*.rst") if p.name != index_rst_path.name])
        for rst_file_stem in generated_files:
            f.write(f"   {rst_file_stem}\n")
    print(f"Generated index: {index_rst_path}")
    print(
        f"\nMake sure to add '{output_dir_name}/{index_rst_path.name}' to your main Sphinx toctree in 'docs/source/index.rst' or other relevant files."
    )

if __name__ == "__main__":
    main()