"""
Conditional Validation for SUEWS Configuration

This module provides conditional validation functionality that can be used
independently of the main data_model package structure.
"""

from typing import List, Dict, Any, Optional, Union, Set
import warnings
import logging

logger = logging.getLogger(__name__)


class ValidationResult:
    """Simple validation result class."""
    def __init__(self, status: str = 'passed', errors: List[str] = None, 
                 warnings: List[str] = None, skipped: List[str] = None,
                 validated_methods: Set[str] = None, skipped_methods: Set[str] = None):
        self.status = status
        self.errors = errors or []
        self.warnings = warnings or []
        self.skipped = skipped or []
        self.validated_methods = validated_methods or set()
        self.skipped_methods = skipped_methods or set()


def analyze_config_methods(config_data: Dict[str, Any]) -> Dict[str, bool]:
    """
    Analyze configuration to determine which methods are enabled.
    
    Args:
        config_data: Configuration dictionary
        
    Returns:
        Dictionary of method flags
    """
    try:
        physics = config_data.get('model', {}).get('physics', {})
        
        # Extract method values
        diagmethod_val = physics.get('diagmethod', 0)
        if isinstance(diagmethod_val, dict):
            diagmethod_val = diagmethod_val.get('value', 0)
        
        roughmethod_val = physics.get('roughlenmommethod', 1)
        if isinstance(roughmethod_val, dict):
            roughmethod_val = roughmethod_val.get('value', 1)
        
        netrad_val = physics.get('netradiationmethod', 0)
        if isinstance(netrad_val, dict):
            netrad_val = netrad_val.get('value', 0)
        
        emissions_val = physics.get('emissionsmethod', 0)
        if isinstance(emissions_val, dict):
            emissions_val = emissions_val.get('value', 0)
        
        storage_val = physics.get('storageheatmethod', 0)
        if isinstance(storage_val, dict):
            storage_val = storage_val.get('value', 0)
        
        # Determine active methods
        methods = {
            'diagmethod_most': diagmethod_val == 0,
            'diagmethod_rst': diagmethod_val == 1,
            'diagmethod_variable': diagmethod_val == 2,
            'roughness_variable': roughmethod_val == 2,
            'netradiation_spartacus': netrad_val >= 1000,
            'emissions_advanced': emissions_val >= 4,
            'storage_estm': storage_val in [4, 5],
        }
        
        # Variable diagmethod uses both MOST and RST
        if methods['diagmethod_variable']:
            methods['diagmethod_most'] = True
            methods['diagmethod_rst'] = True
        
        return methods
        
    except Exception as e:
        logger.warning(f"Error analyzing config methods: {e}")
        return {}


def validate_most_parameters(sites: List[Dict], site_indices: Optional[List[int]] = None) -> List[str]:
    """Validate MOST-specific parameters."""
    errors = []
    indices = site_indices or range(len(sites))
    
    for i in indices:
        if i >= len(sites):
            continue
            
        site_props = sites[i].get('properties', {})
        
        # Extract values (handle RefValue objects)
        z0m = _extract_value(site_props.get('z0m_in'))
        if z0m is not None and (z0m <= 0 or z0m > 10):
            errors.append(f"[MOST] site[{i}].z0m_in={z0m} invalid (must be 0 < z0m ≤ 10)")
        
        zdm = _extract_value(site_props.get('zdm_in'))
        if zdm is not None and (zdm < 0 or zdm > 50):
            errors.append(f"[MOST] site[{i}].zdm_in={zdm} invalid (must be 0 ≤ zdm ≤ 50)")
    
    return errors


def validate_rst_parameters(sites: List[Dict], site_indices: Optional[List[int]] = None) -> List[str]:
    """Validate RST-specific parameters."""
    errors = []
    indices = site_indices or range(len(sites))
    
    for i in indices:
        if i >= len(sites):
            continue
            
        site_props = sites[i].get('properties', {})
        
        bldg_height = _extract_value(site_props.get('building_height'))
        if bldg_height is not None and (bldg_height <= 0 or bldg_height > 200):
            errors.append(f"[RST] site[{i}].building_height={bldg_height} invalid (must be 0 < height ≤ 200)")
        
        bldg_std = _extract_value(site_props.get('building_height_std'))
        if bldg_std is not None and (bldg_std < 0 or bldg_std > 100):
            errors.append(f"[RST] site[{i}].building_height_std={bldg_std} invalid (must be 0 ≤ std ≤ 100)")
    
    return errors


def validate_variable_roughness_parameters(sites: List[Dict], site_indices: Optional[List[int]] = None) -> List[str]:
    """Validate variable roughness parameters."""
    errors = []
    indices = site_indices or range(len(sites))
    
    for i in indices:
        if i >= len(sites):
            continue
            
        site_props = sites[i].get('properties', {})
        
        tree_h_ev = _extract_value(site_props.get('tree_height_evergreen'))
        if tree_h_ev is not None and (tree_h_ev <= 0 or tree_h_ev > 50):
            errors.append(f"[VAR_ROUGH] site[{i}].tree_height_evergreen={tree_h_ev} invalid (must be 0 < height ≤ 50)")
        
        tree_h_dec = _extract_value(site_props.get('tree_height_deciduous'))
        if tree_h_dec is not None and (tree_h_dec <= 0 or tree_h_dec > 50):
            errors.append(f"[VAR_ROUGH] site[{i}].tree_height_deciduous={tree_h_dec} invalid (must be 0 < height ≤ 50)")
    
    return errors


def _extract_value(param: Any) -> Optional[float]:
    """Extract numeric value from parameter (handles RefValue objects)."""
    if param is None:
        return None
    
    if isinstance(param, dict):
        return param.get('value')
    
    if hasattr(param, 'value'):
        return param.value
    
    try:
        return float(param)
    except (ValueError, TypeError):
        return None


def validate_suews_config_conditional(config_data: Dict[str, Any], strict: bool = False, verbose: bool = False) -> ValidationResult:
    """
    Main conditional validation function.
    
    Args:
        config_data: Configuration dictionary
        strict: Whether to treat warnings as errors
        verbose: Whether to print validation progress
        
    Returns:
        ValidationResult object
    """
    if verbose:
        print("=== SUEWS Conditional Validation ===")
    
    # Analyze which methods are enabled
    methods = analyze_config_methods(config_data)
    
    if verbose:
        enabled = [k for k, v in methods.items() if v]
        print(f"Active methods: {', '.join(enabled) if enabled else 'None'}")
    
    errors = []
    skipped = []
    validated_methods = set()
    skipped_methods = set()
    
    sites = config_data.get('site', [])
    
    # Conditional validation based on enabled methods
    if methods.get('diagmethod_most') and not methods.get('diagmethod_variable'):
        if verbose:
            print("DiagMethod MOST: validating MOST parameters")
        errors.extend(validate_most_parameters(sites))
        validated_methods.add('MOST')
        skipped_methods.add('RST')  # RST is skipped when only MOST is enabled
        skipped_methods.add('VARIABLE_ROUGHNESS')
    elif not methods.get('diagmethod_most'):
        skipped.append("MOST parameter validation (DiagMethod != MOST)")
        skipped_methods.add('MOST')
    
    if methods.get('diagmethod_rst'):
        if verbose:
            print("DiagMethod RST: validating RST parameters")
        errors.extend(validate_rst_parameters(sites))
        validated_methods.add('RST')
        if not methods.get('diagmethod_variable'):
            skipped_methods.add('MOST')  # MOST is skipped when only RST is enabled
    else:
        skipped.append("RST parameter validation (DiagMethod != RST)")
        skipped_methods.add('RST')
    
    if methods.get('diagmethod_variable'):
        if verbose:
            print("DiagMethod VARIABLE: validating both MOST and RST")
        # Variable method validates both MOST and RST
        errors.extend(validate_most_parameters(sites))
        errors.extend(validate_rst_parameters(sites))
        validated_methods.add('VARIABLE')
        validated_methods.add('MOST')
        validated_methods.add('RST')
    
    if methods.get('roughness_variable'):
        if verbose:
            print("Variable roughness: validating vegetation parameters")
        errors.extend(validate_variable_roughness_parameters(sites))
        validated_methods.add('VARIABLE_ROUGHNESS')
    else:
        skipped.append("Variable roughness validation (RoughnessMethod != VARIABLE)")
        skipped_methods.add('VARIABLE_ROUGHNESS')
    
    # Determine status
    status = 'failed' if errors else 'passed'
    
    result = ValidationResult(
        status=status,
        errors=errors,
        warnings=[],
        skipped=skipped,
        validated_methods=validated_methods,
        skipped_methods=skipped_methods
    )
    
    if verbose:
        print(f"\n{'❌ FAILED' if errors else '✅ PASSED'}: {len(errors)} errors, {len(skipped)} skipped")
        if errors:
            for error in errors:
                print(f"  • {error}")
        if skipped:
            print(f"Skipped: {', '.join(skipped)}")
    
    return result


# Integration functions for SUEWSConfig
def enhanced_from_yaml_validation(config_data: Dict[str, Any], strict: bool = True) -> Optional[Dict[str, Any]]:
    """Enhanced YAML validation for SUEWSConfig.from_yaml integration."""
    result = validate_suews_config_conditional(config_data, verbose=True)
    
    if result['errors']:
        error_msg = f"SUEWS Configuration Validation Failed: {len(result['errors'])} errors\n"
        error_msg += "\n".join(f"  - {err}" for err in result['errors'])
        
        if strict:
            raise ValueError(error_msg)
        else:
            warnings.warn(f"Validation errors found: {error_msg}")
    
    return result


def enhanced_to_df_state_validation(config_data: Dict[str, Any], strict: bool = False) -> Optional[Dict[str, Any]]:
    """Enhanced validation for to_df_state integration."""
    result = validate_suews_config_conditional(config_data, verbose=False)
    
    if result['errors']:
        error_msg = f"Configuration validation found {len(result['errors'])} issues before to_df_state conversion"
        
        if strict:
            detailed_msg = f"\n{error_msg}:\n" + "\n".join(f"  - {err}" for err in result['errors'])
            raise ValueError(detailed_msg)
        else:
            warnings.warn(f"{error_msg}. Continuing with conversion using available valid data.")
    
    return result