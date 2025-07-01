"""
YAML Annotator for SUEWS Configuration

This module provides functionality to generate annotated YAML files
with inline comments showing validation issues and how to fix them.
"""

import yaml
from typing import Dict, List, Any, Optional, Tuple
from pathlib import Path
import re
from collections import defaultdict
from io import StringIO


class ValidationIssue:
    """Represents a validation issue at a specific location in the config."""
    
    def __init__(self, path: str, param: str, message: str, fix: str, level: str = "WARNING"):
        self.path = path
        self.param = param
        self.message = message
        self.fix = fix
        self.level = level
    
    def __repr__(self):
        return f"ValidationIssue({self.path}, {self.param}, {self.level})"


class YAMLAnnotator:
    """Annotates YAML configuration files with validation feedback."""
    
    def __init__(self):
        self.issues: List[ValidationIssue] = []
        self.issues_by_path: Dict[str, List[ValidationIssue]] = defaultdict(list)
    
    def add_issue(self, path: str, param: str, message: str, fix: str, level: str = "WARNING"):
        """Add a validation issue."""
        issue = ValidationIssue(path, param, message, fix, level)
        self.issues.append(issue)
        # Normalize path for matching
        normalized_path = path.replace("Site ", "sites[").replace("]:", "]")
        self.issues_by_path[normalized_path].append(issue)
    
    def annotate_yaml_string(self, yaml_str: str, original_data: Dict[str, Any]) -> str:
        """
        Annotate a YAML string with validation issues.
        
        Args:
            yaml_str: The original YAML string
            original_data: The parsed YAML data
            
        Returns:
            Annotated YAML string with comments
        """
        lines = yaml_str.split('\n')
        annotated_lines = []
        current_path = []
        indent_stack = [0]
        
        for i, line in enumerate(lines):
            # Skip empty lines and comments
            if not line.strip() or line.strip().startswith('#'):
                annotated_lines.append(line)
                continue
            
            # Calculate indentation
            indent = len(line) - len(line.lstrip())
            
            # Update path based on indentation
            while indent_stack and indent <= indent_stack[-1]:
                indent_stack.pop()
                if current_path:
                    current_path.pop()
            
            # Parse the line to get key
            if ':' in line and not line.strip().startswith('-'):
                key = line.split(':')[0].strip()
                current_path.append(key)
                indent_stack.append(indent)
            elif line.strip().startswith('- '):
                # Handle list items
                if current_path and current_path[-1] == 'sites':
                    # This is a site entry
                    current_path.append('[0]')  # Simplified - would need proper indexing
                    indent_stack.append(indent)
            
            # Add the original line
            annotated_lines.append(line)
            
            # Check if we have issues for this path
            path_str = '.'.join(current_path)
            matching_issues = self._find_issues_for_path(path_str, line)
            
            if matching_issues:
                # Add comment lines for each issue
                for issue in matching_issues:
                    comment_indent = ' ' * (indent + 2)
                    annotated_lines.append(f"{comment_indent}# ⚠️  {issue.level}: {issue.message}")
                    annotated_lines.append(f"{comment_indent}# 💡 FIX: {issue.fix}")
                    if issue.param not in line:
                        # Add example fix
                        annotated_lines.append(f"{comment_indent}# EXAMPLE:")
                        annotated_lines.append(f"{comment_indent}{issue.param}: {{value: <INSERT_VALUE>}}")
        
        return '\n'.join(annotated_lines)
    
    def _find_issues_for_path(self, current_path: str, line: str) -> List[ValidationIssue]:
        """Find issues matching the current path."""
        matching_issues = []
        
        # Direct path matching
        for path, issues in self.issues_by_path.items():
            if path in current_path or current_path in path:
                for issue in issues:
                    # Check if this is the right location for the issue
                    if issue.param in line or self._should_add_param_here(current_path, issue):
                        matching_issues.append(issue)
        
        return matching_issues
    
    def _should_add_param_here(self, current_path: str, issue: ValidationIssue) -> bool:
        """Determine if a missing parameter should be added at this location."""
        # This is a simplified check - in practice would need more sophisticated matching
        issue_path_parts = issue.path.lower().split('/')
        current_path_parts = current_path.lower().split('.')
        
        # Check if we're in the right section
        for part in issue_path_parts:
            if part in current_path_parts:
                return True
        
        return False
    
    def generate_annotated_file(
        self, 
        input_path: Path, 
        output_path: Optional[Path] = None,
        validation_report: Optional[Any] = None
    ) -> Path:
        """
        Generate an annotated YAML file with validation feedback.
        
        Args:
            input_path: Path to the original YAML file
            output_path: Path for the annotated file (default: input_path.stem + '_annotated.yml')
            validation_report: Optional validation report to extract issues from
            
        Returns:
            Path to the generated annotated file
        """
        # Read original file
        with open(input_path, 'r') as f:
            yaml_content = f.read()
            
        # Parse YAML
        data = yaml.safe_load(yaml_content)
        
        # If validation report provided, extract issues
        if validation_report:
            self._extract_issues_from_report(validation_report)
        
        # Create a simpler annotation approach
        annotated_content = self._create_annotated_yaml(yaml_content, data)
        
        # Determine output path
        if output_path is None:
            output_path = input_path.parent / f"{input_path.stem}_annotated.yml"
        
        # Write annotated file
        with open(output_path, 'w') as f:
            f.write(annotated_content)
        
        return output_path
    
    def _create_annotated_yaml(self, yaml_content: str, data: Dict[str, Any]) -> str:
        """Create annotated YAML with a simpler approach."""
        header = """# ANNOTATED SUEWS CONFIGURATION
# ================================
# This file has been annotated with validation feedback.
# Look for:
#   ⚠️  WARNING: Validation issues that should be addressed
#   💡 FIX: Suggested solutions
#   EXAMPLE: Sample values you can use
#
# To use this file:
# 1. Review all WARNING comments
# 2. Add missing parameters as suggested
# 3. Replace <INSERT_VALUE> with appropriate values
# 4. Remove comment lines after fixing
# ================================

"""
        
        # Start with the header
        result = header
        
        # Add the original YAML
        result += yaml_content + "\n\n"
        
        # Add a section with all validation issues organized by site
        result += "\n# ================================\n"
        result += "# VALIDATION ISSUES TO ADDRESS:\n"
        result += "# ================================\n\n"
        
        # Group issues by site
        sites_issues = defaultdict(list)
        for issue in self.issues:
            # Extract site from path
            if 'sites[' in issue.path:
                site_match = re.search(r'sites\[(\d+)\]', issue.path)
                if site_match:
                    site_idx = int(site_match.group(1))
                    sites_issues[site_idx].append(issue)
        
        # Add issues for each site
        for site_idx in sorted(sites_issues.keys()):
            site_name = data['sites'][site_idx].get('name', f'Site {site_idx + 1}')
            result += f"# {site_name}:\n"
            result += f"# {'='*len(site_name)}:\n"
            
            # Group by issue type
            by_type = defaultdict(list)
            for issue in sites_issues[site_idx]:
                # Extract the component (e.g., bldgs, conductance, etc.)
                if '/land_cover/' in issue.path:
                    component = issue.path.split('/land_cover/')[-1].split('/')[0]
                    by_type[f"land_cover.{component}"].append(issue)
                elif '/conductance' in issue.path:
                    by_type['conductance'].append(issue)
                elif '/anthropogenic_emissions/co2' in issue.path:
                    by_type['anthropogenic_emissions.co2'].append(issue)
                else:
                    by_type['other'].append(issue)
            
            # Add issues by type
            for component, issues in sorted(by_type.items()):
                result += f"#\n# {component}:\n"
                for issue in issues:
                    result += f"#   ⚠️  {issue.level}: {issue.message}\n"
                    result += f"#   💡 FIX: {issue.fix}\n"
                    
                    # Add example based on parameter
                    if issue.param == 'bldgh':
                        result += f"#   EXAMPLE: bldgh: {{value: 20.0}}  # typical urban building height\n"
                    elif issue.param == 'faibldg':
                        result += f"#   EXAMPLE: faibldg: {{value: 0.5}}  # frontal area index\n"
                    elif issue.param == 'thermal_layers':
                        result += f"#   EXAMPLE:\n"
                        result += f"#     thermal_layers:\n"
                        result += f"#       dz: {{value: [0.2, 0.1, 0.1, 0.5, 1.6]}}  # layer thickness (m)\n"
                        result += f"#       k: {{value: [1.2, 1.1, 1.1, 1.5, 1.6]}}  # conductivity (W/m/K)\n"
                        result += f"#       rho_cp: {{value: [1.2e6, 1.1e6, 1.1e6, 1.5e6, 1.6e6]}}  # heat capacity (J/m3/K)\n"
                    elif 'g_' in issue.param or issue.param in ['s1', 's2', 'kmax']:
                        result += f"#   EXAMPLE: {issue.param}: {{value: <see documentation>}}\n"
                    elif 'co2' in issue.param.lower():
                        result += f"#   EXAMPLE: {issue.param}: {{value: <see documentation>}}\n"
                    result += "#\n"
            
            result += "\n"
        
        return result
    
    def _extract_issues_from_report(self, report):
        """Extract issues from a validation report."""
        # This would extract issues from the validation report
        # Implementation depends on the report format
        pass


def create_example_annotated_yaml():
    """Create an example annotated YAML file to demonstrate the functionality."""
    
    example_yaml = """sites:
  - name: Urban Test Site
    gridiv: 1
    properties:
      lat: {value: 51.5}
      lng: {value: -0.1}
      land_cover:
        bldgs:
          sfr: {value: 0.4}
"""
    
    # Create annotator and add issues
    annotator = YAMLAnnotator()
    
    # Add building height issue
    annotator.add_issue(
        path="sites[0]/properties/land_cover/bldgs",
        param="bldgh",
        message="Building height is required (building fraction: 40%)",
        fix="Add building height in meters",
        level="WARNING"
    )
    
    # Add frontal area issue
    annotator.add_issue(
        path="sites[0]/properties/land_cover/bldgs",
        param="faibldg",
        message="Frontal area index recommended for accurate aerodynamics",
        fix="Add frontal area index (typical range: 0.1-0.7)",
        level="WARNING"
    )
    
    # Add thermal layers issue
    annotator.add_issue(
        path="sites[0]/properties/land_cover/bldgs",
        param="thermal_layers",
        message="Thermal layer properties required for heat storage calculations",
        fix="Add thermal_layers with dz, k, and rho_cp arrays",
        level="WARNING"
    )
    
    # Add conductance issue
    annotator.add_issue(
        path="sites[0]/properties",
        param="conductance",
        message="Missing surface conductance parameters for evapotranspiration",
        fix="Add conductance block with g_max, g_k, g_sm, s1, s2 parameters",
        level="WARNING"
    )
    
    # Generate annotated YAML
    annotated = annotator.annotate_yaml_string(example_yaml, yaml.safe_load(example_yaml))
    
    # Add header
    header = """# ANNOTATED SUEWS CONFIGURATION
# ================================
# This file shows validation issues found in your configuration.
# Look for:
#   ⚠️  WARNING: Issues that should be addressed
#   💡 FIX: How to fix each issue
#   EXAMPLE: Sample parameter blocks you can use
# ================================

"""
    
    return header + annotated


if __name__ == "__main__":
    # Generate example
    example = create_example_annotated_yaml()
    print(example)