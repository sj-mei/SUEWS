#!/usr/bin/env python
"""
Test inline annotation - much better for users!
"""

import yaml
import os

print("INLINE ANNOTATION TEST")
print("="*60)

# Create a simple test config with missing parameters
yaml_content = """# Test Configuration
sites:
  - name: Test Site
    gridiv: 1
    properties:
      lat: {value: 51.5}
      lng: {value: -0.1}
      
      # Land cover - buildings missing parameters
      land_cover:
        bldgs:
          sfr: {value: 0.40}
          # Missing: bldgh, faibldg
          
        grass:
          sfr: {value: 0.30}
          
      # Missing conductance section entirely
      
      # CO2 parameters incomplete
      anthropogenic_emissions:
        co2: {}  # Empty but valid
"""

# Save test file
test_file = "test_inline.yml"
with open(test_file, 'w') as f:
    f.write(yaml_content)

print("1. Loading config with missing parameters...")
print("-"*60)

from supy.data_model import SUEWSConfig

config = SUEWSConfig.from_yaml(test_file)

print("-"*60)

# Generate inline annotated file
print("\n2. Generating inline annotated YAML...")
annotated_file = config.generate_annotated_yaml(test_file)

print(f"\n✓ Generated: {annotated_file}")

# Show the annotated content
print("\n3. Inline annotated content:")
print("="*60)

with open(annotated_file, 'r') as f:
    content = f.read()
    print(content)

print("\n" + "="*60)
print("BENEFITS OF INLINE ANNOTATION:")
print("✅ Parameters appear exactly where they need to be added")
print("✅ Users can simply uncomment and adjust values")
print("✅ No need to search through separate sections")
print("✅ Preserves original YAML structure")
print("✅ Much easier to fix issues!")

# Cleanup
for f in [test_file, annotated_file]:
    if os.path.exists(f):
        os.remove(f)