#!/usr/bin/env python
"""
Demonstrate how inline annotation helps users fix their config.

Shows a before/after comparison.
"""

print("INLINE ANNOTATION DEMONSTRATION")
print("="*60)

# Show original config with problems
original = """# Urban Config
sites:
  - name: City Centre
    gridiv: 1
    properties:
      lat: {value: 51.5}
      lng: {value: -0.1}
      land_cover:
        bldgs:
          sfr: {value: 0.45}
          # Missing critical parameters!
"""

print("\n1. ORIGINAL CONFIG (with missing parameters):")
print("-"*60)
print(original)

# Show annotated version
annotated = """# Urban Config
sites:
  - name: City Centre
    gridiv: 1
    properties:
      lat: {value: 51.5}
      lng: {value: -0.1}
      land_cover:
        bldgs:
          sfr: {value: 0.45}
          # ⚠️  MISSING: Building height required (fraction: 45.0%)
          # 💡 ADD:
          # bldgh: {value: 20.0}  # building height in meters
          
          # ⚠️  MISSING: Frontal area index needed for wind calculations  
          # 💡 ADD:
          # faibldg: {value: 0.5}  # frontal area index (0.1-0.7)
"""

print("\n2. ANNOTATED VERSION (with inline guidance):")
print("-"*60)
print(annotated)

# Show fixed version
fixed = """# Urban Config
sites:
  - name: City Centre
    gridiv: 1
    properties:
      lat: {value: 51.5}
      lng: {value: -0.1}
      land_cover:
        bldgs:
          sfr: {value: 0.45}
          bldgh: {value: 25.0}  # adjusted for my city
          faibldg: {value: 0.6}  # dense urban area
"""

print("\n3. FIXED VERSION (after uncommenting and adjusting):")
print("-"*60)
print(fixed)

print("\n" + "="*60)
print("BENEFITS:")
print("✅ Users see EXACTLY where to add parameters")
print("✅ No searching through documentation")
print("✅ Ready-to-use values with explanations")
print("✅ Just uncomment, adjust, and remove warnings")
print("✅ Much better than 'Missing parameter X' messages!")

print("\nWORKFLOW:")
print("1. Run validation → Get summary")
print("2. Generate annotated YAML")
print("3. Open file, search for ⚠️")
print("4. Uncomment 💡 ADD blocks")
print("5. Adjust values for your site")
print("6. Remove warning comments")
print("7. Re-validate to confirm")