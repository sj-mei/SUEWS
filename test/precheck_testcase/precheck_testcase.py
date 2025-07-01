from supy.data_model.core import run_precheck, SUEWSConfig
import os

# List of test YAMLs
yaml_cases = [
    "test/precheck_testcase/precheck_testcase1.yml",
    "test/precheck_testcase/precheck_testcase2.yml",
    "test/precheck_testcase/precheck_testcase3.yml",
]

for input_yaml in yaml_cases:
    print("\n" + "=" * 100)
    print(f" Running precheck on: {input_yaml}")

    # Step 1: Run precheck
    try:
        run_precheck(input_yaml)
        print(f" Precheck completed successfully for: {input_yaml}")
    except Exception as e:
        print(f" Precheck failed for {input_yaml} with error:\n{e}")
        continue  # Skip Pydantic step if precheck fails

    # Step 2: Load py0_ updated YAML for Pydantic validation
    output_yaml = os.path.join(
        os.path.dirname(input_yaml),
        f"py0_{os.path.basename(input_yaml)}"
    )

    print(f"\n Loading updated YAML for Pydantic validation: {output_yaml}")
    try:
        config = SUEWSConfig.from_yaml(output_yaml)
        print(f" Pydantic validation passed for: {output_yaml}")
    except Exception as e:
        print(f" Pydantic validation failed for {output_yaml} with error:\n{e}")

    print("=" * 100 + "\n")
