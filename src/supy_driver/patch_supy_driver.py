import subprocess
import sys
import os
from pathlib import Path


def main():
    fn_this_script = sys.argv[0]
    current_source_dir = Path(fn_this_script).parent

    fn_supy_driver = Path(sys.argv[1]).name  # path to f2py or f2py part of numpy
    p_supy_driver =Path.cwd() / fn_supy_driver

    output_dir = sys.argv[2]


    print(f"path to supy_driver: {p_supy_driver}")

    # patch supy_driver.py
    with open(p_supy_driver, "r") as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith("import _supy_driver"):
                lines[i] = """
try:
    from . import _supy_driver
except ImportError:
    try:
        import _supy_driver
    except ImportError:
        raise ImportError("Cannot import _supy_driver")


"""
                break

    # write back to supy_driver.py
    with open(p_supy_driver, "w") as f:
        f.writelines(lines)






    # Move generated files to the output directory
    try:
        subprocess.check_call(
            [
                "python",
                os.path.join(current_source_dir, "move_output_gen.py"),
                fn_supy_driver,
                output_dir,
            ]
        )
        print("Output files moved successfully")
    except subprocess.CalledProcessError as e:
        print(f"Error moving output files: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
