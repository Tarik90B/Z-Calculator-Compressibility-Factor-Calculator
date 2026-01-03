"""
Z – Compressibility Factor Calculator
Copyright (c) 2026 Tarik

Independent implementation of AGA8 according to ISO 12213.
Licensed under the MIT License.
"""

import sys
from pathlib import Path

def resource_path(filename, writable=False):
    if writable:
        if hasattr(sys, "_MEIPASS"):
            # EXE: write next to EXE
            base_path = Path(sys.executable).parent
            print("1:", base_path)
        else:
            # Python: write to cwd
            base_path = Path.cwd()
            print("2:", base_path)
    else:
        if hasattr(sys, "_MEIPASS"):
            # EXE: read from bundled _MEIPASS
            base_path = Path(sys._MEIPASS)
            print("3:", base_path)
        else:
            # Python: read from project folder
            base_path = Path(__file__).resolve().parent
            print("4:", base_path)
    return base_path / filename