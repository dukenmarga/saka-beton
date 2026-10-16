"""
Build the bundled app and zip it for distribution:

    python build.py

Same script on Linux and Windows; the zip is named after the platform it was
built on. Requires the dev dependencies (pip install -r requirements-dev.txt).
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

from version import __version__

root = Path(__file__).parent
name = f"saka-beton-{'windows' if os.name == 'nt' else 'linux'}-{__version__}"

subprocess.run(
    [sys.executable, "-m", "PyInstaller", "--noconfirm", "saka_beton.spec"],
    cwd=root,
    check=True,
)

archive = shutil.make_archive(
    str(root / "dist" / name), "zip", root_dir=root / "dist", base_dir="saka-beton"
)
print(f"\n{archive}  ({Path(archive).stat().st_size / 1024 ** 2:.0f} MB)")
