import sys
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

from pipeline import main
from presets import FULL_LENGTH_BAM_DEFAULTS


if __name__ == "__main__":
    main(FULL_LENGTH_BAM_DEFAULTS)
