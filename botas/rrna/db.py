from importlib.resources import files
from pathlib import Path

def get_default_rrna_db() -> Path:
    """
    Return the path to the bundled default rRNA database.
    """
    return files("botas").joinpath("data/rrna_mini.fa.gz")
