from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("botas-rnaseq")
except PackageNotFoundError:
    __version__ = "0.1.5"
