from setuptools import Extension, setup

ext_modules = [
    Extension(
        "botas.quantify._core_fast",
        sources=["botas/quantify/_core_fast.c"],
        optional=True,
    )
]

setup(ext_modules=ext_modules)