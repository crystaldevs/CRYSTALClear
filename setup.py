import pathlib

import setuptools

# PyPI shows the README as the project page.
long_description = (pathlib.Path(__file__).parent / "README.md").read_text(encoding="utf-8")

setuptools.setup(
    name="CRYSTALClear",
    version="0.2.18",
    author="CRYSTAL devs",
    author_email="crystalunito@gmail.com",
    description="Python framework for the CRYSTAL code.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/crystaldevs/CRYSTALClear",
    project_urls={
        "Documentation": "https://crystaldevs.github.io/CRYSTALClear/",
        "Bug Tracker": "https://github.com/crystaldevs/CRYSTALClear/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(
        include=['CRYSTALClear', 'CRYSTALClear.*']),
    python_requires=">=3.9",
    install_requires=[
        "pymatgen",
        "mendeleev",
        "ase"
    ]
)
