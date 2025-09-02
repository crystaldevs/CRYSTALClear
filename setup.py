import setuptools

long_description = 'Python framework for the <a href="https://www.crystal.unito.it/index.php">CRYSTAL code</a>. Forked from <a href="https://github.com/crystal-code-tools/CRYSTALpytools">CRYSTALpytools</a>.'

setuptools.setup(
    name="CRYSTALClear",
    version="0.2.10",
    author_email="crystalunito@gmail.com",
    description="Python framework for the CRYSTAL code.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/crystaldevs/CRYSTALClear",
    project_urls={
        "Bug Tracker": "https://github.com/crystaldevs/CRYSTALClear/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(
        include=['CRYSTALClear', 'CRYSTALClear.*']),
    # python_requires=">=3.8",
    install_requires=[
        "pymatgen",
        "mendeleev",
        "ase"
    ]
)
