from gogetem import __version__
from setuptools import setup, find_packages

with open("requirements.txt", "r") as f:
    requirements = f.readlines()

setup(
    name="gogetem",
    version=__version__,
    author="Dillon Barker",
    author_email="dillon.barker@phac-aspc.gc.ca",
    description="Download sequence data for Gene Ontology terms",
    url="https://github.com/dorbarker/gogetem",
    packages=find_packages(),
    install_requires=[req for req in requirements if req[:2] != "# "],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={"console_scripts": ["gogetem=gogetem.gogetem:main"]},
)
