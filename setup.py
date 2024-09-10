from setuptools import setup, find_packages

# Load the README file as long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="cnvizard",
    version="0.2",
    description="A tool for visualizing germline copy number variants",
    long_description=long_description,
    long_description_content_type="text/markdown",  # If your README is in Markdown format
    author="Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft",
    author_email="jerkrause@ukaachen.de",
    url="https://github.com/IHGGM-Aachen/CNVizard",
    packages=find_packages(include=["cnvizard"]),
    install_requires=[
        "streamlit",
        "pandas",
        "numpy",
        "pyarrow",
        "fastparquet",
        "plotly",
        "python-dotenv",
        "xlsxwriter",
        "seaborn",
        "cyvcf2",
        "pillow",
    ],
    include_package_data=True,
    package_data={
        'cnvizard': [
            'resources/*.txt',
            'resources/candidate_lists/*.txt',
            'resources/references/*.parquet'
        ]
    },
    entry_points={
        "console_scripts": [
            "cnvizard=cnvizard.run:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.12.4",
)
