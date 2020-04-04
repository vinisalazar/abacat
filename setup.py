import setuptools

setuptools.setup(
    name="abacat",
    version="0.0.2a",
    author="Vini Salazar",
    author_email="viniws@gmail.com",
    description="abacat - A BActerial genome Curation and Annotation Toolkit",
    long_description="Abacat (pronounced 'ABBA-cat') is a toolkit for working with bacterial whole genome sequencing (WGS) data. It provides Python objects to represent elements which are common to WGS analysis workflows, such as coding sequence (CDS) files, containing genes or proteins, alignment methods, or statistics about your sequences.",
    url="https://github.com/vinisalazar/abacat",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.6",
    ],
    packages=setuptools.find_packages(),
    scripts=[
        "abacat/prodigal.py",
        "abacat/pipelines/annotate.py",
        "abacat/pipelines/phenotyping.py",
        "abacat/deprecated/prokka.py",
    ],
    include_package_data=True,
    python_requires=">=3.6",
    install_requires=[
        "biopython",
        "pandas",
        "argparse",
        "matplotlib",
        "scipy",
        "pytest",
    ],
)
