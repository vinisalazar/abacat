import setuptools

setuptools.setup(
    name="abacat",
    version="0.0.1a",
    author="Vini Salazar",
    author_email="viniws@gmail.com",
    description="abacat - A BActerial genome Curation and Annotation Toolkit",
    url="https://github.com/vinisalazar/abacat",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.6"
    ],
    packages=setuptools.find_packages(),
    include_package_data=True,
    python_requires=">=3.6",
    install_requires=[
        "biopython",
        "pandas",
        "argparse",
        "matplotlib",
        "scipy",
        "pytest"
    ]
)
