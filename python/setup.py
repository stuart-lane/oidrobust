from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="oidrobust",
    version="0.1.0",
    author="Stuart Lane",
    author_email="stuart.lane@bristol.ac.uk",
    description="A package for robust overidentification testing in linear IV models.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/stuart-lane/oidrobust",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Statistics",
        "Development Status :: 4 - Beta", 
        "Intended Audience :: Science/Research",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.20.0", 
        "scipy>=1.7.0",
        "pandas>=1.3.0",
        "patsy>=0.5.2",
        "typing-extensions>=4.0.0", 
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=22.0",
            "flake8>=4.0",
            "mypy>=0.900",  
            "isort>=5.0",   
        ],
        "docs": [
            "sphinx>=4.0",
            "jupyter>=1.0",
            "sphinx-rtd-theme>=1.0", 
            "sphinx-autodoc-typehints>=1.12", 
        ],
    },
    package_data={
        "oidrobust": ["py.typed"], 
    },
    zip_safe=False,
)