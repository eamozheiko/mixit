from setuptools import setup, find_packages

setup(
    name="mixit",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[],
    entry_points={
        'console_scripts': [
            'mixit=mixit.main:main',
        ],
    },
    author="Your Name",
    author_email="your.email@example.com",
    description="Mix your DNA Simulator 2000",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/eamozheiko/mixit",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
) 