import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ion_automation",
    version="1.0.2",
    author="Kelsey Zhu",
    author_email="khzhu@@users.noreply.github.com",
    description="Automated workflow for Oncomine Precision Assay NGS sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/khzhu/ion_automation",
    project_urls={
        "Bug Tracker": "https://github.com/khzhu/ion_automation/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)