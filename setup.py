import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="package", 
    version="0.1.1",
    author="Renee Oles",
    author_email="roles@health.ucsd.edu",
    description="Package: description",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Package",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Development Status :: 1 - Beta"],
    python_requires='>=3.5',
    install_requires=[
        "setuptools"],
    entry_points = {
        'console_scripts': ['package = package.main:cli']
    },
    include_package_data=True,
    zip_safe=False,
    package_data = {'package': ['workflow/*']}
)
