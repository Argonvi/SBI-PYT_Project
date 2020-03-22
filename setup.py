from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='complexconstructor',
      version='0.1',
      description='Generate macrocomplexes superimposing \
      paired interacting elements.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      keywords='macrocomplex bioinformatics structural pdb',
      url='https://github.com/Argonvi/SBI-PYT_Project',
      author='Paula Gomis Rosa, Arturo González Vilanova, Marta López Balastegui',
      author_email='marta.lopez01@estudiant.upf.edu, paula.gomis01@estudiant.upf.edu, arturo.gonzalez01@estudiant.upf.edu',
      packages=setuptools.find_packages(),
      install_requires=['biopython'],
      include_package_data=True)
