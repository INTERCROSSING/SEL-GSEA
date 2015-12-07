import distutils.core
from distutils.core import setup

version_number = '0.1'
setup(name = 'gsea',
      version = version_number,
      description = 'calculate the significance of the pathways from SNP data',
      long_description = open('README.md').read(),
      author = 'Alexandra Vatsiou',
      author_email = 'alex.vatsiou@gmail.com',
      packages = ['gea','gsea'],
      requires = ['numpy',
                  'scipy',
                  'pylab'],
      )
