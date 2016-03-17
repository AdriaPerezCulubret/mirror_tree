from setuptools import setup

setup(name='mtree',
      version='1.0',
      description='Predict protein interactions by the mirror tree approach',
      url='http://github.com/storborg/funniest',
      author='S. Castillo-Lara, J. Martí, A. Pérez ',
      author_email='adriperez24@gmail.com',
      license='GPL 3.0',
      packages=['mtree'],
      scripts=['bin/mtree'],
      install_requires=[
          'numpy',
          'biopython',
          'scipy',
          'argparse',
      ],
      zip_safe=False)