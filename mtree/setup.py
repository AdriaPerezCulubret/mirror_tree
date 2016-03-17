from setuptools import setup

def readme():
  with open('README.md') as r:
    return r.read()


setup(name='mtree',
      version='1.0',
      description='Predict protein interactions by the mirror tree approach',
      long_description=readme(),
      url='https://github.com/joanmarticarreras/mirror_tree/tree/master/mtree',
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