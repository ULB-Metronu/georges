"""Installation with setuptools or pip."""
from setuptools import setup, find_packages
import os
import ast


def get_version_from_init():
    """Obtain library version from main init."""
    init_file = os.path.join(
        os.path.dirname(__file__), 'georges', '__init__.py'
    )
    with open(init_file) as fd:
        for line in fd:
            if line.startswith('__version__'):
                return ast.literal_eval(line.split('=', 1)[1].strip())


with open('README.md') as f:
    readme = f.read()

with open('COPYING') as f:
    lic = f.read()


setup(
    name='georges',
    version=get_version_from_init(),
    description='Georges: a Python interface for Zgoubi.',
    long_description=readme,
    author='Cédric Hernaslteens',
    author_email='cedric.hernalsteens@ulb.ac.be',
    url='https://github.com/chernals/georges',
    license=lic,
    packages=find_packages(exclude=('tests', 'docs', 'examples')),
    install_requires=[
        'georges-core',
        'matplotlib',
        'mypy',
        'numba',
        'numpy>=1.14.0',
        'pandas>=0.22.0',
        'pint',
        'plotly',
        'pyarrow',
        'pyyaml',
        'scipy>=1.0.0',
        'uproot',
    ],
    package_data={'georges': []},
)
