from setuptools import setup, find_packages
import os, ast


def get_version_from_init():
    init_file = os.path.join(
        os.path.dirname(__file__), 'georges', '__init__.py'
    )
    with open(init_file, 'r') as fd:
        for line in fd:
            if line.startswith('__version__'):
                return ast.literal_eval(line.split('=', 1)[1].strip())


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    lic = f.read()


setup(
    name='georges',
    version=get_version_from_init(),
    description='Georges',
    long_description=readme,
    author='Cédric Hernaslteens',
    author_email='cedric.hernalsteens@iba-group.com',
    url='https://github.com/chernals/georges',
    license=lic,
    packages=find_packages(exclude=('tests', 'docs', 'examples')),
    install_requires=[
        'numpy>=1.14.0',
        'pandas>=0.22.0',
        'scipy>=1.0.0',
        'matplotlib>=2.1.1',
        'jinja2>=2.9.6',
        'xlrd>=1.1.0',
        'lmfit>=0.9.7',
        'm2r>=0.1.13',
        'deap>=1.2.2',
        'xlrd>=1.1',
    ],
    data_files=[('bin/', ['bin/madx'])],  # Install MAD-X in f"{sys.prefix}/bin/madx"
)
