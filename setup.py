from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    lic = f.read()

setup(
    name='georges',
    version='2018.2',
    description='Georges',
    long_description=readme,
    author='Cédric Hernaslteens',
    author_email='cedric.hernalsteens@iba-group.com',
    url='https://github.com/chernals/georges',
    license=lic,
    packages=find_packages(exclude=('tests', 'docs', 'examples'))
)
