from setuptools import setup, find_packages

setup(
    name='pytransport',
    version='0.1',
    packages=find_packages(exclude=["docs", "tests", "obsolete"]),
    # Not sure how strict these need to be...
    install_requires=["matplotlib",
                      "numpy",
                      "scipy"],
    # Some version of python2.7
    python_requires="==2.7.*",

    author='JAI@RHUL',
    author_email='stewart.boogert@rhul.ac.uk',
    description="Convert TRANSPORT models and load TRANSPORT output.",
    url='https://bitbucket.org/jairhul/pytransport/'
)
