from setuptools import setup, find_packages

try:
    import pypandoc
    long_description = pypandoc.convert_file("README.md", "rst")
except ImportError:
    print ("Warning: pypandoc module not found, could not convert"
           " Markdown to reStructuredText." )
    long_description = ""


setup(
    name='pymadx',
    version='1.0.0',
    packages=find_packages(exclude=["docs", "tests", "obsolete"]),
    # Not sure how strict these need to be...
    install_requires=["matplotlib >= 1.7.1",
                      "numpy >= 1.3.0"],
    # Some version of python2.7
    python_requires="==2.7.*",

    author='JAI@RHUL',
    author_email='laurie.nevay@rhul.ac.uk',
    description="Write MADX models and load MADX output.",
    long_description=long_description,
    url='https://bitbucket.org/jairhul/pymadx/',
    license='GPL3',
    keywords='madx accelerator twiss ptc',
)
