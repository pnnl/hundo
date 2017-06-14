import io
from os.path import dirname, join
from setuptools import setup


long_description = open('README.rst').read()


def get_version(relpath):
    """Read version info from a file without importing it"""
    for line in io.open(join(dirname(__file__), relpath), encoding="cp437"):
        if "__version__" in line:
            if '"' in line:
                # __version__ = "0.9"
                return line.split('"')[1]
            elif "'" in line:
                return line.split("'")[1]


setup(
    name='hundo',
    version=get_version("hundo/__init__.py"),
    url='https://github.com/pnnl/hundo',
    license='MIT',
    author='Joe Brown',
    author_email='joe.brown@pnnl.gov',
    description='Amplicon processing protocol',
    long_description=long_description,
    packages=['hundo'],
    package_data={'': ['hundo/Snakefile',
                       ]},
    include_package_data=True,
    install_requires=[],
    entry_points={
          'console_scripts': [
              'hundo = hundo.hundo:cli'
          ]
    },
)
