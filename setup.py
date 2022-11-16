import os
from setuptools import setup, find_packages, Command


__version__ = None
exec(open('tadtool/version.py').read())

class CleanCommand(Command):
    """
    Custom clean command to tidy up the project root.
    """
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info ./htmlcov')


setup(
    name='tadtool',
    version=__version__,
    description='Assistant to find cutoffs in TAD calling algorithms.',
    packages=find_packages(exclude=["test"]),
    install_requires=[
        'numpy>=1.9.0',
        'matplotlib>=3.6.0',
        'progressbar2',
        'future',
    ],
    author='Vaquerizas lab',
    author_email='kai.kruse@mpi-muenster.mpg.de',
    url='https://github.com/vaquerizaslab/tadtool',
    keywords=['bioinformatics', 'hi-c', 'genomics', 'tad'],
    classifiers=[],
    scripts=['bin/tadtool'],
    cmdclass={
        'clean': CleanCommand
    },
)
