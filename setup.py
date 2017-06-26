from setuptools import setup, find_packages

setup(
    name='tadtool',
    version='0.71',
    description='Assistant to find cutoffs in TAD calling algorithms.',
    packages=find_packages(exclude=["test"]),
    install_requires=[
        'numpy>=1.9.0',
        'matplotlib',
        'progressbar2',
        'future'
    ],
    author='Vaquerizas lab',
    author_email='kai.kruse@mpi-muenster.mpg.de',
    url='https://github.com/vaquerizaslab/tadtool',
    download_url='https://github.com/vaquerizaslab/tadtool/tarball/0.61',
    keywords=['bioinformatics', 'hi-c', 'genomics', 'tad'],
    classifiers=[],
    scripts=['bin/tadtool']
)
