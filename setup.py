from setuptools import setup, find_packages

setup(
    name='tadtool',
    version='0.4',
    description='Assistant to find cutoffs in TAD calling algorithms.',
    packages=find_packages(exclude=["test"]),
    install_requires=[
        'numpy',
        'matplotlib',
        'progressbar2'
    ],
    url='https://github.com/vaquerizaslab/tadtool',
    download_url='https://github.com/vaquerizaslab/tadtool/tarball/0.4',
    keywords=['bioinformatics', 'hi-c', 'genomics', 'tad'],
    classifiers=[],
    scripts=['bin/tadtool']
)
