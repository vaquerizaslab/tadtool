from setuptools import setup, find_packages

setup(
    name='tadtool',
    version='0.4.0',
    description='Assistant to find cutoffs in TAD calling algorithms.',
    packages=find_packages(exclude=["test"]),
    install_requires=[
        'numpy',
        'matplotlib',
        'progressbar2'
    ],

    scripts=['bin/tadtool']
)
