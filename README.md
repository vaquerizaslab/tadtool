# TADtool

TADtool is an interactive tool for the identification of meaningful parameters in TAD-calling algorithms for Hi-C data.

![TADtool main window](docs/images/tadtool.png)



## Quick start

Installation:

```bash
pip install tadtool
```

Run sample data from [GitHub repo](https://github.com/vaquerizaslab/tadtool/):

```bash
tadtool plot examples/chr12_20-35Mb.matrix.txt examples/chr12_20-35Mb_regions.bed chr12:31000000-33000000
```

This should open the interactive plotting window (see above). Start exploring by clicking in plots - the effects should be self-explanatory.

## Installation

You can install TADtool from the command line using PyPI

```bash
pip install tadtool
```

or download the source from our [GitHub repo](https://github.com/vaquerizaslab/tadtool/) and install manually

```bash
python setup.py install
```

This should install both the Python 2 package and a command-line executable called __tadtool__.

Test the installation running

```bash
tadtool -h
```

and you should see a brief help message.
