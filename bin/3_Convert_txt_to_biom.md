# Converting .txt files to .biom format in order to use them in R with *phyloseq* tables

In case you need to edit your metadata after AMPtk pipeline you will want to go from .txt to .biom so they work with the *phyloseq* package in R. The instructions to do this are here. 

## Steps

### Prerequisites 

Install Python 3 if you have other version. The stable version of .biom format requires Python 3:

On macOS

1. Install Homebrew:

> ``ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)``

2. Check Homebrew installed correctly using:

`brew doctor`

3. Now install latest version of Python 3:

`brew install python 3`

In my case since Python 2.7 is already installed when using **pip** refer to **pip3** instead so it can call the command for Python 3.7 and not the older version


### Install biom-format package 

1. You need numpy in order to install biom, so run: 

`pip install numpy`

`pip install biom-format`

To work with BIOM 2.0+ files: 

`pip install h5py`

To see a list of all biom commands, run:

`biom`

If needed: to enable Bash tab completion of biom commands, add the following line to $HOME/.bash_profile (if on Mac OS X):

`eval "$(_BIOM_COMPLETE=source biom)"`

2. Once you have edited what you wanted in your .txt file, convert it back to .biom using:

`biom convert -i taxonomy.txt -o new_taxonomy.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy`


