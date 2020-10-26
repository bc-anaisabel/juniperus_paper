# Converting .txt files back to .biom format for phyloseq tables

The stable version of biom format requires Python 3: 


## 1. Install Python 3 on MacOS:

Install Homebrew:

> ``ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)``

Check Homebrew installed correctly using:

`brew doctor`

Now install latest version of Python 3:

`brew install python 3`

Since Python 2.7 is already installed when using **pip** refer to **pip3** instead so it can call the command for Python 3.7 and not the older version

## 2. Install biom-format package 

`pip install numpy`

`pip install biom-format`

To work with BIOM 2.0+ files: 

`pip install h5py`

To see a list of all biom commands, run:

`biom`

If needed: to enable Bash tab completion of biom commands, add the following line to $HOME/.bash_profile (if on Mac OS X):

`eval "$(_BIOM_COMPLETE=source biom)"`

Once file has been edited in Excel, to convert from .txt to .biom:

`biom convert -i taxonomy.txt -o new_taxonomy.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy`


