Installation process
--------------------

Download the package folder and place it in your prefered location :open_file_folder:. Then, modify the Python Module Search Path in your shell startup file as follows:

#### for Bash shell (sh):

In your startup file (**`~/.bash_profile`** , **`~/.bash_login`** or **`~/.bashrc`** [Linux users], and **`~/.profile`** [Mac OS users]) include:

```bash
export PYTHONPATH=${PYTHONPATH}:/path/to/package
```

#### for C Shell (csh; tcsh):

In your startup file (**`~/.cshrc`** , **`~/.tcshrc`** or **`~/.login`**) include:

```csh
setenv PYTHONPATH $PYTHONPATH\:/path/to/package
```

### External required packages

* [Matplotlib](https://matplotlib.org/users/installing.html)
* [Numpy](https://www.scipy.org/install.html)
* [Pandas](http://pandas.pydata.org/pandas-docs/stable/install.html) :panda_face:
* [Astropy](http://docs.astropy.org/en/stable/install.html) (optional, but recommended)
* [IPython](https://ipython.org/install.html) (optional, but recommended)
