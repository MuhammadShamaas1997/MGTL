# ~/.profile: executed by the command interpreter for login shells.
# This file is not read by bash(1), if ~/.bash_profile or ~/.bash_login
# exists.
# see /usr/share/doc/bash/examples/startup-files for examples.
# the files are located in the bash-doc package.

# the default umask is set in /etc/profile; for setting the umask
# for ssh logins, install and configure the libpam-umask package.
#umask 022

# if running bash
if [ -n "$BASH_VERSION" ]; then
    # include .bashrc if it exists
    if [ -f "$HOME/.bashrc" ]; then
	. "$HOME/.bashrc"
    fi
fi

# set PATH so it includes user's private bin directories
PATH="$HOME/bin:$HOME/.local/bin:$PATH"
export PYTHONPATH="/usr/local/lib/python3.5/site-packages"
export LD_LIBRARY_PATH="$HOME/install/bin:$HOME/install/lib:$HOME/install/share:$HOME/install/include:$HOME/install/h5utils/bin:$HOME/install/h5utils/lib:$HOME/install/h5utils/share:$HOME/install/h5utils/include:$HOME/install/harminv/bin:$HOME/install/harminv/lib:$HOME/install/harminv/share:$HOME/install/harminv/include:$HOME/install/hdf5-1.10.5/bin:$HOME/install/hdf5-1.10.5/lib:$HOME/install/hdf5-1.10.5/include:$HOME/install/hdf5-1.10.5/share:$HOME/install/libctl/bin:$HOME/install/libctl/lib:$HOME/install/libctl/include:$HOME/install/libctl/share:$HOME/install/libGDSII/include:$HOME/install/libGDSII/lib:$HOME/install/libGDSII/bin:$HOME/install/libGDSII/share:$HOME/install/mpb/include:$HOME/install/mpb/share:$HOME/install/mpb/bin:$HOME/install/mpb/lib:$HOME/install/meep/lib:$HOME/install/meep/bin:$HOME/install/meep/share:$HOME/install/meep/include:$LD_LIBRARY_PATH"
