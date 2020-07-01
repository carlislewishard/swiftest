Why Conda for our research group?

>> "Conda is an open-source package management system and environment management system that runs on Windows, macOS, and Linux. Conda quickly installs, runs, and updates packages and their dependencies. Conda easily creates, saves, loads, and switches between environments on your local computer. It was created for Python programs but it can package and distribute software for any language."

How to install miniconda3 ?

>>go to https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html

>>download the miniconda installer for macOS. In terminal do 
~ bash Miniconda3-latest-MacOSX-x86_64.sh
>>follow the prompt instructions

If you have upgraded to MacOS Catalina AND the conda init failed:
the new default shell is zsh. You will instead need to run source <path to conda>/bin/activate followed by conda init zsh.


Now, let's check that our python executable is running from the miniconda folder 
~ cd /Users/username/miniconda3/bin
~ ls 
>> check that python3.7 is listed, now execute python3.7
~ python3.7
>>> import sys
>>> sys.executable 
you should see one path only "Users/uxername/miniconda3/bin/python3.7

ok, now let's configure conda and actually setup our environment 
First, let's clean our path: 
>> from base:
~ conda config --set auto_activate_base False
This sets the auto config to false so when you are opening a new window, you dont start from conda environment, you have to activate the environment yourself, which makes you more aware of what you are doing and where you are working from 

What to do to activate python when you want to use it?
~ conda activate base

How to deactivate at the end? 
~ conda deactivate 

Pretty simple so far. It gets worse. 

Let's create a channel dedicated to our research group. From /Users/filename : 

~ conda create -n swiftest_conda
>> yes (to proceed)
Congrats, you create your own channel 
Activate your channel 
~ conda activate swiftest_conda 
Let's configure it with conda-forge, a conda channel maintained by a github community that has many packages available and work on fixing those packages as needed
~ conda config --env --add channels conda-forge 
Now we want to make sure that whatever package we decide to use, conda won't use anything that has a more recent update by default, which could lead to bugs. We have full control of the packages and the version we want to use (and we want everyone else to use). To do that: 
~ conda config --env --set channel_priority strict 
To see your channels: 
~ conda config --show channels

now, let's install some useful packages: 
MATPLOTLIB : 
~ conda install pandas scikit-learn matplotlib notebook

Now we have access to numpy, pyparsing, jupyter notebooks, etc. 
if we want, we could run a jupyter notebook now
~ jupyter notebook 
TADAM!


Now that you created your own little channel, you want to put the environment necessary to run python script  in this channel.


1) let's download the environment 
where you download the python script, you can create an environment out of the python files: 
~conda create --name <yourenvname> --file requirements.txt

2) let's put that new environment in our swiftest_conda channel
~ conda config --env --add channels swiftest_conda

Now let's say you notice I don't have in requirements.txt one package that you need to plot, you can add it in the environment and replace requirements.txt in the git. 

Now, you should be able to activate the channel swiftest_conda whenever you need python to plot from swiftest outputs. 

