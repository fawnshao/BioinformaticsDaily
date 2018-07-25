# Create the environment using Python 3 and add the IPython and Juypter system with notebook support.
conda create -n py3k python=3 ipython notebook
# now you need to activate the python 3 environment
source activate py3k

pip install bash_kernel
python -m bash_kernel.install
source deactivate py3k

# Remember to activate the python 3 env first!
source activate py3k
# start up the notebook server
jupyter notebook
# when you're finished deactivate the python 3 env so we can default to python 2.7
source deactivate py3k

##### in TACC
conda create -n py3cooler python=3 cooler
# conda install -c conda-forge -c bioconda cooler
source activate py3cooler
source deactivate py3cooler

conda create -n py36 python=3.6 anaconda
source activate py36
# pip install cooler
source deactivate py36