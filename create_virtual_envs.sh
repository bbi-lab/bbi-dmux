#!/bin/bash

#
# Notes:
#   o  run in bbi-dmux directory
#   o  the Shendure cluster has nodes with Intel cpuid_level=11, cpuid_level=13,
#      and cpuid_level=22.
#      The cpuid_level 22 nodes have instructions (vectorized) that are not
#      part of the cpuid_level 11 architecture. So certain software built on
#      cpuid_level 22 nodes may not run on cpuid_level 11 nodes (I am guessing that
#      the cpuid_level 22 instructions are a superset of those on the cpuid_level
#      11 nodes.) This seems to affect at least the numpy module, if I remember
#      correctly. I suggest installing the modules using a cpuid_level=11 node
#      so that this pipeline runs on all Shendure cluster nodes.
#      The Trapnell cluster nodes are all AMD so this note is irrelevant for
#      this pipeline when installed and run on Trapnell cluster nodes.
#      In order to install this pipeline on a cpuid_level 11 node, use a shell
#      obtained with the command
#
#        qlogin -l mfree=16G -l cpuid_level=11
#
#   o  the pypy virtual environment is used only to run 'make_sample_fastqs.py'
#
#   o  there can be conflicts between site packages installed in the 'global'
#      pypy package and the site packages installed in the virtual environment.
#      I have seen this with the numpy and biopython packages.
#        o  the problem appears to happen when the pypy module is loaded - this
#           has the global site packages. The module must be loaded in order to
#           build the virtual environment and again when the virtual environment
#           is activated.
#        o  loading the pypy3 module makes the module site packages available.
#           One can see this using the 'pypy -m pip list' command.
#        o  building the virtual environment is not affected by the pypy
#           module site-packages. One can see this by looking at the
#           packages in src/pypy_env/lib/pypy*/site-packages.
#        o  running packages in the activated pypy virtual environment may
#           require that the pypy module be loaded because the executables
#           in the pypy virtual environment may require access to those
#           libraries. The libraries are not part of the virtual environment.
#        o  in the case that the required libraries are dynamically linked,
#           one can avoid having to load the pypy module by setting
#           the LD_LIBRARY_PATH environment variable to the directory in
#           the module path that contains the libraries. For example,
#             export LD_LIBRARY_PATH="/net/gs/vol3/software/modules-sw/pypy/3.8-7.3.9/Linux/Ubuntu22.04/x86_64/lib:$LD_LIBRARY_PATH"
#        o  in order to make the pypy_requirements.txt file using 'pip freeze',
#           use the LD_LIBRARY_PATH rather than loading the pypy module. This
#           prevents listing the global site packages in the pypy_requirements
#           file.
#
#            


echo "The virtual environment may depend on the CPU architecture."
echo "Clusters with mixed node architectures may fail, possibly"
echo "with Illegal Instruction core dumps when a python script"
echo "runs in a virtual environment. In this case, you may need" 
echo "to restrict the hardware resource to the architecture in"
echo "which the virtual environment was built. If the cluster has" 
echo "nodes of similar architecture but different generations," 
echo "one may be able to build the virtual environment on a node" 
echo "of the earliest generation."

echo
read -r -n1 -p "Press any key to continue: " key 
echo

#
# Prepare pypy virtual environment.
#
module purge
module load modules modules-init modules-gs
export DIR=$(dirname `readlink -f $0`)
source $DIR/load_pypy_env_reqs.sh

echo 'Cleaning cache directories...'
rm -rf ~/.cache/pip

#
# Remove existing pypy virtual environment in the background, if it exists.
#
if [ -d $DIR/bin/pypy_env ]; then
        echo 'Removing existing pypy virtualenv...'
        mv $DIR/bin/pypy_env $DIR/bin/pypy_env.tmp
        rm -rf $DIR/bin/pypy_env.tmp &
fi

#
# Drop PYTHONPATH environment variable to avoid conflicts.
#
export PYTHONPATH=''

# First, the pypy virtualenv
# It's probably safest to set up the virtual environment using
# pypy3 -m venv because there may be no 'virtualenv' command
# defined for pypy3 but there may be 'virtualenv's defined
# elsewhere.
#
echo 'Building pypy virtualenv...'
pypy3 -m venv $DIR/bin/pypy_env

if [ "$?" != 0 ]
then
  echo "Error: the virtualenv command returned an error."
  exit -1
fi

if [ ! -d $DIR/bin/pypy_env ]
then
  echo "Error: failed to make Python virtual environment in $DIR/bin/pypy_env."
  exit -1
fi

source $DIR/bin/pypy_env/bin/activate

pypy3 -m ensurepip
pip3 install -r $DIR/pypy_requirements.txt

deactivate

