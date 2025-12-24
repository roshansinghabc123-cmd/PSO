Installation
============

Quick Start
-----------

To get started with the **π-Stack Optimizer**, follow these steps to set up your environment and dependencies.

1. **Clone the Repository**

   .. code-block:: bash

      git clone https://github.com/roshansinghabc123-cmd/PSO.git
      cd PSO

2. **Set Up a Virtual Environment**

   It is recommended to use a virtual environment to manage dependencies.

   .. code-block:: bash

      python3 -m venv venv
      source venv/bin/activate

3. **Install Python Dependencies**

   Install the required packages using pip.

   .. code-block:: bash

      pip install -r requirements.txt

4. **Activate CLI Wrappers**

   The repository includes a helper script to set up command-line aliases for ease of use.

   .. code-block:: bash

      source ./activate_pi_stack.sh

   This will make commands like ``pi-stack-generator`` and ``pi-hyperopt`` available in your shell.

xTB Backend Configuration
-------------------------

The optimizer relies on **xTB** (Extended Tight-Binding) for fast quantum chemical energy evaluations. You must have xTB installed and properly configured.

1. **Install xTB**
   
   If you rely on conda:

   .. code-block:: bash

      conda install -c conda-forge xtb

2. **Environment Variables**

   Ensure ``xTB`` is in your PATH. You can also configure the number of threads for parallel calculation:

   .. code-block:: bash

      export OMP_NUM_THREADS=4

Verification
------------

To verify your installation, run a simple test optimization:

.. code-block:: bash

   pi-stack-generator example.xyz --workers 1 --optimizer pso

If the command runs without errors and produces output, your installation is complete.
