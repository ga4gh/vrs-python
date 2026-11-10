# Explore with notebooks

The walkthroughs below are rendered examples, so you can read them directly in
your browser without installing anything. They demonstrate how to use VRS-Python
to build variation models, compute identifiers, translate variant descriptions,
and access sequence data. You can also download a notebook from its page and run
it in Jupyter. Pages show saved outputs, which may differ from a fresh run against
your installed release and configured data services.

| Notebook | What you will learn |
| --- | --- |
| [Quick start](1_Quick_Start.md) | Translate a variant using sequence services |
| [SeqRepo data proxy](2_Exploring_the_SeqRepo_DataProxy.md) | Retrieve sequences and resolve aliases |
| [Basic models](3_Basic_Models.md) | Build alleles and compute identifiers |
| [Allele translator](4_Exploring_the_AlleleTranslator.md) | Translate supported variant nomenclatures |
| [CNV translator](5_Exploring_the_CnvTranslator.md) | Represent copy number variation |
| [Upcoming features](6_Upcoming_features.md) | Explore development examples, subject to change |
| [Uncertain HGVS ranges](7_Uncertain_Range_HGVS.md) | Work with uncertain coordinates |

!!! note "Data access"

    Some notebooks use public demonstration endpoints for sequence and transcript
    data. These endpoints are intended for exploring the examples. For heavy
    workloads or production use, we recommend running your own local services.
    Follow [SeqRepo and UTA](../data-services.md) to set them up, then update the
    notebooks to use your endpoints. Inspect setup cells before running them,
    since they may override environment variables you set in the shell.

## Run the notebooks with MyBinder

[Binder](https://mybinder.org) provides a ready-to-use computing environment for running
the notebooks directly in your browser, with no local installation required.

You can launch the notebooks on Binder [here](https://mybinder.org/v2/gh/ga4gh/vrs-python/v3?labpath=notebooks).

Once launched, select `VRS_kernel` as the notebook kernel.

## Run the notebooks locally

Clone the [VRS-Python repository](https://github.com/ga4gh/vrs-python), then follow the
[Installation steps](../installation.md) to create the virtual environment (`venv`)
and activate it.

!!! note "Notebook dependency"

    If you did not install the `notebooks` optional dependency during setup,
    run this command from the root of the cloned project (`vrs-python/`),
    with the virtual environment activated:

    ```shell
    python3 -m pip install "ga4gh.vrs[extras,notebooks]"
    ```

### Open a notebook

We recommend either of the following options for opening a notebook.

#### JupyterLab

1. In the terminal, activate the project's `venv` virtual environment:

         source venv/bin/activate

2. In the activated virtual environment, start [JupyterLab](https://docs.jupyter.org/en/latest/):

         python3 -m jupyterlab

3. JupyterLab usually opens automatically in a new browser window. If it does not, open
   the URL shown in the terminal.
4. Open a notebook from the `notebooks/` directory.

#### VS Code

1. Install the [Python extension](https://marketplace.visualstudio.com/items?itemName=ms-python.python)
   and [Jupyter extension](https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter).
2. Open the project root in VS Code and open a notebook from the
   `notebooks/` directory.
3. Select **Select Kernel** in the upper-right corner and choose the Python
   interpreter from the repository's `venv` virtual environment.
