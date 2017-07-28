# sPEP-project

The bash script (pipeline.sh) will automatically run prediction softwares.
The following are required:
- Python 2.7
- Pyfasta (Python 2.7)
- Anchor
- SignalP
You also need an internet connection to get data predictions from CBS, Biomine and Jpred servers.

In the bash script, you can edit the paths to the Anchor and SignalP programs before running the pipeline.

To run the pipeline, open a terminal in the folder containing the script and use :
```sh main.sh -d {dirname} {filename.fasta}  ```

If you do not mention a directory in the -d option, the results will be saved in a new directory '{filename}_outputs' in your current directory.
The pipeline generates plots and a SQLite database accessed through the Python class SQLite.
