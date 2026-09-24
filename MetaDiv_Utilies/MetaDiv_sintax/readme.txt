
## Docker requirement

MetaDiv_sintax runs VSEARCH/SINTAX through Docker.

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/).
2. Start Docker Desktop before running the notebook.
3. Open the MetaDiv_sintax script and set the mode (16S, ITS or CO1) 
4. Indicate in the script the name of the reference database to use
5. Run the jupyter notebook cells


MetaDiv uses the following PipeCraft VSEARCH image:

```bash
pipecraft/vsearch:2.30.4-pc1.2.0
