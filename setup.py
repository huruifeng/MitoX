import setuptools

setuptools.setup(
	name="mitox",
	version="0.1",
	description="exploring mitochondrial variants and gene expressions from single cell sequencing assays",
	author='Ruifeng Hu',
	author_email="huruifeng.cn@hotmail.com",
	packages=['mitox','mitox.data'],
	install_requires=['sklearn','seaborn','matplotlib','pysam','HTseq','anndata','pandas','numpy', 'scipy','umap-learn']
)