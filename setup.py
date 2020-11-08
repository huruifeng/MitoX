import setuptools

setuptools.setup(
	name="mitox",
	version="0.1",
	description="exploring mitochondrial variants and gene expressions from single cell sequencing assays",
	author='Ruifeng Hu',
	author_email="huruifeng.cn@hotmail.com",
	packages=['mitox','mitox.data'],
	package_data={'mitox': ['data/human_mt.fasta','data/human_mt.fasta.fai','data/MusMT.fasta','data/MusMT.fasta.fai']},
	include_package_data=True,
	install_requires=['sklearn','seaborn','matplotlib','pysam','HTseq','anndata','pandas','numpy', 'scipy','umap-learn']
)