from setuptools import setup, find_packages

setup(
    name             = 'pacmill',
    version          = '0.1.1',
    description      = 'The `pacmill` python package is a bioinformatics '
                       'pipeline that is developed to process microbial '
                       '16S amplicon sequencing data.',
    license          = 'MIT',
    url              = 'http://github.com/xapple/pacmill/',
    author           = 'Lucas Sinclair',
    author_email     = 'lucas.sinclair@me.com',
    packages         = find_packages(),
    install_requires = ['plumbing>=2.8.1', 'autopaths>=1.4.2', 'fasta>=2.0.4',
                        'pymarktex>=1.4.0', 'biopython', 'pandas', 'sh',
                        'bcbio-gff'],
    long_description = open('README.md').read(),
    long_description_content_type = 'text/markdown',
    include_package_data = True,
)