from setuptools import setup, find_packages

setup(
    name             = 'pacmill',
    version          = '0.3.4',
    description      = 'The `pacmill` python package is a bioinformatics '
                       'pipeline that is developed to process microbial '
                       '16S amplicon sequencing data.',
    license          = 'MIT',
    url              = 'http://github.com/xapple/pacmill/',
    author           = 'Lucas Sinclair',
    author_email     = 'lucas.sinclair@me.com',
    packages         = find_packages(),
    install_requires = ['plumbing>=2.8.9', 'autopaths>=1.4.6', 'fasta>=2.1.6',
                        'pymarktex>=1.4.6', 'seqsearch>=1.2.7' 'biopython',
                        'pandas', 'sh', 'tag', 'shell_command', 'tabulate'],
    long_description = open('README.md').read(),
    long_description_content_type = 'text/markdown',
    include_package_data = True,
)