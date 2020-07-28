from setuptools import setup, find_packages

setup(
    name             = 'pacmill',
    version          = '0.2.0',
    description      = 'The `pacmill` python package is a bioinformatics '
                       'pipeline that is developed to process microbial '
                       '16S amplicon sequencing data.',
    license          = 'MIT',
    url              = 'http://github.com/xapple/pacmill/',
    author           = 'Lucas Sinclair',
    author_email     = 'lucas.sinclair@me.com',
    packages         = find_packages(),
    install_requires = ['plumbing>=2.8.7', 'autopaths>=1.4.6', 'fasta>=2.0.8',
                        'pymarktex>=1.4.4', 'seqsearch>=1.2.3' 'biopython',
                        'pandas', 'sh', 'tag', 'shell_command', 'tabulate'],
    long_description = open('README.md').read(),
    long_description_content_type = 'text/markdown',
    include_package_data = True,
)