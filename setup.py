from setuptools import setup, find_namespace_packages

setup(
    name             = 'pacmill',
    version          = '0.6.1',
    description      = 'The `pacmill` python package is a bioinformatics '
                       'pipeline that is developed to process microbial '
                       '16S amplicon sequencing data.',
    license          = 'MIT',
    url              = 'https://github.com/xapple/pacmill/',
    author           = 'Lucas Sinclair',
    author_email     = 'lucas.sinclair@me.com',
    packages         = find_namespace_packages(),
    install_requires = ['plumbing>=2.10.4', 'autopaths>=1.5.2', 'fasta>=2.2.13',
                        'pymarktex>=1.4.6', 'seqsearch>=2.1.6', 'biopython',
                        'pandas', 'sh', 'tag', 'shell_command', 'tabulate',
                        'openpyxl', 'mock', 'matplotlib', 'regex', 'p_tqdm',
                        'wget', 'brewer2mpl', 'crest4', 'ftputil', 'rpy2',
                        'pystache'],
    long_description = open('README.md').read(),
    long_description_content_type = 'text/markdown',
    include_package_data = True,
)