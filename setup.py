from setuptools import setup, find_packages

setup(
    name             = 'pacmill',
    version          = '0.1.0',
    description      = 'pacmill XYZ.',
    license          = 'MIT',
    url              = 'http://github.com/xapple/pacmill/',
    author           = 'Lucas Sinclair',
    author_email     = 'lucas.sinclair@me.com',
    packages         = find_packages(),
    install_requires = ['sh'],
    long_description = open('README.md').read(),
)
