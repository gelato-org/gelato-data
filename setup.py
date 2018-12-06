from setuptools import setup, find_packages


requires = [
    # list required third-party packages here
    'appdirs',
    'requests',
    'clldutils>=1.5.4',
    'pyglottolog',
    'tabulate',
    'tqdm',
    'bagit>=1.5.4',
    'pybtex',
    'openpyxl',
]

setup(
    name='pygelato',
    version='0.0',
    description='python package for the GELATO repository',
    long_description='',
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
    ],
    author='Robert Forkel',
    author_email='forkel@shh.mpg.de',
    url='',
    keywords='data linguistics genetics',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=requires,
    entry_points={
        'console_scripts': ['gelato=pygelato.__main__:main'],
    },
    tests_require=[],
    test_suite="pygelato")
