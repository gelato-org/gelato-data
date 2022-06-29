from setuptools import setup


setup(
    name='cldfbench_gelato',
    py_modules=['cldfbench_gelato'],
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'cldfbench.dataset': [
            'gelato=cldfbench_gelato:Dataset',
        ]
    },
    install_requires=[
        'cldfbench',
        'csvw',
        'pycldf',
    ],
    extras_require={
        'test': [
            'pytest-cldf',
        ],
    },
)
