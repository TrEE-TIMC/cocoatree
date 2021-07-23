from setuptools import setup, find_packages
from akasthesia import __version__

long_description = "Long description of the awesome coevolution package"

setup(
    name='akasthesia',
    version=__version__,
    description='Awesome coevolution stuff',  # Optional
    long_description=long_description,  # Optional
    long_description_content_type='text/markdown',  # Optional (see note above)
    url='https://github.com/FIXME',  # Optional
    author='William Schmitt',
    author_email='FIXME',
    classifiers=[  # Optional
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers, Scientists',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],

    keywords='coevolution, MSA',
    packages=find_packages(),  # Required
    python_requires='>=3.8, <4',
    install_requires=['numpy'],  # Optional
    extras_require={  # Optional
        'dev': ['flake8'],
        'test': ['pytest'],
    },

    project_urls={  # Optional
        'Bug Reports': 'https://github.com/FIXME',
        'Source': 'https://github.com/FIXME/',
    },
)
