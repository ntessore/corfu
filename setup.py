import setuptools
from corfu import __version__, __author__, __email__

with open('README.rst', 'r') as f:
    long_description = f.read()

setuptools.setup(
    name='corfu',
    version=__version__,
    author=__author__,
    author_email=__email__,
    description='projected angular correlation functions',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url='https://github.com/ntessore/corfu',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    install_requires=[
        'numpy',
        'scipy',
    ],
    py_modules=[
        'corfu',
    ],
    extras_require={
        'docs': [
            'sphinxcontrib.bibtex',
        ]
    },
)
