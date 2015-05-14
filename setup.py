from setuptools import setup, find_packages
import sys


def main():
    if 'develop' not in sys.argv:
        raise NotImplementedError("Use python setup.py develop.")
    setup(
        name="vs_utils",
        url='https://github.com/pandegroup/vs-utils',
        description='Virtual screening and cheminformatics utilities',
        packages=find_packages(),
    )

if __name__ == '__main__':
    main()
